
#include <fstream>
#include <map>
#include <iomanip>
#include <mat_ops.H>
#include <SimpleTrimesh.H>

using namespace ASCbase;

const std::string SimpleTrimesh::A_fname = "SimpleTrimesh.C";

SimpleTrimesh::SimpleTrimesh(const std::string fname, bool check_vertex_order)
{
  init();

  // Check if fname ends with .vert or .face
  std::string last_five = fname.substr(fname.length() - 5);
  std::string tmp_fname = fname;
  if(last_five == ".face" || last_five == ".vert"){
    tmp_fname = fname.substr(0, fname.length() - 5);
  }

  if(!read_vert_file(tmp_fname) || 
     !read_face_file(tmp_fname, check_vertex_order)){
    A_fail = true;
    err_msg(A_fname, "SimpleTrimesh()", "reading of surface files failed\n");
    return;
  }
  compute_areas();
}

SimpleTrimesh::SimpleTrimesh(SimpleTrimesh& src, BoundingVolume &bvol)
{
  if(&src == this) return;

  init();
  src.intersect(bvol, &A_num_verts, &A_vertices, &A_normals, &A_triangles);
}

SimpleTrimesh::SimpleTrimesh(SimpleTrimesh& src, const my_float_t* center,
                             const my_float_t radius)
{
  if(&src == this) return;

  init();
  src.clip(center, radius, &A_num_verts, &A_vertices, &A_normals, &A_triangles);
}

SimpleTrimesh::~SimpleTrimesh()
{
  if(A_vertices) delete [] A_vertices;
  if(A_normals) delete [] A_normals;
  if(A_original_vertices) delete [] A_original_vertices;
  if(A_original_normals) delete [] A_original_normals;
  if(A_vert_octree) delete A_vert_octree;
  init();
}

void
SimpleTrimesh::init()
{
  A_vertices = 0;
  A_normals = 0;
  A_original_vertices = 0;
  A_original_normals = 0;
  A_num_verts = 0;
  A_vert_octree = 0;
  A_fail = false;
}

void
SimpleTrimesh::transform(const my_float_t* R, const my_float_t* T)
{
  size_t arr_len = 3*A_num_verts;

  my_float_t *tmp_verts = new my_float_t[arr_len];
  my_float_t *tmp_normals = new my_float_t[arr_len];
  std::copy(A_vertices, A_vertices + arr_len, tmp_verts);
  std::copy(A_normals, A_normals + arr_len, tmp_normals);

  my_gemm(A_num_verts, 3, 3, 1.0, tmp_normals, 3, R, 3, A_normals, 3, 0.0);
  my_float_t *t = A_vertices;
  for(size_t i = 0; i < A_num_verts; ++i, t+=3) std::copy(T, T + 3, t);
  my_gemm(A_num_verts, 3, 3, 1.0, tmp_verts, 3, R, 3, A_vertices, 3, 1.0);

  if(!A_original_vertices){
    A_original_vertices = tmp_verts;
    A_original_normals = tmp_normals;
  }else{
    delete [] tmp_verts;
    delete [] tmp_normals;
  }
}

void
SimpleTrimesh::inverse_transform(const my_float_t* R, const my_float_t* T)
{
  size_t arr_len = 3*A_num_verts;

  my_float_t *tmp_normals = new my_float_t[arr_len];
  std::copy(A_normals, A_normals + arr_len, tmp_normals);

  my_float_t R_transpose[9];
  for(int i = 0; i < 3; ++i)
    for(int j = 0; j < 3; ++j) R_transpose[3*i +j] = R[3*j + i];
  my_gemm(A_num_verts, 3, 3, 1.0, tmp_normals, 3, R_transpose, 3, A_normals, 3,
          0.0);

  if(!A_original_vertices){
    A_original_normals = tmp_normals;
    A_original_vertices = new my_float_t[arr_len];
    std::copy(A_vertices, A_vertices + arr_len, A_original_vertices);
  }else delete [] tmp_normals;

  my_float_t *tmp_verts = new my_float_t[arr_len];
  my_float_t T_inverse[3]; // Additive inverse :)
  for(int i = 0; i < 3; ++i) T_inverse[i] = -1.0*T[i];
  my_float_t *t = tmp_verts;
  for(size_t i = 0; i < A_num_verts; ++i, t+=3) 
    std::copy(T_inverse, T_inverse + 3, t);
  my_axpy(arr_len, 1.0, A_vertices, 1, tmp_verts, 1);

  my_gemm(A_num_verts, 3, 3, 1.0, tmp_verts, 3, R_transpose, 3, A_vertices, 3,
          0.0);

  delete [] tmp_verts;
}
                      

void
SimpleTrimesh::revert()
{
  if(A_original_vertices){
    size_t arr_len = 3*A_num_verts;
    std::copy(A_original_vertices, A_original_vertices + arr_len, A_vertices);
    std::copy(A_original_normals, A_original_normals + arr_len, A_normals);
  }
}

my_float_t
SimpleTrimesh::get_total_SA() const
{
  my_float_t rv = 0;
  std::vector<triangle_t>::const_iterator face;
  for(face = A_triangles.begin(); face < A_triangles.end(); ++face)
    rv += face->area; 
  return rv;
}

#if 0
bool 
SimpleTrimesh::old_compare(const SimpleTrimesh& target, 
                           const my_float_t max_pt_dist, size_t *num_points,
                           my_float_t* RMSE, size_t *num_faces, 
                           my_float_t *area) const
{
  *area = 0.0;
  *num_faces = 0;
  *num_points = 0;
  *RMSE = -1.0;

  const my_float_t neg_max_pt_dist = -1.0 * max_pt_dist;
  corr_map correspondances;
  const my_float_t* query_pt = A_vertices;
  for(size_t i = 0; i < A_num_verts; ++i, query_pt += 3){
    std::vector<triangle_t>::const_iterator targ_face = target.faces_begin();
    for( ; targ_face != target.faces_end(); ++targ_face){
      pt_plus_dist pt_info; 
      my_float_t d;
 
      // Check if corresponding point is on a face
      bool rv = corresponding_point(query_pt, targ_face->vertices[0],
                                    targ_face->vertices[1], 
                                    targ_face->vertices[2],
                                    &d, pt_info.corr_pt);
      if(rv && neg_max_pt_dist <= d && d <= max_pt_dist){
        pt_info.distance = d;
        corr_map::iterator cspond = correspondances.find(query_pt);
        if(cspond == correspondances.end())
          correspondances[query_pt] = pt_info;
        else{
          my_float_t abs_d = d;
          if(abs_d < 0) abs_d *= -1.0;
          my_float_t abs_saved_d = cspond->second.distance;
          if(abs_saved_d < 0) abs_saved_d *= -1.0;
          if(abs_d < abs_saved_d) cspond->second = pt_info;
        }
      }else{
        for(uint j = 0; j < 3; ++j){
          rv = corresponding_point(query_pt, targ_face->vertices[j],
                                   targ_face->vertices[(j+1) % 3],
                                   &d, pt_info.corr_pt);

          if(neg_max_pt_dist <= d && d <= max_pt_dist){
            pt_info.distance = d;
            corr_map::iterator cspond = correspondances.find(query_pt);
            if(cspond == correspondances.end())
              correspondances[query_pt] = pt_info;
            else{
              my_float_t abs_d = d;
              if(abs_d < 0) abs_d *= -1.0;
              my_float_t abs_saved_d = cspond->second.distance;
              if(abs_saved_d < 0) abs_saved_d *= -1.0;
              if(abs_d < abs_saved_d) cspond->second = pt_info;
            }
          }
        }
      }
    }
  }

  // Calculate RMSE of corresponding points
  *num_points = correspondances.size();
  corr_map::const_iterator cpi;
  my_float_t my_sum = 0.0;
  for(cpi = correspondances.begin(); cpi != correspondances.end(); ++cpi)
    my_sum += cpi->second.distance * cpi->second.distance;
  // Adjust by penalizing unmatched points?
  my_sum += 2.0 * max_pt_dist * (A_num_verts - correspondances.size());
  *RMSE = std::sqrt(my_sum / A_num_verts);
  //if(correspondances.size()) *RMSE = std::sqrt(my_sum / correspondances.size());
  
  if(num_faces){
    size_t face_count = 0;
    my_float_t area_sum = 0;

    std::vector<triangle_t>::const_iterator face = A_triangles.begin();
    for( ; face < A_triangles.end(); ++face){
      bool next = false;
      for(size_t i = 0; i < 3 && !next; ++i)
        if(correspondances.find(face->vertices[i]) == correspondances.end())
          next = true;

      if(!next){
        ++face_count;
        if(area) area_sum += face->area;
      }
    }
    *num_faces = face_count;
    if(area) *area = area_sum;
  }
  return true;
}
#endif

void 
SimpleTrimesh::compare(const SimpleTrimesh& target, 
                       const my_float_t max_pt_dist, size_t *num_points,
                       my_float_t* RMSE, correspond_map *cmap,
                       size_t *num_faces, my_float_t *area) const
{
  if(area) *area = 0.0;
  if(num_faces) *num_faces = 0;
  if(cmap) cmap->clear();
  *num_points = 0;
  *RMSE = my_float_max;

  if(!A_vert_octree || !target.A_map_vert_to_faces.size()){
    std::cerr << "Vert octree or target face map is not initialized\n";
    return;
  }

  if(max_pt_dist > 1.5){
    std::cerr << "Max point distance may be no more than 1.5 in "
              << "SimpleTrimesh::compare\n";
    return;
  }

  // The goal is to have a vector such that element ii holds all the 
  // query points the could have a correspoding point in the target face ii
  std::vector< std::map<const my_float_t*, bool> > faces_and_pts;
  faces_and_pts.resize(target.A_triangles.size());
  const my_float_t *t_vert = target.A_vertices;
  for(size_t n = 0; n < target.A_num_verts; ++n, t_vert += 3){
    const std::vector<point_vci> *query_bin;
    query_bin = A_vert_octree->near_by_points(t_vert, 2.5);

    std::vector<point_vci>::const_iterator q_iter = query_bin->begin();
    for( ; q_iter < query_bin->end(); ++q_iter){
      const my_float_t* q_vert = (*q_iter)->pos;
      //if(dist_squared(t_vert, q_vert) <= 4.0){
      if(dist_squared(t_vert, q_vert) <= 6.25){

        vert_face_mci iii = target.A_map_vert_to_faces.find(t_vert);
        const std::vector<face_vci> &t_faces = iii->second;
        std::vector<face_vci>::const_iterator t_faces_iter = t_faces.begin();
        for( ; t_faces_iter < t_faces.end(); ++t_faces_iter){
          size_t offset = *t_faces_iter - target.A_triangles.begin();
          //(*(faces_and_pts.begin() + offset))[q_vert] = true;
          std::map<const my_float_t*, bool> &map_ref = 
            *(faces_and_pts.begin() + offset);
          std::map<const my_float_t*, bool>::iterator map_iter;
          map_iter = map_ref.find(q_vert);
          if(map_iter == map_ref.end())
            map_ref.insert(std::pair<const my_float_t*, bool>(q_vert, true));
        }
      }
    }
  }

  // Intialize the corresponding points vector
  std::vector<pt_plus_dist> corr_pts(A_num_verts);
  std::vector<pt_plus_dist>::iterator pt;
  for(pt = corr_pts.begin(); pt < corr_pts.end(); ++pt)
    pt->distance = my_float_max;

  // Find the "best" corrsponding target point for each query vertex
  std::vector< std::map<const my_float_t*, bool> >::const_iterator fnp_i;
  std::vector<triangle_t>::const_iterator targ_face;
  targ_face = target.A_triangles.begin();
  for(fnp_i = faces_and_pts.begin(); fnp_i < faces_and_pts.end(); 
      ++fnp_i, ++targ_face){
    std::map<const my_float_t*, bool>::const_iterator qpt_i;
    for(qpt_i = fnp_i->begin(); qpt_i != fnp_i->end(); ++qpt_i){
      const my_float_t* query_pt = qpt_i->first;
      const size_t offset = (query_pt - A_vertices) / 3;
      std::vector<pt_plus_dist>::iterator corr_pt = corr_pts.begin() + offset;

      my_float_t d;
      my_float_t my_pt[3]; 

      // Check if corresponding point is on a face
      bool rv = corresponding_point(query_pt, targ_face->vertices[0],
                                    targ_face->vertices[1], 
                                    targ_face->vertices[2], 
				    corr_pt->distance, &d, my_pt);
      //if(rv && neg_max_pt_dist <= d && d <= max_pt_dist){
      if(rv){
        my_float_t abs_d = d;
        if(abs_d < 0) abs_d *= -1.0;
        //if(abs_d > max_pt_dist) continue;
        my_float_t abs_saved_d = corr_pt->distance;
        if(abs_saved_d < 0) abs_saved_d *= -1.0;
        if(abs_d < abs_saved_d){
          corr_pt->distance = d;
          std::copy(my_pt, my_pt + 3, corr_pt->corr_pt);
        }
      }else{
        for(uint j = 0; j < 3; ++j){
          rv = corresponding_point(query_pt, targ_face->vertices[j],
                                   targ_face->vertices[(j+1) % 3], &d, my_pt);
          my_float_t abs_d = d;
          if(abs_d < 0) abs_d *= -1.0;
          //if(abs_d > max_pt_dist) continue;
          my_float_t abs_saved_d = corr_pt->distance;
          if(abs_saved_d < 0) abs_saved_d *= -1.0;
          if(abs_d < abs_saved_d){
            corr_pt->distance = d;
            std::copy(my_pt, my_pt + 3, corr_pt->corr_pt);
          }
        }
      }
    } 
  }

  // compute num_points & RMSE
  const my_float_t neg_max_pt_dist = -1.0 * max_pt_dist;
  size_t &num_pts = *num_points;
  my_float_t SSD = 0.0;
  std::map<const my_float_t*, bool> kept_query_verts;
  const my_float_t *q_vert = A_vertices;
  for(pt = corr_pts.begin(); pt < corr_pts.end(); ++pt, q_vert += 3){
    // Moved this back up since most of the points do not actually have 
    // any assigned correspondences if they are too far away
    if(pt->distance < neg_max_pt_dist || max_pt_dist < pt->distance) continue;

    if(cmap){
      correspond_t c_entry;
      std::copy(q_vert, q_vert + 3, c_entry.model_pt);
      std::copy(pt->corr_pt, pt->corr_pt + 3, c_entry.target_pt);
      cmap->insert(std::pair<my_float_t, correspond_t>(pt->distance, c_entry));
    }

    SSD += pt->distance * pt->distance;
    ++num_pts;
    kept_query_verts.insert(kept_query_verts.end(),
                            std::pair<const my_float_t*, bool>(q_vert, true));
  }
  SSD += max_pt_dist*max_pt_dist*(A_num_verts - num_pts);
  //*RMSE = std::sqrt(SSD / num_pts);
  // we want total number of points
  *RMSE = std::sqrt(SSD / A_num_verts);
//  std::cout << "Actual RMSE: " << *RMSE << std::endl;

  // compute area, num_faces of query surface
  if(num_faces){
    face_vci q_face;
    for(q_face = A_triangles.begin(); q_face < A_triangles.end(); ++q_face){
      bool include_face = true;
      for(size_t i = 0; i < 3 && include_face; ++i){
        std::map<const my_float_t*, bool>::const_iterator kept_q_vert;
        kept_q_vert = kept_query_verts.find(q_face->vertices[i]);
        if(kept_q_vert == kept_query_verts.end()) include_face = false;
      }
      if(area && include_face){
        *area += q_face->area;
        ++(*num_faces); 
      }
    }
  }
}

void 
SimpleTrimesh::clip(const my_float_t* center, const my_float_t radius)
{
  size_t num_new_verts = 0;
  my_float_t *new_verts = 0;
  my_float_t *new_normals = 0;
  std::vector<triangle_t> new_faces; 
  clip(center, radius, &num_new_verts, &new_verts, &new_normals, &new_faces);

  A_num_verts = num_new_verts;
  if(A_vertices) delete [] A_vertices;
  A_vertices = new_verts;
  if(A_normals) delete [] A_normals;
  A_normals = new_normals;
  // this could be very expensive
  A_triangles = new_faces;
}

bool
SimpleTrimesh::write(const std::string ofname) const
{
  if(A_fail) return false;

  std::ofstream face_file, vert_file;
  if(!open_ofstream(vert_file, ofname + ".vert")) return false;
  if(!open_ofstream(face_file, ofname + ".face")) return false;

  //vert_file.setf(std::ios_base::fixed, std::ios_base::floatfield);
  vert_file << std::fixed << std::setprecision(3);
  vert_file << "# MSMS solvent excluded surface vertices ... \n"
            << "#faces  #sphere density probe_r -- #sphere is unknown \n"
            << std::setw(7) << A_num_verts
            << std::setw(8) << 0 << "  0.00  0.00";
  const my_float_t* v = A_vertices; 
  const my_float_t* n = A_normals; 
  for(size_t vert_idx = 0; vert_idx < A_num_verts; ++vert_idx){
    vert_file << "\n";
    for(size_t i = 0; i < 3; ++i, ++v) vert_file << std::setw(9) << *v << " ";
    for(size_t i = 0; i < 3; ++i, ++n) vert_file << std::setw(9) << *n << " ";
    vert_file << "      0       0  0";
  }

  face_file << "# MSMS solvent excluded surface faces ... \n"
            << "#faces  #sphere density probe_r -- #sphere is unknown \n"
            << std::setw(7) << A_triangles.size()
            << std::setw(8) << 0 << "  0.00  0.00";
  std::vector<triangle_t>::const_iterator face;
  for(face = A_triangles.begin(); face != A_triangles.end(); ++face){
    face_file << "\n";
    for(size_t i = 0; i < 3; ++i)
      face_file << std::setw(7) << 1 + (face->vertices[i] - A_vertices)/3 
                << " ";
    face_file << "0    0";
  }

  return true;
}

size_t 
SimpleTrimesh::intersect(BoundingVolume &bvol, size_t *new_num_verts, 
                         my_float_t **new_vertices, my_float_t **new_normals, 
		         std::vector<triangle_t> *new_faces)
{
  // Init a map with (old addr to keep, 0)
  const my_float_t* v = A_vertices; 
  const my_float_t* n = A_normals; 
  typedef std::pair<const my_float_t*, locations> loc_pair;
  locations locs = {0, 0, 0};
  std::map<const my_float_t*, locations> my_map;
  for(size_t i = 0; i < A_num_verts; ++i, v += 3, n += 3){
    if(bvol.contains(v)){
      locs.old_n = n;
      my_map.insert(my_map.end(), loc_pair(v, locs));
    }
  }

  // Add to the map the points corresponding to the triangles with at least
  // one point falling inside the bounding volume
  std::vector<triangle_t>::const_iterator old_face_i = A_triangles.begin();
  std::vector<bool> inside(3, false);
  std::map<const my_float_t*, locations> add_to_map;
  for( ; old_face_i < A_triangles.end(); ++old_face_i){
    for(size_t i = 0; i < 3; ++i)
      if(my_map.find(old_face_i->vertices[i]) != my_map.end()) inside[i] = true;
    if(!inside[0] && !inside[1] && !inside[2]) continue;

    for(size_t i = 0; i < 3; ++i)
      if(!inside[i]){
        locs.old_n = old_face_i->normals[i];
        add_to_map.insert(loc_pair(old_face_i->vertices[i], locs));
      }
    std::fill(inside.begin(), inside.end(), false);
  }
  my_map.insert(add_to_map.begin(), add_to_map.end());
  copy_kept(my_map, new_vertices, new_normals, new_faces);

  *new_num_verts = my_map.size();
  return *new_num_verts;
}


size_t 
SimpleTrimesh::clip(const my_float_t* center, const my_float_t radius,
                    size_t *new_num_verts, my_float_t **new_vertices, 
                    my_float_t **new_normals, 
		    std::vector<triangle_t> *new_faces)
{
  // Init a map with (old addr to keep, 0)
  const my_float_t r_squared = radius*radius;
  const my_float_t* v = A_vertices; 
  const my_float_t* n = A_normals; 
  typedef std::pair<const my_float_t*, locations> loc_pair;
  locations locs = {0, 0, 0};
  std::map<const my_float_t*, locations> my_map;
  for(size_t i = 0; i < A_num_verts; ++i, v += 3, n += 3){
    if(dist_squared(center, v) >= r_squared){
      locs.old_n = n;
      my_map.insert(my_map.end(), loc_pair(v, locs));
    }
  }
  copy_kept(my_map, new_vertices, new_normals, new_faces);

  *new_num_verts = my_map.size();
  return *new_num_verts;
}

bool
SimpleTrimesh::read_vert_file(const std::string fname)
{
  std::ifstream vert_file;
  if(!open_ifstream(vert_file, fname + ".vert")) return false;

  // Eat the first two lines
  std::string line;
  std::getline(vert_file, line);
  std::getline(vert_file, line);

  // We only use the A_num_verts field and use getline to eat the rest of the
  // line (as part of the loop below).
  vert_file >> A_num_verts;

  A_vertices = new my_float_t[3*A_num_verts];
  A_normals = new my_float_t[3*A_num_verts];
  my_float_t *vert = A_vertices;
  my_float_t *normal = A_normals;
  // Notice that the getline is used only to eat up the remaining portion 
  // of the line and that it is ok to eat it here since we have a part of 
  // the third line left at the beginning of the loop
  size_t i;
  for(i = 0; i < A_num_verts && std::getline(vert_file, line); ++i){
    vert_file >> *vert >> *(vert + 1) >> *(vert + 2)
              >> *normal >> *(normal + 1) >> *(normal + 2);
    vert += 3;
    normal += 3;
  }
  if(i != A_num_verts){
    std::cerr << "Unable to read the entire vertex file: " 
              << fname << std::endl;
    return false;
  }
  
  return true;
}

bool
SimpleTrimesh::read_face_file(const std::string fname, bool check_vertex_order)
{
  std::ifstream face_file;
  if(!open_ifstream(face_file, fname + ".face")) return false;

  // Eat the first two lines
  std::string line;
  std::getline(face_file, line);
  std::getline(face_file, line);

  size_t num_faces, num_spheres;
  my_float_t density, probe_radius;
  face_file >> num_faces >> num_spheres >> density >> probe_radius;

//  my_float_t min_dist = my_float_max, max_dist = 0.0;
//  std::multimap<my_float_t, bool> edge_lengths;

  size_t v_idz[3];
  triangle_t tmp;
  A_triangles.reserve(num_faces);
  for(size_t i = 0; i < num_faces && std::getline(face_file, line); ++i){
    face_file >> v_idz[0] >> v_idz[1] >> v_idz[2];
    for(size_t j = 0; j < 3; ++j){
      size_t idx = 3*(v_idz[j] - 1);
      tmp.vertices[j] = A_vertices + idx;
      tmp.normals[j] = A_normals + idx;
    }
#if 0
    for(int j = 0; j < 3; ++j){
      my_float_t d = dist_squared(tmp.vertices[j], tmp.vertices[(j+1)%3]);
      min_dist = (min_dist > d ? d : min_dist);
      max_dist = (max_dist < d ? d : max_dist); 
      edge_lengths.insert(std::pair<my_float_t, bool>(d, true));
    }
#endif

    // Since this check is expected to be used only at template creation time,
    // in the interest of having less code we use the corresponding point 
    // function although the centroid is already in the plane.
    if(check_vertex_order){
      my_float_t centroid[3];
      for(size_t j = 0; j < 3; ++j)
        centroid[j] = (tmp.vertices[0][j] + tmp.vertices[1][j] + 
                       tmp.vertices[2][j]) / 3.0;  

      my_float_t d;
      my_float_t pt[3];
      if(!corresponding_point(centroid, tmp.vertices[0], tmp.vertices[1],
                              tmp.vertices[2], my_float_max, &d, pt)){
        std::swap(tmp.vertices[1], tmp.vertices[2]);
        std::swap(tmp.normals[1], tmp.normals[2]);
      }
    }
    A_triangles.push_back(tmp);
  }
#if 0
  std::multimap<my_float_t, bool>::iterator e_iter;
  for(e_iter = edge_lengths.begin(); e_iter != edge_lengths.end(); ++e_iter){
    std::cout << e_iter->first << ",";
  }
  std::cout << "\n";


  std::cout << "min triangle side length " << std::sqrt(min_dist) << "\n"
            << "max triangle side length " << std::sqrt(max_dist) << "\n";
#endif

  return true;
}

void
SimpleTrimesh::copy_kept(std::map<const my_float_t*, locations> &V_keep,
                         my_float_t **new_vertices, my_float_t **new_normals, 
                         std::vector<triangle_t> *new_faces)
{
  // Copy old addr[0:3] to new addr[0:3] 
  // and set (old addr to keep, [new addr of vertex, new addr of normal])
  size_t len = 3*V_keep.size();
  *new_vertices = new my_float_t[len];
  *new_normals = new my_float_t[len];
  std::map<const my_float_t*, locations>::iterator map_i;
  my_float_t *new_v = *new_vertices;
  my_float_t *new_n = *new_normals;
  for(map_i = V_keep.begin(); map_i != V_keep.end(); ++map_i){
    std::copy(map_i->first, map_i->first + 3, new_v); 
    std::copy(map_i->second.old_n, map_i->second.old_n + 3, new_n); 
    map_i->second.v = new_v;
    map_i->second.n = new_n;
    new_v += 3;
    new_n += 3;       
  }

  // Copy all triangles with 3 vertices in map by pointing to new addr
  new_faces->clear();
  std::vector<triangle_t>::const_iterator old_face_i = A_triangles.begin();
  triangle_t new_face;
  for( ; old_face_i < A_triangles.end(); ++old_face_i){
    *(new_face.vertices + 2) = 0;
    for(size_t i = 0; i < 3; ++i){
      std::map<const my_float_t*, locations>::const_iterator m_ci;
      m_ci = V_keep.find(old_face_i->vertices[i]);
      if(m_ci == V_keep.end()) break;
      new_face.vertices[i] = m_ci->second.v;
      new_face.normals[i] = m_ci->second.n;
    }
    if(*(new_face.vertices + 2)){
      new_face.area = old_face_i->area;
      new_faces->push_back(new_face);
    }
  }
}

void
SimpleTrimesh::compute_areas()
{
  my_float_t U[3], V[3], U_x_V[3];
  std::vector<triangle_t>::iterator face;
  for(face = A_triangles.begin(); face < A_triangles.end(); ++face){
    std::copy(face->vertices[0], face->vertices[0] + 3, U); 
    std::copy(face->vertices[2], face->vertices[2] + 3, V); 
    my_axpy(3, -1.0, face->vertices[1], 1, U, 1);
    my_axpy(3, -1.0, face->vertices[1], 1, V, 1);
    cross(U, V, U_x_V); 
    face->area = 0.5 * normalize(U_x_V);
  }
}


void
SimpleTrimesh::build_vert_face_map()
{
  if(A_map_vert_to_faces.size()) return;

  face_vci face_iter = A_triangles.begin();
  for( ; face_iter < A_triangles.end(); ++face_iter){
    const triangle_t &face = *face_iter;
    for(size_t i = 0; i < 3; ++i){
      vert_face_mi map_iter = A_map_vert_to_faces.find(face.vertices[i]);

      if(map_iter == A_map_vert_to_faces.end()){
        std::vector<face_vci> tmp;
        tmp.push_back(face_iter);
        A_map_vert_to_faces[face.vertices[i]] = tmp;
      }else map_iter->second.push_back(face_iter);
    }
  }
}

void
SimpleTrimesh::build_vert_octree()
{
  if(A_vert_octree) return;

  A_vert_vector.reserve(A_num_verts);  
  my_float_t* v = A_vertices;
  point_t pt(NO_ALLOCATION);
  for(size_t i = 0; i < A_num_verts; ++i, v += 3){
    pt.pos = v; 
    A_vert_vector.push_back(pt); 
  }
  //A_vert_octree = new SimpleOctree<point_vci>(2.0, 10);
  A_vert_octree = new SimpleOctree<point_vci>(2.5, 10);
  A_vert_octree->build(A_vert_vector.begin(), A_vert_vector.end());
}
