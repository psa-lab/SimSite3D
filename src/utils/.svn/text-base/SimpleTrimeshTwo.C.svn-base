
#include <iomanip>
#include <stream_basics.H>
#include <SimpleTrimeshTwo.H>


using namespace SimSite3D::geometry;


SimpleTrimeshTwo::SimpleTrimeshTwo(const std::string fname,
                                   const bool load_mutable_faces_in)
{
  init();
  if(!read_vert_file(fname)) A_fail = true;
  if(load_mutable_faces_in){
    if(!load_mutable_faces(fname)) A_fail = true;
  }
}

SimpleTrimeshTwo::~SimpleTrimeshTwo()
{
  if(A_vertices) delete [] A_vertices;
  if(A_normals) delete [] A_normals;
  if(A_original_vertices) delete [] A_original_vertices;
  if(A_original_normals) delete [] A_original_normals;
  init();
}
  
void
SimpleTrimeshTwo::transform(const my_float_t* R, const my_float_t* T)
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
SimpleTrimeshTwo::inverse_transform(const my_float_t* R, const my_float_t* T)
{
  size_t arr_len = 3*A_num_verts;

  my_float_t *tmp_normals = new my_float_t[arr_len];
  std::copy(A_normals, A_normals + arr_len, tmp_normals);

  my_float_t R_transpose[9];
  for(int i = 0; i < 3; ++i)
    for(int j = 0; j < 3; ++j) R_transpose[3*i +j] = R[3*j + i];
  my_gemm(A_num_verts, 3, 3, 1.0, tmp_normals, 3, R_transpose, 3, 
          A_normals, 3, 0.0);

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

bool
SimpleTrimeshTwo::compare(const my_float_t *other_vert_begin, 
                          const uint num_other_vert, my_float_t *closest_pts, 
                          my_float_t *dists,
                          const my_float_t *other_normals_begin,
                          my_float_t *best_Ns, my_float_t *N_dists,
                          const my_float_t max_dist) const
{
  std::fill(dists, dists + num_other_vert, my_float_max);
  std::fill(closest_pts, closest_pts + 3*num_other_vert, my_float_max);

  my_float_t *cp = closest_pts;
  my_float_t *d = dists;
  const my_float_t *other_vert_end = other_vert_begin + 3 * num_other_vert;
  int cnt = 0;
  for(const my_float_t *v = other_vert_begin; v < other_vert_end; v += 3){
    face_vci face;
    for(face = A_triangles.begin(); face < A_triangles.end(); ++face){
      compute_close_point(v, face->vertices[0], face->vertices[1], 
                          face->vertices[2], *d, d, cp);
    }
    ++cnt;
    d += 1;
    cp += 3; 
  } 
  return true;
}


void
SimpleTrimeshTwo::compute_features(const my_float_t *dists, 
                                   const my_float_t max_dist, 
                                   const my_float_t *norm_dots,
                                   my_float_t *vertex_count, 
                                   my_float_t *face_count, my_float_t *RMSE,
                                   my_float_t *RMS_norm_err,
                                   my_float_t *area, 
                                   my_float_t *max_q_pt_dist_p) const
{
  *vertex_count = 0.0;
  *area = 0.0;
  *RMSE = 0.0;
  *RMS_norm_err = 0.0;
  *face_count = 0.0;    
  my_float_t &max_q_pt_dist = *max_q_pt_dist_p;
  max_q_pt_dist = 0.0;

  // We may wish to "throw away" triangle for which one or more vertices
  // have negative dot product with the db surface

  // Spin through each vertex
  std::vector<bool> close_enough(A_num_verts);
  std::fill(close_enough.begin(), close_enough.end(), false);
  std::vector<bool>::iterator ce = close_enough.begin();
  const my_float_t *vert_end = A_vertices + 3*A_num_verts;
  const my_float_t *curr_dist = dists;
  const my_float_t *curr_ndot = 0;
  if(norm_dots) curr_ndot = norm_dots;
  for(const my_float_t *v = A_vertices; v < vert_end; v += 3){
    if(-1.0*max_dist <= *curr_dist && *curr_dist <= max_dist){
      *ce = true;
      *RMSE += (*curr_dist) * (*curr_dist);
      *vertex_count += 1.0;
      if(*curr_dist < -1.0*max_q_pt_dist || max_q_pt_dist < *curr_dist)
        max_q_pt_dist = *curr_dist;
      if(norm_dots) *RMS_norm_err += (*curr_ndot) * (*curr_ndot);
    }
    ++ce;
    ++curr_dist;
    if(norm_dots) ++curr_ndot;
  }
  size_t num_missed = A_num_verts - static_cast<int>(*vertex_count);
  if(num_missed > 0) max_q_pt_dist = max_dist;
  *RMSE += (num_missed)*max_dist*max_dist;
  *RMSE = std::sqrt(*RMSE / A_num_verts);
  *RMS_norm_err += (A_num_verts - static_cast<int>(*vertex_count))*4;
  *RMS_norm_err = std::sqrt(*RMS_norm_err / A_num_verts);
   
  // Spin through each triangle
  face_vci face;
  for(face = A_triangles.begin(); face < A_triangles.end(); ++face){
    int idz[] = { (face->vertices[0] - A_vertices)/3,
                  (face->vertices[1] - A_vertices)/3,
                  (face->vertices[2] - A_vertices)/3};
    if(close_enough[idz[0]] && close_enough[idz[1]] && close_enough[idz[2]]){
      *area += face->area;
      *face_count += 1.0;
    }
  }   
} 

void
SimpleTrimeshTwo::compute_close_point(const my_float_t *pt, 
                                      const my_float_t *v0, 
                                      const my_float_t *v1, 
                                      const my_float_t *v2, 
				      const my_float_t prev_best_d,
                                      my_float_t *best_d, 
                                      my_float_t *closest_pt) const
{
  my_float_t d = my_float_max;
  my_float_t cp[3];
  std::fill(cp, cp+3, my_float_max);

  // Check if projection of pt to face plane is contained in the face
  bool rv = corresponding_point(pt, v0, v1, v2, prev_best_d, &d, cp);
  if(rv) update_closest_pt(best_d, closest_pt, d, cp);
  else{
    const my_float_t *V[] = { v0, v1, v2 };
    for(uint j = 0; j < 3; ++j){
      rv = corresponding_point(pt, V[j], V[(j+1) % 3], &d, cp);
      update_closest_pt(best_d, closest_pt, d, cp);
    }
  }
}

bool
SimpleTrimeshTwo::write(const std::string ofname) const
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
  for(face_vci face = A_triangles.begin(); face != A_triangles.end(); ++face){
    face_file << "\n";
    for(size_t i = 0; i < 3; ++i)
      // Stored faces are zero based indexed to vertices -- MSMS writes out
      // 1 based -- stay consistent.
      face_file << std::setw(7) << 1 + (face->vertices[i] - A_vertices)/3 
                << " ";
    face_file << "0    0";
  }

  return true;
}

bool
SimpleTrimeshTwo::read_face_file(const std::string fname, 
                                 bool check_vertex_order, 
                                 bool compute_max_dist_to_centroid)
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
  init_triangle_data_structure(A_num_verts, num_faces);
    
  uint v_idz[3];
  for(size_t i = 0; i < num_faces && std::getline(face_file, line); ++i){
    face_file >> v_idz[0] >> v_idz[1] >> v_idz[2];

    const my_float_t *tmp_vertices[3], *tmp_normals[3];
    for(size_t j = 0; j < 3; ++j){
      // MSMS vertex indices are 1 based -- we need zero based
      size_t idx = 3*(v_idz[j] - 1);
      tmp_vertices[j] = A_vertices + idx;
      tmp_normals[j] = A_normals + idx;
    }

    // Since this check is expected to be used only at template creation time,
    // in the interest of having less code we use the corresponding point 
    // function although the centroid is already in the plane.
    if(check_vertex_order){
      my_float_t centroid[3];
      for(size_t j = 0; j < 3; ++j)
        centroid[j] = (tmp_vertices[0][j] + tmp_vertices[1][j] +
                       tmp_vertices[2][j]) / 3.0;

      my_float_t d;
      my_float_t pt[3];
      if(!corresponding_point(centroid, tmp_vertices[0], tmp_vertices[1],
                              tmp_vertices[2], my_float_max, &d, pt)){
        std::swap(tmp_vertices[1], tmp_vertices[2]);
        std::swap(tmp_normals[1], tmp_normals[2]);
      }
    }

    // This should be computed at template generation time and stored in the
    // .csv file rather than computing at screening time -- we want
    // to make screening as fast as possible
    if(compute_max_dist_to_centroid){
      my_float_t centroid[3];
      for(size_t j = 0; j < 3; ++j)
        centroid[j] = (tmp_vertices[0][j] + tmp_vertices[1][j] +
                       tmp_vertices[2][j]) / 3.0;
      
      for(size_t j = 0; j < 3; ++j){
         my_float_t d2 = dist_squared(centroid, tmp_vertices[j]);
         A_max_dist_to_cent = 
           (A_max_dist_to_cent < d2 ? d2 : A_max_dist_to_cent);
      }
      A_max_dist_to_cent = sqrt(A_max_dist_to_cent);
    }

    update_triangle_data_structure(i, v_idz, tmp_vertices, tmp_normals);
  }

  compute_areas(); 
  return true;
}

void
SimpleTrimeshTwo::init()
{
  A_vertices = 0;
  A_normals = 0;
  A_original_vertices = 0;
  A_original_normals = 0;
  A_max_dist_to_cent = -1.0;
  A_num_verts = 0;
  A_fail = false;
}

bool
SimpleTrimeshTwo::read_vert_file(const std::string fname)
{
  std::ifstream vert_file;
  if(!open_ifstream(vert_file, fname + ".vert")) return false;
  
  // Eat the first two lines
  std::string line;
  std::getline(vert_file, line);
  std::getline(vert_file, line);
  
  size_t num_spheres;
  my_float_t density, probe_radius;
  vert_file >> A_num_verts >> num_spheres >> density >> probe_radius;
  
  A_vertices = new my_float_t[3*A_num_verts];
  A_normals = new my_float_t[3*A_num_verts];
  my_float_t *vert = A_vertices;
  my_float_t *normal = A_normals;
  // Notice that the getline is used only to eat up the remaining portion 
  // of the line and that it is ok to eat it here since we have a part of 
  // the third line left at the beginning of the loop
  for(size_t i = 0; i < A_num_verts && std::getline(vert_file, line); ++i){
    vert_file >> *vert >> *(vert + 1) >> *(vert + 2)
              >> *normal >> *(normal + 1) >> *(normal + 2);
    vert += 3;
    normal += 3;
  }

  return true;
}

bool
SimpleTrimeshTwo::remove_small_triangles()
{
  if(A_mutable_faces.size() == 0){
    std::cerr << "Cannot adjust the mesh if mutable faces were not loaded\n";
    return false;
  }

  const my_float_t min_edge_len = 0.01;
  const my_float_t min_len_squared = min_edge_len * min_edge_len;

  // Get a vector of pointers to the beginning of the vertices for each
  // edge with length less than min_edge_len
  std::vector< local_edge_ptrs_t > short_edges_vec;

  // Get a map of each vertex to a vector of triangles that contain that
  // vertex
  mutable_vert_face_map V2F;

  mutable_face_vi f;
  for(f = A_mutable_faces.begin(); f != A_mutable_faces.end(); ++f){
    for(int i = 0; i < 3; ++i){
      mutable_vert_face_mi v2f_iter = V2F.find(f->vertices[i]);
      if(v2f_iter == V2F.end()){
        V2F[f->vertices[i]] = std::vector<mutable_face_vi>();
        v2f_iter = V2F.find(f->vertices[i]);
      }
      v2f_iter->second.push_back(f);

      // For each edge, if the edg length squared is < min_edge_len squared,
      // put the two endpoints in a vector
      my_float_t edge_len_squared = 
        dist_squared(f->vertices[i], f->vertices[(i+1)%3]);
      if(edge_len_squared < min_len_squared){
        local_edge_ptrs_t tmp;
        tmp.vert0_p = &(f->vertices[i]);
        tmp.vert1_p = &(f->vertices[(i+1)%3]);
        tmp.normal0_p = &(f->normals[i]);
        tmp.normal1_p = &(f->normals[(i+1)%3]);
        short_edges_vec.push_back(tmp);
      }
    }
  }

  std::vector<bool> remove_verts(A_num_verts);
  std::fill(remove_verts.begin(), remove_verts.end(), false);
  std::vector<bool> remove_faces(A_mutable_faces.size());
  std::fill(remove_faces.begin(), remove_faces.end(), false);

  std::vector< local_edge_ptrs_t >::iterator SE;
  for(SE = short_edges_vec.begin(); SE < short_edges_vec.end(); ++SE){
    // If at any time we run (and we should) across edges with V0 and V1 
    // pointing to the same vertex we know that the edge was already processed
    // and we can safely ignore this edge
    if(SE->vert0_p == SE->vert1_p) continue;

    // Average the two vertices and their corresponding normals and store 
    // in the first vertex and mark second as vertex to remove
    my_axpy(3, 1.0, *(SE->vert1_p), 1, *(SE->vert0_p), 1); 
    for(int i = 0; i < 3; ++i) (*(SE->vert0_p))[i] *= 0.5;
    my_axpy(3, 1.0, *(SE->normal1_p), 1, *(SE->normal0_p), 1); 
    normalize(*(SE->normal0_p)); 
    remove_verts[ (*(SE->vert1_p) - A_vertices) / 3 ] = true;

    // Update all triangles that used the second vertex and normal to 
    // point to the first vertex and normal
    mutable_vert_face_mi v2f_iter = V2F.find(*(SE->vert1_p));
    if(v2f_iter == V2F.end()){
      std::cerr << "The second vertex was not found in the vertex to faces "
                << "map\n";
      return false;
    }
  
    std::vector<mutable_face_vi>::iterator face_ii; 
    for(face_ii = v2f_iter->second.begin(); face_ii < v2f_iter->second.end();
        ++face_ii){
      mutable_face_vi &face = *face_ii;

      for(int i = 0; i < 3; ++i){
        // Remove up to 2 triangles that share the current edge (actually need 
        // to mark for deletion -- cannot remove them as they will invalidate 
        // the triangle iterators stored in V2F.).  If the edge is on the 
        // boundary, then there will be only one triangle.
        if(face->vertices[i] == *(SE->vert0_p)){
          remove_faces[ face - A_mutable_faces.begin() ] = true;
          break;

        // Change the vertex and normal that current point to the second
        // pair to point to the first vertex and normal
        }else if(face->vertices[i] == *(SE->vert1_p)){
          face->vertices[i] = *(SE->vert0_p);
          face->normals[i] = *(SE->normal0_p);
        }
      }
    }
  }
  return true;
}

void
SimpleTrimeshTwo::compute_areas()
{
  my_float_t U[3], V[3], U_x_V[3];
  for(face_vi face = A_triangles.begin(); face < A_triangles.end(); ++face){
    std::copy(face->vertices[0], face->vertices[0] + 3, U);
    std::copy(face->vertices[2], face->vertices[2] + 3, V);
    my_axpy(3, -1.0, face->vertices[1], 1, U, 1);
    my_axpy(3, -1.0, face->vertices[1], 1, V, 1);
    cross(U, V, U_x_V);
    face->area = 0.5 * normalize(U_x_V);
  }
}

bool
SimpleTrimeshTwo::load_mutable_faces(const std::string faces_fname)
{
  std::ifstream face_file;
  if(!open_ifstream(face_file, faces_fname + ".face")) return false;

  // Eat the first two lines
  std::string line;
  std::getline(face_file, line);
  std::getline(face_file, line);

  size_t num_faces, num_spheres;
  my_float_t density, probe_radius;
  face_file >> num_faces >> num_spheres >> density >> probe_radius;
  A_mutable_faces.reserve(num_faces);

  uint v_idz[3];
  for(size_t i = 0; i < num_faces && std::getline(face_file, line); ++i){
    face_file >> v_idz[0] >> v_idz[1] >> v_idz[2];

    mutable_face_t my_face;
    for(size_t j = 0; j < 3; ++j){
      // MSMS vertex indices are 1 based -- we need zero based
      size_t idx = 3*(v_idz[j] - 1);
      my_face.vertices[j] = A_vertices + idx;
      my_face.normals[j] = A_normals + idx;
    }
    my_face.area = -1.0;
    A_mutable_faces.push_back(my_face);
  }

  return true;
}
