#ifndef IMMOVABLE_TRIMESH_HEADER_INCLUDED
#define IMMOVABLE_TRIMESH_HEADER_INCLUDED

#include <SimpleTrimeshTwo.H>
#include <DistanceArray.H>

namespace SimSite3D{
namespace geometry{

//! Simple derived class that gives us a unique number for 2/3rds of triangle
class two_thirds_triangle : public dir_point_t{
public:

  two_thirds_triangle(alloc_t a = ALLOC_POSITION) : dir_point_t(a)
  { ; }

  two_thirds_triangle(const two_thirds_triangle& p) : dir_point_t(p)
  { do_copy(p); }

  const two_thirds_triangle&
  operator=(const two_thirds_triangle& p)
  {
    if(&p != this){
      do_copy(p);
      dir_point_t::operator= (p);
    }
    return *this;
  }

  ~two_thirds_triangle()
  { ; }

  int triangle_index;  //!< Index of the triangle in the input file

private:
  void
  do_copy(const two_thirds_triangle& p)
  { triangle_index = p.triangle_index; }

};

class ImmovableTrimesh : public SimpleTrimeshTwo{
public:

  ImmovableTrimesh(const std::string fname, 
                   const my_float_t max_pt2surf_dist = 1.5) 
    : SimpleTrimeshTwo(fname)
  {
    if(!read_face_file(fname, true, true)){
      set_fail_flag();
      return;
    }

    A_max_pt2surf_dist = max_pt2surf_dist;
    const my_float_t bin_width = 2.0;
    const my_float_t max_pt2pt_dist = 
      sqrt(max_dist_to_centroid() * max_dist_to_centroid() + 
           max_pt2surf_dist*max_pt2surf_dist);
    A_grid.setup_grid(vertices_begin(), number_of_vertices(), bin_width, 
                      max_pt2pt_dist);
  }

  ~ImmovableTrimesh()
  {
  
  }

  //! Find closest point on surface mesh for the given point (pt)
  bool
  closest_point(const my_float_t *pt, my_float_t *dist, my_float_t *cp,
                const my_float_t max_dist = 1.5) const
  {
    *dist = my_float_max;
    std::fill(cp, cp + 3, my_float_max);

    if(max_dist > A_max_pt2surf_dist){
      std::cerr << "ImmovableTrimesh::closest_point requires a max point to "
                << "surface distance <= " << A_max_pt2surf_dist << " (A)\n";
      return false;
    }
    const my_float_t max_pt2pt_sq_dist = 
      max_dist_to_centroid() * max_dist_to_centroid() + max_dist*max_dist;

    DistanceArray::bin_vci bin_iter = A_grid.get_bin(pt); 
    if(bin_iter == A_grid.bins_end()) return false;

    // For each point (vertex) in the bin:
    //   if(dist(vertex, pt) <= max_pt2pt_dist)
    //     keep all triangles containing vertex
    //     uniquify using triangle index

    // Find all of the triangles we need to test
    std::map<int, bool> have_triangle;
    std::vector<const my_float_t*>::const_iterator vert = bin_iter->begin();
//    for( ; vert < bin_iter->end() && !close_enough; ++vert){
    for( ; vert < bin_iter->end(); ++vert){
      if(dist_squared(pt, *vert) > max_pt2pt_sq_dist) continue;

      // Look up the triangles for current vertex using the vertex's index
      int vert_idx = (*vert - vertices_begin()) / 3;
      const dir_point_storage<two_thirds_triangle>& vert_tris = 
        get_vert_triangles(vert_idx);

      const my_float_t *c_vert = *vert;
      dir_point_storage<two_thirds_triangle>::const_iterator delta;
      for(delta = vert_tris.begin(); delta != vert_tris.end(); delta += 2){
        if(have_triangle.find(delta->triangle_index) != have_triangle.end())
          continue;

        // No need to save triangles -- just compute closest point
        have_triangle[delta->triangle_index] = true;
        SimpleTrimeshTwo::compute_close_point(pt, c_vert, delta->pos, 
                                              (delta + 1)->pos, *dist, dist, cp);
//        if(-1.0*dist_tol < *dist && *dist < dist_tol){
//          close_enough = true;
//          break;
//        }
      }
    }
    return true;
  }

  //! Find closest point on surface mesh for the given point (pt)
  /*! This method is reasonably efficient since it uses a 3D bin to determine
   *  which triangles need to be considered for each vertex
   */
  bool
  compare(const my_float_t *other_vert_begin, const uint num_other_vert,
          my_float_t *closest_pts, my_float_t *dists,
          const my_float_t *other_normals_begin = 0, 
          my_float_t *best_Ns = 0, my_float_t *N_dists = 0,
          const my_float_t max_dist = 1.5) const
  {
//    std::cout << "\n\nTime to compare surfaces\n";
    std::fill(dists, dists + num_other_vert, my_float_max);
    std::fill(closest_pts, closest_pts + 3*num_other_vert, my_float_max);
    if(other_normals_begin && best_Ns && N_dists){
      std::fill(best_Ns, best_Ns + 3*num_other_vert,
                my_float_max);
      std::fill(N_dists, N_dists + num_other_vert, my_float_max);
    }

    if(max_dist > A_max_pt2surf_dist){
      std::cerr << "ImmovableTrimesh::compare requires a max point to "
                << "surface distance <= " << A_max_pt2surf_dist << " (A)\n";
      return false;
    }

    // It is relatively easy to show that an upper bound on a vertex from the
    // other surface to this mesh is given by the sqrt of the square of the 
    // maximum allowed distance + the square of the maximum distance from
    // the centroid of any triangle to its vertex farthest away from the
    // centroid
    const my_float_t max_pt2pt_sq_dist =
      max_dist_to_centroid() * max_dist_to_centroid() + max_dist*max_dist;

    // Get all the bins that contain at least one vertex from the other surface
    DistanceArray::bin2vert_mmap b2v_map;
    A_grid.get_bins(other_vert_begin, num_other_vert, &b2v_map);

    DistanceArray::bin2vert_mmci b2v_iter;
    for(b2v_iter = b2v_map.begin(); b2v_iter != b2v_map.end(); ++b2v_iter){
      if(b2v_iter->first == A_grid.bins_end()) break;
 
      const my_float_t *other_vert = b2v_iter->second; 
      const size_t other_vert_idx = (other_vert - other_vert_begin) / 3;
      my_float_t *curr_dist = &dists[other_vert_idx];
      my_float_t *curr_cp = &closest_pts[3*other_vert_idx];

#if 0
      std::cout.precision(6);
      std::cout << "vert idx: " << other_vert_idx << "\n";
      std::cout << "finding closest point for the vertex " 
                << other_vert[0] << " "
                << other_vert[1] << " "
                << other_vert[2] << "\n";
#endif

      // Find all of the triangles we need to test
      std::map<int, bool> have_triangle;
      std::vector<const my_float_t*>::const_iterator self_vert = 
        b2v_iter->first->begin();
      for( ; self_vert < b2v_iter->first->end(); ++self_vert){
        if(dist_squared(other_vert, *self_vert) > max_pt2pt_sq_dist) continue;

        // Look up the triangles for current vertex using the vertex's index
        int self_vert_idx = (*self_vert - vertices_begin()) / 3;
        const dir_point_storage<two_thirds_triangle>& vert_tris =
          get_vert_triangles(self_vert_idx);
    
        const my_float_t *c_vert = *self_vert;
        dir_point_storage<two_thirds_triangle>::const_iterator delta;
        for(delta = vert_tris.begin(); delta != vert_tris.end(); delta += 2){
          if(have_triangle.find(delta->triangle_index) != have_triangle.end())
            continue;

          // No need to save triangles -- just compute closest point
          have_triangle[delta->triangle_index] = true;

          if(other_normals_begin && best_Ns && N_dists){
            compute_close_point(other_vert, 
                                other_normals_begin + 3*other_vert_idx,
                                c_vert, 
                                normals_begin() + (c_vert - vertices_begin()),
                                *delta, *(delta + 1), curr_dist, curr_cp,
                                N_dists + other_vert_idx,
                                best_Ns + 3*other_vert_idx);


          }else
            SimpleTrimeshTwo::compute_close_point(other_vert, c_vert, 
                                                  delta->pos, (delta + 1)->pos,
                                                  *curr_dist, curr_dist, curr_cp);
#if 0
          std::cout << "current best cp: " 
                    << curr_cp[0] << " "
                    << curr_cp[1] << " "
                    << curr_cp[2] << "\n";
          std::cout << "current best dist: " << *curr_dist << "\n\n";
#endif
//        if(-1.0*dist_tol < *dist && *dist < dist_tol){
//          close_enough = true;
//          break;
//        }
        }
      }
//      ++curr_dist;
//      curr_cp += 3;
    }
    return true;
  }

  const dir_point_storage<two_thirds_triangle>&
  get_vert_triangles(const int vert_idx) const
  { return A_triangles_mess[vert_idx]; }

  //! Since this is an Immovable mesh, transform is not allowed
  void
  transform(const my_float_t* R, const my_float_t* T)
  {
    std::cerr << "Transformations of ImmovableTrimesh are not allowed\n";
    return;
  }

  //! Since this is an Immovable mesh, transform is not allowed
  void
  inverse_transform(const my_float_t* R, const my_float_t* T)
  {
    std::cerr << "Transformations of ImmovableTrimesh are not allowed\n";
    return;
  }

  //! No op because transformations are not allowed
  void
  revert()
  { return; }

  //! Get transform to move current vertices to original position of
  //! the vertices
  inline void
  get_current_inverse_3D_transform(Quaternion *Q, my_float_t *T) const
  {
    std::cerr << "Transformations of ImmovableTrimesh are not allowed\n";
    return;
  }


  inline void
  get_bins(const my_float_t *vert_begin, const uint num_vert,
           DistanceArray::bin2vert_mmap *pts_in_bins) const
  {
    A_grid.get_bins(vert_begin, num_vert, pts_in_bins);
  }

  //! Constant iterator to one past the last bin in the vector
  inline DistanceArray::bin_vci
  bins_end() const
  { return A_grid.bins_end(); }

private:

  void 
  init()
  { ; }

  void
  compute_close_point(const my_float_t *pt, const my_float_t *N_in, 
                      const my_float_t *v0, const my_float_t *n0, 
                      const two_thirds_triangle &delta1, 
                      const two_thirds_triangle &delta2, my_float_t *best_d, 
                      my_float_t *closest_pt, my_float_t *best_N_dist,
                      my_float_t *computed_N) const
  {
    my_float_t d = my_float_max;
//    my_float_t N_d = my_float_max;
    my_float_t cp[3];
    std::fill(cp, cp+3, my_float_max);
  
    // Check if projection of pt to face plane is contained in the face
    // NOTE: this assumes that best_d holds the previous best distance
    bool rv = corresponding_point(pt, v0, 
                                  static_cast<const my_float_t*>(delta1.pos), 
                                  static_cast<const my_float_t*>(delta2.pos), 
				  *best_d, &d, cp);
    if(rv){
      bool pos_was_updated = update_closest_pt(best_d, closest_pt, d, cp);

      // Need to think about how we wish to approach this and with what
      // accuracy
      // Compute normals -- try average for first 

      // if normals are more than 90 degrees off we may wish to drop
      // corresponding point -- for now we just assume closest point is best
      // and sum of normals will tell us otherwise
      if(pos_was_updated && best_N_dist && computed_N){
        my_float_t N_avg[3];
        std::copy(n0, n0 + 3, N_avg);
        my_axpy(3, 1.0, delta1.dir, 1, N_avg, 1);
        my_axpy(3, 1.0, delta2.dir, 1, N_avg, 1);
        normalize(N_avg);
        
        // compute 1/2 L2 norm
        // It is easy to see that if the two vectors U & V have unit length
        // the L2 norm is 2(1 - dot(U,V))
        *best_N_dist = 1 - dot(N_avg, N_in);
         
        std::copy(N_avg, N_avg + 3, computed_N); 
      }
      
    }else{ 
      const my_float_t *V[] = { v0, delta1.pos, delta2.pos};
      const my_float_t *N[] = { n0, delta1.dir, delta2.dir};
      for(uint j = 0; j < 3; ++j){
        rv = corresponding_point(pt, V[j], V[(j+1) % 3], &d, cp);
        bool pos_was_updated = update_closest_pt(best_d, closest_pt, d, cp);

        // Need to think about how we wish to approach this and with what
        // accuracy
        // Compute normals -- try average for first 
  
        // if normals are more than 90 degrees off we may wish to drop
        // corresponding point -- for now we just assume closest point is best
        // and sum of normals will tell us otherwise
        if(pos_was_updated && best_N_dist && computed_N){
          my_float_t N_avg[3];
          std::copy(N[j], N[j] + 3, N_avg);
          my_axpy(3, 1.0, N[(j+1) % 3], 1, N_avg, 1);
          normalize(N_avg);
  
          // compute 1/2 L2 norm 
          *best_N_dist = 1 - dot(N_avg, N_in);
  
          std::copy(N_avg, N_avg + 3, computed_N);
        }

      }
    }
  }


#if 0
  void
  compute_close_point(const my_float_t *pt, const my_float_t *v0,
                      const my_float_t *v1, const my_float_t *v2,
                      my_float_t *best_d, my_float_t *closest_pt) const
  {
    my_float_t d = my_float_max;
    my_float_t cp[3];
    std::fill(cp, cp+3, my_float_max);

    // Check if projection of pt to face plane is contained in the face
    bool rv = corresponding_point(pt, v0, v1, v2, &d, cp);
    if(rv) update_closest_pt(best_d, closest_pt, d, cp);
    else{
      const my_float_t *V[] = { v0, v1, v2 };
      for(uint j = 0; j < 3; ++j){
        rv = corresponding_point(pt, V[j], V[(j+1) % 3], &d, cp);
        update_closest_pt(best_d, closest_pt, d, cp);
      }
    }
  }

  inline void
  update_closest_pt(my_float_t *best_d, my_float_t *closest_pt,
                  const my_float_t d, const my_float_t *cp) const
  {
    my_float_t abs_d = (d < 0.0 ? -1.0 * d : d);
    my_float_t abs_saved_d = (*best_d < 0.0 ? -1.0 * (*best_d) : (*best_d));
    if(abs_d < abs_saved_d){
      *best_d = d;
      std::copy(cp, cp + 3, closest_pt);
    }
  }
#endif

  void
  init_triangle_data_structure(const uint num_vert, const uint num_tri)
  {
    A_triangles_mess.reserve(num_vert);
    // Be generous and assume a soft max of 10 triangles per vertex
    dir_point_storage<two_thirds_triangle> tmp_store(20,3);
    for(uint i = 0; i < num_vert; ++i) A_triangles_mess.push_back(tmp_store);
    SimpleTrimeshTwo::init_triangle_data_structure(num_vert, num_tri);
  }

  void
  update_triangle_data_structure(const uint tri_idx, 
                                 const uint *vert_idz,
                                 const my_float_t **vertices, 
                                 const my_float_t **normals)
  {
    for(int j = 0; j < 3; ++j){
      two_thirds_triangle tmp;
      tmp.triangle_index = tri_idx;

      int next_j = (j + 1) % 3;
      std::copy(vertices[next_j], vertices[next_j] + 3, tmp.pos);
      std::copy(normals[next_j], normals[next_j] + 3, tmp.dir);
      A_triangles_mess[vert_idz[j] - 1].push_back(tmp);

      next_j = (j + 2) % 3;
      std::copy(vertices[next_j], vertices[next_j] + 3, tmp.pos);
      std::copy(normals[next_j], normals[next_j] + 3, tmp.dir);
      A_triangles_mess[vert_idz[j] - 1].push_back(tmp);
    }
    SimpleTrimeshTwo::update_triangle_data_structure(tri_idx, vert_idz,
                                                     vertices, normals);
  }


  DistanceArray A_grid;
  my_float_t A_max_pt2surf_dist;  //!< Maximum distance for pt 2 surface distance to be reliable
  std::vector< dir_point_storage<two_thirds_triangle> > A_triangles_mess; //!< Triangle mess

};

} // end namespace geometry
} // end namespace SimSite3D
#endif
