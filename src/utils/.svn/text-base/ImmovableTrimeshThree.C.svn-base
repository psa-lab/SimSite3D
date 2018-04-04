#include <ImmovableTrimeshThree.H>
#include <iomanip>


using namespace SimSite3D::geometry;

ImmovableTrimeshThree::ImmovableTrimeshThree(const std::string fname,
                                             const my_float_t max_pt2surf_dist)
 : SimpleTrimeshThree()
{  
  if(read_vert_file(fname, false) == false || 
     // testing method -- so set these to true for now
     read_face_file(fname, 0, true, true) == false){
    set_fail_flag();
    return;
  }

  A_max_pt2surf_dist = max_pt2surf_dist;
  const my_float_t bin_width = 2.0;
  const my_float_t max_pt2pt_dist =
    sqrt(max_dist_to_centroid() * max_dist_to_centroid() +
         max_pt2surf_dist*max_pt2surf_dist);
  std::cout << "Max pt 2 pt distance: " << max_pt2pt_dist << "\n";
  A_grid.setup_grid(vertices_begin(), vertices_end(), faces_begin(),
                    faces_end(), bin_width, max_pt2pt_dist);
}

bool
ImmovableTrimeshThree::compare(const my_float_t *other_vert_begin, 
                               const size_t num_other_vert,
                               my_float_t *closest_pts, my_float_t *dists,
                               my_float_t *Ns_for_cp,
                               const my_float_t max_dist) const
{
  // Initialize the correspondences arrays
  std::fill(dists, dists + num_other_vert, my_float_max);
  std::fill(closest_pts, closest_pts + 3*num_other_vert, my_float_max);
  if(Ns_for_cp)
    std::fill(Ns_for_cp, Ns_for_cp + 3*num_other_vert, my_float_max);
  
  if(max_dist > A_max_pt2surf_dist){
    std::cerr << "ImmovableTrimeshThree::closest_point requires a max point to "
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

  // Get the bin for each vertex
  FaceBins::bin2vert_mmap verts_N_bins;
  A_grid.get_bins(other_vert_begin, num_other_vert, &verts_N_bins);
  if(verts_N_bins.size() < 1) return false;

  // For each bin, check the vertices in that bin against all of the 
  // triangles
  FaceBins::bin2vert_mmci lwr_bnd = verts_N_bins.begin();
  FaceBins::bin2vert_mmci upr_bnd = verts_N_bins.begin();
  for( ; lwr_bnd != verts_N_bins.end(); lwr_bnd = upr_bnd){
    upr_bnd = verts_N_bins.upper_bound(lwr_bnd->first);

    // For each face in the bin, check all vertices in the bin
    std::vector<FaceAttrib::vci>::const_iterator bin_iter;
    bin_iter = lwr_bnd->first->begin(); 
    for( ; bin_iter < lwr_bnd->first->end(); ++bin_iter){
      FaceAttrib::vci F = *bin_iter;

      // Check each vertex in [lwr_bnd, upr_bnd[  w.r.t. the current face F
      FaceBins::bin2vert_mmci c_pair;
      for(c_pair = lwr_bnd; c_pair != upr_bnd; ++c_pair){
        // debugging check 
        if(c_pair->first != lwr_bnd->first) 
          std::cerr << "vertex range has different bins" << std::endl;

        // Hmm cannot win -- either have to cycle in inner loop over 
	// faces in the bin or the 4 correspondence arrays :(
	// It is quite possible that looping over the faces in the bin
	// could be the better choice -- or it could depend on hardware as 
	// well
	const my_float_t *V = c_pair->second;
	const size_t three_X_V_idx = V - other_vert_begin;
	const size_t V_idx = three_X_V_idx / 3;
        my_float_t *best_d = dists + V_idx;

	my_float_t cp[3], d, N_for_cp[3];
        if(Ns_for_cp) F->closest_point(V, *best_d, &d, cp, N_for_cp);
        else F->closest_point(V, *best_d, &d, cp, 0);
	
	// We may wish to check vectors here and if the dot product of
	// N_for_cp and N is < 0.0, then ignore this cp
	my_float_t abs_best_d = (*best_d > 0 ? *best_d : -1.0 * (*best_d));
	my_float_t abs_d = (d > 0 ? d : -1.0 * d);
	if(abs_d < abs_best_d){
          my_float_t *best_cp = closest_pts + three_X_V_idx;
          *best_d = d;
	  std::copy(cp, cp + 3, best_cp);
	  if(Ns_for_cp){
	    my_float_t *N_for_best_cp = Ns_for_cp + three_X_V_idx;
	    std::copy(N_for_cp, N_for_cp + 3, N_for_best_cp);
	  }
	}
      }
    }
  }

  return true;
}
