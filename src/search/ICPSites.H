#ifndef ICP_SITES_HEADER_FILE_INCLUDED
#define ICP_SITES_HEADER_FILE_INCLUDED

#include <ModelSitemap.H>
#include <DbaseSitemap.H>
#include <AlignDataWithPreIKData.H>

namespace SimSite3D{
namespace ICP{

template<typename Score_T, typename Align_T> bool
fine_tune_surfaces(ModelSitemap *model, const DbaseSitemap *dset,
                   Align_T *align, Score_T *score_method, 
                   const my_float_t min_delta = 1E-06,
                   const size_t max_iter=100, const bool verbose=false)
{
  if(verbose) std::cout << "\n\n\n\n Inside fine tune align\n";
  const my_float_t RMSE_TOL = 0.05;
  my_float_t prev_RMSE = my_float_max;
  my_float_t* query_pts = 0, *db_pts = 0;
  size_t npts;

  // Save the initial alignment in the event we need to revert to it
  Align_T initial_align(*align);

  for(size_t n = 0; n < max_iter; ++n){
    // Get current correspondences
    align->score = score_method->score(*model, *dset, align);
    my_float_t RMSE = 
      score_method->correspondences(*model, &query_pts, &db_pts, &npts);

    if(RMSE <= RMSE_TOL) break;    
    if(prev_RMSE - RMSE < min_delta) break;
    prev_RMSE = RMSE;

    // We require at least 3 correspondences to fully describe a 3D alignment
    if(npts < 3){
      // Clean up mem
      if(query_pts) delete [] query_pts;
      if(db_pts) delete [] db_pts;
      query_pts = 0;
      db_pts = 0;

      return false;
    }

    // Get the transformation 
    Quaternion Q;
    my_float_t T[3];
    lse_3D_fit(query_pts, db_pts, npts, &Q, T);
  
    // Clean up mem
    if(query_pts) delete [] query_pts;
    if(db_pts) delete [] db_pts;
    query_pts = 0;
    db_pts = 0;
  
    // Apply the most recent transformation
    my_float_t R[9];
    Q.get_ortho_rot_mat(R);
  
    // This is a very inelegant solution, but we need to concentrate our
    // resources elsewhere at this point (rather than on numerics).
    // The main issue we have been facing is with a small number of points
    // 3  (the minimum required).  When we have a small number of points, 
    // sometimes too many of the corresponding points are the same;  in such 
    // a cases, the solution is degenerate and the roots of the characteristic
    // equation contain a complex component.
    if(std::isnan(T[0]) || std::isnan(T[1]) || std::isnan(T[2])){
       std::cout << "Fine tuning failed.\n"; 
       std::cout << "Number of pairs of corresponding points: "
                 << npts << "\n\n";
       *align = initial_align;
       return false;
    } 

    // Apply the transform
    model->inverse_transform(R, T, score_method->uses_surface_mesh(),
                             score_method->uses_hbond_surfaces());

  }

  // Get a single transformation that moves the dbase site to the query
//  model->get_current_inverse_3D_transform(&(align->Q), align->T);
//  align->Q.get_ortho_rot_mat(align->R);
//  align->score = (align->score - A_mu) / A_sigma;

  return true;
}

template<typename Score_T, typename Align_T> bool
fine_tune_caps_and_surf(ModelSitemap *model, const DbaseSitemap *dset,
                        Align_T *align, Score_T *score_method, 
                        const my_float_t surf_pt_W, const my_float_t hb_cap_W,
                        const my_float_t min_delta = 1E-06,
                        const size_t max_iter=100, const bool verbose=false)
{
  if(verbose) std::cout << "\n\n\n\n Inside fine tune align\n";
  const my_float_t RMSE_TOL = 0.05;
  my_float_t prev_RMSE = my_float_max;
  my_float_t* query_pts = 0, *db_pts = 0, *q_polar_pts = 0, *db_polar_pts = 0;
  size_t npts, n_polar_pts;

  // Save the initial alignment in the event we need to revert to it
  Align_T initial_align(*align);

  for(size_t n = 0; n < max_iter; ++n){

    // Get current correspondences
    align->score = score_method->score(*model, *dset, align);
    score_method->correspondences(*model, &query_pts, &db_pts, &npts);
    score_method->polar_correspondences(*model, &q_polar_pts, &db_polar_pts,
                                        &n_polar_pts);

    size_t tot_npts = npts + n_polar_pts;
    my_float_t *q_pts = new my_float_t[6*tot_npts];
    my_float_t *d_pts = q_pts + 3*tot_npts;

    // Need to be careful with "zero" weights
    tot_npts = 0;
    if(surf_pt_W > 0){
      tot_npts = npts;
      std::copy(query_pts, query_pts + 3*npts, q_pts);
      std::copy(db_pts, db_pts + 3*npts, d_pts);
      if(hb_cap_W > 0){
        tot_npts += n_polar_pts;
        std::copy(q_polar_pts, q_polar_pts + 3*n_polar_pts, q_pts + 3*npts);
        std::copy(db_polar_pts, db_polar_pts + 3*n_polar_pts, d_pts + 3*npts);
      }
    }else if(hb_cap_W > 0){
      tot_npts = n_polar_pts;
      std::copy(q_polar_pts, q_polar_pts + 3*n_polar_pts, q_pts);
      std::copy(db_polar_pts, db_polar_pts + 3*n_polar_pts, d_pts);
    }

    Quaternion Q;
    my_float_t T[3];
    my_float_t RMSE = 0.0;

    // We require at least 3 correspondences to fully describe a 3D alignment
    if(tot_npts > 2){
      // Fill the weight vector if needed
      my_float_t *W = 0;
      if(surf_pt_W != hb_cap_W){
        W = new my_float_t[tot_npts];
        if(surf_pt_W > 0){
          std::fill(W, W + npts, surf_pt_W);
          if(hb_cap_W > 0) std::fill(W + npts, W + tot_npts, hb_cap_W);
        }else if(hb_cap_W > 0) std::fill(W, W + n_polar_pts, hb_cap_W);
      }

      // Compute the RMSE before computing the transform
      my_float_t *tmp = new my_float_t[3*tot_npts];
      vector(3*tot_npts, q_pts, db_pts, tmp);
      my_float_t N = 0.0;
      if(surf_pt_W != hb_cap_W){
        int cnt = 0;
        my_float_t *w = W;
        for(size_t i = 0; i < 3*tot_npts; ++i, ++cnt){
          RMSE += (*w) * tmp[i]*tmp[i];
          if(cnt >= 3){
            N += (*w);
            cnt = 0;
            ++w;    
          } 
        }
      }else{
        N = tot_npts;
        for(size_t i = 0; i < 3*tot_npts; ++i) RMSE += tmp[i]*tmp[i];
      }
      RMSE = std::sqrt(RMSE / N);

      if(tmp) delete[] tmp;
      tmp = 0;

      // see if this makes sense
      if(verbose) std::cout << "round " << n << " RMSD: " << RMSE << "\n";
      if(RMSE <= RMSE_TOL) break;    
      if(prev_RMSE - RMSE < min_delta) break;
      prev_RMSE = RMSE;

  
      // Get most recent transformation 
      lse_3D_fit(q_pts, d_pts, tot_npts, &Q, T, W);
      if(q_pts) delete [] q_pts;
      if(W) delete [] W;
      W = d_pts = q_pts = 0;
    }

    // Clean up mem
    if(query_pts) delete [] query_pts;
    if(db_pts) delete [] db_pts;
    if(q_polar_pts) delete [] q_polar_pts;
    if(db_polar_pts) delete [] db_polar_pts;
    query_pts = 0;
    db_pts = 0; 
    q_polar_pts = 0;
    db_polar_pts = 0;

    // We require at least 3 correspondences to fully describe a 3D alignment
    if(tot_npts < 3) return false;

    // Apply the most recent transformation
    my_float_t R[9];
    Q.get_ortho_rot_mat(R);

    // This is a very inelegant solution, but we need to concentrate our
    // resources elsewhere at this point (rather than on numerics).
    // The main issue we have been facing is with a small number of points
    // 3  (the minimum required).  When we have a small number of points, 
    // sometimes too many of the corresponding points are the same;  in such 
    // a cases, the solution is degenerate and the roots of the characteristic
    // equation contain a complex component.
    if(std::isnan(T[0]) || std::isnan(T[1]) || std::isnan(T[2])){
       std::cout << "Fine tuning failed.\n";
       std::cout << "Number of pairs of corresponding points: "
                 << npts << "\n\n";
       *align = initial_align;
       return false;
    }

    // Apply the transform
    model->inverse_transform(R, T, score_method->uses_surface_mesh(),
                             score_method->uses_hbond_surfaces());
      
  }

  // Get a single transformation that moves the dbase site to the query
//  model->get_current_inverse_3D_transform(&(align->Q), align->T);
//  align->Q.get_ortho_rot_mat(align->R);
//  align->score = (align->score - A_mu) / A_sigma;
    
  return true;
}



} // end namespace ICP
} // end namespace SimSite3D
#endif
