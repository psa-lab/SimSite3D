/******************************************************************************
 * Copyright (c) 2006,2007, Michigan State University (MSU) Board of Trustees.
 *   All rights reserved.
 *
 * This file is part of the SimSite3D Software project.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 * Authors: Jeffrey Van Voorst, vanvoor4@msu.edu
 *          Leslie Kuhn, Ph.D., KuhnL@msu.edu 
 *****************************************************************************/

#ifndef POINT_AND_SURFACE_SCORING_HEADER_FILE
#define POINT_AND_SURFACE_SCORING_HEADER_FILE

#include <WeightedSumsScore.H>

namespace SimSite3D{

//! Test class for scoring rigid alignments of dbase sitemaps to a given model
//! using both point and surface features
class point_and_surf_score : public WeightedSumsScore{
public:
  typedef std::less<my_float_t>		score_cmp;  //! A lesser score is more favorable
  //! Default constructor for site map point scoring and surface mesh scoring
  //! of (aligned) site maps
  point_and_surf_score();

  //! basic destruction
  virtual
  ~point_and_surf_score();

  virtual my_float_t
  correspondences(const ModelSitemap &model, my_float_t **query_pts_ptr, 
                  my_float_t **db_pts_ptr, size_t *npts);

  virtual void
  correspondences(const ModelSitemap &model, 
                  std::vector<const my_float_t *> *q_pts,
                  std::vector<const my_float_t *> *corr_pts);
#if 0
  // No idea why the compiler won't inherit this one in point_and_surf_score
  //! Used for IK where not all of the caps may be mobile
  virtual void 
  polar_correspondences(const ModelSitemap &model,
                        m_cap_iters_vci cap_iters_beg,
                        m_cap_iters_vci cap_iters_end,
                        std::vector<const my_float_t *> *q_pts,
                        std::vector<const my_float_t *> *corr_pts) const
  {
    std::cerr << "point_and_surf_score polar correspondences are NULL\n";
  }

  virtual my_float_t
  polar_correspondences(const ModelSitemap &model, my_float_t **vertices_ptr,
                        my_float_t **closest_pts_ptr, size_t *npts)
  {
    std::cerr << "point_and_surf_score polar correspondences are NULL\n";
  }
#endif

  //! Given an alignment of dbase site to the query, score the two sites 
  //! to the query relative using the site map points
  /*!
   * @param model Const ref to the model site
   * @param dbase Const ref to the sitemap aligned to the model sitemap
   * @param scores iterator to the score data class
   * @return The score of the alignment
   */
  //template< class align_T_iterator > my_float_t 
  my_float_t
  score(const ModelSitemap &model, const DbaseSitemap& dbase, 
        rigid_align_t *scores)
        //align_T_iterator scores)
  {
    WeightedSumsScore::score(model, dbase, scores);
    //my_float_t pt_score = WeightedSumsScore::score(model, dbase, scores);
  
    std::vector<my_float_t> &terms = scores->terms;
    uint start_idx = terms.size();
    terms.resize(terms.size() + A_num_terms);
    std::fill(terms.end() - A_num_terms, terms.end(), 0.0);
  
    my_float_t &closest_polar_sum = terms[0];
    my_float_t &polar_mismatch_sum = terms[1];
    my_float_t &AA_DD_sum = terms[2];
    my_float_t &doneptor_sum = terms[3];
    my_float_t &phobic_count = terms[4];
    my_float_t &unsat_polar_count = terms[5]; 
  
    // Surface complementarity
    my_float_t &num_surf_pts = terms[start_idx];
    my_float_t &num_faces = terms[start_idx + 1];
    my_float_t &RMSE = terms[start_idx + 2];
    my_float_t &area = terms[start_idx + 3];
    my_float_t &RMS_norm_err = terms[start_idx + 4];
    my_float_t &max_q_pt_dist = terms[start_idx + 5];
  
    const geometry::TransformableTrimesh& model_surf = 
      model.binding_site_mesh_handle();
    const geometry::SimpleTrimeshTwo& dbase_surf = 
      dbase.binding_site_mesh_handle();
    A_max_surf_pt_dist = dbase.max_corr_surf_pt_dist();
    init_storage(model_surf);

//    std::cout << "number of model surf verts: "  
//              << model_surf.number_of_vertices() << "\n";
//    std::cout << "number of dset surf verts: "  
//              << dbase_surf.number_of_vertices() << "\n";
  
    dbase_surf.compare(model_surf.vertices_begin(), 
                       model_surf.number_of_vertices(), A_surf_closest_pts, 
                       A_dists, model_surf.normals_begin(), A_best_normals, 
                       A_norm_dists, A_max_surf_pt_dist);
    model_surf.compute_features(A_dists, A_max_surf_pt_dist, A_norm_dists, 
  //  model_surf.compute_features(A_dists, 1.5, A_norm_dists, 
                                &num_surf_pts, &num_faces, &RMSE, &RMS_norm_err,
                                &area, &max_q_pt_dist);
    if(scale_terms()){
      // The terms before start_idx were scaled in the base class
      std::vector<my_float_t>::const_iterator max_val =
        model.get_max_feature_vals().begin() + start_idx;
      std::vector<my_float_t>::iterator t = terms.begin() + start_idx;
      for( ; t < terms.end(); ++t, ++max_val) *t /= *max_val;
  
      return SCALED_CONSTANT_TERM + 
             SCALED_POLAR_SUM_W * closest_polar_sum + 
             SCALED_POLAR_MISMATCH_W * polar_mismatch_sum + 
             SCALED_DD_OR_AA_POLAR_SUM_W * AA_DD_sum +
             SCALED_DONEPTOR_POLAR_SUM_W * doneptor_sum +
             SCALED_HPHOBIC_COUNT_W * phobic_count + 
             SCALED_UNSAT_POLAR_W * unsat_polar_count +
             SCALED_SURFACE_RMSE_W * RMSE;
  
    }else{
      return CONSTANT_TERM + 
             POLAR_SUM_W * closest_polar_sum + 
             POLAR_MISMATCH_W * polar_mismatch_sum + 
             DD_OR_AA_POLAR_SUM_W * AA_DD_sum +
             DONEPTOR_POLAR_SUM_W * doneptor_sum +
             HPHOBIC_COUNT_W * phobic_count + 
             UNSAT_POLAR_W * unsat_polar_count +
             SURFACE_RMSE_W * RMSE;
    }
  }
  
  virtual inline bool
  uses_surface_mesh() 
  { return true; }

  virtual inline bool
  uses_hbond_surfaces() 
  { return false; }

private:

  //! Initialize the closest point and distance variables
  inline void
  init_storage(const geometry::TransformableTrimesh& model_surf)
  {
    const int N = model_surf.number_of_vertices();
    if(N <= A_max_num_verts) return;
  
    if(A_surf_closest_pts){
      delete[] (A_surf_closest_pts);
     A_surf_closest_pts = 0;
     A_dists = 0;
    }
  
    A_max_num_verts = N;
    A_surf_closest_pts = new my_float_t[8 * N];
    A_dists = &A_surf_closest_pts[3*N];
    A_best_normals = &A_surf_closest_pts[4*N];
    A_norm_dists = &A_surf_closest_pts[7*N];
  }

  int A_max_num_verts;
  my_float_t A_max_surf_pt_dist;
  my_float_t *A_surf_closest_pts;
  my_float_t *A_dists; 
  my_float_t *A_best_normals;
  my_float_t *A_norm_dists;

  static const my_float_t CONSTANT_TERM;
  static const my_float_t POLAR_SUM_W;
  static const my_float_t POLAR_MISMATCH_W;
  static const my_float_t DD_OR_AA_POLAR_SUM_W;
  static const my_float_t DONEPTOR_POLAR_SUM_W;
  static const my_float_t HPHOBIC_COUNT_W;
  static const my_float_t UNSAT_POLAR_W;
  static const my_float_t SURFACE_RMSE_W;
  static const my_float_t SCALED_CONSTANT_TERM;
  static const my_float_t SCALED_POLAR_SUM_W;
  static const my_float_t SCALED_POLAR_MISMATCH_W;
  static const my_float_t SCALED_DD_OR_AA_POLAR_SUM_W;
  static const my_float_t SCALED_DONEPTOR_POLAR_SUM_W;
  static const my_float_t SCALED_HPHOBIC_COUNT_W;
  static const my_float_t SCALED_UNSAT_POLAR_W;
  static const my_float_t SCALED_SURFACE_RMSE_W;
  static const int A_num_terms;
  static const std::string A_fname;   //!< Source file name
};
}
#endif
