
#ifndef HBOND_SURFACES_SCORE_HEADER_INCLUDED
#define HBOND_SURFACES_SCORE_HEADER_INCLUDED

#include <point_and_surf_score.H> 

namespace SimSite3D{

class HbondSurfacesScore : public point_and_surf_score{
public:
  typedef std::less<my_float_t>         score_cmp;  //! A lesser score is more favorable

  typedef std::vector<ModelHbondSurfaces::m_surf_vci>::const_iterator
    m_cap_iters_vci;

  //! Do nothing constructor
  HbondSurfacesScore() : point_and_surf_score() 
  { A_num_polar_pts = 0; }

  //! basic destruction
  virtual
  ~HbondSurfacesScore()
  { ; }

  //! Given an alignment of dbase site to the query, score the two sites 
  //! to the query relative using the site map points
  /*!
   * @param model Const ref to the model site
   * @param dbase Const ref to the sitemap aligned to the model sitemap
   * @param scores iterator to the score data class
   * @return The score of the alignment
   */
  //template< class align_T_iterator >my_float_t 
  my_float_t
  score(const ModelSitemap &model, const DbaseSitemap& dbase, 
        rigid_align_t *scores)
        //align_T_iterator scores)
  {
    //my_float_t old_score = point_and_surf_score::score(model, dbase, scores);
    point_and_surf_score::score(model, dbase, scores);
  
    std::vector<my_float_t> &sums = scores->terms;
    uint start_idx = sums.size();
    sums.resize(sums.size() + A_num_terms);
    std::fill(sums.end() - A_num_terms, sums.end(), 0.0);
    A_num_polar_pts = 0;

  
  #if 0  
   0) (17) Number of query donor/acceptor group matches
   1) (18) sum of dot product weighted point matches AA,DD
   2) (19) sum of linear error weighted matches (tol - dist) / tol -- AA, DD
   3) (20) sum of squared error weighted matches (tol - dist) / tol -- AA, DD
   4) (21) sum of linear error * dot product weighted matches -- AA, DD
   5) (22) sum of squared error * dot product weighted matches -- AA, DD
   6) (23) Number of query doneptor group matches
   7) (24) sum of dot product weighted point matches (doneptors)
   8) (25) sum of linear error weighted matches (tol - dist) / tol -- doneptors
   9) (26) sum of squared error weighted matches (tol - dist) / tol -- doneptors
   10) (27) sum of linear error * dot product weighted matches -- doneptors
   11) (28) sum of squared error * dot product weighted matches -- doneptors
   12) (29) underestimate of complementary surface area
  #endif
  
    HbondSurfaces<hbond_surface_t>& db_surfs = dbase.hbond_surfaces();
    const ModelHbondSurfaces* model_surfs = model.model_hbond_surfaces();
  
    // Save correlated points for ICP in the future
  
    // Use brute force search for now -- just get this thing coded
    std::vector<model_hbond_surf_t>::const_iterator m_surf_iter;
    for(m_surf_iter = model_surfs->surf_caps_begin();
        m_surf_iter != model_surfs->surf_caps_end(); ++m_surf_iter){
      const model_hbond_surf_t& m_surf = *m_surf_iter;
      A_num_polar_pts += m_surf.num_verts();
  
  /*
        std::cout << "Query HB cap for residue (" << m_surf.atom->chainID << ") "
	          << PDB_residues::residue_to_string(m_surf.atom->res)
		  << m_surf.atom->res_num << " "
		  << PDB_residues::atom_to_string(m_surf.atom->name) << std::endl;
  */
  
      uint N = m_surf.num_terms() * m_surf.num_verts();
      my_float_t* AA_DD_terms = new my_float_t[N];
      my_float_t* doneptor_terms = new my_float_t[N];
      my_float_t surface_area;
      m_surf.compute_best_terms(db_surfs.surf_caps_begin(),
                                db_surfs.surf_caps_end(), AA_DD_terms, 
                                doneptor_terms, &surface_area, A_max_corr_dist);
   
      // Sum the point contributions
      const uint doneptor_idx = start_idx + m_surf.num_terms();
      uint AA_DD_pt_cnt = 0;
      uint doneptor_pt_cnt = 0;
      for(uint j = 0; j < m_surf.num_verts(); ++j){
        if(AA_DD_terms[j*m_surf.num_terms()]) ++AA_DD_pt_cnt;
        else if(doneptor_terms[j*m_surf.num_terms()]) ++doneptor_pt_cnt;
  
        // We want count of groups not points for terms 0 & 6
        for(uint i = 1; i < m_surf.num_terms(); ++i){ 
          sums[start_idx + i] += AA_DD_terms[j*m_surf.num_terms() + i];
          sums[doneptor_idx + i] += doneptor_terms[j*m_surf.num_terms() + i];
        }
      }
  
      if(AA_DD_pt_cnt > 0 && AA_DD_pt_cnt > doneptor_pt_cnt) ++sums[start_idx];
      else if(doneptor_pt_cnt > 0) ++sums[doneptor_idx];
  
      sums[start_idx + 12] += surface_area;
  
      delete [] AA_DD_terms;
      delete [] doneptor_terms;
    }
  
    if(scale_terms()){
      // The terms before start_idx were scaled in the base class
      std::vector<my_float_t>::const_iterator max_val =
        model.get_max_feature_vals().begin() + start_idx;
      std::vector<my_float_t>::iterator t = sums.begin() + start_idx;
      for( ; t < sums.end(); ++t, ++max_val) *t /= *max_val;
  
      return SCALED_CONSTANT_TERM + 
             SCALED_SURFACE_RMSE_W * sums[12] + 
             SCALED_QUAD_DIST_N_DOT_PROD_DD_OR_AA_SUM_W * sums[start_idx + 4] +
             SCALED_QUAD_DIST_N_DOT_PROD_DNEPTOR_SUM_W * sums[start_idx + 10];
    }else{
      return CONSTANT_TERM + 
             SURFACE_RMSE_W * sums[12] + 
             QUAD_DIST_N_DOT_PROD_DD_OR_AA_SUM_W * sums[start_idx + 4] +
             QUAD_DIST_N_DOT_PROD_DNEPTOR_SUM_W * sums[start_idx + 10];
    }
  
  }

  virtual my_float_t
  polar_correspondences(const ModelSitemap &model, my_float_t **vertices_ptr,
                        my_float_t **closest_pts_ptr, size_t *npts) const;

  virtual void
  polar_correspondences(const ModelSitemap &model,
                        m_cap_iters_vci cap_iters_beg,
                        m_cap_iters_vci cap_iters_end,
                        std::vector<const my_float_t *> *q_pts,
                        std::vector<const my_float_t *> *corr_pts) const;

  inline bool
  uses_hbond_surfaces()
  { return true; }

private:
  size_t A_num_polar_pts;  //!< number of points in all polar caps
  static const std::string A_fname;   //!< Source file name
  static const uint A_num_terms;
  static const my_float_t A_max_corr_dist;
  static const my_float_t CONSTANT_TERM;
  static const my_float_t SURFACE_RMSE_W;
  static const my_float_t QUAD_DIST_N_DOT_PROD_DD_OR_AA_SUM_W;
  static const my_float_t QUAD_DIST_N_DOT_PROD_DNEPTOR_SUM_W;
  static const my_float_t SCALED_CONSTANT_TERM;
  static const my_float_t SCALED_SURFACE_RMSE_W;
  static const my_float_t SCALED_QUAD_DIST_N_DOT_PROD_DD_OR_AA_SUM_W;
  static const my_float_t SCALED_QUAD_DIST_N_DOT_PROD_DNEPTOR_SUM_W;
};

} // end namespace SimSite3D
#endif
