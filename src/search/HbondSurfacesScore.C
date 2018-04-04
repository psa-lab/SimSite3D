
#include <HbondSurfacesScore.H>

using namespace ASCbase;

const std::string HbondSurfacesScore::A_fname = "HbondSurfacesScore.C";
//const my_float_t HbondSurfacesScore::A_max_corr_dist = 1.0;
const my_float_t HbondSurfacesScore::A_max_corr_dist = 1.5;
const uint HbondSurfacesScore::A_num_terms = 13;

#if 0
const my_float_t HbondSurfacesScore::CONSTANT_TERM = -1.05127;
const my_float_t HbondSurfacesScore::SURFACE_RMSE_W = 0.707929;
const my_float_t HbondSurfacesScore::QUAD_DIST_N_DOT_PROD_DD_OR_AA_SUM_W = 
  -0.0122393;
const my_float_t HbondSurfacesScore::QUAD_DIST_N_DOT_PROD_DNEPTOR_SUM_W =
  0.00790193;
#endif

#if 1
// "New" weights for terms -- June 2010
// the naming is fubar at the present -- more than likely this is linear 
// distance and dot product weighted and not quadratic
const my_float_t HbondSurfacesScore::CONSTANT_TERM = -2.42522;
const my_float_t HbondSurfacesScore::SURFACE_RMSE_W = 1.94214;
const my_float_t HbondSurfacesScore::QUAD_DIST_N_DOT_PROD_DD_OR_AA_SUM_W = 
  -0.0100155;
const my_float_t HbondSurfacesScore::QUAD_DIST_N_DOT_PROD_DNEPTOR_SUM_W =
  0.00583479;

#endif
const my_float_t HbondSurfacesScore::SCALED_CONSTANT_TERM = -1.57122;
const my_float_t HbondSurfacesScore::SCALED_SURFACE_RMSE_W = 1.92629;
const my_float_t HbondSurfacesScore::SCALED_QUAD_DIST_N_DOT_PROD_DD_OR_AA_SUM_W 
 = -1.91508;
const my_float_t HbondSurfacesScore::SCALED_QUAD_DIST_N_DOT_PROD_DNEPTOR_SUM_W =
  -0.029968;

#if 0
my_float_t
HbondSurfacesScore::score(const ModelSitemap &model, const DbaseSitemap& dbase,
                          rigid_align_vi scores)
{
  my_float_t old_score = point_and_surf_score::score(model, dbase, scores);

  std::vector<my_float_t> &sums = scores->terms;
  uint start_idx = sums.size();
  sums.resize(sums.size() + A_num_terms);
  std::fill(sums.end() - A_num_terms, sums.end(), 0.0);
  A_num_polar_pts = 0;

#if 0  
 0) (16) Number of query donor/acceptor group matches
 1) (17) sum of dot product weighted point matches AA,DD
 2) (18) sum of linear error weighted matches (tol - dist) / tol -- AA, DD
 3) (19) sum of squared error weighted matches (tol - dist) / tol -- AA, DD
 4) (20) sum of linear error * dot product weighted matches -- AA, DD
 5) (21) sum of squared error * dot product weighted matches -- AA, DD
 6) (22) Number of query doneptor group matches
 7) (23) sum of dot product weighted point matches (doneptors)
 8) (24) sum of linear error weighted matches (tol - dist) / tol -- doneptors
 9) (25) sum of squared error weighted matches (tol - dist) / tol -- doneptors
 10) (26) sum of linear error * dot product weighted matches -- doneptors
 11) (27) sum of squared error * dot product weighted matches -- doneptors
 12) (28) underestimate of complementary surface area
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
#endif


my_float_t HbondSurfacesScore::
polar_correspondences(const ModelSitemap &model, my_float_t **vertices_ptr,
                      my_float_t **closest_pts_ptr, size_t *npts) const
{
  const ModelHbondSurfaces* hb_caps = model.model_hbond_surfaces();
  
  // Count number of vertices
  size_t N = 0;
  ModelHbondSurfaces::m_surf_vci cur_cap_iter = hb_caps->surf_caps_begin();
  for( ; cur_cap_iter < hb_caps->surf_caps_end(); ++cur_cap_iter)
    N += cur_cap_iter->num_verts();

  // Grab mem and init vars
  *vertices_ptr = new my_float_t[3*N];
  *closest_pts_ptr = new my_float_t[3*N];
  *npts = 0;
  my_float_t SE = 0.0;

  my_float_t *vert_out = *vertices_ptr;
  my_float_t *cp_out = *closest_pts_ptr;
  cur_cap_iter = hb_caps->surf_caps_begin();
  for( ; cur_cap_iter < hb_caps->surf_caps_end(); ++cur_cap_iter){ 
    const model_hbond_surf_t& m_surf = *cur_cap_iter;
    
    const my_float_t *saved_cp = m_surf.closest_pts_begin();
    const my_float_t *saved_dist = m_surf.closest_pts_dists_begin();
    geometry::vertex_vci V;
    for(V = m_surf.verts_begin(); V != m_surf.verts_end(); ++V){
      if(*saved_dist <= A_max_corr_dist){
        std::copy(V->pos, V->pos + 3, vert_out);
        std::copy(saved_cp, saved_cp + 3, cp_out);
        vert_out += 3;
        cp_out += 3;
        SE += (*saved_dist)*(*saved_dist);
        *npts += 1;
      }else{
        SE += A_max_corr_dist*A_max_corr_dist;
      }
      saved_cp += 3;
      ++saved_dist; 
    }
  }

  return std::sqrt(SE/N);
}

// NEED TO FIX THIS????
// Probably move to ModelHbondSurfaces at some point?
void HbondSurfacesScore::
polar_correspondences(const ModelSitemap &model,
                      m_cap_iters_vci cap_iters_beg,
                      m_cap_iters_vci cap_iters_end,
                      std::vector<const my_float_t *> *q_pts,
                      std::vector<const my_float_t *> *corr_pts) const
{
  //const ModelHbondSurfaces* model_surfs = model.model_hbond_surfaces();
  q_pts->resize(A_num_polar_pts);
  corr_pts->resize(A_num_polar_pts);

  //std::vector<model_hbond_surf_t>::const_iterator m_surf_iter;
  //m_surf_vci m_surf_iter;
  std::vector<const my_float_t *>::iterator q_out = q_pts->begin();
  std::vector<const my_float_t *>::iterator corr_pt_out = corr_pts->begin();
//  for(m_surf_iter = model_surfs->surf_caps_begin();
//      m_surf_iter != model_surfs->surf_caps_end(); ++m_surf_iter){
  //std::map<m_surf_vci, bool>::const_iterator cap_map_iter = cap_iters.begin();
  m_cap_iters_vci cur_cap_iter= cap_iters_beg;
  for( ; cur_cap_iter < cap_iters_end; ++cur_cap_iter){ 
    const model_hbond_surf_t& m_surf = **cur_cap_iter;
    
    const my_float_t *saved_cp = m_surf.closest_pts_begin();
    const my_float_t *saved_dist = m_surf.closest_pts_dists_begin();
    geometry::vertex_vci V;
    for(V = m_surf.verts_begin(); V != m_surf.verts_end(); ++V){
      if(*saved_dist <= A_max_corr_dist){
        *q_out = V->pos;
        *corr_pt_out = saved_cp; 
      }else{
        *q_out = 0;
        *corr_pt_out = 0;
      }
      ++q_out;
      ++corr_pt_out;
      ++saved_dist;
      saved_cp += 3;
    }
  }
}
