/******************************************************************************
 * Copyright (c) 2006-2010, Michigan State University (MSU) Board of Trustees.
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
 * Authors: Jeffrey Van Voorst, jeff.vanvoorst@gmail.com
 *          Leslie Kuhn, Ph.D., KuhnL@msu.edu 
 *****************************************************************************/

#include <WeightedSumsScore.H>
#include <math_basics.H>
#include <iostream>

using namespace SimSite3D;

const std::string WeightedSumsScore::A_fname = "WeightedSumsScore.C";

#if 0
// "new" weights -- Jan, 2010
const my_float_t WeightedSumsScore::CONSTANT_TERM = -0.0346633;
const my_float_t WeightedSumsScore::POLAR_SUM_W = 0.0;
const my_float_t WeightedSumsScore::POLAR_MISMATCH_W = 0.0;
const my_float_t WeightedSumsScore::DD_OR_AA_POLAR_SUM_W = -0.0399575;
const my_float_t WeightedSumsScore::DONEPTOR_POLAR_SUM_W = 0.00100798;
const my_float_t WeightedSumsScore::HPHOBIC_COUNT_W = -0.0252981;
const my_float_t WeightedSumsScore::UNSAT_POLAR_W = 0.0;
#endif

#if 1
// "new" weights -- June, 2010
const my_float_t WeightedSumsScore::CONSTANT_TERM = 0.0054731; 
const my_float_t WeightedSumsScore::POLAR_SUM_W = 0.0;
const my_float_t WeightedSumsScore::POLAR_MISMATCH_W = 0.0;
const my_float_t WeightedSumsScore::DD_OR_AA_POLAR_SUM_W = -0.0456776;
const my_float_t WeightedSumsScore::DONEPTOR_POLAR_SUM_W = -0.012702;
const my_float_t WeightedSumsScore::HPHOBIC_COUNT_W = -0.0323914;
const my_float_t WeightedSumsScore::UNSAT_POLAR_W = 0.0;
#endif

const my_float_t WeightedSumsScore::SCALED_CONSTANT_TERM = 0.0373266;
const my_float_t WeightedSumsScore::SCALED_POLAR_SUM_W = 0.0;
const my_float_t WeightedSumsScore::SCALED_POLAR_MISMATCH_W = 0.0;
const my_float_t WeightedSumsScore::SCALED_DD_OR_AA_POLAR_SUM_W = -1.39892;
const my_float_t WeightedSumsScore::SCALED_DONEPTOR_POLAR_SUM_W = 0.0329575;
const my_float_t WeightedSumsScore::SCALED_HPHOBIC_COUNT_W = -0.367527;
const my_float_t WeightedSumsScore::SCALED_UNSAT_POLAR_W = 0.0;

const my_float_t WeightedSumsScore::A_max_dist = 1.5;

#if 0
my_float_t 
WeightedSumsScore::score(const ModelSitemap &model, const DbaseSitemap &dbase,
                         rigid_align_vi scores)
//                         rigid_align_vi scores, bool save_correspondences)
{
  init_storage(model);
  std::vector<my_float_t> &terms = scores->terms;
  terms.resize(10);
  std::fill(terms.begin(), terms.end(), 0.0);

  my_float_t &closest_polar_sum = terms[0];
  my_float_t &polar_mismatch_sum = terms[1];
  my_float_t &AA_DD_sum = terms[2];
  my_float_t &doneptor_sum = terms[3];
  my_float_t &phobic_count = terms[4];
  my_float_t &unsat_polar_count = terms[5];

  my_float_t &best_polar_sum = terms[6];
  my_float_t &best_AA_DD_sum = terms[7];
  my_float_t &best_doneptor_sum = terms[8];
  my_float_t &best_polar_mismatch_sum = terms[9];

  std::fill(scores->match_print.begin(), scores->match_print.end(), false);
  std::vector<bool>::iterator M_hit = scores->match_print.begin();
  std::fill(A_closest_pts, A_closest_pts + 4*A_max_num_pts, my_float_max);

  const HbondPoints& M_hbonds = model.hbond_points();
  const HphobPoints& M_hphobs = model.hphob_points();
  const HbondPoints& S_hbonds = dbase.hbond_points();
  const HphobPoints& S_hphobs = dbase.hphob_points();

  my_float_t *cp = A_closest_pts;
  my_float_t *cp_dist = A_pt_dists;

  // Hbond counts
  for(hbond_fit_pt_vci M = M_hbonds.fit_pts_beg(); M < M_hbonds.fit_pts_end(); 
      ++M, ++M_hit, cp += 3, ++cp_dist){

    hbond_fit_pt_vec::float_const_iter_map close_S_hbonds;
    S_hbonds.close_fit_pts(M->pos, A_max_dist, &close_S_hbonds);
    hbond_fit_pt_vci best_pt = S_hbonds.fit_pts_end();
    hbond_fit_pt_vci closest_pt = S_hbonds.fit_pts_end();
    // Let the first idx (0) be for complementary AA or DD matches 
    // and the second idx (1) be for complementary doneptor - polar matches
    // and the third idx (2) be for acceptor-donor matches
    my_float_t best_scores[] = {-2.0, -2.0, -2.0};
    my_float_t closest_scores[] = {-2.0, -2.0, -2.0};
    my_float_t shortest_sq_dists[] = 
      {A_max_dist*A_max_dist, A_max_dist*A_max_dist, A_max_dist*A_max_dist};
    hbond_fit_pt_vec::float_const_iter_map::const_iterator close_pt;
    for(close_pt = close_S_hbonds.begin(); close_pt != close_S_hbonds.end();
        ++close_pt){
      my_float_t dotprod = dot(M->dir, close_pt->second->dir);
      my_float_t sq_dist = dist_squared(M->pos, close_pt->second->pos);

      int idx = 0;
      // Acceptor-donor mismatch
      if((close_pt->second->act_type == DONOR && M->act_type == ACCEPTOR) ||
         (close_pt->second->act_type == ACCEPTOR && M->act_type == DONOR))
        idx = 2;
      // Doneptor match
      else if(close_pt->second->act_type == DONEPTOR || M->act_type == DONEPTOR)
        idx = 1;
 
      if(dotprod > best_scores[idx]) best_scores[idx] = dotprod;
      if(sq_dist < shortest_sq_dists[idx]){
        shortest_sq_dists[idx] = sq_dist;
        closest_scores[idx] = dotprod;
        if(idx < 2) closest_pt = close_pt->second;
      }
    }

    // Best scoring match in the neighborhood
    if(best_scores[0] > 0.0 || best_scores[1] > 0.0){
      if(best_scores[0] > best_scores[1]){
        best_AA_DD_sum += best_scores[0];
        best_polar_sum += best_scores[0];
      }else{
        best_doneptor_sum += best_scores[1];
        best_polar_sum += best_scores[1];
      }
    }else if(best_scores[2] > 0.0) best_polar_mismatch_sum += best_scores[2];
    else ++unsat_polar_count;
         
    // Closest polar point to neighborhood center
    if(shortest_sq_dists[0] < A_max_dist*A_max_dist || 
       shortest_sq_dists[1] < A_max_dist*A_max_dist){
      if(shortest_sq_dists[0] < shortest_sq_dists[1]){
        if(closest_scores[0] > 0.0){
          AA_DD_sum += closest_scores[0]; 
          closest_polar_sum += closest_scores[0]; 
          *M_hit = true;
          std::copy(closest_pt->pos, closest_pt->pos + 3, cp);
          *cp_dist = shortest_sq_dists[0];
        }
      }else if(closest_scores[1] > 0.0){
        doneptor_sum += closest_scores[1]; 
        closest_polar_sum += closest_scores[1]; 
        *M_hit = true;
        std::copy(closest_pt->pos, closest_pt->pos + 3, cp);
        *cp_dist = shortest_sq_dists[1];
      }
    }else if(shortest_sq_dists[2] < A_max_dist*A_max_dist &&
             closest_scores[2] > 0.0) polar_mismatch_sum += closest_scores[2];
  }

  // Hydrophobic Counts
  int my_cnt = 0;
  for(hphob_point_lci H = M_hphobs.begin(); H != M_hphobs.end(); 
      ++H, ++M_hit, cp += 3, ++cp_dist){
    my_float_t d_phil, d_phob;
    S_hbonds.closest_fit_pt(H->pos, &d_phil);
    hphob_point_lci hphob_cp = S_hphobs.closest_point(H->pos, &d_phob);

    if(d_phob <= A_max_dist && d_phob <= d_phil){
      ++phobic_count;
      *M_hit = true;
      std::copy(hphob_cp->pos, hphob_cp->pos + 3, cp);
      *cp_dist = d_phob;
    }
    ++ my_cnt;
  }

  // Scale the terms if desired
  if(scale_terms()){
    std::vector<my_float_t>::const_iterator max_val =
      model.get_max_feature_vals().begin();
    std::vector<my_float_t>::iterator t; 
    for(t = terms.begin(); t < terms.end(); ++t, ++max_val) *t /= *max_val;

    return SCALED_CONSTANT_TERM + 
           SCALED_POLAR_SUM_W * closest_polar_sum + 
           SCALED_POLAR_MISMATCH_W * polar_mismatch_sum + 
           SCALED_DD_OR_AA_POLAR_SUM_W * AA_DD_sum +
           SCALED_DONEPTOR_POLAR_SUM_W * doneptor_sum +
           SCALED_HPHOBIC_COUNT_W * phobic_count + 
           SCALED_UNSAT_POLAR_W * unsat_polar_count;
  }else{
    return CONSTANT_TERM + 
           POLAR_SUM_W * closest_polar_sum + 
           POLAR_MISMATCH_W * polar_mismatch_sum + 
           DD_OR_AA_POLAR_SUM_W * AA_DD_sum +
           DONEPTOR_POLAR_SUM_W * doneptor_sum +
           HPHOBIC_COUNT_W * phobic_count + 
           UNSAT_POLAR_W * unsat_polar_count;
  }
}
#endif

void
WeightedSumsScore::init_storage(const ModelSitemap &model)
{
  const HbondPoints& M_hbonds = model.hbond_points();
  const HphobPoints& M_hphobs = model.hphob_points();
  const size_t N = M_hbonds.fit_pts_size() + M_hphobs.size();
  if(N <= A_max_num_pts) return;

  clear_mem();
  A_closest_pts = 0;
  A_pt_dists = 0;
  A_max_num_pts = N;

  A_closest_pts = new my_float_t[4*N];
  A_pt_dists = &A_closest_pts[3*N];
}

my_float_t WeightedSumsScore::
correspondences(const ModelSitemap &model, my_float_t **Q_pts_ptr, 
                my_float_t **D_pts_ptr, size_t *npts)
{
  const HbondPoints& M_hbonds = model.hbond_points();
  const HphobPoints& M_hphobs = model.hphob_points();
  const size_t N = M_hbonds.fit_pts_size() + M_hphobs.size();
  *Q_pts_ptr = new my_float_t[3*N];
  *D_pts_ptr = new my_float_t[3*N];
  *npts = 0;

  my_float_t SE = 0.0;
  my_float_t *Q_out = *Q_pts_ptr;
  my_float_t *D_out = *D_pts_ptr;
  const my_float_t *cp = A_closest_pts;
  const my_float_t *cp_dist = A_pt_dists;

  hbond_fit_pt_vci M;
  for(M = M_hbonds.fit_pts_beg(); M < M_hbonds.fit_pts_end(); ++M){
    if(*cp_dist <= A_max_dist){
      std::copy(M->pos, M->pos + 3, Q_out);
      std::copy(cp, cp + 3, D_out);
      Q_out += 3;
      D_out += 3;
      SE += (*cp_dist) * (*cp_dist);
      *npts += 1;
    }else SE += A_max_dist * A_max_dist;
   
    cp += 3;
    ++cp_dist;
  }

  for(hphob_point_lci H = M_hphobs.begin(); H != M_hphobs.end(); ++H){
    if(*cp_dist <= A_max_dist){
      std::copy(H->pos, H->pos + 3, Q_out);
      std::copy(cp, cp + 3, D_out);
      Q_out += 3;
      D_out += 3;
      SE += (*cp_dist) * (*cp_dist);
      *npts += 1;
    }else SE += A_max_dist * A_max_dist;
    
    cp += 3;
    ++cp_dist;
  }

  return std::sqrt(SE/N);
}

void WeightedSumsScore::
polar_correspondences(const ModelSitemap &model,
                      std::vector<const my_float_t *> *q_pts,
                      std::vector<const my_float_t *> *corr_pts) const
{
  const HbondPoints& M_hbonds = model.hbond_points();
  const HphobPoints& M_hphobs = model.hphob_points();
  const size_t N = M_hbonds.fit_pts_size() + M_hphobs.size();
  q_pts->resize(N);
  corr_pts->resize(N);

  const my_float_t *cp = A_closest_pts;
  const my_float_t *cp_dist = A_pt_dists;
  std::vector<const my_float_t *>::iterator q_pts_iter = q_pts->begin();
  std::vector<const my_float_t *>::iterator corr_pts_iter = corr_pts->begin();

  hbond_fit_pt_vci M;
  for(M = M_hbonds.fit_pts_beg(); M < M_hbonds.fit_pts_end(); ++M){
    if(*cp_dist <= A_max_dist){
      *q_pts_iter = M->pos;
      *corr_pts_iter = cp;
    }else{
      *q_pts_iter = 0;
      *corr_pts_iter = 0;
    }
    cp += 3;
    ++cp_dist;
    ++q_pts_iter;
    ++corr_pts_iter;
  }

  for(hphob_point_lci H = M_hphobs.begin(); H != M_hphobs.end(); ++H){
    if(*cp_dist <= A_max_dist){
      *q_pts_iter = H->pos;
      *corr_pts_iter = cp;
    }else{
      *q_pts_iter = 0;
      *corr_pts_iter = 0;
    }
    cp += 3;
    ++cp_dist;
    ++q_pts_iter;
    ++corr_pts_iter;
  }
}
