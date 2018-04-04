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

/*
 * $Source: /psa/share/repository/pfizer_proj/src/search/EllipsoidalHbondsScore.C,v $
 * $Revision: 1.3 $
 * $Author: vanvoor4 $
 * $Date: 2008-01-04 15:50:04 $
 * 
 * $Log: not supported by cvs2svn $
 * Revision 1.2  2007/11/01 16:36:20  vanvoor4
 * Now pass in the parameters obj instead of those needed
 *
 * Revision 1.1  2007/10/11 16:21:25  vanvoor4
 * Initial checkin
 *
 * 
 * 
 */

#include <EllipsoidalHbondsScore.H>
#include <math_basics.H>
#include <iostream>

using namespace SimSite3D;

const std::string EllipsoidalHbondsScore::_fname = "EllipsoidalHbondsScore.C";
// Weights are not set in stone at this point.


EllipsoidalHbondsScore::EllipsoidalHbondsScore(ModelSitemap* model_in, 
                                               const SearchParameters& 
                                               params)
 : ScoreRigidAlignments<std::greater<my_float_t> >(model_in, params)
{
 
}

EllipsoidalHbondsScore::~EllipsoidalHbondsScore()
{

}
  

my_float_t 
EllipsoidalHbondsScore::score(const DbaseSitemap& search, rigid_align_vi scores)
{
  scores->point_count = 0;
  std::fill(scores->match_print.begin(), scores->match_print.end(), false);
  //std::vector<bool>::iterator M_hit = scores->match_print.begin();

  const HbondPoints& M_hbonds = model->hbond_points();
  //const HphobPoints& M_hphobs = model->hphob_points();
  const HbondPoints& S_hbonds = search.hbond_points();
  //const HphobPoints& S_hphobs = search.hphob_points();

  scores->terms.clear();

  // Compute the simple approximation to "overlap" of hbonding volumes
  // Strategy is highest scoring pair wins all, and remove both the model and
  // search point from consideration -- rinse and repeat

  // Note: not sure on efficiency here
  std::multimap<my_float_t, S_and_M> interact_pairs;
  S_and_M pt_pair;
  hbond_ideal_pt_vci M;
  for(M = M_hbonds.ideal_pts_beg(); M < M_hbonds.ideal_pts_end(); ++M){
    pt_pair.M = M;
    hbond_ideal_pt_vci S;
    for(S = S_hbonds.ideal_pts_beg(); S < S_hbonds.ideal_pts_end(); ++S){
      pt_pair.S = S;
      my_float_t W = ellipsoidal_hbonds(M, S);
      if(W < my_float_max)
        interact_pairs.insert(std::pair<my_float_t, S_and_M>(W, pt_pair));
    }
  }
  
  my_float_t AA_sum = 0;
  my_float_t DD_sum = 0;
  my_float_t AD_sum = 0;
  my_float_t Nstar_sum = 0;  // Doneptor to any 
  size_t num_M_pts = (M_hbonds.ideal_pts_end() - M_hbonds.ideal_pts_beg());
  size_t num_S_pts = (S_hbonds.ideal_pts_end() - S_hbonds.ideal_pts_beg());
  std::map<hbond_ideal_pt_vci, bool> used_M_pts;
  std::map<hbond_ideal_pt_vci, bool> used_S_pts;

  std::multimap<my_float_t, S_and_M>::reverse_iterator P;
  for(P = interact_pairs.rbegin(); P != interact_pairs.rend(); ++P){
    hbond_ideal_pt_vci M_pt = P->second.M;
    hbond_ideal_pt_vci S_pt = P->second.S;
    if(used_M_pts.find(M_pt) != used_M_pts.end()  
       || used_S_pts.find(S_pt) != used_S_pts.end()) continue;

    used_M_pts[M_pt] = true; 
    used_S_pts[S_pt] = true; 
    if(M_pt->act_type == DONEPTOR || S_pt->act_type == DONEPTOR) 
      Nstar_sum += P->first;
    else if(M_pt->act_type == ACCEPTOR && S_pt->act_type == ACCEPTOR) 
      AA_sum += P->first;
    else if(M_pt->act_type == DONOR && S_pt->act_type == DONOR) 
      DD_sum += P->first;
    else AD_sum += P->first;

    if(used_M_pts.size() >= num_M_pts || used_S_pts.size() >= num_S_pts) break;
  }
  scores->AA_sum = AA_sum;
  scores->DD_sum = DD_sum;
  scores->AD_sum = AD_sum;
  scores->Nstar_sum = Nstar_sum; 
  scores->score = AA_sum + DD_sum + Nstar_sum - AD_sum;
  return scores->score;
}

my_float_t
EllipsoidalHbondsScore::ellipsoidal_hbonds(hbond_ideal_pt_vci M,  
                                           hbond_ideal_pt_vci S)
{
  // Vector from the model to source point
  my_float_t V_ms[3];
  for(uint i = 0; i < 3; ++i) V_ms[i] = S->pos[i] - M->pos[i];
  // square of the norm of V_ms
  my_float_t norm_2 = 0;
  for(uint i = 0; i < 3; ++i) norm_2 += V_ms[i] * V_ms[i];
  if(norm_2 > 4.0) return my_float_max;
  // Out of plane distance
  my_float_t oopd = 0;
  for(uint i = 0; i < 3; ++i) oopd += V_ms[i] * M->dir[i];
  oopd = (oopd >= 0 ? oopd : -1 * oopd);
  if(oopd > 1) return my_float_max;
  // In plane distance
  my_float_t ipd = std::sqrt(-1 * oopd*oopd + norm_2);
  if(ipd > 2) return my_float_max;
 
  // Check oopd and ipd so that we do not have divison by small numbers or 0 
  my_float_t ratio = 0;
  if(norm_2 < 0.0625) ratio = 0;
  else if(oopd < 10E-08) ratio = ipd / 2;
  else if(ipd < 10E-08) ratio = oopd / 1;
  else{
    my_float_t frac = ipd/oopd;
    ratio = oopd / std::sqrt(1 / (1 + 0.25 * frac*frac));
  }

  // Block function to have "full" strength if the distance M to S is less
  // than half of the distance to the ellipse 0.25 x^2 + y^2 = 1 which is 
  // moved to M and has the minor axis parallel to the normal M.dir and the 
  // major axis perpendicular to M.dir and perpendicular to M.dir cross S.dir.
  my_float_t rv = 0;
  if(ratio <= 0.5) rv = 1;
  else rv = 1 - 2 * (ratio - 0.5);
  my_float_t M_dot_S = 0;
  for(uint i = 0; i < 3; ++i) M_dot_S += M->dir[i] * S->dir[i];
  rv *= M_dot_S;
  return (rv > 0 ? rv : 0);
}

/*
  size_t num_M_pts = (M_hbonds.ideal_pts_end() - M_hbonds.ideal_pts_beg());
  size_t num_S_pts = (S_hbonds.ideal_pts_end() - S_hbonds.ideal_pts_beg());
  std::map<hbond_ideal_pt_vci, bool> used_M_pts;
  std::map<hbond_ideal_pt_vci, bool> used_S_pts;

  std::multimap<my_float_t, S_and_M>::const_iterator P;
  for(P = dist_pairs.begin(); P != dist_pairs.end(); ++P){
    if(used_M_pts.find(P->second.M) != used_M_pts.end()
       || used_S_pts.find(P->second.S) != used_S_pts.end()) continue;

    used_M_pts[P->second.M] = true;
    used_S_pts[P->second.S] = true;
    ellipsoidal_hbonds(P->second.M, P->second.S);

    if(used_M_pts.size() >= num_M_pts || used_S_pts.size() >= num_S_pts) break;
  }

*/
