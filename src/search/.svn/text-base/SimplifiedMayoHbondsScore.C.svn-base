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
 * $Source: /psa/share/repository/pfizer_proj/src/search/SimplifiedMayoHbondsScore.C,v $
 * $Revision: 1.1 $
 * $Author: vanvoor4 $
 * $Date: 2008-05-13 16:48:25 $
 * 
 * $Log: not supported by cvs2svn $
 * 
 */

#include <SimplifiedMayoHbondsScore.H>
#include <math_basics.H>
#include <iostream>
#include <mat_ops.H>

using namespace SimSite3D;

const std::string SimplifiedMayoHbondsScore::A_fname = "SimplifiedMayoHbondsScore.C";
const my_float_t SimplifiedMayoHbondsScore::V_zero = 8.0;
const my_float_t SimplifiedMayoHbondsScore::d_zero = 2.8;
const my_float_t SimplifiedMayoHbondsScore::d_zero_squared = d_zero*d_zero;

SimplifiedMayoHbondsScore::SimplifiedMayoHbondsScore(ModelSitemap* model_in, 
                                               const SearchParameters& params)
 : ScoreRigidAlignments<std::less<my_float_t> >(model_in, params)
{
 
}

SimplifiedMayoHbondsScore::~SimplifiedMayoHbondsScore()
{

}
  

my_float_t 
SimplifiedMayoHbondsScore::score(const DbaseSitemap& search, rigid_align_vi scores)
{
  scores->point_count = 0;
  std::fill(scores->match_print.begin(), scores->match_print.end(), false);

  const HbondPoints& M_hbonds = model->hbond_points();
  const HbondPoints& S_hbonds = search.hbond_points();

  scores->terms.clear();

  // Note: not sure on efficiency here
  std::multimap<my_float_t, S_and_M> interact_pairs;
  S_and_M pt_pair;
  hbond_ideal_pt_vci M;
  for(M = M_hbonds.ideal_pts_beg(); M < M_hbonds.ideal_pts_end(); ++M){
    pt_pair.M = M;
    hbond_ideal_pt_vci S;
    for(S = S_hbonds.ideal_pts_beg(); S < S_hbonds.ideal_pts_end(); ++S){
      pt_pair.S = S;
      my_float_t W = simp_mayo_fun(M, S);
      // Assumption is an "hbond" with energy less than -0.01 is too weak
      if(W < -0.01)
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
  scores->terms.resize(num_M_pts, 0);

  // Strategy is highest scoring pair wins all, and remove both the model and
  // search point from consideration -- rinse and repeat

  //std::multimap<my_float_t, S_and_M>::reverse_iterator P;
  std::multimap<my_float_t, S_and_M>::iterator P;
  for(P = interact_pairs.begin(); P != interact_pairs.end(); ++P){
    hbond_ideal_pt_vci M_pt = P->second.M;
    hbond_ideal_pt_vci S_pt = P->second.S;
    if(used_M_pts.find(M_pt) != used_M_pts.end()  
       || used_S_pts.find(S_pt) != used_S_pts.end()) continue;
 
    used_M_pts[M_pt] = true; 
    used_S_pts[S_pt] = true; 

    scores->terms[M_pt - M_hbonds.ideal_pts_beg()] = P->first;
    if(M_pt->act_type == DONEPTOR || S_pt->act_type == DONEPTOR) 
      Nstar_sum += P->first;
    else if(M_pt->act_type == ACCEPTOR && S_pt->act_type == ACCEPTOR) 
      AA_sum += P->first;
    else if(M_pt->act_type == DONOR && S_pt->act_type == DONOR) 
      DD_sum += P->first;
    else{
      AD_sum += P->first;
      scores->terms[M_pt - M_hbonds.ideal_pts_beg()] *= -1;
    }

    if(used_M_pts.size() >= num_M_pts || used_S_pts.size() >= num_S_pts) break;
  }
  scores->AA_sum = AA_sum;
  scores->DD_sum = DD_sum;
  scores->AD_sum = AD_sum;
  scores->Nstar_sum = Nstar_sum; 
  scores->score = AA_sum + DD_sum + Nstar_sum - AD_sum;
  return scores->score;
}

// This is getting to be a pain, besides we don't have a feeling for how
// much the preacceptor-acceptor-H angle should affect the strength of the
// hbond.  It seems unlikely that 180 degrees would be the best.
// We also want a more smooth function that is a bit more soft on 
// theta (DHA angle), so we throw out the exponential part of the term.
my_float_t
SimplifiedMayoHbondsScore::simp_mayo_fun(hbond_ideal_pt_vci M,  
                                         hbond_ideal_pt_vci S)
{
  const atom_vci& M_atom = M->atom;

/*
    std::cout << PDB_residues::residue_to_string(M->atom->res) 
              << M->atom->res_num << " "
              << PDB_residues::atom_to_string(M->atom->name)
              << " (" << M->atom->atom_num << ") -----  ";
    std::cout << PDB_residues::residue_to_string(S->atom->res) 
              << S->atom->res_num << " "
              << PDB_residues::atom_to_string(S->atom->name)
              << " (" << S->atom->atom_num << ") ";
*/

  // Query prot donor atom OR
  // Search prot donor atom OR
  // Query and search prot doneptor atoms
  if(M->act_type == ACCEPTOR || S->act_type == ACCEPTOR ||
     (M->act_type == DONEPTOR && S->act_type == DONEPTOR)){
    // dir is M->dir
    my_float_t my_O[3];
    std::copy(M_atom->pos, M_atom->pos + 3, my_O);
    my_axpy(3, 2.8, M->dir, 1, my_O, 1);
    
    // If the distance is greater than 4.0 (A) ignore it
    my_float_t d_squared = dist_squared(S->atom->pos, my_O);
    if(d_squared > 16.0) return my_float_max;

    // Compute theta
    my_float_t S_H_pos[3];
    std::copy(S->atom->pos, S->atom->pos + 3, S_H_pos);
    my_axpy(3, 1.0, S->dir, 1, S_H_pos, 1);
    my_float_t O_to_H_dir[3];
    for(size_t i = 0; i < 3; ++i) O_to_H_dir[i] = S_H_pos[i] - my_O[i];
    normalize(O_to_H_dir);
    // If cos_theta > 0, ignore this match
    my_float_t cos_theta = dot(O_to_H_dir, S->dir); 
    if(cos_theta >= 0) return my_float_max;

    // Because of numerical issues, cos(theta) can be slightly less than -1.0
    // which will cause std::acos to spit out NaN
    if(cos_theta < -1) cos_theta = -1;

    // If cos_phi > 0, ignore this match
    my_float_t H_to_O_dir[3];   // Dir from search H to fictitious carbonyl O
    my_float_t C_to_O_dir[3];   // Direction of C to carbonyl O bond
    for(size_t i = 0; i < 3; ++i){
      H_to_O_dir[i] = -1*O_to_H_dir[i];
      C_to_O_dir[i] = -1*M->dir[i];
    } 
    my_float_t cos_phi = dot(H_to_O_dir, C_to_O_dir);
    if(cos_phi >= 0) return my_float_max;

    // Attempt to make a "flat well" 
    if(d_squared < d_zero_squared){
      my_float_t d = std::sqrt(d_squared);
      d += d_zero - d;
      d_squared = d*d;
    }

    my_float_t my_ratio = d_zero_squared/d_squared;
    my_float_t my_ratio_5th = std::pow(my_ratio, 5);
    my_float_t my_ratio_6th = my_ratio_5th * my_ratio;
    return V_zero * (5*my_ratio_6th - 6*my_ratio_5th) * cos_theta*cos_theta;
  }

  // Query prot acceptor atom OR
  // Search prot acceptor atom
  if(M->act_type == DONOR || S->act_type == DONOR){
    // We are laying down a fictitious amino group
    // N to H is the opposite direction as M->dir
    
    my_float_t my_H[3], my_N[3];
    std::copy(M->atom->pos, M->atom->pos + 3, my_H);
    std::copy(M->atom->pos, M->atom->pos + 3, my_N);
    my_axpy(3, 1.8, M->dir, 1, my_H, 1); 
    my_axpy(3, 2.8, M->dir, 1, my_N, 1); 

    // If the distance is greater than 4.0 (A) ignore it
    my_float_t d_squared = dist_squared(S->atom->pos, my_N);
    if(d_squared > 16.0) return my_float_max;

    // H to S acceptor dir
    my_float_t H_to_S_Acc_dir[3];
    for(size_t i = 0; i < 3; ++i) H_to_S_Acc_dir[i] = S->atom->pos[i] - my_H[i];
    normalize(H_to_S_Acc_dir);

    my_float_t cos_theta = dot(M->dir, H_to_S_Acc_dir);
    if(cos_theta >= 0) return my_float_max;
    // Because of numerical issues, cos(theta) can be slightly less than -1.0
    // which will cause std::acos to spit out NaN
    if(cos_theta < -1) cos_theta = -1;

    // Compute the F() value
    my_float_t F = cos_theta*cos_theta;

    // SP2 or CO2 (CO2 is given to ASP and GLU sidechain Os)
    if(S->atom->orbit == SP2 || S->atom->orbit == CO2){
      my_float_t Pre_to_S_Acc_dir[3];
      std::copy(S->atom->pos, S->atom->pos + 3, Pre_to_S_Acc_dir);
      my_axpy(3, -1.0, S->carbon_nbr->pos, 1, Pre_to_S_Acc_dir, 1);
      normalize(Pre_to_S_Acc_dir);
      my_float_t cos_phi = dot(H_to_S_Acc_dir, Pre_to_S_Acc_dir);
      if(cos_phi >= 0) return my_float_max;

      // the estimate of the out of plane angle is given by the angle
      // between the normal of the plane given by the search neighbor C,
      // search Acceptor, and my_N to some other normal -- 

    // A do nothing if statement from the "full" function
    }else if(S->atom->orbit == SP3);

    // Attempt to make a "flat well" 
    if(d_squared < d_zero_squared){
      my_float_t d = std::sqrt(d_squared);
      d += d_zero - d;
      d_squared = d*d;
    }

    my_float_t my_ratio = d_zero_squared/d_squared;
    my_float_t my_ratio_5th = std::pow(my_ratio, 5);
    my_float_t my_ratio_6th = my_ratio_5th * my_ratio;
    return V_zero * (5*my_ratio_6th - 6*my_ratio_5th) * F;
  }
  return my_float_max;
}

/*
my_float_t
SimplifiedMayoHbondsScore::simp_mayo_fun(hbond_ideal_pt_vci M,  
                                         hbond_ideal_pt_vci S)
{
  // For now we keep things very simple and use
  // F(theta, phi, gamma) = cos^2(theta)exp(-(pi - theta)^6) cos^2 (phi)
  // in the case of M atom is donor
  // If the query protein atom is a donor, then set down O-C inline 
  // If the query protein atom is an acceptor, then set down H-N
  // If doneptor, then check what the search atom type is first
  const atom_vci& M_atom = M->atom;

    std::cout << PDB_residues::residue_to_string(M->atom->res) 
              << M->atom->res_num << " "
              << PDB_residues::atom_to_string(M->atom->name)
              << " (" << M->atom->atom_num << ") -----  ";
    std::cout << PDB_residues::residue_to_string(S->atom->res) 
              << S->atom->res_num << " "
              << PDB_residues::atom_to_string(S->atom->name)
              << " (" << S->atom->atom_num << ") ";

  // Query prot donor atom OR
  // Search prot donor atom OR
  // Query and search prot doneptor atoms
  if(M->act_type == ACCEPTOR || S->act_type == ACCEPTOR ||
     (M->act_type == DONEPTOR && S->act_type == DONEPTOR)){
    // dir is M->dir
    my_float_t my_O[3];
    std::copy(M_atom->pos, M_atom->pos + 3, my_O);
    my_axpy(3, 2.8, M->dir, 1, my_O, 1);
    
    // If the distance is greater than 4.0 (A) ignore it
    my_float_t d_squared = dist_squared(S->atom->pos, my_O);
    if(d_squared > 16.0) return my_float_max;

    // Compute theta
    my_float_t S_H_pos[3];
    std::copy(S->atom->pos, S->atom->pos + 3, S_H_pos);
    my_axpy(3, 1.0, S->dir, 1, S_H_pos, 1);
    my_float_t O_to_H_dir[3];
    for(size_t i = 0; i < 3; ++i) O_to_H_dir[i] = S_H_pos[i] - my_O[i];
    normalize(O_to_H_dir);
    // If cos_theta > 0, ignore this match
    my_float_t cos_theta = dot(O_to_H_dir, S->dir); 
    if(cos_theta >= 0) return my_float_max;

    // Because of numerical issues, cos(theta) can be slightly less than -1.0
    // which will cause std::acos to spit out NaN
    if(cos_theta < -1) cos_theta = -1;

    // If cos_phi > 0, ignore this match
    my_float_t H_to_O_dir[3];   // Dir from search H to fictitious carbonyl O
    my_float_t C_to_O_dir[3];   // Direction of C to carbonyl O bond
    for(size_t i = 0; i < 3; ++i){
      H_to_O_dir[i] = -1*O_to_H_dir[i];
      C_to_O_dir[i] = -1*M->dir[i];
    } 
    my_float_t cos_phi = dot(H_to_O_dir, C_to_O_dir);
    if(cos_phi >= 0) return my_float_max;

    my_float_t my_ratio = d_zero_squared/d_squared;
    return V_zero * (5*std::pow(my_ratio, 6) - 6*std::pow(my_ratio, 5)) 
      * cos_theta*cos_theta  
      * std::exp(-1*std::pow(M_PI - std::acos(cos_theta), 6)) 
      * cos_phi*cos_phi;
  }

  // Query prot acceptor atom OR
  // Search prot acceptor atom
  if(M->act_type == DONOR || S->act_type == DONOR){
    // We are laying down a fictitious amide nitrogen
    // N to H is the opposite direction as M->dir
    
    my_float_t my_H[3], my_N[3];
    std::copy(M->atom->pos, M->atom->pos + 3, my_H);
    std::copy(M->atom->pos, M->atom->pos + 3, my_N);
    my_axpy(3, 1.8, M->dir, 1, my_H, 1); 
    my_axpy(3, 2.8, M->dir, 1, my_N, 1); 

    // If the distance is greater than 4.0 (A) ignore it
    my_float_t d_squared = dist_squared(S->atom->pos, my_N);
    if(d_squared > 16.0) return my_float_max;

    // H to S acceptor dir
    my_float_t H_to_S_Acc_dir[3];
    for(size_t i = 0; i < 3; ++i) H_to_S_Acc_dir[i] = S->atom->pos[i] - my_H[i];
    normalize(H_to_S_Acc_dir);

    my_float_t cos_theta = dot(M->dir, H_to_S_Acc_dir);
    if(cos_theta >= 0) return my_float_max;
    // Because of numerical issues, cos(theta) can be slightly less than -1.0
    // which will cause std::acos to spit out NaN
    if(cos_theta < -1) cos_theta = -1;

    // Compute the F() value
    my_float_t F = cos_theta*cos_theta 
      * std::exp(-1*std::pow(M_PI - acos(cos_theta), 6));

    // SP2 or CO2 (CO2 is given to ASP and GLU sidechain Os)
    if(S->atom->orbit == SP2 || S->atom->orbit == CO2){
      my_float_t Pre_to_S_Acc_dir[3];
      std::copy(S->atom->pos, S->atom->pos + 3, Pre_to_S_Acc_dir);
      my_axpy(3, -1.0, S->carbon_nbr->pos, 1, Pre_to_S_Acc_dir, 1);
      normalize(Pre_to_S_Acc_dir);
      my_float_t cos_phi = dot(H_to_S_Acc_dir, Pre_to_S_Acc_dir);
      if(cos_phi >= 0) return my_float_max;

      // the estimate of the out of plane angle is given by the angle
      // between the normal of the plane given by the search neighbor C,
      // search Acceptor, and my_N to some other normal


      F *= cos_phi*cos_phi;
    }else if(S->atom->orbit == SP3) F *= F;

    my_float_t my_ratio = d_zero_squared/d_squared;
    return V_zero * (5*std::pow(my_ratio, 6) - 6*std::pow(my_ratio, 5)) * F;
  }
  return my_float_max;
}
*/
