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
 * Authors: Jeffrey Van Voorst, jeff.vanvoorst@gmail.com
 *          Leslie Kuhn, Ph.D., KuhnL@msu.edu 
 *****************************************************************************/

/*
 * $Source: /psa/share/repository/pfizer_proj/src/gen_points/HbondPoints.C,v $
 * $Revision: 1.16 $
 * $Author: vanvoor4 $
 * $Date: 2008/07/23 14:22:19 $
 * 
 * $Log: HbondPoints.C,v $
 * Revision 1.16  2008/07/23 14:22:19  vanvoor4
 * Changed int to size_t in a loop and added some comments to
 * the header file.
 *
 * Revision 1.15  2008/05/15 17:27:54  vanvoor4
 * Changed stuff around and made somethings more modular.
 *
 * Revision 1.14  2008/03/31 18:00:22  vanvoor4
 * Updated to partially support metal sitemap points.  Also,
 * some of the code and names were changed to better reflect their current
 * use or meaning.
 *
 * Revision 1.13  2008/01/04 21:40:03  vanvoor4
 * Had forgotten to remove print statement.
 *
 * Revision 1.12  2008/01/04 21:24:55  vanvoor4
 * Added check for metals to allow template points to be closer
 * to metals than prot atoms.
 *
 * Revision 1.11  2008/01/04 18:21:49  vanvoor4
 * Removed extra print statements
 *
 * Revision 1.10  2008/01/04 18:12:46  vanvoor4
 * Changed handling of water molecules to be better suited for
 * users.
 *
 * Revision 1.9  2007/12/17 21:28:51  vanvoor4
 * Added support for water and metal points
 *
 * Revision 1.8  2007/11/01 16:29:56  vanvoor4
 * Paths now are not expected to end with a trailing slash '/'.
 *
 * Revision 1.7  2007/09/24 15:51:00  vanvoor4
 * Changed from generating several rotations to the BKP Horn
 * 3pt alignment using 2 rotations.
 *
 * Revision 1.6  2007/08/29 20:32:30  vanvoor4
 * Updated to respect the handling of hbond points and standarize
 * (somewhat) the accessor methods
 *
 * Revision 1.5  2007/08/21 18:12:57  vanvoor4
 * A large number of changes
 *
 * Revision 1.4  2007/02/07 15:25:17  vanvoor4
 * Had to rename all the point_t and derived typedefs to the hbond_point_t
 * type to avoid conflict with the point_t already in ../basics/types.H.
 * In the future it may be worth the time to look at this more closely
 * and merge the two types if possible.
 * Also removed the PDB namespace and use SimSite3D instead.
 *
 * Revision 1.3  2006/11/16 20:24:25  vanvoor4
 * Removed debug statements
 *
 * Revision 1.2  2006/10/20 13:19:27  vanvoor4
 * Many diff mods
 *
 * Revision 1.1  2006/08/25 17:21:11  vanvoor4
 * Initial checkin
 *
 *
 */

#include <sstream>
#include <cmath>
#include <HbondPoints.H>
#include <basics.H>
#include <mat_ops.H>
#include <algorithm>
#include <iomanip>
#include <PDB_metals.H>

using namespace SimSite3D;

const my_float_t HbondPoints::MAXRADDIST = 9.0;
const my_float_t HbondPoints::MAXBINDDIST = 5.0;
const my_float_t HbondPoints::MIN_VDW_DIST = 2.5;
const my_float_t HbondPoints::MIN_VDW_DIST_METAL = 2.0;
const my_float_t HbondPoints::MIN_VDW_DIST_2 = 
  HbondPoints::MIN_VDW_DIST * HbondPoints::MIN_VDW_DIST;
const my_float_t HbondPoints::MIN_VDW_DIST_METAL_2 = 
  HbondPoints::MIN_VDW_DIST_METAL * HbondPoints::MIN_VDW_DIST_METAL;
const atom_vci HbondPoints::NO_CARBON_NEIGHBOR = atom_t::NULL_ATOM_VCI;
const std::string HbondPoints::_fname = "HbondPoints.C";

HbondPoints::HbondPoints(std::istream& in, const uint num_lines,
                         PDBBase& rad_atoms, PDBBase& points)
{
  std::string line;
  uint prev_atom_num = 0;
  uint curr_atom_num = 1000000000;
  uint prev_obj_num = 0;
  uint curr_obj_num = 1000000000;
  A_fail = true;
  
  hbond_point_t fit_pt;
  std::vector<int> num_fit_pts;
  size_t i;
  for(i = 0; i < num_lines && std::getline(in, line); ++i){
    if(line.substr(0,23) == "</hydrogen_bond_points>"){
      std::stringstream my_sstr;
      my_sstr << "Only found " << i << " hbond points; "
              << "expected " << num_lines;
      warn(_fname, "cstr(site map)", my_sstr.str());
      break;
    }
    std::vector<std::string> toks;
    string_tok(line, &toks, '|');

    uint pt_num;
    my_strtoui(toks[0], &pt_num);        // template point number
    my_strtoui(toks[6], &curr_atom_num); // associated atom number/serial
    my_strtoui(toks[7], &curr_obj_num);  // H, lone pair, etc number
    
    if(prev_atom_num != curr_atom_num || prev_obj_num != curr_obj_num){
      prev_atom_num = curr_atom_num;
      prev_obj_num = curr_obj_num;

      hbond_ideal_point_t ideal_pt;
      std::vector<std::string> coords_str;
      string_tok(toks[8], &coords_str, ' ');
      if(coords_str.size() != 3){ /*error*/}
      for(uint j = 0; j < coords_str.size(); ++j)
        my_strtof(coords_str[j], &(ideal_pt.pos[j]));

      ideal_pt.pt_num = curr_obj_num;
      if(pt_num > 1000 && pt_num < 2001) ideal_pt.act_type = ACCEPTOR;
      else if(pt_num > 2000 && pt_num < 3001) ideal_pt.act_type = DONOR;
      else if(pt_num > 3000 && pt_num < 4001) ideal_pt.act_type = DONEPTOR;
      else{
        // bad hbond point type
      }

      // Prior to 2.7.3 labels files have only 9 columns
      if(toks.size() > 9){
        uint atom_num;
        my_strtoui(toks[9], &atom_num);
        if(atom_num){
          ideal_pt.carbon_nbr = rad_atoms.get_atom(atom_num);
          if(ideal_pt.carbon_nbr == rad_atoms.atoms_end())
            std::cerr << "unable to find a match for atom number " 
                      << atom_num << " in " << rad_atoms.name () << std::endl;
        }else ideal_pt.carbon_nbr = HbondPoints::NO_CARBON_NEIGHBOR;
      }else{
        std::string msg = "\"Old-Style\" sitemaps -- simplifiedmayohbonds ";
        warn(_fname, "HbondPoints(existing_sitemap)", msg + "will not work");
        ideal_pt.carbon_nbr = HbondPoints::NO_CARBON_NEIGHBOR;
      }

      ideal_pt.atom = rad_atoms.get_atom(curr_atom_num);
      if(ideal_pt.atom == rad_atoms.atoms_end()){
        std::cerr << "unable to find a match for atom number " 
                  << curr_atom_num << " in " << rad_atoms.name () << std::endl;
        return;
      }
      unit_vector(ideal_pt.dir, ideal_pt.pos, ideal_pt.atom->pos);
      ideal_pts.push_back(ideal_pt);
      fit_pt.act_type = ideal_pt.act_type;
      num_fit_pts.push_back(0);
    }

    atom_vci point = points.get_atom(pt_num);
    if(point == rad_atoms.atoms_end()){
      std::cerr << "unable to find a match for atom number " 
                << pt_num << " in " << rad_atoms.name () << std::endl;
      return;
    }
    std::copy(point->pos, point->pos + 3, fit_pt.pos);
    fit_pt.atom = ideal_pts.back().atom;
    unit_vector(fit_pt.dir, fit_pt.pos, ideal_pts.back().atom->pos);
    fit_pts.push_back(fit_pt);
    ++(num_fit_pts.back());
  }

  init_iterators(&ideal_pts, &fit_pts, num_fit_pts);
  A_fail = false;

  // pointer is not used
  bounding_volume = 0;
}

HbondPoints::HbondPoints(const GenPointsParameters::hbond_method_t hbond_density, 
                         const std::string param_path, 
                         BoundingVolume* vol_in,
                         const bool compute_volume)
{
  A_fail = true;
  bounding_volume = vol_in;
  load_hbond_angles(hbond_density, param_path);
  build_residue_table();

  my_float_t tmp[3];
  A_compute_volume = compute_volume;
  if(A_compute_volume){
    std::ifstream in;
    if(!open_ifstream(in, param_path + "/spherical_cap.txt")) return;
  
    for(std::string line; std::getline(in, line); ){
      std::vector<std::string> toks;
      string_tok(line, &toks, ' ');
      for(uint i = 0; i < 3; ++i) my_strtof(toks[i], &tmp[i]);
      A_hbond_vol.push_back(tmp);
    }
  }
  A_fail = false;
}

void 
HbondPoints::find_binding_site_hbonds(const PDBStructure* prot,
                                      const std::vector<std::string>& H2O_res,
                                      const bool include_metals, 
                                      const interact_atoms_vec& hbond_atoms,
                                      std::vector<atom_vci>* bind_atoms,
                                      std::vector<atom_vci>* rad_atoms)
{
  // Check if metal atoms are near the binding site
  std::vector<atom_vci> nearby_metals;
  if(include_metals){
    std::vector<atom_vci>::const_iterator metal;
    for(metal = prot->metals_begin(); metal < prot->metals_end(); ++metal){
      if(bounding_volume->BIND_vol_contains((*metal)->pos)){
        bind_atoms->push_back(*metal);
        nearby_metals.push_back(*metal);
      }
      if(bounding_volume->RAD_vol_contains((*metal)->pos))
        rad_atoms->push_back(*metal);
    }
  }

  // Check if the given waters are near the binding site
  std::vector<atom_vci> nearby_waters;
  if(H2O_res.size())
    get_nearby_waters(prot, H2O_res, &nearby_waters, bind_atoms, rad_atoms);

  // Generate the template points for the possible hbonding atoms
  std::vector<int> num_fit_pts;
  interact_atoms_vci a;
  for(a = hbond_atoms.begin(); a < hbond_atoms.end(); ++a)
    get_atom_hbonds(a->chain, a->residue, a->atom, *rad_atoms, &ideal_pts, 
                    &fit_pts, &num_fit_pts);

  if(A_compute_volume) return;

  // Generate the template points for the explicitly specified bound waters
  if(H2O_res.size()){
    WaterPoints water_hbond_pts;
    std::vector<atom_vci>::const_iterator water;
    for(water = nearby_waters.begin(); water < nearby_waters.end(); ++water)
      add_water_points(*water, *rad_atoms, &water_hbond_pts, &ideal_pts, 
                       &fit_pts, &num_fit_pts);
  }

  // Add in template points for metals near the binding site
  // NOTE: do NOT move this before calculating hbonds since the initialization
  // of the iterators for the ideal points depends on the current ordering
  // of the fit points.  
  std::vector<atom_vci>::const_iterator metal;
  for(metal = nearby_metals.begin(); metal < nearby_metals.end(); ++metal)
    add_metal_points(*metal, *rad_atoms, &ideal_pts, &fit_pts, &num_fit_pts);

  init_iterators(&ideal_pts, &fit_pts, num_fit_pts);
}

void
HbondPoints::get_nearby_waters(const PDBStructure* prot,
                               const std::vector<std::string>& H2O_res,
                               std::vector<atom_vci>* nearby_waters,
                               std::vector<atom_vci>* bind_atoms,
                               std::vector<atom_vci>* rad_atoms)
{
  // Build a map of the protein structure's water molecules keyed by residue 
  // name (string)
  std::map<std::string, atom_vci> prot_H2Oz;
  std::vector<atom_vci>::const_iterator w;
  for(w = prot->waters_beg(); w < prot->waters_end(); ++w){
    // Skip hydrogen atoms (only interested at the molecule level here)
    if((*w)->name != O) continue;

    std::ostringstream Wstr;
    if((*w)->chainID != ' ') Wstr << (*w)->chainID;
    Wstr << (*w)->res_num;
    if((*w)->iCode != ' ') Wstr << (*w)->iCode;
    prot_H2Oz[Wstr.str()] = *w;
  }

  for(size_t i = 0; i < H2O_res.size(); ++i){
    std::map<std::string, atom_vci>::const_iterator prot_H2O_iter;
    prot_H2O_iter = prot_H2Oz.find(H2O_res[i]);
    if(prot_H2O_iter == prot_H2Oz.end()){
      std::ostringstream msg;
      msg << "The water molecule " << H2O_res[i] 
          << " does not exist in the protein structure file " << prot->name();
      warn(_fname, "find_binding_site_hbonds", msg.str());
      continue;
    }

    atom_vci water = prot_H2O_iter->second;
    if(bounding_volume->BIND_vol_contains(water->pos)){
      bind_atoms->push_back(water);
      nearby_waters->push_back(water);
    }
    if(bounding_volume->RAD_vol_contains(water->pos))
      rad_atoms->push_back(water);
  }
}

void
HbondPoints::init_iterators(hbond_ideal_pt_vec* ideal_pts_p,
                            hbond_fit_pt_vec* fit_pts_p,
                            std::vector<int> num_fit_pts)
{
  hbond_fit_pt_vi fit_pt = fit_pts_p->begin();
  hbond_ideal_pt_vi ideal_pt = ideal_pts_p->begin(); 
  std::vector<int>::const_iterator num = num_fit_pts.begin();
  for( ; ideal_pt < ideal_pts_p->end() && num < num_fit_pts.end(); 
      ++ideal_pt, ++num){
    ideal_pt->fit_pts_beg = fit_pt;
    for(int j = 0; j < *num; ++j, ++fit_pt) fit_pt->ideal_pt = ideal_pt;
    ideal_pt->fit_pts_end = fit_pt;
  }
}

bool
HbondPoints::get_atom_hbonds(const chain_const_iter chain, 
                             const residue_vci residue,
                             const atom_vci B_atom, 
                             const std::vector<atom_vci>& rad_atoms,
                             hbond_ideal_pt_vec* ideal_pts_p,
                             hbond_fit_pt_vec* fit_pts_p, 
                             std::vector<int>* num_fit_pts)
{
  hbond_data_t tmp_data = {B_atom->name, 0, 0};
  const my_float_t *positions [] = {0, B_atom->pos, 0};
  atom_vci carbon_nbr = residue->atoms_begin;

  // carbonyl oxygen
  if(B_atom->name == O){
    residue_vci next_residue = residue + 1;
    if(next_residue < chain->residues_end){

      // Check for chain break -- kludge at this point since we are not 
      // Keeping track of chain breaks as specified in the coordinate file
      carbon_nbr = residue->get_atom(C);
      atom_vci N_atom = next_residue->get_atom(N);

      // N-C bond length is 1.32 (A) -- allow up to 1.5 (A) to be nice
      if(carbon_nbr != atom_t::NULL_ATOM_VCI &&
         N_atom != atom_t::NULL_ATOM_VCI &&
         2.25 > dist_squared(carbon_nbr->pos, N_atom->pos)){

        positions[0] = carbon_nbr->pos;
        positions[2] = N_atom->pos;
        tmp_data.angles = &(A_O_angles[MAIN_CHAIN]);
        tmp_data.triad = &(hbond_triads[0]);
      }else
        warn(_fname, "get_atom_hbonds()",
             "Break in the chain does not allow computation of hbonds points "
             "for the carbonyl oxygen");
    }else
      warn(_fname, "get_atom_hbonds()",
           "Cannot compute the hbonds points for the carbonyl oxygen; "
           "it is the last residue in the chain");

  // main chain nitrogen
  }else if(B_atom->name == N){
    if(residue != chain->residues_begin){
      residue_vci previous_residue = residue - 1;

      carbon_nbr = previous_residue->get_atom(C);
      atom_vci O_atom = previous_residue->get_atom(O);

      // N-C bond length is 1.32 (A) -- allow up to 1.5 (A) to be nice
      if(carbon_nbr != atom_t::NULL_ATOM_VCI &&
         O_atom != atom_t::NULL_ATOM_VCI &&
         2.25 > dist_squared(carbon_nbr->pos, B_atom->pos)){

        positions[0] = carbon_nbr->pos;
        positions[2] = O_atom->pos;
        tmp_data.angles = &(A_N_angles[MAIN_CHAIN]);
        tmp_data.triad = &(hbond_triads[1]);
      }else
        warn("HbondPoints.C", "get_atom_hbonds",
             "Break in the chain does not allow computation of hbonds points "
             "for the main chain nitrogen");
    }else
      warn("HbondPoints.C", "get_atom_hbonds",
           "Cannot compute the hbonds points for the main chain nitrogen; "
           "it is the first residue in the chain");

  // polar side chains
  }else{
    residue_table_t::const_iterator tbl_iter;
    tbl_iter = hbond_residue_table.find(residue->name);
    if(tbl_iter != hbond_residue_table.end()){
      const std::vector<hbond_data_t> &tmp_vec = tbl_iter->second;
      for(uint i = 0; i < tmp_vec.size(); ++i){
        if(B_atom->name == tmp_vec[i].atom){
          tmp_data.triad = tmp_vec[i].triad;         
          tmp_data.angles = tmp_vec[i].angles; 
          sidechain_nhbr_pos(residue, *(tmp_data.triad), positions, 
                             &carbon_nbr);
          break;
        }
      }
    }else{
      // warn about no hbonds for this residue
      std::cerr << "could not find residue " 
                << PDB_residues::residue_to_string(residue->name) 
                << " in table\n";
    }
  }

  if(tmp_data.angles == 0 || tmp_data.triad == 0 
     || positions[0] == 0 || positions[2] == 0) return false;

  if(A_compute_volume)
    calc_volume(B_atom, carbon_nbr, tmp_data, positions, rad_atoms);
  else
    calc_positions(B_atom, carbon_nbr, tmp_data, positions, rad_atoms, 
                   ideal_pts_p, fit_pts_p, num_fit_pts);
  return true;
}

bool 
HbondPoints::calc_positions(const atom_vci atom, const atom_vci carbon_nbr,
                            const hbond_data_t& hbond_info, 
                            const my_float_t **positions,
                            const std::vector<atom_vci>& rad_atoms,
                            hbond_ideal_pt_vec* ideal_pts_p,
                            hbond_fit_pt_vec* fit_pts_p,
                            std::vector<int>* num_fit_pts)
{
/*
  std::cout << "Calculating the hbond template points for " 
            << PDB_residues::residue_to_string(atom->res)
            << atom->res_num << " " << PDB_residues::atom_to_string(atom->name)
            << "\n";
*/

  // Get transformation to move from origin to B + neighbors A & C
  my_float_t trans_A[3], trans_C[3];
  for(uint i = 0; i < 3; i++){
    trans_A[i] = positions[0][i] - positions[1][i];
    trans_C[i] = positions[2][i] - positions[1][i];
  }
  my_float_t R[9];
  get_local_orientation(trans_A, trans_C, R); 

  if(hbond_info.angles->size() > 19){
    std::cerr << "The number of angles is greater than the static array size\n";
    exit(-1);
  }
  my_float_t offset[60]; // enough for 20 pts
  my_float_t pts[60];    // enough for 20 pts
  angles_vec::const_iterator grp;
  uint pt_num = 1;
  for(grp = hbond_info.angles->begin(); grp < hbond_info.angles->end(); ++grp){
    my_float_t* cur_pt = pts;
    my_float_t* cur_offset = offset;
    
    // Fit points
    uint num_pts;
    for(num_pts = 0; num_pts < grp->alphas.size(); ++num_pts){
      std::copy(positions[1], positions[1] + 3, cur_pt);

    // Reduced the distance of hbond points to 2.5
    //compute_offset(2.5, grp->alphas[num_pts],
                   //grp->betas[num_pts], cur_offset);
    compute_offset(hbond_info.triad->bond_len, grp->alphas[num_pts],
                   grp->betas[num_pts], atom, cur_offset);

      cur_offset += 3;
      cur_pt += 3;
    }
   
    // "Ideal point" 
    ++num_pts;
    std::copy(positions[1], positions[1] + 3, cur_pt);

    // Reduced the distance of hbond points to 2.5
    //compute_offset(2.5, grp->ideal_alpha, grp->ideal_beta, cur_offset);
    compute_offset(hbond_info.triad->bond_len, grp->ideal_alpha,
                   grp->ideal_beta, atom, cur_offset);

    my_gemm(num_pts, 3, 3, 1.0, offset, 3, R, 3, pts, 3, 1.0);
     
    // Add fit points if they are inside and do not bump too much
    uint cnt = 0;
    cur_pt = pts;
    for(uint i = 0; i < num_pts - 1; ++i, cur_pt += 3){
      //std::cout << "Point " << i+1 << " ";
      if(add_fit_point(atom, cur_pt, hbond_info.triad->type, rad_atoms, 
                       fit_pts_p)) ++cnt;
    }

    // if at least 1 fit point was added, add the "ideal point"
    if(cnt) add_ideal_point(atom, carbon_nbr, cur_pt, hbond_info.triad->type, 
                            pt_num, cnt, ideal_pts_p, num_fit_pts);
    ++pt_num;
  }

  return true;
}

void 
HbondPoints::add_water_points(const atom_vci water_O, 
                              const std::vector<atom_vci>& rad_atoms,
                              WaterPoints* water_hbond_pts,
                              hbond_ideal_pt_vec* ideal_pts_p,
                              hbond_fit_pt_vec* fit_pts_p,
                              std::vector<int>* num_fit_pts)
{
  // Get transformation to move from origin to center of water_O and place
  // the vector bisecting the H-O-H angle on the X-axis.
  // Assumes H2 and H1 immediately follow the iterator pointing to the oxygen
  // of the water atom
  my_float_t B[3], trans_H2[3];
  atom_vci water_H2 = water_O + 1;
  atom_vci water_H1 = water_O + 2;
  for(size_t i = 0; i < 3; ++i)
    B[i] = 0.5 * (water_H1->pos[i] + water_H2->pos[i]);
  std::copy(water_H2->pos, water_H2->pos + 3, trans_H2);
  my_axpy(3, -1.0, water_O->pos, 1, B, 1);
  my_axpy(3, -1.0, water_O->pos, 1, trans_H2, 1);
  my_float_t R[9];
  get_local_orientation(B, trans_H2, R); 

  // Transform the water hbond points using the computed rotation and 
  // transformation
  my_float_t trans_R[9];
  for(size_t i = 0; i < 3; ++i)
    for(size_t j = 0; j < 3; ++j) trans_R[3*i + j] = R[3*j + i];
  water_hbond_pts->transform(trans_R, water_O->pos);

  // For each ideal point (2 hydrogen atoms and 2 lone pairs of electrons)
  // if at least one fit point is kept, keep that ideal point.
  uint ideal_pt_num = 1;
  hbond_ideal_pt_vec::const_iterator Wi = water_hbond_pts->ideal_pts_begin(); 
  for( ; Wi < water_hbond_pts->ideal_pts_end(); ++Wi){
    uint cnt = 0;
    hbond_fit_pt_vec::const_iterator fit_pt;
    for(fit_pt = Wi->fit_pts_beg ; fit_pt < Wi->fit_pts_end; ++fit_pt)
      if(add_fit_point(water_O, fit_pt->pos, fit_pt->act_type, rad_atoms, 
                       fit_pts_p)) ++cnt;

    // if at least 1 fit point was added, add the "ideal point"
    if(cnt) add_ideal_point(water_O, NO_CARBON_NEIGHBOR, Wi->pos, Wi->act_type, 
                            ideal_pt_num, cnt, ideal_pts_p, num_fit_pts);
    ++ideal_pt_num;
  }
 
  water_hbond_pts->revert();
}

void 
HbondPoints::add_metal_points(const atom_vci metal, 
                              const std::vector<atom_vci>& rad_atoms,
                              hbond_ideal_pt_vec* ideal_pts_p,
                              hbond_fit_pt_vec* fit_pts_p,
                              std::vector<int>* num_fit_pts)
{
//  std::cout << "Adding metal points for " << PDB_metals::lookup(metal->name) << " " << metal->atom_num << std::endl;
  DiscreteSphere sphere(DISCRETE_SPHERE_LEVEL_TWO);
  sphere.scale(PDB_metals::lookup(metal->name)->interact_pt_rad);
  const my_float_t* orig_pt;
  int cnt = 0;
  for(orig_pt = sphere.pts_begin(); orig_pt < sphere.pts_end(); orig_pt += 3){
    my_float_t metal_pt[3];
    for(size_t i = 0; i < 3; ++i) metal_pt[i] = metal->pos[i] + orig_pt[i]; 
    if(add_fit_point(metal, metal_pt, ACCEPTOR, rad_atoms, fit_pts_p)) ++cnt;
  }

  // if at least 1 fit point was added, add the "ideal point" even though
  // this "ideal point will not be used"
  if(cnt){
    my_float_t no_pos[] = {0, 0, 0};
    add_ideal_point(metal, NO_CARBON_NEIGHBOR, no_pos, ACCEPTOR, 
                    0, cnt, ideal_pts_p, num_fit_pts);
  }
}

void
HbondPoints::compute_offset(const my_float_t bond_len, const my_float_t alpha,
                            const my_float_t beta, const atom_vci atom,
                            my_float_t *offset)
{
  // probably should store cos and sin of angles rather than computing them
  // each time the are run across or it may make even more sense to just
  // store the relative positions 

  // One potential area of confusion is the alignment code places the U and
  // V vectors in the XY plane and aligns U with the X axis.  This means that
  // the "out-of-plane" angle is with respect to the XY plane or Z=0 plane.

  my_float_t cos_alpha = std::cos(alpha); 
  my_float_t cos_beta = std::cos(beta); 
  my_float_t sin_alpha = std::sin(alpha); 
  my_float_t sin_beta = std::sin(beta); 

#if 0
  offset[0] = bond_len * cos_alpha;
  offset[1] = bond_len * cos_beta * sin_alpha, 
  offset[2] = bond_len * sin_beta * sin_alpha; 
#endif

  // We have given Lys NZ an orbital of SP4
  // We treat the OH of Tyr as SP2 since it is more likely to donate than
  // accept and tends to donate in the plane of the phenyl ring
  if(atom->res == TYR || (atom->orbit != SP3 && atom->orbit != SP4)){
    offset[0] = bond_len * cos_beta * cos_alpha;
    offset[1] = bond_len * cos_beta * sin_alpha, 
    offset[2] = -1.0 * bond_len * sin_beta;
  }else{
    offset[0] = bond_len * cos_alpha;
    offset[1] = bond_len * cos_beta * sin_alpha;
    offset[2] = bond_len * sin_beta * sin_alpha;
  }
}

bool
HbondPoints::sidechain_nhbr_pos(residue_vci residue,
                                const hbond_triad_t& triad,
                                const my_float_t **positions,
                                atom_vci* carbon_nbr)
{
  uint nfound = 0;
  atom_vci atom = residue->atoms_begin;
  for( ; nfound < 2 && atom != residue->atoms_end; atom++){
    if(atom->name == triad.nbr_A){
      // The carbon covalently bonded to the hbonder is listed first in the 
      // triad
      *carbon_nbr = atom;
      positions[0] = atom->pos; 
      nfound++;
    }else if(atom->name == triad.nbr_C){
      positions[2] = atom->pos; 
      nfound++;
    }
  }

  if(nfound < 2){
    std::ostringstream msg;
    msg << "Unable to find the neighbors for "
        << PDB_residues::residue_to_string(residue->name) 
        << " " << residue->number << "(" << residue->chainID << ") "
        << PDB_residues::atom_to_string(triad.B);
    warn(_fname, "sidechain_nhbr_pos", msg.str());
    return false;
  }
  return true;
}

void
HbondPoints::build_residue_table()
{
  atom_t atom; 
  for(uint i = 0; i < num_hbond_triads; ++i){
    atom.res = hbond_triads[i].residue;
    atom.name = hbond_triads[i].B;

    hbond_data_t my_data;
    my_data.atom = atom.name;
    my_data.triad = &(hbond_triads[i]);
    if(atom.is_oxygen()) my_data.angles = &(A_O_angles[atom.res]);
    else if(atom.is_nitrogen()){
      // Arginine is special -- its epsilon nitrogen is similar to the 
      // nitrogen in the tryptophan ring.
      if(atom.res == ARG && atom.name == NE) 
        my_data.angles = &(A_N_angles[TRP]); 
      else my_data.angles = &(A_N_angles[atom.res]);
    }else{
      // should never get here
      std::ostringstream ostr;
      ostr << "I do not know anything about the atom type specified in the " 
           << i+1 << "th row\n of hbond_triads";
      warn(_fname, "build_residue_table()", ostr.str());
      continue;
    }

    hbond_residue_table[atom.res].push_back(my_data);
  }
}

bool 
HbondPoints::load_hbond_angles(GenPointsParameters::hbond_method_t hb_density, 
                               std::string param_path) 
{
  std::string hbonds_path = param_path;
  switch(hb_density){
    case GenPointsParameters::OPTIMUM_HBONDS:
    hbonds_path += "/optimum_hbonds.dat";
    break;
    case GenPointsParameters::MIN_HBONDS:
    hbonds_path += "/minimal_hbonds.dat";
    break;
    case GenPointsParameters::SPARSE_HBONDS:
    hbonds_path += "/sparse_hbonds.dat";
    break;
    case GenPointsParameters::DENSE_HBONDS:
    err_msg(_fname, "load_hbond_angles()", 
            "Dense hbonds data file does not exist ");
    return false;
    break;
  default:
    warn(_fname, "load_hbond_angles()", "Bad hbonds density was given");
    return false;
    break;
  }

  if(!read_angles_file(hbonds_path, add_favored_angles)) return false;
  return read_angles_file(param_path + "/optimum_hbonds.dat", add_ideal_angles);
}

bool
HbondPoints::read_angles_file(std::string fname, add_angles_t add_angles)
{
  std::ifstream in;
  if(!open_ifstream(in, fname)) return false;

  for(std::string line; std::getline(in, line); ){
    if(strip_trailing_comments(&line, "%#") == 0) continue;

    std::vector<std::string> tokens;
    string_tok(line, &tokens, ',');
    if(tokens.size() != 6){
      std::ostringstream ostr;
      ostr << "Expected a csv line with 6 columns, but found " << tokens.size()
           << " columns.  Line is: \n\t" << line << "\n"
           << "Please repair the file " << fname << "\n";
      err_msg(_fname, "load_hbond_angles()", ostr.str());
      return false; 
    }

    if(tokens[1] == "O") add_angles(&A_O_angles, tokens);
    else if(tokens[1] == "N") add_angles(&A_N_angles, tokens);
    else{
      std::ostringstream ostr;
      ostr << "Expected an O or N, but saw " << tokens[1] << " in column 2\n";
      ostr << line << "\n\twill be ignored.\n";
      warn(_fname, "load_hbond_angles()", ostr.str());
      continue;
    }
  }
  in.close();

  return true;
}

bool
HbondPoints::add_favored_angles(std::map<residue_type, angles_vec>* angles,
                                const std::vector<std::string>& tokens)
{
  residue_type residue = PDB_residues::string_to_residue(tokens[0]);
  if(residue == UNKNOWN_RESIDUE){
    std::ostringstream ostr;
    ostr << "Unknown residue identifier (" << tokens[0] << ")\n"
         << "\tCannot compute hydrogen bond positions for this residue type\n";
    warn(_fname, "add_favored_angles()", ostr.str());
    return false;
  }
  std::map<residue_type, angles_vec>::iterator mi = angles->find(residue);
  
  uint idx;
  my_strtoui(tokens[3], &idx);

  // The residue has not been seen before
  if(mi == angles->end()){
    angles_vec tmp_vec(idx);
    (*angles)[residue] = tmp_vec;
    mi = angles->find(residue);

  // The residue is in the map but doesn't have enough angle objects
  }else if(mi->second.size() < idx) mi->second.resize(idx); 

  my_float_t angle;
  my_strtof(tokens[4], &angle);
  angles_vec& A = mi->second;
  A[idx - 1].alphas.push_back(angle / 180 * M_PI); 
  my_strtof(tokens[5], &angle);
  A[idx - 1].betas.push_back(angle / 180 * M_PI); 
  A[idx - 1].ideal_alpha = my_float_max;
  A[idx - 1].ideal_beta = my_float_max;

  return true;
}

bool
HbondPoints::add_ideal_angles(std::map<residue_type, angles_vec>* angles,
                              const std::vector<std::string>& tokens)
{
  // Find the residue in the angles table and make sure there is a entry for
  // the residue and that there are enough sets of favored angles
  residue_type residue = PDB_residues::string_to_residue(tokens[0]);
  if(residue == UNKNOWN_RESIDUE){
    std::ostringstream ostr;
    ostr << "Unknown residue identifier (" << tokens[0] << ")\n"
         << "\tCannot compute hydrogen bond positions for this residue type\n";
    warn(_fname, "add_ideal_angles()", ostr.str());
    return false;
  }
  angles_vec& my_angles = angles->find(residue)->second;
  uint idx;
  my_strtoui(tokens[3], &idx);
  if(my_angles.size() < idx){
    std::ostringstream ostr;
    ostr << "Too many optimum angles specified for " << tokens[0] << " " ;
    ostr << tokens[1] << " " << tokens[2] << " " << tokens[3] << " ";
    ostr << "; additional angles are ignored.\n";
    warn(_fname, "add_ideal_angles()", ostr.str());
    return false;
  }

  // Get an iterator to the element in the list which corresponds to the 
  // number in column 4, and update the optimum angles
  angles_vec::iterator blah = my_angles.begin() + idx - 1;
  my_float_t tmp;
  my_strtof(tokens[4], &tmp);
  blah->ideal_alpha = tmp / 180 * M_PI; 
  my_strtof(tokens[5], &tmp);
  blah->ideal_beta = tmp / 180 * M_PI; 
 
  return true; 
}

bool
HbondPoints::add_fit_point(const atom_vci atom, const my_float_t* pt_pos, 
                           const interactionType pt_act_type, 
                           const std::vector<atom_vci>& rad_atoms, 
                           hbond_fit_pt_vec* fit_pts_p)
{
  if(!bounding_volume->contains(pt_pos)){
    //std::cout << " was outside the bounding volume\n";
    return false;
  }

  // Check for steric clash between a dummy atom centered at the point
  // and any protein atom
  std::vector<atom_vci>::const_iterator rad_atom;
  for(rad_atom = rad_atoms.begin() ; rad_atom != rad_atoms.end(); ++rad_atom){
    if(*rad_atom == atom) continue;

    my_float_t d_squared = dist_squared((*rad_atom)->pos, pt_pos);
    if(d_squared < MIN_VDW_DIST_METAL_2){
        //std::cout << " was too close to the atom " 
                  //<< PDB_residues::residue_to_string((*rad_atom)->res)
                  //<< (*rad_atom)->res_num << " " 
                  //<< PDB_residues::atom_to_string((*rad_atom)->name)
                  //<< " ( dist. " << std::sqrt(d_squared) << ")\n";
      return false;
    }if((*rad_atom)->res != PDB_METAL && d_squared < MIN_VDW_DIST_2){
      //std::cout << " was too close to the atom " 
                //<< PDB_residues::residue_to_string((*rad_atom)->res)
                //<< (*rad_atom)->res_num << " " 
                //<< PDB_residues::atom_to_string((*rad_atom)->name)
                //<< " ( dist. " << std::sqrt(d_squared) << ")\n";
      return false;
    }
  }

  //std::cout << " will be added to sitemap\n";
  // Add point to fit points
  hbond_point_t tmp_pt;
  tmp_pt.act_type = pt_act_type;
  tmp_pt.atom = atom;
  std::copy(pt_pos, pt_pos + 3, tmp_pt.pos);
  fit_pts_p->push_back(tmp_pt);
  return true;
}

void
HbondPoints::add_ideal_point(const atom_vci atom, const atom_vci carbon_nbr,
                             const my_float_t* pt_pos, 
                             const interactionType pt_act_type, 
                             const uint pt_num, const uint num_pts,
                             hbond_ideal_pt_vec* ideal_pts_p,
                             std::vector<int>* num_fit_pts)
{
  num_fit_pts->push_back(num_pts);
  hbond_ideal_point_t tmp_pt;
  tmp_pt.carbon_nbr = carbon_nbr;
  tmp_pt.pt_num = pt_num;
  tmp_pt.act_type = pt_act_type;
  tmp_pt.atom = atom;
  std::copy(pt_pos, pt_pos + 3, tmp_pt.pos);
  ideal_pts_p->push_back(tmp_pt);
}

bool 
HbondPoints::calc_volume(const atom_vci atom, const atom_vci carbon_nbr,
                         const hbond_data_t& hbond_info, 
                         const my_float_t **positions,
                         const std::vector<atom_vci>& rad_atoms)
{
  if(atom->res_num != 43 || atom->name != OG1 || atom->res != THR) 
    return true;


  std::cout << "<tmp residue=" << PDB_residues::residue_to_string(atom->res)
            << " res_num=" << atom->res_num 
            << " atom=\"" << PDB_residues::atom_to_string(atom->name) 
            << "\">\n";

  // Get transformation to move from origin to B + neighbors A & C
  my_float_t trans_A[3], trans_C[3];
  for(uint i = 0; i < 3; i++){
    trans_A[i] = positions[0][i] - positions[1][i];
    trans_C[i] = positions[2][i] - positions[1][i];
  }
  my_float_t R[9];
  get_local_orientation(trans_A, trans_C, R); 

  //my_float_t offset[] = {0.0, 0.0, 0.0};
  angles_vec::const_iterator grp;
  my_float_t* pts = new my_float_t[3*A_hbond_vol.size()];
  for(grp = hbond_info.angles->begin(); grp < hbond_info.angles->end(); ++grp){
    my_float_t* pt = pts;
    std::cout << "B: " << positions[1][0] << " " << positions[1][1] << " "
              << positions[1][2] << "\n";
    for(size_t i = 0; i < A_hbond_vol.size(); ++i, pt += 3)
      std::copy(positions[1], positions[1] + 3, pt);
    std::cout << "B: " << pts[0] << " " << pts[1] << " "
              << pts[2] << "\n";


    // 1) orient the point cloud with respect to the "ideal" positions
    my_float_t cos_alpha = std::cos(grp->ideal_alpha); 
    my_float_t cos_beta = std::cos(grp->ideal_beta); 
    my_float_t sin_alpha = std::sin(grp->ideal_alpha); 
    my_float_t sin_beta = std::sin(grp->ideal_beta); 

    // 
    my_float_t R_Z[] = {     cos_alpha, sin_alpha, 0,
                        -1.0*sin_alpha, cos_alpha, 0,
                                     0,         0, 1};
    my_float_t R_X[] = {  1,             0,        0,    
                          0,      cos_beta, sin_beta,
                          0, -1.0*sin_beta, cos_beta};
    my_float_t tmp[9];
    my_gemm(3, 3, 3, 1.0, R_X, 3, R, 3, tmp, 3, 0.0);
    my_float_t RR[9];
    my_gemm(3, 3, 3, 1.0, R_Z, 3, tmp, 3, RR, 3, 0.0);

    // 2) Do the transformation
    my_gemm(A_hbond_vol.size(), 3, 3, 1.0, A_hbond_vol.begin(), 3, RR, 3, pts, 3, 1.0);

    // 3) Clip the transformed points
    pt = pts;
    for(size_t i = 0; i < A_hbond_vol.size(); ++i, pt += 3){
      if(!bounding_volume->contains(pt)) continue;
    
      // Check for steric clash between a dummy atom centered at the point
      // and any protein atom
      std::vector<atom_vci>::const_iterator rad_atom;
      bool print_pt = true;
      for(rad_atom = rad_atoms.begin() ; rad_atom != rad_atoms.end(); 
          ++rad_atom){
        if(*rad_atom == atom) continue;
    
        my_float_t d_squared = dist_squared((*rad_atom)->pos, pt);
        if(d_squared < MIN_VDW_DIST_METAL_2){
            //std::cout << " was too close to the atom " 
                      //<< PDB_residues::residue_to_string((*rad_atom)->res)
                      //<< (*rad_atom)->res_num << " " 
                      //<< PDB_residues::atom_to_string((*rad_atom)->name)
                      //<< " ( dist. " << std::sqrt(d_squared) << ")\n";
          print_pt = false;
          break;
        }if((*rad_atom)->res != PDB_METAL && d_squared < MIN_VDW_DIST_2){
          //std::cout << " was too close to the atom " 
                    //<< PDB_residues::residue_to_string((*rad_atom)->res)
                    //<< (*rad_atom)->res_num << " " 
                    //<< PDB_residues::atom_to_string((*rad_atom)->name)
                    //<< " ( dist. " << std::sqrt(d_squared) << ")\n";
          print_pt = false;
          break;
        }
      }
      if(print_pt) std::cout << pt[0] << " " << pt[1] << " " << pt[2] << "\n";
    }

    // 4) Print the remaining points 
    //pt = pts;
    //for(int i = 0; i < A_hbond_vol.size(); ++i, pt += 3)
      //std::cout << pt[0] << " " << pt[1] << " " << pt[2] << "\n";
  }
  delete [] pts;

  std::cout << "</tmp>\n";
  return true;
}
