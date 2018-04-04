/******************************************************************************
 * Copyright (c) 2006,2007, Michigan State University (MSU) Board of Trustees.
 *   All rights reserved.
 *
 * This file is part of the ASCbase Software project.
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
 * $Source: /psa/share/repository/pfizer_proj/src/gen_points/HphobPoints.C,v $
 * $Revision: 1.11 $
 * $Author: vanvoor4 $
 * $Date: 2009/01/12 21:09:00 $
 * 
 * $Log: HphobPoints.C,v $
 * Revision 1.11  2009/01/12 21:09:00  vanvoor4
 * Didn't realize that the list erase returns an iter to the previous element.
 * This makes the code more clear and less of a chance to reference an invalid
 * iterator.
 *
 * Revision 1.10  2008/05/15 17:37:07  vanvoor4
 * Added and corrected the support for a 2.5 (A) pseudo surf.
 *
 * Revision 1.9  2007/12/17 21:26:24  vanvoor4
 * Added waters to be considered as hphill pts.
 *
 * Revision 1.8  2007/11/01 20:38:22  vanvoor4
 * Added 2 more levels of spheres.
 * Only hphob atoms are considered neighbors of hphob pts.
 *
 * Revision 1.7  2007/09/26 14:49:53  vanvoor4
 * Added support for hphob 3.0 (A) shell method.
 *
 * Revision 1.6  2007/08/29 20:31:22  vanvoor4
 * updated to respect the handling of hbond points
 *
 * Revision 1.5  2007/08/21 18:23:41  vanvoor4
 * Changes required to support the point_storage scheme
 *
 * Revision 1.4  2007/03/06 19:39:34  vanvoor4
 * Added support for changing the cluster radius
 *
 * Revision 1.3  2007/02/07 15:29:47  vanvoor4
 * Reorganized to allow the testing any number of methods to specify
 * hydrophobicity.
 *
 * Revision 1.2  2006/11/16 20:24:45  vanvoor4
 * Added CVS header
 *
 *
 */

#include <cstring>
#include <iomanip>
#include <sstream>
#include <hydro.H>
#include <complete_link_clustering.H>
#include <mat_ops.H>
#include <Quaternion.H>
#include <HphobPoints.H>

using namespace ASCbase;

const std::string HphobPoints::A_fname = "HphobPoints.C";
const my_float_t HphobPoints::MIN_HYDRO_DIST = 3.0;
const my_float_t HphobPoints::MAX_HYDRO_DIST = 4.5;
const my_float_t HphobPoints::HPHIL_CUTOFF = 100;
hphob_triad_map_t HphobPoints::hphob_nbrs_tbl;

HphobPoints::HphobPoints(std::istream& in, const uint num_lines,
                         PDBBase& atoms, PDBBase& points_in)
{
  hphob_point_t pt;
  std::string line;
  for(uint i = 1; i <= num_lines && std::getline(in, line); ++i){
    atom_vci point = points_in.get_atom(i);
    std::copy(point->pos, point->pos + 3, pt.pos);
    points.push_back(pt); 
  }  
}

HphobPoints::HphobPoints(const GenPointsParameters::hphob_method_t method_in, 
                         const double cluster_diameter_in,
                         const sphere_sample_level_t sample_level_in)
{
  method = method_in;
  cluster_diameter = cluster_diameter_in;
  sphere_sample_level = sample_level_in;
}

bool 
HphobPoints::gen_points(const interact_atoms_vec& hphob_atoms,
                        const hphob_triad_vec& hphob_triads,
                        const std::vector<atom_vci>& rad_atoms,
                        BoundingVolume* site_vol)
{
  if(method == GenPointsParameters::THREE_PROTEIN_ATOMS || 
     method == GenPointsParameters::THREE_MORE_HPHOB)
    return grid_point_methods(rad_atoms, site_vol);
  else if(method == GenPointsParameters::ATOM_CENTERS)
    return get_hphobic_atoms(rad_atoms);
  else if(method == GenPointsParameters::PSEUDO_SURFACE)
    return pseudo_surface_points(hphob_atoms, hphob_triads, rad_atoms, 
                                 site_vol);
  return false;
}

bool 
HphobPoints::grid_point_methods(const std::vector<atom_vci>& rad_atoms,
                                BoundingVolume* site_vol)
{
  // Cull the points based on hueristical methods
  hphob_point_vec kept_pts;
  const my_float_t *cur_pt = site_vol->grid_points_begin();
  for( ; cur_pt < site_vol->grid_points_end(); cur_pt += 3){
    uint hphob_sum = 0;
    uint hphil_sum = 0;
    std::vector<atom_vci>::const_iterator rad_atom;
    for(rad_atom = rad_atoms.begin(); rad_atom != rad_atoms.end(); rad_atom++){
      my_float_t my_dist = dist((*rad_atom)->pos, cur_pt);
      if(my_dist < MIN_HYDRO_DIST){
        hphob_sum = 0;
        hphil_sum = 0;
        break;
      }else if(my_dist < MAX_HYDRO_DIST){
        // May need to fix this if we get problems with array overruns -- most
        // likely everything will get mapped to UNKNOWN_RESIDUE and atom 
        // giving zero as value and mapping it to HPHOB
        if((*rad_atom)->res != HOH && 
           hydro[(*rad_atom)->res][(*rad_atom)->name] < HPHIL_CUTOFF) 
	  hphob_sum++;
        else hphil_sum++;
      }
    }    

    bool add = false;
    switch(method){
    case GenPointsParameters::THREE_PROTEIN_ATOMS:
      if((hphob_sum + hphil_sum) > 2) add = true;
      break;
    case GenPointsParameters::THREE_MORE_HPHOB:
      if((hphob_sum - hphil_sum) > 2) add = true;
      break;
    //case THREE_MORE_HPHOB_AND_AT_MOST_ONE_HPIL:
      //if((hphob_sum - hphil_sum) > 2 && hphil_sum < 2) add = true;
      //break;
    default:
      std::string msg = "The set hphob method is of incompatible type with";
      warn(A_fname, "grid_point_methods", msg + " this routine", std::cerr);
      return false;
    }

    if(add) add_point(cur_pt, &kept_pts);
  }

  // Cluster the points
  complete_link_clustering(kept_pts, &points, cluster_diameter);
 
  // Assign all protein atoms within MAX_HYDRO_DIST as the "atoms making 
  // the interaction"
  hphob_point_li hphob_pt = points.begin();
  for(; hphob_pt != points.end(); ++hphob_pt){
    std::vector<atom_vci>::const_iterator rad_atom;
    for(rad_atom = rad_atoms.begin(); rad_atom != rad_atoms.end(); rad_atom++)
      if((*rad_atom)->is_hydrophobic() && 
         dist(hphob_pt->pos, (*rad_atom)->pos) <= MAX_HYDRO_DIST)
        hphob_pt->atoms.push_back(*rad_atom);
 
    // If no protein atom within MAX_HYDRO_DIST is hydrophobic, remove this
    // hydrophobic point
    if(hphob_pt->atoms.empty()){
      hphob_pt = points.erase(hphob_pt); 
      --hphob_pt;
    } 
  }

  return true;
}

void 
HphobPoints::cull_too_close_to_polar(hbond_fit_pt_vci hbonds_begin,
                                     hbond_fit_pt_vci hbonds_end)
{
  my_float_t threshold = cluster_diameter/2.0;
  hbond_fit_pt_vci hbond_pt;
  for(hbond_pt = hbonds_begin; hbond_pt < hbonds_end; ++hbond_pt){
    hphob_point_li hphob_pt;
    for(hphob_pt = points.begin(); hphob_pt != points.end(); ++hphob_pt)
      if(dist(hphob_pt->pos, hbond_pt->pos) < threshold){
        hphob_pt = points.erase(hphob_pt);
        // The assumption here is that at most one hydrophobic point will be
        // within the ball centered at hphob_points[i].pos with radius
        // threshold.  Whether this holds depends on the threshold AND
        // the hydrophobic cluster radius.
        break;
//        --hphob_pt;
      }
  }
}

bool 
HphobPoints::get_hphobic_atoms(const std::vector<atom_vci>& rad_atoms)
{
  hphob_point_vec tmp_pts;
  std::vector<atom_vci>::const_iterator rad_atom;
  for(rad_atom = rad_atoms.begin(); rad_atom != rad_atoms.end(); ++rad_atom){
    atom_vci atom = *rad_atom;
    bool add = false;
    switch(atom->res){
    case ALA:
      if(atom->name == CB) add = true;
      break;
    case ILE:
      if(atom->name == CB || atom->name == CG1 || atom->name == CG2 
         || atom->name == CD) add = true;
      break;
    case LEU:
      if(atom->name == CB || atom->name == CG || atom->name == CD1 
         || atom->name == CD2) add = true;
      break;
    case MET:
      if(atom->name == CB || atom->name == CG || atom->name == SD 
         || atom->name == CE) add = true;
      break;
    case PHE:
      if(atom->name == CB || atom->name == CG || atom->name == CG
         || atom->name == CD || atom->name == CD || atom->name == CZ) 
        add = true;
      break;
    case PRO:
      if(atom->name == CB || atom->name == CG || atom->name == CD) add = true;
      break;
    case THR:
      if(atom->name == CG2) add = true;
      break;
    case TRP:
      if(atom->name == CB || atom->name == CG || atom->name == CD2
         || atom->name == CE3 || atom->name == CZ2 || atom->name == CZ3
         || atom->name == CH2) add = true;
    case TYR:
      if(atom->name == CB || atom->name == CD1 || atom->name == CD2 
         || atom->name == CE1 || atom->name == CE2) add = true;
      break;
    case VAL:
      if(atom->name == CB || atom->name == CG1 || atom->name == CG2)
        add = true;
      break;
    default:
      break;
    } 
    if(add) add_point(atom->pos, &points);
  }

  return true;
}

void 
HphobPoints::add_point(const my_float_t* pos, hphob_point_vec* pts_vec)
{
  hphob_point_t my_point;
  std::memcpy(my_point.pos, pos, 3*my_float_size);
  pts_vec->push_back(my_point);
}

void 
HphobPoints::add_point(const my_float_t* pos, hphob_point_list* pts_list)
{
  hphob_point_t my_point;
  std::memcpy(my_point.pos, pos, 3*my_float_size);
  pts_list->push_back(my_point);
}

bool 
HphobPoints::pseudo_surface_points(const interact_atoms_vec& hphob_atoms,
                                   const hphob_triad_vec& hphob_triads,
                                   const std::vector<atom_vci>& rad_atoms,
                                   BoundingVolume* site_vol)
{
  DiscreteSphere orig_sphere(sphere_sample_level);
  orig_sphere.scale(2.5);

  my_float_t MIN_DIST_SQUARED = BoundingVolume::MIN_VDW_DIST;
  MIN_DIST_SQUARED *= BoundingVolume::MIN_VDW_DIST;
  my_float_t MAX_DIST_SQUARED = MAX_HYDRO_DIST * MAX_HYDRO_DIST;

  interact_atoms_vci atom_i;
  hphob_triad_vci triad_i = hphob_triads.begin();
  for(atom_i = hphob_atoms.begin(); atom_i < hphob_atoms.end(); ++atom_i){
    const atom_vci atom = atom_i->atom;

    // 1) Search through the atom's residue to find the corresponding atoms
    //    in the triad
    const residue_vci residue = atom_i->residue;
    const hphob_triad_t* triad = *triad_i;
    atom_vci nbr_A = atom_t::NULL_ATOM_VCI;
    atom_vci nbr_C = atom_t::NULL_ATOM_VCI;
    for(atom_vci a = residue->atoms_begin; a < residue->atoms_end; ++a){
      if(a->name == triad->nbr_A) nbr_A = a; 
      if(a->name == triad->nbr_C) nbr_C = a; 
    }
    ++triad_i; 

    if(nbr_A == atom_t::NULL_ATOM_VCI || nbr_C == atom_t::NULL_ATOM_VCI){
      std::ostringstream ostr; 
      ostr << "Unable to compute hydrophobic points for: "
           << PDB_residues::residue_to_string(atom->res) << " "
           << atom->res_num 
           << PDB_residues::atom_to_string(atom->name);
      warn(A_fname, "psuedo_surface_points", ostr.str());
      continue;
    }

    // 2) Get the transformation to move the "ball" to the center of B
    my_float_t U[6];
    my_float_t* V = U + 3;
    std::copy(nbr_A->pos, nbr_A->pos + 3, U);
    std::copy(nbr_C->pos, nbr_C->pos + 3, V);
    my_axpy(3, -1.0, atom->pos, 1, U, 1);
    my_axpy(3, -1.0, atom->pos, 1, V, 1);
    my_float_t R[9];
    get_local_orientation(U, V, R);
    DiscreteSphere aligned_sphere(orig_sphere, R, atom->pos);

    const my_float_t* pt;
    for(pt = aligned_sphere.pts_begin(); pt < aligned_sphere.pts_end(); pt+=3){
      if(!site_vol->contains(pt))  continue;

      hphob_point_t my_point;
      std::copy(pt, pt + 3, my_point.pos);
      bool discard_point = false;

      std::vector<atom_vci>::const_iterator rad_iter = rad_atoms.begin();
      for( ; rad_iter != rad_atoms.end(); rad_iter++){
        atom_vci rad_atom = *rad_iter;

        // Point cannot be too close to any atom
        my_float_t d_squared = MIN_DIST_SQUARED;
        if(rad_atom != atom){
          d_squared = dist_squared(rad_atom->pos, pt);
          if(d_squared < MIN_DIST_SQUARED){
            discard_point = true;
            break;
          }
        }

        // Assign all hphob protein atoms within MAX_HYDRO_DIST as the 
        // "atoms making the interaction"
        const hphob_triad_t* not_used;
        if(atom_is_hphob(rad_atom->res, rad_atom->name, &not_used) &&
           d_squared <= MAX_DIST_SQUARED)
          my_point.atoms.push_back(rad_atom);
      }
      if(!discard_point) points.push_back(my_point);
    }
  }
  return true;
}

void 
HphobPoints::init_hphob_nbrs_tbl()
{
  if(hphob_nbrs_tbl.size()) return;

  const hphob_triad_t* triad = hphob_triads;
  for(uint i = 0; i < num_hphob_triads; ++i, ++triad)
    hphob_nbrs_tbl[triad->residue][triad->B] = triad;
}

bool
HphobPoints::atom_is_hphob(const residue_type res, const atom_type atom,
                           const hphob_triad_t** triad_ptr)
{   
  if(hphob_nbrs_tbl.size() == 0) init_hphob_nbrs_tbl();

  hphob_triad_map_t::const_iterator res_iter = hphob_nbrs_tbl.find(res);
  if(res_iter == hphob_nbrs_tbl.end()) return false;
  std::map<atom_type, const hphob_triad_t*>::const_iterator atom_iter
    = res_iter->second.find(atom);
  if(atom_iter == res_iter->second.end()) return false;
  *triad_ptr = atom_iter->second;
  return true;
}
