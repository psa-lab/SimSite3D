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
#include <SitemapAtomsFile.H>
#include <algorithm>
#include <mol2File.H>

using namespace SimSite3D;
const std::string SitemapAtomsFile::A_fname = "SitemapAtomsFile.C";
const my_float_t SitemapAtomsFile::A_max_interact_dist = 12.0;

SitemapAtomsFile::SitemapAtomsFile(const std::string filename, 
                                   const verbose_level_t verbosity)
 : PDBBase(filename, verbosity)
{
  if(fail()) return;
  std::ifstream pdb_file;
  if(open_ifstream(pdb_file, filename)){
    if(!read_data(pdb_file)) red_light();
    pdb_file.close();
  }else red_light();
}

SitemapAtomsFile::SitemapAtomsFile(const atom_map_t& a_tbl)
{
  for(atom_map_t::const_iterator m = a_tbl.begin(); m != a_tbl.end(); ++m)
    append_atom(*(m->second));
}

SitemapAtomsFile::SitemapAtomsFile(PDBStructure *prot, const atom_map_t& a_tbl,
                                   const std::string lig_fname, 
                                   const int min_chain_size,
                                   const my_float_t max_interact_dist)
{
  CoordFile* lig = 0;

  // Mol2 and PDB files are supported
  std::string::size_type pos = lig_fname.rfind(".");
  if(pos == std::string::npos){
    std::string msg = "The ligand file name must include either the .pdb ";
    msg += "or .mol2 file extension as part of its name\n";
    err_msg(A_fname, "SitemapAtomsFile()", msg);
    red_light();
    return;
  }
  std::string ext = lig_fname.substr(pos + 1);

  if(ext == "pdb") lig = new PDBBase(lig_fname);
  else if(ext == "mol2") lig = new mol2File(lig_fname);
  else{
    std::string msg = "The ligand file name must include either the .pdb ";
    msg += "or .mol2 file extension as part of its name\n";
    err_msg(A_fname, "SitemapAtomsFile()", msg);
    red_light();
    return;
  }
  if(lig == 0 || lig->fail()){
    err_msg(A_fname, "SitemapAtomsFile()", "Unable to read the ligand file");
    red_light();
    if(lig) delete lig;
    lig = 0;
    return;
  }

  chain_const_iter chain;
  for(chain = prot->chains_begin(); chain < prot->chains_end(); chain++){
    if(chain->residues_end - chain->residues_begin < min_chain_size)
      continue;

    residue_vci res = chain->residues_begin;
    for( ; res < chain->residues_end; ++res){
      bool inc_residue = false;
      for(atom_vci atom = res->atoms_begin; 
          atom < res->atoms_end && !inc_residue; ++atom){
        // Hydrogen atoms should be ignored as a_tbl should not include
        // any hydrogen atoms
        if(a_tbl.find(atom->atom_num) != a_tbl.end()){
          inc_residue = true;
          break;
        }else{
          my_float_t d;
          lig->closest_atom(atom->pos, &d);
          if(d <= max_interact_dist){
            inc_residue = true;
            break;
          }
        }        
      } 
  
      if(inc_residue)
        for(atom_vci atom = res->atoms_begin; atom < res->atoms_end; ++atom) 
          // Skip hydrogen atoms
          if(atom->name != H && atom->name != D) append_atom(*atom);
    }
  }

  // Assumption is a_tbl will only include metals if they are desired
  std::vector<atom_vci>::const_iterator metal_iter = prot->metals_begin();
  for( ; metal_iter != prot->metals_end(); ++metal_iter){
    if(a_tbl.find((*metal_iter)->atom_num) != a_tbl.end())
      append_atom(**metal_iter);
    else{
      my_float_t d;
      lig->closest_atom((*metal_iter)->pos, &d);
      if(d <= max_interact_dist) append_atom(**metal_iter);
    }
  }

  // Assumption is a_tbl will only include waters if they are asked for
  // We also don't want to add waters that are nearby, but not included in
  // a_tbl since we don't want to clutter the rad file with waters we are
  // not interested in from an interaction viewpoint.
  std::vector<atom_vci>::const_iterator w_iter = prot->waters_beg();
  for( ; w_iter != prot->waters_end(); ++w_iter){
    if(a_tbl.find((*w_iter)->atom_num) != a_tbl.end())
      append_atom(**w_iter);
  }

  green_light();

  delete lig;
}

SitemapAtomsFile::SitemapAtomsFile(PDBStructure *prot, 
                                   BoundingVolume* site_vol,
                                   const int min_chain_size,
                                   const atom_map_t& a_tbl)
{
  chain_const_iter chain;
  for(chain = prot->chains_begin(); chain < prot->chains_end(); chain++){
    if(chain->residues_end - chain->residues_begin < min_chain_size)
      continue;

    residue_vci res = chain->residues_begin;
    for( ; res < chain->residues_end; ++res){
      for(atom_vci atom = res->atoms_begin; atom < res->atoms_end; ++atom){
        // Ignore hydrogen atoms
        if(atom->name == H || atom->name == D) continue;
        if(site_vol->BIND_vol_contains(atom->pos)){
          for(atom_vci a = res->atoms_begin; a < res->atoms_end; ++a){
            // Skip hydrogen atoms
            if(a->name != H && a->name != D) append_atom(*a);
          }
          break;
        }
      }
    }
  }

  // Assumption is a_tbl will only include metals if they are desired
  std::vector<atom_vci>::const_iterator m_iter = prot->metals_begin(); 
  for( ; m_iter != prot->metals_end(); ++m_iter){
    atom_vci metal = *m_iter;
    if(a_tbl.find(metal->atom_num) != a_tbl.end()) append_atom(*metal);
    else{
      if(site_vol->BIND_vol_contains(metal->pos)) append_atom(*metal);
    }
  }

  // Assumption is a_tbl will only include waters if they are asked for
  // We also don't want to add waters that are nearby, but not included in
  // a_tbl since we don't want to clutter the rad file with waters we are
  // not interested in from an interaction viewpoint.
  std::vector<atom_vci>::const_iterator w_iter = prot->waters_beg();
  for( ; w_iter != prot->waters_end(); ++w_iter){
    if(a_tbl.find((*w_iter)->atom_num) != a_tbl.end())
      append_atom(**w_iter);
  }

  green_light();
}

SitemapAtomsFile::~SitemapAtomsFile()
{
}

void
SitemapAtomsFile::write_xyzr(std::ostream &out, const atom_vci close_atom, 
                             uint* close_atom_idx)
{
  uint atom_idx = 0;
  for(atom_vci a = atoms_begin(); a != atoms_end(); ++a){
    // Skip hydrogen atoms
    if(a->name != H && a->name != D){
      write_xyzr_line(out, a->pos, a->vdw_radius);
      if(close_atom == a) *close_atom_idx = atom_idx;
      ++atom_idx;
    }
  }
}
