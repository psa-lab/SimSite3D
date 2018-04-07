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
#include <iostream>
#include <iomanip>
#include <sstream>
#include <basics.H>
#include <PDBStructure.H>

using namespace SimSite3D;


const std::string PDBStructure::_fname = "PDBStructure.C";
const residue_t PDBStructure::NULL_RESIDUE;
const residue_vec PDBStructure::NULL_RESIDUE_VECTOR;
const residue_vci PDBStructure::NULL_RESIDUE_VCI = 
  PDBStructure::NULL_RESIDUE_VECTOR.begin();

PDBStructure::PDBStructure(const std::string filename,
                           const bool ignore_altLocs,
                           const verbose_level_t verbosity)
 : PDBBase(filename, verbosity)
{
  if(fail()) std::cout << "PBDBase::PDBBase() failed!!!!" << std::endl;
  if(fail() || !read_file(filename, ignore_altLocs)) red_light();
  else{
    green_light();
    form_chains_and_residues();
    report_stats(std::cout, verbosity);
  }
}

PDBStructure::~PDBStructure()
{
}

bool
PDBStructure::read_file(std::string filename, const bool ignore_altLocs)
{
  std::ifstream pdb_file;
  bool rv = open_ifstream(pdb_file, filename);
  if(rv){
    rv = read_data(pdb_file, ignore_altLocs);
    pdb_file.close();
  }
  return rv;
}

void
PDBStructure::form_chains_and_residues()
{
  atom_vci atom = atoms_begin();
  residue_vi current_residue = residues.end();
  residue_t tmp(atom);
  bool residue_ended = true;
  for(atom = atoms_begin(); atom < atoms_end(); atom++){
    if(atom->is_water()){
      A_waters.push_back(atom);
      if(!residue_ended){
        residue_ended = true;
        current_residue->atoms_end = atom;
      }
    }else if(atom->is_metal()){
      metals.push_back(atom);
      if(!residue_ended){
        residue_ended = true;
        current_residue->atoms_end = atom;
      }
    }else if(current_residue == residues.end()
             || atom->res_num != current_residue->number 
             || atom->iCode != current_residue->icode 
             || atom->chainID != current_residue->chainID){
      if(!residue_ended) current_residue->atoms_end = atom;
      if(add_residue(atom, &current_residue)) residue_ended = false;
      // If there is an error, we want to end the residue
      else residue_ended = true;
    }
  }
  // Need to check -- the input file might not contain HETATMs 
  if(residue_ended == false) current_residue->atoms_end = atom;

  // Form the chains -- assumes no break in the chains (i.e. no insertion
  // of ATOM or HETATM records with a different chain ID in between atoms
  // belonging to the same chain
  chains.resize(1);
  chain_iter current_chain = chains.begin();
  current_chain->chainID = residues.begin()->chainID;
  current_chain->residues_begin = residues.begin();
  current_chain->residues_end =  residues.end();
  residue_vci residue;
  for(residue = residues.begin(); residue < residues.end(); residue++){
    if(residue->chainID != current_chain->chainID){
      current_chain->residues_end = residue;
      chain_t chain;
      chains.push_back(chain);
      current_chain = chains.end() - 1;
      current_chain->chainID = residue->chainID;
      current_chain->residues_begin = residue;
      current_chain->residues_end = residues.end();
    }
  }
}

void
PDBStructure::report_stats(std::ostream& out, const verbose_level_t verbosity)
{
  if(verbosity){
    out << "Stats of PDB file " << CoordFile::name() << "\n";
    out << "Number of chains:   " << chains.size() << "\n";
    for(chain_iter chain = chains.begin(); chain < chains.end(); chain++)
      out << "Number of residues in chain " << chain->chainID << " :" 
          << chain->residues_end - chain->residues_begin << "\n";
  }
}

// Might make sense to speed this up or provide each atom with a 
// pointer to its residue
residue_vci
PDBStructure::get_residue(atom_vci atom) const
{
  if(atom < atoms_begin() or atoms_end() <= atom){
    std::cerr << "protein atom is not from this PDB structure\n";
    return NULL_RESIDUE_VCI;
  }

  chain_const_iter chain;
  for(chain = chains.begin(); chain < chains.end(); ++chain)
    if(chain->chainID == atom->chainID) break;
  if(chain->chainID != atom->chainID) return NULL_RESIDUE_VCI;

  residue_vci residue;
  for(residue = chain->residues_begin; residue < chain->residues_end; ++residue)
    if(residue->number == atom->res_num && residue->icode == atom->iCode)
      return residue;

  // shouldn't get here
  std::cout << "in get_residue " << PDB_residues::residue_to_string(atom->res)
            << atom->res_num << " " 
            << PDB_residues::atom_to_string(atom->name) << std::endl;
  std::cerr << "protein atom is not from this PDB structure\n";
  return NULL_RESIDUE_VCI; 
}

bool
PDBStructure::add_residue(atom_vci initial_atom, residue_vi* current_residue)
{
  // If the residue is specified using ATOM records consider it as part of
  // a polypeptide chain
  residue_t tmp(initial_atom);
  if(!tmp.is_HET){
    residues.push_back(tmp);
    *current_residue = residues.end() - 1;
    return true;
  }


  // If the residue is specified using HETATM records, check if it has a
  // corresponding MODRES entry.  If it is a MODRES, consider it as part of
  // a polypeptide chain
  bool added = false;
  modres_record_vci mr = modres_records_begin();
  for( ; mr != modres_records_end(); ++mr)
    if(mr->contains(initial_atom)){
      tmp.is_MODRES = true;
      residues.push_back(tmp);
      *current_residue = residues.end() - 1;
      added = true;
      break;
    }

  // If the residue is neither one of the 20 amino acids nor has a MODRES
  // entry, consider this residue as a separate entity (ligand, cofactor, etc).
  if(!added){
    het_record_vci hr = het_records_begin();
    for( ; hr != het_records_end(); ++hr)
      if(hr->contains(initial_atom)){
        A_hetgroups.push_back(tmp);
        *current_residue = A_hetgroups.end() - 1;
        added = true;
        break;
      }
  }
  if(!added){
    std::ostringstream ostr;
    ostr << "Unknown residue/hetgroup for atom to initialize residue\n"
         << "\tResidue: " << initial_atom->res_str
         << "\tAtom: " << initial_atom->name_str << "\n"
         << "Determined the HET entry from the first atom for this HET group\n";
    warn(_fname, "add_residue", ostr.str());
    add_het_record(initial_atom);
    A_hetgroups.push_back(tmp);
    *current_residue = A_hetgroups.end() - 1;
  }
  return true;
}

