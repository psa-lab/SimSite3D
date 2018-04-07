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
#include <sstream>
#include <iomanip>
#include <basics.H>
#include <PDBBase.H>
#include <PDBHetatmTypes.H>
#include <PDB_metals.H>

using namespace SimSite3D;

const std::string PDBBase::A_fname = "PDBBase.C";

PDBBase::PDBBase()
{
}

PDBBase::PDBBase(const std::string filename, const verbose_level_t verbosity)
 : CoordFile(filename, verbosity)
{
}

PDBBase::~PDBBase()
{
}

bool 
PDBBase::read_data(std::ifstream &fin, const bool ignore_altLocs)
{
  // If we want to ignore alternate locations, we need an easy to 
  // implement method that checks if alternate locations exist.
  // Using a map is the easiest method.  One issue is alternate locations
  // for a given atom are not required (by the PDB specifications) to be
  // on consecutive lines.  In particular if there is a large loop, it is
  // perfectly acceptible to have the first specify the N resdiues for location
  // 'A', then the N residues for location 'B', etc.
  std::map<std::string, bool> atom_ids;
  bool altLoc_warning_was_printed = false;

  bool has_multiple_models = false;
  unsigned int model_num = 0;
  for(std::string line; std::getline(fin, line); ){
    pdb_record_type r_type = get_record_type(line);
    if(r_type == PDB_MODEL_RECORD){
      has_multiple_models = true;
      my_strtoui(line.substr(10,4).c_str(), &(model_num));
      if(model_num > 1) break;
      continue;
    }else if(r_type == PDB_MODRES_RECORD){
      read_modres_record(line); 
      continue;
    }else if(r_type == PDB_HET_RECORD){
      read_het_record(line);
      continue;
    }else if(r_type == PDB_NOT_ATOM_RECORD) continue; 

    atom_t atom;
    atom.is_hetero = false;
    if(r_type == PDB_HEAVY_ATOM){
      if(!get_res_info(line.substr(12,4), line.substr(17,3), &atom)){
        std::ostringstream msg;
        msg << "Unknown atom or residue found in the file:\n\t"
            << name() << "\nLine is: \n" << line;
        warn(A_fname, "read_data()", msg.str());
      }
    }else{ 
      if(r_type == PDB_HYDROGEN_ATOM){
        atom.res = PDB_residues::string_to_residue(line.substr(17,3));
        if(atom.res == UNKNOWN_RESIDUE){
          std::ostringstream msg;
          msg << "Unknown atom or residue found in the file:\n\t"
              << name() << "\nLine is: \n" << line;
          warn(A_fname, "read_data()", msg.str());
        }
      }else{
        atom.res_str = line.substr(17,3);
        atom.is_hetero = true;
      }

      atom.name_str = line.substr(12,4);
      atom.name = PDBHetatmTypes::string_to_atom(atom.name_str);
      const atom_info_t* het_info = PDBHetatmTypes::lookup_vdw_rad(atom.name);
      if(het_info){
        atom.vdw_radius = het_info->vdw_radius; 
        atom.orbit = het_info->orbit;
      }else{
        std::ostringstream msg;
        msg << "Unknown atom or residue found in the file:\n\t"
            << name() << "\nLine is: \n" << line;
        warn(A_fname, "read_data()", msg.str());
      }

      // If the atom is a heavy HETATM, check if it is a metal atom.
      // Silently ignore for now since we need to handle SO4, PO4, etc 
      // We want to only designate a metal atom as a "METAL" in terms
      // of site map points and scoring if that atom is not covalently 
      // attached to a ligand (i.e. ligands, such as MAP, contain metal
      // atoms that should be considered as part of the ligand and not
      // as part of the binding site).  I need to think how this might
      // affect protein-ligand scoring.
      if(r_type == PDB_HEAVY_HET && 
         PDB_metals::lookup(atom.name_str, atom.res_str,
                            SILENTLY_IGNORE_UKNOWN_METALS))
        atom.res = PDB_METAL;

      // Special case for waters -- vdw radius is taken from the Li & 
      // Nussinov paper referenced in PDB_residues
      if(atom.res_str == "HOH"){
        if(atom.name == O) atom.vdw_radius = 1.36;
	else atom.vdw_radius = 1.0;
	atom.res = HOH;
        atom.act_type = DONEPTOR;
        atom.hydro = 345;
      }
    }

    atom.altLoc = line[16];
    atom.chainID = line[21];
    atom.iCode = line[26];
    my_strtoui(line.substr(6,5).c_str(), &(atom.atom_num));
    my_strtoui(line.substr(22,4).c_str(), &(atom.res_num));
    // Coords for ATOM lines are 30-37,38-45,46-53 (when counting from 0)
    for(int i = 0; i < 3; i++)
      my_strtof(line.substr(30 + 8*i, 38 + 8*i), &(atom.pos[i]));
    my_strtof(line.substr(54,6), &(atom.occupancy));
    my_strtof(line.substr(60,6), &(atom.tempFactor));

    if(ignore_altLocs){
      std::stringstream atom_key;
      atom_key << line.substr(12,4) << atom.chainID 
               << std::setw(4) << atom.res_num << atom.iCode;
      if(atom_ids.find(atom_key.str()) != atom_ids.end()){
        if(!altLoc_warning_was_printed){
          std::string msg;
          msg = "SimSite3D does not support alternate locations for atoms.\n";
          msg += "SimSite3D uses the first alternate location for each atom.\n";
          msg += "If you wish to use additional alternate locations, ";
          msg += "please create a protein pdb\n";
          msg += "file for each desired combination of alternate locations\n";
          warn(A_fname, "read_data()", msg);
          altLoc_warning_was_printed = true;
        }
        continue;
      }else 
        atom_ids.insert(atom_ids.begin(), 
                        std::pair<std::string, bool>(atom_key.str(), true));
    }    

    append_atom(atom);
  }

  if(has_multiple_models){
    std::cout << name() << " has multiple models\n"
              << "\tThe first model will be used\n"
              << "If you wish to use multiple models, please split the models "
              << "into separate files\n";
  }

  return true;
}

pdb_record_type
PDBBase::get_record_type(const std::string line)
{
  pdb_record_type rv = PDB_NOT_ATOM_RECORD;
  std::string label = line.substr(0,6); 
  if(label == "HET   ") rv = PDB_HET_RECORD;
  else if(label == "ATOM  ") rv = PDB_HEAVY_ATOM;
  else if(label == "HETATM") rv = PDB_HEAVY_HET;
  else if(label == "MODEL ") return PDB_MODEL_RECORD;
  else return PDB_NOT_ATOM_RECORD; 

  // Warn if we have a "short" atom line
  if(line.length() < 66){
    std::string msg = "Short atom line found: \"";
    msg += std::string("\"\n") + line + "\nSkipping ...";
    warn(A_fname, "get_record_type()", msg);
    rv = PDB_NOT_ATOM_RECORD;
  // Silently ignore hydrogens -- for now
  }else if(is_hydrogen(line.substr(12,4))){
    if(rv == PDB_HEAVY_ATOM) rv = PDB_HYDROGEN_ATOM;
    else rv = PDB_HYDROGEN_HET;
  }
  return rv;
}

bool
PDBBase::write(const std::string filename)
{
  std::ofstream out;
  if(!open_ofstream(out, filename)) return false;
  out << A_header;

  out.setf(std::ios_base::fixed, std::ios_base::floatfield);
  for(atom_vci a = atoms_begin(); a < atoms_end(); ++a){
    const atom_t& atom = *a;
    out.precision(3);
    if(atom.is_hetero){
      out << "HETATM" << std::setw(5) << atom.atom_num << " " << atom.name_str
          << atom.altLoc << atom.res_str << " ";
    }else{
      out << "ATOM  " << std::setw(5) << atom.atom_num << " ";
      if(atom.name_str.length()) out << atom.name_str;
      else out << PDB_residues::atom_to_string(atom.name);
      out << atom.altLoc << PDB_residues::residue_to_string(atom.res) << " ";
    }
    out << atom.chainID << std::setw(4) << atom.res_num << atom.iCode << "   ";
    for(uint j = 0; j < 3; ++j) out << std::setw(8) << atom.pos[j];
    out.precision(2);
    out << std::setw(6) << atom.occupancy
        << std::setw(6) << atom.tempFactor << "\n";
  }
  return true;
}

const bool
PDBBase::const_write(const std::string filename) const
{
  std::ofstream out;
  if(!open_ofstream(out, filename)) return false;
  out << A_header;

  out.setf(std::ios_base::fixed, std::ios_base::floatfield);
  for(atom_vci a = atoms_begin(); a < atoms_end(); ++a){
    const atom_t& atom = *a;
    out.precision(3);
    if(atom.is_hetero){
      out << "HETATM" << std::setw(5) << atom.atom_num << " " << atom.name_str
          << atom.altLoc << atom.res_str << " ";
    }else{
      out << "ATOM  " << std::setw(5) << atom.atom_num << " ";
      if(atom.name_str.length()) out << atom.name_str;
      else out << PDB_residues::atom_to_string(atom.name);
      out << atom.altLoc << PDB_residues::residue_to_string(atom.res) << " ";
    }
    out << atom.chainID << std::setw(4) << atom.res_num << atom.iCode << "   ";
    for(uint j = 0; j < 3; ++j) out << std::setw(8) << atom.pos[j];
    out.precision(2);
    out << std::setw(6) << atom.occupancy
        << std::setw(6) << atom.tempFactor << "\n";
  }
  return true;
}

bool 
PDBBase::write_xyzr(std::ostream &out, bool include_metals,
                    const std::vector<std::string> &waters)
{
  if(out.fail()){
    err_msg(A_fname, "write_xyzr()", "ostream.fail() is true");
    return false;
  }

  out.setf(std::ios_base::fixed, std::ios_base::floatfield);
  for(atom_vci a = atoms_begin(); a < atoms_end(); ++a){
    // Skip hydrogen atoms
    if(a->name == H || a->name == D) continue;
    my_float_t radius = -1.0;

    // Keep only those water molecules in the waters vector
    if(a->res == HOH){
      if(a->name != O) continue;   
      if(!waters.size()) continue;
     
      std::ostringstream Wstr;
      if(a->chainID != ' ') Wstr << a->chainID;
      Wstr << a->res_num;
      if(a->iCode != ' ') Wstr << a->iCode;
      std::vector<std::string>::const_iterator w;
      bool not_found = true;
      for(w = waters.begin(); w < waters.end() && not_found; ++w)
        if(*w == Wstr.str()) not_found = false;
      if(not_found) continue; 
      // Reduced the radius to that of Li & Nussinov 
      radius = 1.36;

    // Keep the SimSite3D aware metals if include_metals is true
    }else if(a->res == PDB_METAL){
      if(!include_metals) continue;

      const pdb_metal_info_t* m_info = PDB_metals::lookup(a->name);
      if(m_info){
        const atom_info_t* a_info = AtomTypes::lookup_vdw_rad(a->name);
        if(a_info) radius = a_info->vdw_radius;
        else{
          err_msg(A_fname, "write_xyzr()", "Metal atom radius lookup failed");
          return false;
        }
      }else{
        err_msg(A_fname, "write_xyzr()", "Metal lookup failed");
        return false;
      }
    // Use the radii for all atoms in PDB_residues and for C, N, O, S, and P
    // for atoms in non standard residues.  Note: that we handled metals
    // and waters before getting to this point.  Also, the code is getting to
    // the point that prepared structures should have any irrelavent peptides,
    // olgionucleotides, ligands, and junk removed.
    }else{ 
      const pdb_atom_info_t* atom_info;
      atom_info = PDB_residues::get_atom_info(a->name, a->res);
      if(atom_info) radius = atom_info->vdw_radius;
      else{
        const atom_info_t* tmp_info = AtomTypes::lookup_vdw_rad(a->name);
        // If the atom is one of C, N, O, S, or P use the vdw radii from 
        // AtomTypes, otherwise skip the atom.
        if(tmp_info && (a->name == C || a->name == N || a->name == O || 
                        a->name == S || a->name == P)){
          radius = tmp_info->vdw_radius;
        }else{
          std::string msg = "Unknown atom type -- unable to find a radius for ";
          msg += "the atom type: ";
          warn(A_fname, "write_xyzr()", msg + a->name_str);
          continue;
        }
      }
    }

    out.precision(3);
    out << std::setw(8) << a->pos[0] << " " << std::setw(8) << a->pos[1] << " "
        << std::setw(8) << a->pos[2] << " ";
    out.precision(2);
    out << radius << "\n";
  }
  return true;
}

bool
PDBBase::is_hydrogen(const std::string name)
{
  char c13 = name[0];
  char c14 = toupper(name[1]);
  if((c13 == ' ' || isdigit(c13)) && (c14 == 'H' || c14 == 'D')) return true;
  // The remediated PDB files have a different scheme
  if(c13 == 'H' && isalpha(c14) && isdigit(name[2])) return true;
  return false;
}

void 
PDBBase::set_header(const std::string header_in)
{
  A_header = header_in; 
}

bool
PDBBase::get_res_info(const std::string atom_name, const std::string res_name,
                      atom_t* atom)
{
  const pdb_atom_info_t* atom_info;
  atom_info = PDB_residues::get_atom_info(atom_name, res_name);
  if(!atom_info) return false;

  atom->name = atom_info->atom;
  // Kludge for MAIN_CHAIN
  if(atom_info->res != MAIN_CHAIN) atom->res = atom_info->res;
  else atom->res = PDB_residues::string_to_residue(res_name);
  atom->vdw_radius = atom_info->vdw_radius;
  atom->orbit = atom_info->orbit;
  atom->act_type = atom_info->act_type;
  atom->charge = atom_info->charge;
  atom->hydro = PDB_residues::get_atom_hydrophobicity(atom->res, atom->name);

  return true;
}

void
PDBBase::read_het_record(std::string& het_line)
{
  het_record_type tmp;

  tmp.hetID = het_line.substr(7, 3); 
  tmp.chainID = het_line[12];
  my_strtoui(het_line.substr(13,4), &tmp.seqNum);
  tmp.iCode = het_line[17];
  my_strtoui(het_line.substr(20,5), &tmp.numHetAtoms);
  tmp.text = het_line.substr(30);
  A_het_records.push_back(tmp);
}

void
PDBBase::add_het_record(atom_vci atom)
{
  het_record_type tmp;
  tmp.hetID = atom->res_str;
  tmp.chainID = atom->chainID;
  tmp.seqNum =  atom->res_num;
  tmp.iCode = atom->iCode;
  tmp.numHetAtoms = 0;
  tmp.text = "Added by SimSite3D -- no HET entry";
  A_het_records.push_back(tmp);
}

void
PDBBase::read_modres_record(std::string& modres_line)
{
  modres_record_type tmp;
  tmp.idCode = modres_line.substr(7, 4);
  tmp.resName = modres_line.substr(12, 3);
  tmp.chainID = modres_line[16];
  my_strtoui(modres_line.substr(18, 4), &tmp.seqNum);
  tmp.iCode = modres_line[22];
  tmp.stdRes = PDB_residues::string_to_residue(modres_line.substr(24, 3));
  tmp.comment = modres_line.substr(29);
  A_modres_records.push_back(tmp);
}
