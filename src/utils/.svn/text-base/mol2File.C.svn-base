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
#include <sstream>
#include <iomanip>
#include <mol2File.H>
#include <basics.H>
#include <Timer.H>
#include <list>
#include <set>
#include <stack>

using namespace SimSite3D;

const my_float_t mol2File::ATOMS_TOO_CLOSE = 1.75;
const my_float_t mol2File::ATOMS_TOO_FAR = 4.5;
const std::string mol2File::_fname = "mol2File.C";
const my_float_t mol2File::MINIMAL_CHARGE = 0.5;

void
mol2File::report_stats(std::ostream& out)
{
//blah
}

mol2File::mol2File(const std::string filename, const verbose_level_t verbosity)
 : CoordFile(filename, verbosity)
{
  a_num_heavy_atoms = 0;
  if(fail()) return;

  std::ifstream mol2_file;
  if(!open_ifstream(mol2_file, CoordFile::name())) red_light();
  else{  
    atom_ids_are_consecutive = true;
    is_frag_file = false;
    init_atom_table_by_string();
    init_bond_table_by_string();
    if(!read_data(mol2_file)) red_light();
  }
  if(!fail()){
    frag_atoms.resize(num_atoms(), true);
    assign_act_types();
  }
}

/* -- defunct
mol2File::mol2File(const mol2File& src, const CoordFile& prot_atoms,
                   const my_float_t* R, const my_float_t* T)
 : CoordFile(src)
{
  if(&src == this) return;

  a_num_heavy_atoms = 0;
  is_frag_file = true;

  // This is a copy cstr so we need to make a copy of the source atoms 
  // first and then transform them.
  point_storage<atom_t, std::list<atom_t> > tmp_atoms;
  for(atom_vci src_atom = src.atoms_begin(); src_atom != src.atoms_end();             ++src_atom) tmp_atoms.push_back(*src_atom);
  
  // transform the copy's atoms
  tmp_atoms.transform(R, T); 

  // Note the advance by one -- will not use atom_nums[0];
  std::vector<int> atom_nums(tmp_atoms.size() + 1, -1);
  std::vector<int>::iterator num = atom_nums.begin() + 1;

  uint cnt = 1;
  point_storage<atom_t, std::list<atom_t> >::const_iterator a;
  a = tmp_atoms.begin();
  for(uint i = 0; i < tmp_atoms.size(); ++i, ++num, ++a){
    if(!prot_atoms.overlapping_atoms(*a)){
      append_atom(*a);
      rest_of_lines.push_back(src.rest_of_lines[i]);
      *num = cnt;
      ++cnt;
    }
  }

  // remove the bonds for those atoms that were removed
  num = atom_nums.begin();
  for(mol2_bond_vci b = src.bonds.begin(); b < src.bonds.end(); ++b){
    mol2_bond_t tmp;
    tmp.atom_num1 = *(num + b->atom_num1);
    tmp.atom_num2 = *(num + b->atom_num2);
    if(tmp.atom_num1 > 0 && tmp.atom_num2 > 0){
      tmp.type = b->type;
      tmp.status_bits = b->status_bits;
      bonds.push_back(tmp);
    }
  }
}
*/

mol2File::~mol2File()
{

}

mol2File::mol2File(const mol2File& src, const BoundingVolume &site_vol,
//mol2File::mol2File(const mol2File& src, const CoordFile &site_pts,
                   const CoordFile& prot_atoms, const uint min_num_atoms,
                   const my_float_t* R, const my_float_t* T)
 : CoordFile(src)
{
  if(&src == this) return;

  a_num_heavy_atoms = 0;
  is_frag_file = true;

  // This is a copy cstr so we need to make a copy of the source atoms 
  // first and then transform them.
  atom_vec_t tmp_atoms(src.num_atoms());
  uint max_atom_num = 0;
  atom_vci src_atom;
  for(src_atom = src.atoms_begin(); src_atom != src.atoms_end(); ++src_atom){
    tmp_atoms.push_back(*src_atom);
    max_atom_num = 
      (max_atom_num < src_atom->atom_num ? src_atom->atom_num : max_atom_num);
  }

  // transform the copy's atoms
  tmp_atoms.transform(R, T); 

  // Get all heavy atoms which are within 2.0(A) of any sitemap point and
  // do not significantly bump with the protein.  All atom numbers for those
  // we will keep are first marked as 1 since the numbers may need to be
  // adjusted if there are fragments which have too few atoms
  std::map<uint, atom_vci> H_atoms;
  std::vector<uint> atom_nums(max_atom_num + 1, 0);
  std::vector<uint>::iterator atom_nums_beg = atom_nums.begin();
  for(atom_vci a = tmp_atoms.begin(); a < tmp_atoms.end(); ++a){
    if(a->name == H) H_atoms[a->atom_num] = a;
#if 0
    // Use sitemap points to determine the vol
    else{
      my_float_t d;
      site_pts.closest_point(a->pos, &d);
      if(d < 2.5 && !prot_atoms.overlapping_atoms(*a))
        *(atom_nums_beg + a->atom_num) = 1;
    }
#endif
#if 1
    else if(site_vol.contains(a->pos) && !prot_atoms.overlapping_atoms(*a)) 
      *(atom_nums_beg + a->atom_num) = 1;
#endif
  }

  std::list<std::map<uint, uint> > frags;
  for(mol2_bond_vci b = src.bonds.begin(); b < src.bonds.end(); ++b){
    uint a1 = b->atom_num1;
    uint a2 = b->atom_num2;
    uint keep_a1 = *(atom_nums_beg + b->atom_num1);
    uint keep_a2 = *(atom_nums_beg + b->atom_num2);
    // If we are not keeping either atom, skip this bond
    if(!keep_a1 && !keep_a2) continue;

    // Check if the atoms a1 & a2 already belong to fragments
    std::list<std::map<uint, uint> >::iterator frag;
    std::list<std::map<uint, uint> >::iterator f1 = frags.end();
    std::list<std::map<uint, uint> >::iterator f2 = frags.end();
    for(frag = frags.begin(); frag != frags.end(); ++frag){
      if(keep_a1 && frag->find(a1) != frag->end()) f1 = frag;
      else if(keep_a2 && frag->find(a2) != frag->end()) f2 = frag;
    }

    // We are keeping both atoms
    if(keep_a1 && keep_a2){
      // Atoms are both already marked -- check if they link 2 fragments
      if(f1 != frags.end() && f2 != frags.end()){
        // Atoms are both in the same fragment, nothing to do here
        if(f1 == f2) continue;
        // Copy all atoms from fragment f2 into fragment f1 and delete f2
        std::map<uint, uint>::iterator m_iter;
        for(m_iter = f2->begin(); m_iter != f2->end(); ++m_iter)
          f1->insert(std::pair<uint, uint>(m_iter->first, 0));
        frags.erase(f2);
      // else if only 1 map was found, add the other atom to the found map
      // (they belong to the same fragment as they are covalently bonded 
      //  together)
      }else if(f1 != frags.end()) f1->insert(std::pair<uint, uint>(a2, 0));
      else if(f2 != frags.end()) f2->insert(std::pair<uint, uint>(a1, 0));
      // else both atoms are kept 
      // add the two atoms to their own map (fragment) append the new map 
      // (fragment) to list 
      else{
        std::map<uint, uint> tmp_map;
        tmp_map.insert(std::pair<uint, uint>(a1, 0));
        tmp_map.insert(std::pair<uint, uint>(a2, 0));
        frags.push_back(tmp_map);
      }
    // We are only keeping atom a1 and it is not in a fragment, give it its
    // own fragment
    }else if(keep_a1 && f1 == frags.end()){
     std::map<uint, uint> tmp_map;
     tmp_map.insert(std::pair<uint, uint>(a1, 0));
     frags.push_back(tmp_map);
    // We are only keeping atom a2 and it is not in a fragment, give it its
    // own fragment
    }else if(keep_a2 && f2 == frags.end()){
     std::map<uint, uint> tmp_map;
     tmp_map.insert(std::pair<uint, uint>(a2, 0));
     frags.push_back(tmp_map);
    }
  }

  // Zero out the atoms which belong to fragments with too few heavy atoms
  std::list<std::map<uint, uint> >::iterator frag;
  for(frag = frags.begin(); frag != frags.end(); ++frag){
    if(frag->size() < min_num_atoms){
      std::map<uint, uint>::iterator atom;
      for(atom = frag->begin(); atom != frag->end(); ++atom)
        *(atom_nums_beg + atom->first) = 0;
    }
  }

  // Assign (counting from 1) consecutive atom numbers to the kept atoms
  uint cnt = 0;
  for(std::vector<uint>::iterator n = atom_nums_beg; n < atom_nums.end(); ++n)
    *n = (*n ? ++cnt : 0);

  // store bonds
  std::map<uint, atom_vci>::const_iterator H;
  for(mol2_bond_vci b = src.bonds.begin(); b < src.bonds.end(); ++b){
    uint a1 = *(atom_nums_beg + b->atom_num1);
    uint a2 = *(atom_nums_beg + b->atom_num2);
    // both atoms in this bond are heavy atoms we want to keep
    if(a1 && a2) add_bond(a1, a2, b->type, b->status_bits);
    // Only the first atom in this bond is kept, and the second could be a 
    // hydrogen atom
    else if(a1){
      H = H_atoms.find(b->atom_num2);
      if(H != H_atoms.end()){
        *(atom_nums_beg + b->atom_num2) = ++cnt;  
        add_bond(a1, cnt, b->type, b->status_bits);
      } 
    // Only the second atom in this bond is kept, and the first could be a 
    // hydrogen atom
    }else if(a2){
      H = H_atoms.find(b->atom_num1);
      if(H != H_atoms.end()){
        *(atom_nums_beg + b->atom_num1) = ++cnt;  
        add_bond(cnt, a2, b->type, b->status_bits);
      } 
    }
  }

  // store atoms -- the write routine blindly assigns consecutive numbers
  // to the atoms bases on their order in the vector ==> need to place the
  // atoms in the vector ordered by atom number.  Do the usual cheese and
  // make use of a map
  std::map<uint, atom_vi> atom_map; 
  frag_atoms.resize(tmp_atoms.size(), false);
  for(atom_vi a = tmp_atoms.begin(); a < tmp_atoms.end(); ++a){
    uint atom_num = *(atom_nums_beg + a->atom_num);
    if(atom_num) atom_map[atom_num] = a;
  }
  
  std::map<uint, atom_vi>::iterator A_iter;
  for(A_iter = atom_map.begin(); A_iter != atom_map.end(); ++A_iter){
    atom_vi atom = A_iter->second;    
    uint idx = atom->atom_num - 1;
    atom->atom_num = A_iter->first; 
    append_atom(*atom);
    frag_atoms[idx] = true; 
  }

  generate_nbrs_vec();
}

const uint 
mol2File::num_heavy_atoms()
{
  if(a_num_heavy_atoms) return a_num_heavy_atoms;

  for(atom_vci a = atoms_begin(); a != atoms_end(); ++a)
    if(a->name != H) ++a_num_heavy_atoms; 

  return a_num_heavy_atoms;
}

void
mol2File::add_bond(const uint a1, const uint a2, const bond_type type, 
                   const std::string status_bits)
{
  mol2_bond_t bond;
  bond.atom_num1 = a1;
  bond.atom_num2 = a2;
  bond.type = type;
  bond.status_bits = status_bits;
  bonds.push_back(bond);
}

bool
mol2File::read_data(std::ifstream &mol2_file)
{  
  // For now assume only one of each section
  bool have_molecule = false;
  bool have_atom = false;
  bool have_bond = false;
  uint num_bonds;
  uint atoms_cnt;

  std::string line;
  std::string tag = "@<TRIPOS>";
  while(std::getline(mol2_file, line)){
    if(line.substr(0,9) == tag){
      std::string section = line.substr(9);
      // We must use substr here since it is likely that files generated on 
      // other operating systems (besides *nix) will get pushed through
      // this function.
      if(section.substr(0,8) == "MOLECULE" && !have_molecule){
        if(!read_molecule_section(mol2_file, &atoms_cnt, &num_bonds))
          return false;
        have_molecule = true; 
      }else if(section.substr(0,4) == "ATOM" && have_molecule && !have_atom){
        if(!read_atom_section(mol2_file, atoms_cnt)) return false;
        have_atom = true;
      }else if(section.substr(0,4) == "BOND" && have_atom && !have_bond){
        if(!read_bond_section(mol2_file, num_bonds)) return false;
        have_bond = true; 
        green_light();
      }
    }
  }   
  
  if(have_molecule && have_atom && have_bond) return true;
  else return false;
}

bool
mol2File::write(const std::string ofname)
{
  std::ofstream out;
  if(!open_ofstream(out, ofname)) return false;
  if(atom_types_to_str.size() == 0) init_atom_types_to_str();

  if(is_frag_file){
    char tmp[80];
    out << "# SimSite3D fragment of " << name() << "\n"
        << "# Created: " << Timer::get_local_time(tmp, 80) << "\n";
  }

  std::string mol_name, path;
  if(!get_path_and_struct_id(ofname, &path, &mol_name)) return false;
  out << "@<TRIPOS>MOLECULE\n" << mol_name << "\n"
      << std::setw(5) << atoms_end() - atoms_begin() << " "
      << std::setw(5) << bonds.size() << " "
      << std::setw(5) << 0 << " "
      << std::setw(5) << 0 << " "
      << std::setw(5) << 0 << "\n"
      << "SMALL\n";
  // we don't recompute charges for fragments
  out << A_charge_type << "\n\n";

  out << "@<TRIPOS>ATOM\n";
  uint cnt = 1;
  out << std::fixed << std::setprecision(4);
  for(atom_vci a = atoms_begin(); a < atoms_end(); ++a, ++cnt){
    out << std::setw(7) << cnt << " "
        << std::left << std::setw(8) << a->name_str;
    for(uint i = 0; i < 3; ++i) out << std::right << std::setw(10) << a->pos[i];
    out << " " << std::left << std::setw(6) 
        << atom_types_to_str[a->name][a->orbit] << std::right;
    if(a->subst_id){
      out << " " << a->subst_id;
      if(a->subst_name.length()){
        out << " " << a->subst_name;
        // We don't want to write out the charge sums, although it could be
        // interesting or possibly useful for debugging
        if(A_charge_type != "NO_CHARGES") 
          out << " " << std::setw(6) << a->orig_charge;
      }      
    }
    out << "\n";
  } 

  init_bond_types_to_str();
  out << "@<TRIPOS>BOND\n";
  cnt = 1;
  for(mol2_bond_vci b = bonds.begin(); b < bonds.end(); ++b, ++cnt)
    out << std::setw(7) << cnt << " "
        << std::setw(7) << b->atom_num1 << " "
        << std::setw(7) << b->atom_num2 << " "
        << std::setw(3) << bond_types_to_str[b->type]
        << b->status_bits << "\n";

  return true;
}

bool 
mol2File::read_molecule_section(std::ifstream &mol2_file, 
                                uint* num_atoms_out,
                                uint* num_bonds_out)
{
  *num_bonds_out = 0;
  *num_atoms_out = 0;

  // Assume we have a "correct" mol2 header
  std::getline(mol2_file, A_mol_name);
  mol2_file >> *num_atoms_out >> *num_bonds_out;  
  std::string line;
  std::getline(mol2_file, line);
  std::getline(mol2_file, A_mol_type);
  std::getline(mol2_file, A_charge_type);

  if(!mol2_file){
    err_msg(_fname, "read_molecule_section", 
            "Unable to read the molecule section");
    return false;
  }
  return true;
}
 
bool 
mol2File::read_atom_section(std::ifstream &mol2_file, const uint atoms_cnt)
{
  if(atom_types_to_str.size() == 0) init_atom_types_to_str();
  atom_t atom;
  my_float_t pos[3];
  std::string line;
  uint i;
  for(i = 0; i < atoms_cnt && std::getline(mol2_file, line); ++i){
    if(line.substr(0,9) == "@<TRIPOS>") break;

    std::istringstream lstr(line);
    lstr >> atom.atom_num >> atom.name_str >> pos[0] >> pos[1] >> pos[2];
    std::copy(pos, pos + 3, atom.pos);
    std::string sybil_atom_type;
    lstr >> sybil_atom_type;

    // Work around for fixed format mol2 issues (Pfizer)
    if(sybil_atom_type.length() > 5){
      my_strtoi(sybil_atom_type.substr(5), &(atom.subst_id));
      sybil_atom_type = sybil_atom_type.substr(0, 5);
    }else lstr >> atom.subst_id;
    lstr >> atom.subst_name >> atom.charge;
    // Save as orig charge as well since we might call this.charge_sums()
    atom.orig_charge = atom.charge; 

    // Work around for H_ADD 
    if(sybil_atom_type == "H_ADD") sybil_atom_type = "H";

    // Work around for Du entry in sybil atom type field
    if(sybil_atom_type == "Du"){
      //if(atom.name_str.length() > 2){
       // bad_file_line(_fname, "read_atom_section", name(), line);
        //return false;
     // }
      sybil_atom_type = atom.name_str[0];
      if(atom.name_str.length() >= 2)
        sybil_atom_type += tolower(atom.name_str[1]);
    }

    if(!look_up_atom_by_string(sybil_atom_type, &atom)){
      bad_file_line(_fname, "read_atom_section", name(), line);
      return false;
    }
    
    // Add in the (adjusted) hydrophobicity values
    if(atom.name == O) atom.hydro = 295;
    else if(atom.name == N) atom.hydro = 115;
    else if(atom.name == C || atom.name == S) atom.hydro = -155;
    else atom.hydro = 82;

    // Check for premature eof or other failure
    if(!mol2_file){
      std::ostringstream msg;
      msg << "Incomplete atom line in " << name();
      err_msg(_fname, "read_atom_section", msg.str());
      return false;
    }

    append_atom(atom);    
    if(atom.atom_num - 1 != i) atom_ids_are_consecutive = false;
  }
  
  if(i < atoms_cnt){
    std::ostringstream msg;
    msg << name() << "\nAtom section: Expected to read " << atoms_cnt 
        << " lines.\nOnly " << i << " lines were found.\n";
    err_msg(_fname, "read_atom_section", msg.str());
    return false;
  }
  return true;
}

bool
mol2File::read_bond_section(std::ifstream& mol2_file, const uint num_bonds)
{
  bonds.clear();
  bonds.reserve(num_bonds);

  std::string line;
  uint i;
  for(i = 0; i < num_bonds && std::getline(mol2_file, line); ++i){
    if(line.substr(0,9) == "@<TRIPOS>") break;

    mol2_bond_t bond;
    std::string bond_str;
    std::istringstream lstr(line);
    uint bond_id;
    lstr >> bond_id >> bond.atom_num1 >> bond.atom_num2 >> bond_str 
         >> bond.status_bits;
    look_up_bond_by_string(bond_str, &bond);

    // Check for premature eof or other failure
    if(!mol2_file){
      std::ostringstream msg;
      msg << "Incomplete bond line in " << name();
      err_msg(_fname, "read_bond_section", msg.str());
      return false;
    }
    bonds.push_back(bond);
  }

  if(i < num_bonds){
    std::ostringstream msg;
    msg << name() << "\nBond section: Expected to read " << num_bonds 
        << " lines.\nOnly " << i << " lines were found.\n";
    err_msg(_fname, "read_atom_section", msg.str());
    return false;
  }
  return true;
}

void
mol2File::init_atom_table_by_string()
{
  for(uint i = 0; i < mol2_atom_table_size; ++i)
    atoms_by_type_str[mol2_atom_table[i].atom_type_str] = &(mol2_atom_table[i]);
}

void
mol2File::init_bond_table_by_string()
{
  for(uint i = 0; i < mol2_bond_table_size; ++i)
    bonds_by_type_str[mol2_bond_table[i].bond_type_str] = 
      mol2_bond_table[i].type;
}

void
mol2File::init_bond_types_to_str()
{
  for(uint i = 0; i < mol2_bond_table_size; ++i)
    bond_types_to_str[mol2_bond_table[i].type] = 
      mol2_bond_table[i].bond_type_str;
}

void
mol2File::init_atom_types_to_str()
{
  for(uint i = 0; i < mol2_atom_table_size; ++i)
    atom_types_to_str[mol2_atom_table[i].atom][mol2_atom_table[i].orbit] =
      mol2_atom_table[i].atom_type_str;
}


bool
mol2File::look_up_atom_by_string(const std::string atom_str, atom_t* atom)
{
  std::map<std::string, const mol2_atom_info_t*>::const_iterator iter;
  iter = atoms_by_type_str.find(atom_str);
  if(iter != atoms_by_type_str.end()){
    atom->name = iter->second->atom;
    atom->orbit = iter->second->orbit;
    atom->vdw_radius = iter->second->vdw_radii;
    return true;
  }
  std::ostringstream msg;
  msg << "Unknown mol2 atom type and orbit (" << atom_str << ")";
  err_msg(_fname, "look_up_atom_by_string", msg.str());
  return false;
}

bool
mol2File::look_up_bond_by_string(const std::string bond_str, mol2_bond_t* bond)
{
  std::map<std::string, bond_type>::const_iterator iter;
  iter = bonds_by_type_str.find(bond_str);
  if(iter != bonds_by_type_str.end()){
    bond->type = iter->second;
    return true;
  }
  return false;
}

bool
mol2File::assign_act_types()
{
// not sure if it is OK to comment this out -- need to go through all of the
// posibilities again
#if 0
  if(is_frag_file){
    std::string msg = "This function cannot be used on ligand fragments,";
    msg += " it expects a chemically\n correct molecule";
    err_msg(_fname, "assign_act_types", msg);
    return false;
  }
#endif

  // Don't do this twice
  if(A_nbrs_vec.size()) return false;

  generate_nbrs_vec();
#if 0
  // Build a vector of the neighbors for each atom
  A_nbrs_vec.resize(num_atoms() + 1);
  for(mol2_bond_vci b_iter = bonds.begin(); b_iter != bonds.end(); ++b_iter){
    size_t idx1 = b_iter->atom_num1 - 1;
    size_t idx2 = b_iter->atom_num2 - 1;
    (A_nbrs_vec.begin() + idx1)->push_back(atoms_begin() + idx2);
    (A_nbrs_vec.begin() + idx2)->push_back(atoms_begin() + idx1);
  }
#endif

  // Explicity coding of some Dock rules for hbonding atoms -- 
  // because we need to modify the stored atoms
  // we must use the protected method of accessing the atoms
  for(atom_vi a = atoms_storage_begin(); a != atoms_storage_end(); ++a){
    if(a->name == O){
      if(a->orbit == SP2){
        a->act_type = ACCEPTOR;
        continue;
      }

      bool N_nbr = false;
      bool H_nbr = false;
      bool C_nbr = false;
      const std::vector<atom_vci> &nbrs = get_nbrs(a);
      std::vector<atom_vci>::const_iterator nbr_iter;
      for(nbr_iter = nbrs.begin(); nbr_iter < nbrs.end(); ++nbr_iter){
        if((*nbr_iter)->name == H) H_nbr = true;
        else if((*nbr_iter)->name == N) N_nbr = true;
        else if((*nbr_iter)->name == C) C_nbr = true;
      }
      if(H_nbr && C_nbr && nbrs.size() == 2)a->act_type = DONEPTOR;
// Try this
      else if(!N_nbr) a->act_type = ACCEPTOR;
      
#if 0
      else if(!H_nbr && nbrs.size() == 2){
        a->act_type = ACCEPTOR;
      }
      else if(!H_nbr && !N_nbr && nbrs.size() == 1){
        a->act_type = ACCEPTOR;
      }
#endif

    }else if(a->name == N){
      bool O_nbr = false;
      bool H_nbr = false;
      const std::vector<atom_vci> &nbrs = get_nbrs(a);
      std::vector<atom_vci>::const_iterator nbr_iter;
      for(nbr_iter = nbrs.begin(); nbr_iter < nbrs.end(); ++nbr_iter){
        if((*nbr_iter)->name == H) H_nbr = true;
        else if((*nbr_iter)->name == O) O_nbr = true;
      }
      if(H_nbr) a->act_type = DONOR;
      else if(nbrs.size() < 3 && O_nbr == false) a->act_type = ACCEPTOR;
      else a->act_type = NOTHING; 
    }else if(a->name == F) a->act_type = ACCEPTOR;
    else if(a->name == CL) a->act_type = ACCEPTOR;
  }  

  return true;
}

bool
mol2File::calc_charge_sums()
{
  if(A_charge_type == "NO_CHARGES"){
    std::string msg = "This function cannot be used on ligands without ";
    msg += " charges.";
    err_msg(_fname, "calc_charge_sums", msg);
    return false;
  }

#if 0
  if(is_frag_file){
    std::string msg = "This function cannot be used on ligand fragments,";
    msg += " it expects a chemically\n correct molecule";
    err_msg(_fname, "calc_charge_sums", msg);
    return false;
  }
#endif

  // 0. Zero out the charges
  for(atom_vi A = atoms_storage_begin(); A < atoms_storage_end(); ++A){
    A->charge = 0.0;
  }

  /* 1. Sum charges of all H neighbors to heavy atoms (non-H atoms) and record
   * the new charges in "charge_sum" entry.
   */
  std::vector<my_float_t> temp_charge_sums(num_atoms());
  std::vector<my_float_t>::iterator charge_sum_iter = temp_charge_sums.begin();
  for(atom_vci A = atoms_begin(); A < atoms_end(); ++A, ++charge_sum_iter){
    *charge_sum_iter = A->orig_charge;
    if(A->name == H) continue;

    const std::vector<atom_vci> &nbrs = get_nbrs(A);
    for(size_t i = 0; i < nbrs.size(); ++i)
      if(nbrs[i]->name == H) *charge_sum_iter += nbrs[i]->orig_charge;
  }

  /* 2. Sum charges of all non-H neighbors to atom N or O. If they only have C
   * neighbors and these C neighbors have other N or O neighbors, sum also 
   * the charges of their N or O neighbors.
   */
  charge_sum_iter = temp_charge_sums.begin();
  for(atom_vi A = atoms_storage_begin(); A < atoms_storage_end(); ++A){
    A->charge = *charge_sum_iter;
    if(A->name != N && A->name != O){
      ++charge_sum_iter; 
      continue;
    }


    int num_of_non_H_nbrs = 0;
    bool C_nbr_only = true;   
    const std::vector<atom_vci> &nbrs = get_nbrs(A);
    atom_vci hvy_nbr = atom_t::NULL_ATOM_VCI;;
    for(size_t i = 0; i < nbrs.size(); ++i){
      if(nbrs[i]->name == H) continue;
    
      hvy_nbr = nbrs[i]; 
      A->charge += temp_charge_sums[nbrs[i] - atoms_begin()]; 
      ++num_of_non_H_nbrs;
      if(nbrs[i]->name != C) C_nbr_only = false;
    }

    /* If there is only 1 direct non_H neighbor, we need to sum up charges on
     * the 2nd layer neighbors 
     */
    if(num_of_non_H_nbrs == 1){
      const std::vector<atom_vci> &second_nbrs = get_nbrs(hvy_nbr);
      for(size_t i = 0; i < second_nbrs.size(); ++i){
        if(second_nbrs[i]->name == H || second_nbrs[i] == A) continue;

        A->charge += temp_charge_sums[second_nbrs[i] - atoms_begin()];
      }

    /* If all neighbors to atom i are carbons, we need to sum charges of any 
     * other N/O neighbor to these C neighbors 
     *          l----Carbon 
     * atom i---l---Carbon 
     *          l---Carbon------N/O 
     */
    }else if(C_nbr_only){
      for(size_t i = 0; i < nbrs.size(); ++i){
        if(nbrs[i]->name != C) continue;

        const std::vector<atom_vci> &C_nbrs = get_nbrs(nbrs[i]);
        for(size_t j = 0; j < C_nbrs.size(); ++j)
          if(C_nbrs[j] != A && (C_nbrs[j]->name == O || C_nbrs[j]->name == N))
            A->charge += temp_charge_sums[C_nbrs[j] - atoms_begin()];
      }
    }
    ++charge_sum_iter; 
  }

  /* Now we want to remove those small charges and charges for non-heavy atoms 
   */
  for(atom_vi A = atoms_storage_begin(); A < atoms_storage_end(); ++A){
//    std::cout << A->name_str << " charge : " << A->charge << "    ";
    if(A->name == H || A->name == C || A->name == S) A->charge = 0.0;
    else if(-1.0*MINIMAL_CHARGE < A->charge &&
            A->charge < MINIMAL_CHARGE) A->charge = 0.0;
//    std::cout << "\"new\" charge : " << A->charge << "\n";
  }

  return true;
}

void 
mol2File::generate_nbrs_vec()
{
  // Build a vector of the neighbors for each atom
  A_nbrs_vec.resize(num_atoms() + 1);
  for(mol2_bond_vci b_iter = bonds.begin(); b_iter != bonds.end(); ++b_iter){
    size_t idx1 = b_iter->atom_num1 - 1;
    size_t idx2 = b_iter->atom_num2 - 1;
    (A_nbrs_vec.begin() + idx1)->push_back(atoms_begin() + idx2);
    (A_nbrs_vec.begin() + idx2)->push_back(atoms_begin() + idx1);
  }
}


/*
mol2File::mol2File(std::string fname)
{
  _fail = true;
  num_atoms = 0;
  orig_positions = NULL;
  rings = NULL;
  atom_ids_are_consecutive = true;
  mol2_file_name = fname;

  std::ifstream mol2_file;
  if(!open_ifstream(mol2_file, fname)) return;
  
  // For now assume only one of each section
  bool have_molecule = false;
  bool have_atom = false;
  bool have_bond = false;
  uint num_bonds = 0;

  std::string line;
  std::string tag = "@<TRIPOS>";
  while(std::getline(mol2_file, line)){
    if(line.substr(0,9) == tag){
      std::string section = line.substr(9);
      if(section == "MOLECULE" && !have_molecule){
        if(!read_molecule_section(mol2_file, &molecule, num_atoms, num_bonds))
          return;
        have_molecule = true; 
      }else if(section == "ATOM" && have_molecule && !have_atom){
	positions = new my_float_t[3 * num_atoms];
        if(!read_atom_section(mol2_file, &molecule, positions, num_atoms)) 
          return;
        have_atom = true;
      }else if(section == "BOND" && have_atom && !have_bond){
        if(!read_bond_section(mol2_file, &molecule, &bonds, num_bonds)) return;
        have_bond = true; 
        _fail = false;
      }
    }
  }   
  orig_positions = positions;
  enumerate_rings();
}

mol2File::~mol2File()
{
  if(orig_positions != positions && orig_positions) delete [] orig_positions;
  orig_positions = 0;
  molecule.clear();
  bonds.clear();
  stuff.clear();
  if(positions) delete [] positions;
  positions = 0;
  if(rings) delete rings;
  rings = 0;
  _fail = true;
}

void
mol2File::revert_positions()
{
  if(positions != orig_positions)
    std::copy(orig_positions, orig_positions + 3*num_atoms, positions);
}

void
mol2File::transform_positions(const my_float_t* R, const my_float_t* T)
{
  if(!orig_positions || orig_positions == positions){
    orig_positions = new my_float_t[3*num_atoms];
    memcpy(orig_positions, positions, 3*num_atoms*my_float_size);
    move_positions(num_atoms, positions, orig_positions, R, T);
  }else{
    my_float_t* prev_pos = new my_float_t[3*num_atoms];
    memcpy(prev_pos, positions, 3*num_atoms*my_float_size);
    move_positions(num_atoms, positions, prev_pos, R, T);
    delete[] prev_pos;
  }
}

void 
mol2File::setup_stuff()
{
  // make this map a "static" so that it is only run once per program.
  
  // SP1 as orig_atom
  stuff[SP1].push_back(SP3);

  // SP2 as orig_atom
  stuff[SP2].push_back(SP2);
  stuff[SP2].push_back(SP3);
  stuff[SP2].push_back(SP4);
  stuff[SP2].push_back(AR);
  stuff[SP2].push_back(PL3);
  stuff[SP2].push_back(AM);

  // SP3 as orig_atom
  stuff[SP3].push_back(SP1);
  stuff[SP3].push_back(SP2);
  stuff[SP3].push_back(SP3);
  stuff[SP3].push_back(SP4);
  stuff[SP3].push_back(AR);
  stuff[SP3].push_back(PL3);
  stuff[SP3].push_back(AM);

  // SP4 as orig_atom
  stuff[SP4].push_back(SP2);
  stuff[SP4].push_back(SP3);
  
  // AROMATIC as orig_atom
  stuff[AR].push_back(SP2);
  stuff[AR].push_back(SP3);
  stuff[AR].push_back(PL3);
  stuff[AR].push_back(AR);
  
  // TRIGONAL_PLANAR as orig_atom
  stuff[PL3].push_back(SP3);
  stuff[PL3].push_back(SP2);
  stuff[PL3].push_back(AR);

  // AMIDE as orig_atom
  stuff[AM].push_back(AM);
  stuff[AM].push_back(AR);
  stuff[AM].push_back(SP2);
  stuff[AM].push_back(SP3);

}

void
mol2File::mark_kept_heavy_atoms(const mol2_atom_vec& atoms, 
                                const my_float_t* trans_pos,
                                const pts_vci pts_begin, const pts_vci pts_end,
                                my_float_t radius, 
                                std::vector<bool>* kept_atoms)
{
  // Get centroid of template points
  my_float_t centroid[3];
  memset(centroid, 0, 3*my_float_size);
  for(pts_vci pt = pts_begin; pt < pts_end; ++pt)
    for(uint i = 0; i < 3; i++) centroid[i] += pt->pos[i];
  uint sz = pts_end - pts_begin;
  for(uint i = 0; i < 3; i++) centroid[i] /= sz;

  // Find farthest point from centroid
  my_float_t max_dist = 0;
  for(pts_vci pt = pts_begin; pt < pts_end; ++pt){
    my_float_t tmp = dist(centroid, pt->pos);
    if(tmp > max_dist) max_dist = tmp;
  }

  my_float_t rad = max_dist + 0.5;
  std::vector<bool>::iterator keep_atom = kept_atoms->begin(); 
  for(mol2_atom_vci atom = atoms.begin(); atom < atoms.end(); 
      ++atom, ++keep_atom, trans_pos += 3)
    // Skip hydrogens since their inclusion is based on the heavy atoms to
    // which they are covalently bonded
    if(atom->name != H && dist(trans_pos, centroid) <= rad) *keep_atom = true;
}

void
mol2File::keep_bond(mol2_atom_vci a, mol2_atom_vci b, 
                    bond_marks* my_bonds)
{
  if(a < b) (*my_bonds)[a][b] = true;
  else (*my_bonds)[b][a] = true;
}

void
mol2File::cleave_bonds(const mol2_atom_vec& full_atom_set, 
                       const bond_map& full_bond_set,
                       std::vector<bool>* kept_atoms,
                       bond_marks* kept_bonds)
{

  int count = 0;
  for(uint i = 0; i < kept_atoms->size(); ++i)
    if((*kept_atoms)[i]) ++count;
  std::cout << "\nbefore ring check\n";
  std::cout << "total num atoms in lig: " << kept_atoms->size() << "\n";
  std::cout << "num atoms in lig frag: " << count << "\n";


  // Find all the rings which contain two or more atoms that are to be kept
  uint num_rings = rings->num_cycles();
  std::vector<uint> kept_in_ring(num_rings, 0);
  mol2_atom_vci atom = full_atom_set.begin();
  for(uint i = 0; i < kept_atoms->size(); ++i, ++atom){
    if((*kept_atoms)[i] == false) continue;

    std::vector<uint> rnums = rings->cycles_that_contain(atom);
    std::vector<uint>::const_iterator r;
    for(r = rnums.begin(); r < rnums.end(); ++r) ++(kept_in_ring[*r]);
  }

  // Get the rings which have two or more atoms to be kept and mark all as kept
  for(uint i = 0; i < kept_in_ring.size(); ++i){
    if(kept_in_ring[i] < 2) continue;

    // Inefficient but we typically have fewer than 5 rings
    std::vector<mol2_atom_vci> atoms_in_ring;
    rings->get_cycle(i, &atoms_in_ring); 
    atom = full_atom_set.begin();
    for(uint j = 0; j < kept_atoms->size(); ++j, ++atom){
      if((*kept_atoms)[j] == true) continue;

      std::vector<mol2_atom_vci>::iterator ring_atom = atoms_in_ring.begin(); 
      for( ; ring_atom < atoms_in_ring.end(); ++ring_atom)
        if(*ring_atom == atom){
          (*kept_atoms)[j] == true;
          break;
        }
    }
  }

  // Keep all bonds for which we wish to keep both atoms
  // Inefficient at the present
  const mol2_atom_vci atoms_beg = full_atom_set.begin();
  atom = full_atom_set.begin();
  for(uint i = 0; i < kept_atoms->size(); ++i, ++atom){
    if((*kept_atoms)[i] == false) continue;

    std::vector<mol2_atom_vi>::const_iterator nbr;
    for(nbr = atom->neighbors.begin(); nbr < atom->neighbors.end(); ++nbr)
      if((*kept_atoms)[*nbr - atoms_beg]) keep_bond(atom, *nbr, kept_bonds);
      else if((*nbr)->name == H){
        keep_bond(atom, *nbr, kept_bonds);
        (*kept_atoms)[*nbr - atoms_beg] = true;
      }
  }

count = 0;
  for(uint i = 0; i < kept_atoms->size(); ++i)
    if((*kept_atoms)[i]) ++count;
  std::cout << "\nafter ring check\n";
  std::cout << "total num atoms in lig: " << kept_atoms->size() << "\n";
  std::cout << "num atoms in lig frag: " << count << "\n";
}

mol2File::mol2File(const mol2File& src, const my_float_t* rotmat, 
                   const my_float_t *tvec, const pts_vci pts_begin, 
                   const pts_vci pts_end, my_float_t radius)
{
  _fail = true;
  num_atoms = 0;
  orig_positions = NULL;
  positions = NULL;
  atom_ids_are_consecutive = true;
  rings = 0;
  if(src.fail()) return;

  // Get a local copy of the current atom positions from src and transform the 
  // positions using the given transformation.  
  my_float_t* trans_src_pos = new my_float_t[3 * src.num_atoms];
  memcpy(trans_src_pos, src.positions, 3 * src.num_atoms * my_float_size);
  move_positions(src.num_atoms, trans_src_pos, src.positions, rotmat, tvec);
  // Mark all those atoms which fall inside the region(s) (volume)
  std::vector<bool> kept_atoms(src.num_atoms, false);
  mark_kept_heavy_atoms(src.molecule, trans_src_pos, pts_begin, pts_end, 
                        radius, &kept_atoms);
  // Find and adjust for all the non-cleavable bonds which are set to be
  // cleaved.
  setup_stuff();
  bond_marks kept_bonds;
  // poor hack to keep the copy from replicating the ring information of the
  // source 
  rings = src.rings;
  cleave_bonds(src.molecule, src.bonds, &kept_atoms, &kept_bonds);
  rings = 0;

  ////////
  // Copy kept atoms from src.molecule to this.molecule
  ////////
  uint cnt = 0;
  for(std::vector<bool>::iterator i = kept_atoms.begin(); i < kept_atoms.end();
      ++i)  if(*i) cnt++;
  positions = new my_float_t[3 * cnt];
  std::vector<int> atom_nums(src.molecule.size(), -1);
  int atom_num = -1;
  cnt = 0;
  for(mol2_atom_vci atom = src.molecule.begin(); atom < src.molecule.end(); 
      ++atom, ++cnt)
    if(kept_atoms[cnt]){
      atom_nums[cnt] = ++atom_num;
      molecule.push_back(*atom);   
      memcpy(&positions[3*atom_num], &trans_src_pos[3*cnt], 3*my_float_size);  
      molecule.rbegin()->pos = &positions[3*atom_num];
      // Clear the neighbors lists as they will be created later
      molecule.rbegin()->neighbors.clear();
    }
  orig_positions = positions;
    
  ////////
  // Create the bonds and populate neighbor lists.
  ////////
  const mol2_atom_vci first_src_atom = src.molecule.begin();
  const mol2_atom_vi first_dest_atom = molecule.begin();
  bond_marks::const_iterator orig_pair = kept_bonds.begin();
  for( ; orig_pair != kept_bonds.end(); ++orig_pair){
    std::map<mol2_atom_vci, bool>::const_iterator targ_pair = 
      orig_pair->second.begin();
    for( ; targ_pair != orig_pair->second.end(); ++targ_pair){
      if(targ_pair->second == false){
        std::cout << "found a false one" << std::endl;
	continue;
      }

      mol2_atom_vi a = first_dest_atom;
      mol2_atom_vi b = first_dest_atom;
      a += atom_nums[orig_pair->first - first_src_atom];
      b += atom_nums[targ_pair->first - first_src_atom];
      if(b < a) iter_swap(a,b);

      // stupid const causes [] operator to fail
      // bonds[a][b] = src.bonds[orig_pair->first][targ_pair->first];
      bond_map::const_iterator ii; 
      for(ii = src.bonds.begin(); ii != src.bonds.end(); ++ii)
        if(ii->first->name == a->name){
	  //std::cout << ii->second.size() << std::endl;
          std::map<mol2_atom_vci, bond_t>::const_iterator jj;
          for(jj = ii->second.begin(); jj != ii->second.end(); ++jj)
            if(jj->first->name == b->name){
              bonds[a][b] = jj->second;
              break;
	    }
          break;
        }

      a->neighbors.push_back(b);
      b->neighbors.push_back(a);
    }
  }

  if(trans_src_pos) delete[] trans_src_pos;
  _fail = false;
}

bool 
mol2File::read_molecule_section(std::ifstream &mol2_file, 
                                std::vector<mol2_atom_t>* atoms, uint& num_atoms_out,
                                uint& num_bonds)
{
  // Ignore mol_name, mol_type, charge_type fields for now
  num_bonds = 0;
  num_atoms_out = 0;
  std::string line;
  if(std::getline(mol2_file, line) && line.substr(0,9) != "@<TRIPOS>"){
    // Get number of atoms and bonds
    mol2_file >> num_atoms_out >> num_bonds;  
    if(!mol2_file){
      std::cerr << "Error: Unable to read the molecule section "
                << "(mol2File::read_molecule_section)\n";
      return false;
    }
    atoms->clear();
    if(atoms->max_size() < num_atoms_out){
      std::cerr << "Error: Unable to allocate enough memory for vector "
                << "(mol2File::read_molecule_section)\n";
      return false;
    }
    atoms->reserve(num_atoms_out);
  }
  return true;
}
 
bool 
mol2File::read_atom_section(std::ifstream &mol2_file, 
                            mol2_atom_vec* atoms, my_float_t* atom_pos,
                            uint num_atoms_in)
{
  std::string line;
  my_float_t* cpos = atom_pos;
  for(uint i = 0; i < num_atoms_in && std::getline(mol2_file, line); 
      ++i, cpos += 3){
    if(line.substr(0,9) == "@<TRIPOS>"){
      std::cerr << "Expected atom section to have " << num_atoms_in 
                << " lines.\n" << "Only " << i << " lines were found.\n";
      return false;
    }

    mol2_atom_t atom;
    std::istringstream lstr(line);
    lstr >> atom.atom_num >> atom.atom_label >> cpos[0] >> cpos[1] >> cpos[2];
    atom.pos = cpos;
    std::string type;
    lstr >> type;
    atom.rest_of_line = line.substr(lstr.tellg());
 
    std::string::size_type pos = type.find(".");
    if(pos == std::string::npos){
      atom.orbit = DEFAULT_ORBIT;
      pos = type.length();
    }else atom.orbit = str_to_orbit_type(type.substr(pos));
    atom.name = str_to_atom_type(type.substr(0, pos));    

    // Check for premature eof or other failure
    if(!mol2_file){
      std::cerr << "Error: bad atom line in mol2 file "
                << "(mol2File::read_atom_section)\n";
      return false;
    }
    atoms->push_back(atom);    
   
    if(atom.atom_num - 1 != i) atom_ids_are_consecutive = false;
  }

  return true;
}

bool 
mol2File::read_bond_section(std::ifstream &mol2_file, mol2_atom_vec* atoms, 
                            bond_map* bonds_map, uint num_bonds)
{
  // If we need to handle this do so when we encounter such a situation. 
  if(!atom_ids_are_consecutive) return false;

  std::string line;
  for(uint i = 0; i < num_bonds && std::getline(mol2_file, line); ++i){
    if(line.substr(0,9) == "@<TRIPOS>"){
      std::cerr << "Expected bond section to have " << num_bonds << " lines.\n"
                << "Only " << i << " lines were found.\n";
      return false;
    }
   
    bond_t bond;
    uint orig_atom_id;
    uint targ_atom_id;
    std::string bond_str;
    std::istringstream lstr(line);
    lstr >> bond.id >> orig_atom_id >> targ_atom_id >> bond_str;

    // Check for premature eof or other failure
    if(!mol2_file){
      std::cerr << "Error: bad atom line in mol2 file "
                << "(mol2File::read_bond_section)\n";
      return false;
    }

    // Assume consecutive ids for bonds and atoms
    bond.status_bits = line.substr(lstr.tellg());
    bond.type = str_to_bond_type(bond_str);
    mol2_atom_vi orig_atom = atoms->begin() + orig_atom_id - 1;
    mol2_atom_vi targ_atom = atoms->begin() + targ_atom_id - 1;
    if(orig_atom < targ_atom) (*bonds_map)[orig_atom][targ_atom] = bond;
    else (*bonds_map)[targ_atom][orig_atom] = bond;
  
    // Append to the neighbor vectors
    orig_atom->neighbors.push_back(targ_atom);
    targ_atom->neighbors.push_back(orig_atom);
  }

  return true;
}

bool
mol2File::write(std::string fname)
{
  std::ofstream mol2_file;
  if(!open_ofstream(mol2_file, fname)) return false;

  mol2_file << "# Created by SimSite3D::site_dock\n"
            << "# Date: " << Timer::get_local_time() << "\n\n";

  // Count total number of bonds;
  uint num_bonds = 0;
  for(bond_map::iterator p1 = bonds.begin(); p1 != bonds.end(); ++p1)
    num_bonds += p1->second.size();

  mol2_file << "@<TRIPOS>MOLECULE\n"
            << "*****\n"
            << std::setw(5) << molecule.size() << " "
            << std::setw(5) << num_bonds << "      0      0     0\n"
            << "SMALL\n"
            << "USER_CHARGES\n"
            << "\n";

  mol2_file << "@<TRIPOS>ATOM\n";
  mol2_file << std::fixed << std::setprecision(4);
  int cnt = 1;
  for(mol2_atom_vi atom = molecule.begin(); atom < molecule.end(); ++atom, ++cnt){
    mol2_file << std::setw(6) << cnt << " "
              << std::setw(6) << std::left << atom->atom_label;
    std::right(mol2_file);
    for(uint i = 0; i < 3; ++i) mol2_file << std::setw(10) << atom->pos[i];
    mol2_file << std::setw(3) << atom_type_to_string(atom->name)
              << orbit_type_to_string[atom->orbit]
	      << atom->rest_of_line << "\n";
  }

  mol2_file << "@<TRIPOS>BOND\n";
  cnt = 1;
  for(bond_map::const_iterator p1 = bonds.begin(); p1 != bonds.end(); ++p1){
    std::map<mol2_atom_vci, bond_t>::const_iterator p2;// = p1->second.begin();
    for(p2 = p1->second.begin(); p2 != p1->second.end(); ++p2, ++cnt){
      mol2_file << std::setw(6) << cnt
                << std::setw(6) << p1->first - molecule.begin() + 1
                << std::setw(6) << p2->first - molecule.begin() + 1
	        << " " << bond_type_to_string[p2->second.type] << " " 
                << p2->second.status_bits << "\n";
    }
  }
  return true;
}

void
mol2File::enumerate_rings()
{
  // Build an adjacency list for the bonds as if the molecule was a directed
  // graph (digraph) with each bond representing an edge in both directions
  std::vector<std::vector<mol2_atom_vci> > edges(molecule.size());
  mol2_atom_vci mol_beg = molecule.begin();
  for(bond_map::const_iterator p1 = bonds.begin(); p1 != bonds.end(); ++p1){
    std::map<mol2_atom_vci, bond_t>::const_iterator p2;
    for(p2 = p1->second.begin(); p2 != p1->second.end(); ++p2){
      edges[p1->first - mol_beg].push_back(p2->first);
      edges[p2->first - mol_beg].push_back(p1->first);
    }
  }

  rings = new simple_graph<mol2_atom_t>(molecule, edges);
}

bool
mol2File::bond_is_rotatable(mol2_atom_vci orig, mol2_atom_vci targ, 
                            const bond_map& full_bonds_set)
{
  if(targ < orig){
    mol2_atom_vci tmp = targ;
    targ = orig;
    orig = tmp;
  }

  // Assume only bonds listed as single (1) are rotatable
  //if(full_bonds_set[orig][targ].type != SINGLE_BOND) return false;
  bond_map::const_iterator i; 
  for(i = full_bonds_set.begin(); i != full_bonds_set.end(); ++i)
    if(i->first == orig){
      std::map<mol2_atom_vci, bond_t>::const_iterator j;
      for(j = i->second.begin(); j != i->second.end(); ++j)
        if(j->first == targ){
	  if(j->second.type != SINGLE_BOND) return false;
          break;
	}
      break;
    }
  
  // We do not want to chop off any hydrogens if we keep the heavy atom to
  // which they are attached.
  else if(orig->name == H || targ->name == H) return false; 

  std::vector<orbit_type>::iterator my_hybrid = 
    stuff[orig->orbit].begin();
  for( ; my_hybrid < stuff[orig->orbit].end(); ++my_hybrid)
    if(*my_hybrid == targ->orbit){
      // Need to check for the trigonal_planar case here

      if(can_flex(*(targ)) || can_flex(*(orig))) return true;
    }
   
  return false;
}

bool
mol2File::can_flex(const mol2_atom_t& atom)
{
  if(atom.orbit == SP1 || atom.orbit == AR || atom.orbit == AM) 
    return true;
  if(atom.orbit == SP2){
    if(too_many_H_neighbors(atom, 1)) return false;
    else return true;
  }else if(atom.orbit == SP3){ 
    if(too_many_H_neighbors(atom, 2)) return false;
    else return true;
  }else if(atom.orbit == SP4 && atom.name == N){ 
    if(too_many_H_neighbors(atom, 2)) return false;
    else return true;
  }

  return false;
}

bool 
mol2File::too_many_H_neighbors(const mol2_atom_t& atom, uint max)
{
  // Apparently all isotopes of H are labeled as H in .mol2 ?
  uint cnt = 0;
  std::vector<mol2_atom_vi>::const_iterator neighbor;
  for(neighbor = atom.neighbors.begin(); neighbor < atom.neighbors.end();
      ++neighbor)
    if((*neighbor)->name == H) ++cnt;

  if(cnt > max) return true;
  return false;
}

  */
