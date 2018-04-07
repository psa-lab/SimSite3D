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
#ifndef _MOL2FILE_HEADER_FILE
#define _MOL2FILE_HEADER_FILE

#include <map>
#include <iostream>
#include <simple_graph.H>
#include <mol2_atoms.H>
#include <CoordFile.H>
#include <BoundingVolume.H>

namespace SimSite3D{

  //!< Fields of a mol2 bond record
  typedef struct{
    bond_type type;   //!< Type of bond (single, double, etc.)
    int atom_num1;    //!< First atom in the bond (as listed in the file)
    int atom_num2;    //!< Second atom in the bond (as listed in the file)
    std::string status_bits; //!< Everything else on the line after the 2 atoms
  }mol2_bond_t;

  typedef std::vector<mol2_bond_t> 		mol2_bond_vec;
  typedef mol2_bond_vec::iterator		mol2_bond_vi;
  typedef mol2_bond_vec::const_iterator		mol2_bond_vci;

  //! Read and write mol2 files and perform various operations.
  /*!
   * This class has 2 constructors.  The first is to read the molecule(s) from
   * a given file.  The second constructor is a copy constructor and requires
   * a protein and sitemap.  This copy constructor discards all ligand
   * atoms (in the copy) that are not within 2.0 (A) of any template point.
   * These atoms are then checked to ensure that they do not significantly
   * overlap with the protein.  All of the ligand fragments are checked to 
   * ensure that they each have at least N connected heavy atoms.
   */
  class mol2File : public CoordFile{
  public:

    //! Constructor -- reads in the file and sets up the structure
    /*!
     * @param filename Path of the mol2 file to read
     */
    mol2File(const std::string filename,
             const verbose_level_t verbosity = VERBOSE_SILENT);

/*
    mol2File(const mol2File& src, const CoordFile& prot_atoms, 
             const my_float_t* R, const my_float_t* T);
*/

    //! Copy constructor -- fragments the ligand
    /*!
     * @param src Reference to the full ligand to copy (search ligand)
     * @param site_vol Reference to the sitemap volume estimate
     * @param prot_atoms Reference to the query protein atoms
     * @param min_num_atoms Minimum number of atoms in a ligand fragment
     * @param R Rotation matrix to move ligand to the model ref. frame
     * @param T Translation vector to move ligand to the model ref. frame
     */
//    mol2File(const mol2File& src, const CoordFile &site_pts,
    mol2File(const mol2File& src, const BoundingVolume &site_vol,
             const CoordFile& prot_atoms, const uint min_num_atoms,
             const my_float_t* R, const my_float_t* T);

    //! Do nothin dstor
    ~mol2File();

    void report_stats(std::ostream& out);

    //! Set the current positions to those read from file (original positions).
    void revert_positions();

    //! Write the molecule(s) to the specified file
    /*!
     * @param ofname Name of the file to write molecule(s)
     */
    bool 
    write(std::string ofname);

    //! Get the number of heavy atoms in the molecule(s)
    const uint 
    num_heavy_atoms();

    //! Get a const iterator to the beginning of the boolean vector denoting 
    //! which atoms of the full ligand were kept.
    /*! The vector referred to is of the same length as the number of atoms
     * in the "original" ligand.  The nth bool is true or false if the nth
     * ligand atom is in a ligand fragment or is not in a ligand fragment
     * respectively.
     */
    std::vector<bool>::const_iterator
    frag_atoms_flags_beg()
    { return frag_atoms.begin(); }

    //! Get a const iterator to the end of the boolean vector denoting 
    //! which atoms of the full ligand were kept.
    std::vector<bool>::const_iterator
    frag_atoms_flags_end()
    { return frag_atoms.end(); }

    //! Get the neighbors of atom in this molecule.  If someone wants to
    //! pass in an atom from another molecule, they will get rewarded with
    //! a const ref to an empty vector.
    const std::vector<atom_vci>&
    get_nbrs(const atom_vci atom) const 
    {
      if(atom < atoms_begin() || atom >= atoms_end()){
        warn(_fname, "get_nbrs", "Atom is not in the molecule");
        return A_nbrs_vec.back();
      }
      return *(A_nbrs_vec.begin() + (atom - atoms_begin()));
    }


    std::string
    atom_type_to_str(const atom_vci atom)
    { 
      typedef std::map<orbit_type, std::string> orbit_lvl;
      if(atom_types_to_str.size() == 0) init_atom_types_to_str();
      std::map<atom_type, orbit_lvl>::const_iterator atom_iter;
      atom_iter = atom_types_to_str.find(atom->name);
      if(atom_iter == atom_types_to_str.end()) return "Du.UNK";
      orbit_lvl::const_iterator orbit_iter;
      orbit_iter = atom_iter->second.find(atom->orbit);
      if(orbit_iter == atom_iter->second.end()) return "Du.UNK";
      return orbit_iter->second;
    }

    bool
    calc_charge_sums();

  private:
    //! Simple temporary type to store bond info in a vector
    typedef struct{
      atom_vci nbr;
      bond_type bond;
    }nbr_bond_type;

    //! Top level read function -- read in the molecule(s) from the file
    //! passed to the constructor
    bool read_data(std::ifstream &mol2_file);

    //! Read the first molecule section
    bool read_molecule_section(std::ifstream &mol2_file,
                               uint* num_atoms_out,
                               uint* num_bonds_out);

    //! Read the first atom section after the first molecule section
    bool read_atom_section(std::ifstream &mol2_file, const uint num_atoms);

    //! Read the first bond section after the first atom section
    bool read_bond_section(std::ifstream &mol2_file, const uint num_bonds);

    //! Given mol2 atom and orbital string, get internal info on atom
    bool look_up_atom_by_string(const std::string atom_str, atom_t* atom);
   
    //! Given mol2 bond string, get internal bond type 
    bool look_up_bond_by_string(const std::string bond_str, mol2_bond_t* bond);

    //! Setup the map used to look up atom info based on mol2 atom and orbital
    void init_atom_table_by_string();

    //! Setup the map to convert from mol2 bond string to internal bond type
    void init_bond_table_by_string();
 
    //! Setup the map to convert from interal bond type to mol2 bond string 
    void init_bond_types_to_str();
 
    //! Setup the map to convert from internal atom and orbital to mol2 atom 
    //! string
    void init_atom_types_to_str();

    void add_bond(const uint a1, const uint a2, const bond_type type,
                  const std::string status_bits);

    bool assign_act_types();

    void generate_nbrs_vec();

    std::string A_mol_name;
    std::string A_mol_type;
    std::string A_charge_type;

    std::map<std::string, const mol2_atom_info_t*> atoms_by_type_str;
    std::map<std::string, bond_type> bonds_by_type_str;
    std::map<bond_type, std::string> bond_types_to_str;
    std::map<atom_type, std::map<orbit_type, std::string> > atom_types_to_str;
    std::vector<std::string> rest_of_lines; //! Rest of the atom lines
    mol2_bond_vec bonds;
    bool atom_ids_are_consecutive;
    bool is_frag_file;     //!< True if the molecule is a fragment, else false
    uint a_num_heavy_atoms;
    std::vector<bool> frag_atoms; //!< Bool string of atoms present in fragment
    //std::multimap<atom_vci, nbr_bond_type> A_nbrs_map;
    std::vector<std::vector<atom_vci> > A_nbrs_vec;

    static const my_float_t ATOMS_TOO_FAR;
    static const my_float_t ATOMS_TOO_CLOSE;
    static const std::string _fname;  //!< Name of source file
    static const my_float_t MINIMAL_CHARGE;
  };
}
#endif
