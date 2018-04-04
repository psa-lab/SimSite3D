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

#ifndef PDB_STRUCTURE_HEADER_INCLUDED
#define PDB_STRUCTURE_HEADER_INCLUDED
#include <PDBBase.H>

namespace SimSite3D{

//! Simple package to hold information about a PDB residue
class residue_t{
public:
  residue_t()
  {
    name = NULL_RESIDUE;
    number = 0; 
    chainID = '@';
    icode = '@';
    atoms_begin = atom_t::NULL_ATOM_VCI;
    atoms_end = atom_t::NULL_ATOM_VCI;
    is_MODRES = false;
    is_HET = false;
  }

  //! Initialize residue based on the given atom
  residue_t(atom_vci atom)
  {
    name = atom->res;
    number = atom->res_num;
    chainID = atom->chainID;
    icode = atom->iCode;
    atoms_begin = atom;
    atoms_end = atom;
    is_MODRES = false;
    is_HET = atom->is_hetero;
  }

  ~residue_t() {;}

  atom_vci 
  get_atom(atom_type atom_name) const
  {
    for(atom_vci AA = atoms_begin; AA < atoms_end; ++AA)
      if(AA->name == atom_name) return AA;
    return atom_t::NULL_ATOM_VCI;
  }

  residue_type name;     //!< 3 character name of the residue
  uint number;           //!< Number assigned to residue in the PDB file 
  char icode;            //!< Insertion code assigned to the residue
  char chainID;          //!< Chain ID assigned to the residue
  atom_vci atoms_begin;  //!< Iterator to first atom in the residue
  atom_vci atoms_end;    //!< Iterator to 1 past the last atom in the residue
  bool is_MODRES;        //!< True implies this residue corresponds to a MODRES
  bool is_HET;           //!< True implies this residue corresponds to a HET record
};

typedef residue_t* residue_ptr;
typedef std::vector<residue_t> residue_vec;
typedef std::vector<residue_t>::iterator residue_vi;
typedef std::vector<residue_t>::const_iterator residue_vci;

typedef struct {
  residue_vci residues_begin; //!< Iter to first residue in the chain
  residue_vci residues_end;   //!< Iter to 1 past the last residue
  char chainID;                      //!< Character assigned to the chain
}chain_t;

typedef chain_t* chain_ptr;
typedef std::vector<chain_t> chain_vec;
typedef std::vector<chain_t>::iterator chain_iter;
typedef std::vector<chain_t>::const_iterator chain_const_iter;

//! Basic class to impose the chain/residue/atom hierarchy of a PDB structure
//! file.  Does not contain a default constructor as we currently have no
//! need to write structure files
class PDBStructure : public PDBBase{
public:
  //! Read the heavy atoms stored in the PDB file given by filename
  /*! Uses PDBBase to read the heavy atoms and then imposes/creates the
   * residues and chains from the heavy atoms (in CoordFile).
   *
   * @param filenam Path to the PDB structure file to read
   * @param verbosity whether or not to write some stuff to stdout
   */
  PDBStructure(const std::string filename, const bool ignore_altLocs=false, 
               const verbose_level_t verbosity = VERBOSE_SILENT);

  //! Does nothing 
  ~PDBStructure();

#if 0
  bool
  write_xyzr(std::ostream &out, bool include_metals, 
             const std::vector<std::string> &included_waters, 
             const atom_vci close_atom, uint* close_atom_idx);
#endif

  void report_stats(std::ostream& out, const verbose_level_t verbosity);
  void report_hetatms(std::ostream& out, std::string filename);

  chain_const_iter chains_begin() const { return chains.begin(); }
  chain_const_iter chains_end() const { return chains.end(); }

  residue_vci residues_begin() const { return residues.begin(); }
  residue_vci residues_end() const { return residues.end(); }

  //! Get an iterator to the oxygen of the first water molecule
  std::vector<atom_vci>::const_iterator
  waters_beg() const
  { return A_waters.begin(); }

  //! Get an iterator to 1 past the last HOH
  std::vector<atom_vci>::const_iterator
  waters_end() const
  { return A_waters.end(); }

  //! Get an iterator to the first metal atom
  std::vector<atom_vci>::const_iterator
  metals_begin() const
  { return metals.begin(); }

  //! Get an iterator to 1 past the last metal atom
  std::vector<atom_vci>::const_iterator
  metals_end() const
  { return metals.end(); }

  residue_vci 
  get_residue(atom_vci atom) const;

  static const residue_t NULL_RESIDUE;  //!< the "standard" null residue
  static const residue_vec NULL_RESIDUE_VECTOR; //!< The "standard" null vector of residues -- for the null_residue_vci
  static const residue_vci NULL_RESIDUE_VCI;  //!< the "standard" null residue vector constant iterator
private:

  bool 
  read_file(std::string filename, const bool ignore_altLocs);

  //! Create chains and residues from the atom information
  void 
  form_chains_and_residues();

  bool
  add_residue(atom_vci initial_atom, residue_vi* current_residue);

  residue_vec residues;   //!< All residues in the PDB file 
  residue_vec A_hetgroups;  //!< HET entries that are not also MODRES entries
  chain_vec chains;       //!< All chains in PDB file
  std::vector<atom_vci> A_waters;  //!< All the waters in PDB file
  std::vector<atom_vci> metals;  //!< All the metals in PDB file

  static const std::string _fname;   //!< Name of the source file
};

} // End of namespace SimSite3D

#endif
