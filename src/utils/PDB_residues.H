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
 * $Source: /psa/share/repository/pfizer_proj/src/basics/PDB_residues.H,v $
 * $Revision: 1.6 $
 * $Author: vanvoor4 $
 * $Date: 2008-07-28 15:21:16 $
 * 
 * $Log: not supported by cvs2svn $
 * Revision 1.5  2008/04/17 18:41:34  vanvoor4
 * Added another residue to the enum to distinguish between
 * errors/unknown residues and the NULL residue.
 *
 * Revision 1.4  2008/02/26 19:07:58  vanvoor4
 * Reworked so that the PDB specific items are used and included only
 * by the PDB classes and the mol2 specific items are used and
 * included only by the mol2 classes.
 *
 * Revision 1.3  2007/10/11 16:05:13  vanvoor4
 * Added SimSite3D namespace
 *
 * Revision 1.2  2007/08/24 15:26:49  vanvoor4
 * Updated to reflect the vdw radii in Tripos FF and Li & Nussinov (1998)
 *
 * Revision 1.1  2007/08/21 15:53:56  vanvoor4
 * Initial checkin
 *
 *
 * 
 * 
 */

#ifndef PDB_RESIDUES_HEADER_INCLUDED
#define PDB_RESIDUES_HEADER_INCLUDED

#include <AtomTypes.H>
#include <types.H>
#include <point.H>

namespace SimSite3D{

  //! 20 residues + nonstandard PDB residues
  typedef enum residue_type_enum{
    NULL_RESIDUE = 0,
    ALA,  //!< Alanine
    ARG,  //!< Arginine
    ASN,  //!< Asparagine
    ASP,  //!< Aspartic acid
    CYS,  //!< Cysteine
    GLN,  //!< Glutamine
    GLU,  //!< Glutamic acid
    GLY,  //!< Glycine
    HIS,  //!< Histidine
    ILE,  //!< Isoleucine
    LEU,  //!< Leucine
    LYS,  //!< Lysine
    MET,  //!< Methionine
    PCA,  //!< pyrrolidone carboxylic acid
    PHE,  //!< Phenylalanine
    PRO,  //!< Proline
    SER,  //!< Serine
    THR,  //!< Threonine
    TRP,  //!< Tryptophan
    TYR,  //!< Tyrosine
    VAL,  //!< Valine
    HOH,  //!< Water 
    TPO,  //!< Phosphorylated threonine
    PTR,  //!< Phosphorylated tyrosine
    ACE,
    MAIN_CHAIN,        //!< A wild card used by SimSite3D::HbondPoints
    PDB_METAL,         //!< Designates a metal atom
    UNKNOWN_RESIDUE,
    RESIDUE_ENUM_END   //!< End of enum marker
  }residue_type;

  //! Silly little struct to help convert from string to residue name
  typedef struct{
    residue_type residue;
    std::string name;
  }residue_conv_type;

  //! Unlike hbonding potential, etc the hydrophobicity is defined differently
  //! for the main chain atoms of each residue type
  typedef struct{
    residue_type residue;
    atom_type atom;
    int hydrophobicity;  //!< This is the hydro used by SLIDE that is already adjusted downward by 235
  }pdb_hphob_type;
  typedef std::map<atom_type, int> hphob_tbl_atom_lvl;
  typedef std::map<residue_type, hphob_tbl_atom_lvl> hphob_tbl_res_lvl;

  typedef enum{
    WARN_UNKNOWN_RESIDUES,
    SILENTLY_IGNORE_UNKNOWN_RESIDUES
  }residue_warn_type;

  //! Entry for a PDB atom for a known residue -- at this time the 20 
  //! amino acids and PCA are the only recognized residues.
  typedef struct{
    std::string pdb_atom_str;   //!< 4 char PDB atom name
    std::string pdb_res_str;    //!< 3 char PDB residue name
    atom_type atom;             //!< SimSite3D atom identifier (enum)
    residue_type res;           //!< SimSite3D atom identifier (enum)
    my_float_t vdw_radius;      //!< Van der Waals radius
    my_float_t ms_atm_radius;   //!< MS, AMS, and MSMS unified atomic radius
    orbit_type orbit;           //!< SimSite3D atom orbital identifier (enum)
    atom_level_type level;      //!< Atom level -- not used at present
    interactionType act_type;   //!< Hbond donor/acceptor/doneptor/other
    my_float_t charge;          //!< Formal charge
    my_float_t hydro;
  }pdb_atom_info_t;
  typedef std::map<std::string, const pdb_atom_info_t*> str_tbl_res_lvl;
  typedef std::map<std::string, str_tbl_res_lvl> atom_table_t;

  typedef std::map<residue_type, const pdb_atom_info_t*> enum_tbl_res_lvl;
  typedef std::map<atom_type, enum_tbl_res_lvl> atom_enum_table_t;
  
  //! Simple class to facilitate the lookup of protein residue atoms using the
  //! PDB 3 char residue name and the PDB 4 char atom name
  /*! PDB file van der Waals Radii Definitions and atom oribital types
   *
   * Taken from the Tripos force field [Clark et al 1989], with the changes 
   * listed in Table XI: The Minimum Radii of Atoms, of Li & Nussinov (1998)
   * "A Set of van der Waals and Coulombic Radii of Protein Atoms for Molecular
   * and Solvent-Accessible Surface Calculation, Packing Evaluation, and 
   * Docking", Proteins 32, 111-127.  
   * 
   * This set of radii is used since it is a set common to the Kuhn lab software
   * and was/is used by DrugScore -- Gohlke, Hendlich and Klebe (2000) 
   * "Knowledge-based scoring function to predict protein-ligand interactions",
   * J. Mol. Biol., 295, 337-356.
   */
  class PDB_residues{
  public:
    //! Get the atom info based on residue and atom names 
    static const pdb_atom_info_t*
    get_atom_info(const std::string atom_str, const std::string res_str)
    {
      if(!residue_table.size()) PDB_residues::build_res_tables();
      atom_table_t::const_iterator a_iter = residue_table.find(atom_str);
      if(a_iter == residue_table.end()) return 0;
      str_tbl_res_lvl::const_iterator res_iter = a_iter->second.find(res_str);
      if(res_iter == a_iter->second.end()) return 0;
      return res_iter->second;
    }

    //! Get the atom info based on residue and atom types 
    static const pdb_atom_info_t*
    get_atom_info(const atom_type atom, residue_type residue)
    {
      if(!residue_type_table.size()) PDB_residues::build_res_tables();
      atom_enum_table_t::const_iterator a_iter = residue_type_table.find(atom);
      if(a_iter == residue_type_table.end()) return 0;
      enum_tbl_res_lvl::const_iterator res_iter = a_iter->second.find(residue);
      if(res_iter == a_iter->second.end()) return 0;
      return res_iter->second;
    }

    static int
    get_atom_hydrophobicity(const residue_type res, const atom_type atom)
    {
      if(A_hydrophobic_value.size() == 0) build_atom_hphob_map();

      // Default value for oxygen is 530 - 235
      if(atom == OXT) return 295; 
      hphob_tbl_res_lvl::const_iterator res_iter;
      res_iter = A_hydrophobic_value.find(res);
      if(res_iter == A_hydrophobic_value.end()) 
        return std::numeric_limits<int>::max();
      hphob_tbl_atom_lvl::const_iterator a_iter;
      a_iter = res_iter->second.find(atom);
      if(a_iter == res_iter->second.end())
        return std::numeric_limits<int>::max();
      return a_iter->second;
    }

    //! Convert from 3 char PDB residue string to SimSite3D residue identifier
    static residue_type string_to_residue(const std::string res_name);
 
    //! Convert from SimSite3D residue identifier to 3 char PDB residue string
    static std::string residue_to_string(const residue_type res);

    //! Convert from SimSite3D PDB ATOM atom identifier to 4 char PDB atom string
    static std::string atom_to_string(const atom_type atom);

  private:
    PDB_residues()
    { ; }

    PDB_residues(const PDB_residues &)
    { ; }

    //! Populate the 2D map holding the atoms keyed by atom and residue id
    static void build_res_tables();

    //! Populate the maps to convert between 3 char PDB string and residue id
    static void build_res_conv_maps();

    //! Populate the map to convert between atom id and 4 char PDB atom string
    static void build_atom_conv_map();
  
    static void build_atom_hphob_map();

    static atom_table_t residue_table;   //!< Array holding the PDB atom info for standard residues
    static atom_enum_table_t residue_type_table;
    static std::map<residue_type, std::string> A_residue_to_string;
    static std::map<std::string, residue_type> A_string_to_residue;
    static std::map<atom_type, std::string> A_atom_to_string;
    static hphob_tbl_res_lvl A_hydrophobic_value;

    static const pdb_atom_info_t Pro_main_chain_N;
    static const uint A_res_conv_array_size;
    static const residue_conv_type A_res_conv_array[]; 
    static const uint A_res_array_size;   //!< Number of els in residue array
    static const pdb_atom_info_t A_res_array[]; //!< Array holding residue info
    static const uint A_atom_conv_array_size;
    static const atom_conv_type A_atom_conv_array[];
    static const size_t A_hphob_val_array_size;
    static const pdb_hphob_type A_hphob_val_array[];
    static const std::string A_fname;   //!< Name of the source file
  };
}
#endif
