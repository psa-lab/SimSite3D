
#ifndef PROTEIN_LIGAND_SCORE_HEADER_INCLUDED
#define PROTEIN_LIGAND_SCORE_HEADER_INCLUDED

#include <PDBStructure.H>
#include <mol2File.H>
#include <HbondGeometry.H>

namespace ASCbase{

//! Compute orientation and affinity scores between protein atoms and ligand
//! atoms.  Note: because the protein and ligand are required for scoring but
//! also used by other functions, this class is really only valid if both the
//! ligand and protein objects are still in scope (exist).
class ProtLigScore{
public:

  typedef enum{
    INITIAL = 0,
    INTERFACIAL_ATOM,
    METAL_DIRECT_HBOND,
    SALT_BRIDGE,
    DIRECT_HBOND,
    UNSAT_CHARGE,
    UNSAT_POLAR,
    INTRA_LIGAND_HBOND,
    INTRA_LIGAND_SALT_BRIDGE,
    INTRA_TARGET_HBOND,
    INTRA_TARGET_SALT_BRIDGE
  }atom_act_class_t;
  typedef std::pair<atom_vci, atom_act_class_t> atom_act_pair;
  typedef std::map<atom_vci, atom_act_class_t> atom_act_map;
  typedef atom_act_map::iterator atom_act_mi;
  typedef atom_act_map::const_iterator atom_act_mci;

  typedef struct{
    atom_vci A;
    atom_vci B;
    my_float_t sq_dist;
  }salt_bridge_t;

  typedef struct{
    atom_vci A;
    atom_vci B;
    my_float_t sq_dist;
    my_float_t cos_DHA_angle;
    my_float_t cos_preacc_angle;
  }hbond_t;

  typedef struct{
    int num_metal_hbonds; // needs to be checked
    int num_salt_bridges; // metal sb need checking 
    int num_hbonds;        
    int num_unsat_polar;
    int num_unsat_charge;
    my_float_t total_target_hydro;
    int contact_hphob_hphob;            
    my_float_t total_hphob_hphob_comp;  
  }raw_terms_type;

  ProtLigScore(PDBStructure &prot, mol2File &lig, 
               bool save_prot_lig_acts = true, 
               bool compute_charge_sums = true);

  ~ProtLigScore();

  const my_float_t
  affiscore() const
  { return A_affiscore; }

  const my_float_t
  orientscore() const
  { return A_orientscore; }

  const my_float_t
  ligand_efficiency() const
  { return A_affiefficientscore; }

  void
  print_hbonds(std::ostream &out);

  static void
  score_header_lines(std::ostream &out, bool classic = false);

  void
  classic_report(std::ostream &out);

  void
  report_scores(std::ostream &out, char delim = '|');

  //! Generate the ligand atoms' interaction strings
  /*! This function is public to support the ASCbasePy interface
   */
  void
  gen_lig_act_strings(const char field_delim = '|', 
                      const char prot_delim = ',');

  //! Write the ligand atoms with corresponding interactions from the protein
  void
  write_lig_acts(std::ostream &out = std::cout, char field_delim = '|',
                 char prot_delim = ',');

  std::vector<salt_bridge_t>::const_iterator
  salt_bridges_begin() const
  { return A_prot_lig_SB.begin(); } 

  std::vector<salt_bridge_t>::const_iterator
  salt_bridges_end() const
  { return A_prot_lig_SB.end(); } 

  std::vector<hbond_t>::const_iterator
  hbonds_begin() const
  { return A_prot_lig_hbonds.begin(); }

  std::vector<hbond_t>::const_iterator
  hbonds_end() const
  { return A_prot_lig_hbonds.end(); }

  std::map<residue_vci, bool>::const_iterator
  hphob_SC_begin() const
  { return A_hphob_prot_SC.begin(); }

  std::map<residue_vci, bool>::const_iterator
  hphob_SC_end() const
  { return A_hphob_prot_SC.end(); }

  const std::string
  lig_name() const
  {
    return A_lig->name();
  }

  //! Added for ASCbasePy interface -- why compute it in python when it is 
  //! already computed and debugged in C++
  std::vector<std::string>::const_iterator
  lig_act_strings_begin()
  { 
    if(!A_lig_act_strings.size()) gen_lig_act_strings();
    return A_lig_act_strings.begin(); 
  }  

  //! Added for ASCbasePy interface -- why compute it in python when it is 
  //! already computed and debugged in C++
  std::vector<std::string>::const_iterator
  lig_act_strings_end()
  { 
    // Need to have the check here as well as *_begin() since someone could
    // call *_end() first
    if(!A_lig_act_strings.size()) gen_lig_act_strings();
    return A_lig_act_strings.end(); 
  }  

private:
  void
  compute_scores();

  void
  compute_terms();

  void
  metal_interactions(atom_vci lig_atom, int *num_lig_nbrs);

  void
  prot_lig_atom_interactions(atom_vci lig_atom, int *num_lig_nbrs,
                             int *sum_hphob_hphob,
                             int *num_hphob_hphob_contact_1_lig_atm);
  void
  polar_interaction(atom_vci prot_atom, atom_vci lig_atom, my_float_t sq_dist);

  void
  flag_intra_ligand_hbonds();

  void
  flag_intra_protein_hbonds();

  void
  sum_target_hydrophobicity();

  void
  flag_unsat_polar_atoms(atom_act_map *my_map, int *num_unsat_charge,
                         int *num_unsat_polar);

  void
  print_atom_act_classes(atom_act_map *my_map);

  PDBStructure *A_prot; //!< Pointer to target protein data structures -- only valid while constructor is active
  mol2File *A_lig;      //!< Pointer to ligand to score -- only valid while constructor is active

  HbondGeometry A_hbond_geometry;  //!< Instance of hbond geometry class
  std::map<atom_vci, atom_act_class_t> A_ligand_flag;
  std::map<atom_vci, atom_act_class_t> A_target_flag;
  std::vector<salt_bridge_t> A_prot_lig_SB;  //!< Protein-ligand salt bridges
  std::vector<hbond_t> A_prot_lig_hbonds;        //!< Protein-ligand hbonds
  std::map<residue_vci, bool> A_hphob_prot_SC;  //!< Keep note of hphob sidechains interating with ligand hphob atoms
  raw_terms_type A_raw_terms;
  int A_num_lig_carbons;         //!< Count of carbons in the ligand
  int A_num_exposed_lig_carbons; //!< Count of ligand carbon atoms with no protein neighbors (within 4.5 (A))
  std::vector<my_float_t> A_affi_terms;
  my_float_t A_affiscore;
  my_float_t A_orientscore;
  my_float_t A_affiefficientscore;
  bool A_keep_prot_lig_acts;
  std::vector<std::vector<residue_vci> > A_hphob_prot_lig_SC;
  std::vector<bool> A_lig_exposed_Cz;
  std::vector<std::string> A_lig_act_strings; //!< A string for each ligand atom that denoting interactions with the protein 

  static const my_float_t MIN_BURY_SQUARED_DIST;
  static const my_float_t MAX_BURY_SQUARED_DIST;
  static const my_float_t HYDRO_SQUARED_DIST;
  static const my_float_t MAX_INTERACT_SQUARED_DIST;
  static const my_float_t AFFI_64_WEIGHTS[];
  static const my_float_t ORIENT_99_WEIGHTS[];
};
} // end of namespace ASCbase
#endif
