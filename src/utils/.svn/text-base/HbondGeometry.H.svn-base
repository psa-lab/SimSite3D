
#ifndef HBOND_GEOMETRY_HEADER_FILE_INCLUDED
#define HBOND_GEOMETRY_HEADER_FILE_INCLUDED
#include <mol2File.H>

namespace SimSite3D{
//! A helper class to compute various geometrical features of hydrogen bonds
class HbondGeometry{
public:

  typedef enum{
    ACCEPTOR_ACCEPTOR = ACCEPTOR + ACCEPTOR,
    ACCEPTOR_DONOR = ACCEPTOR + DONOR,
    ACCEPTOR_DONEPTOR = ACCEPTOR + DONEPTOR,
    DONOR_DONOR = DONOR + DONOR,
    DONOR_DONEPTOR = DONOR + DONEPTOR,
    DONEPTOR_DONEPTOR = DONEPTOR + DONEPTOR,
    ACCEPTOR_METAL_1 = ACCEPTOR + METAL_1,
    ACCEPTOR_METAL_2 = ACCEPTOR + METAL_2,
    DONEPTOR_METAL_1 = DONEPTOR + METAL_1,
    DONEPTOR_METAL_2 = DONEPTOR + METAL_2
  }interactionSum;
  
  //! Used to define the hbond coordinate symmetry
  typedef enum{
    NO_FLIP = 0,
    FLIP_NBR_Y, 
    FLIP_ACC_Y,
    FLIP_Z
  }local_coord_flip;

  //! Simple struct to hold protein sidechain donor atom information
  typedef struct{
    residue_type residue;           //!< PDB name of the residue
    atom_type atoms[3];             //!< PDB names of 3 protein atoms to define local coordinate system
    const my_float_t* local_H_pos;  //!< Local H position in positive Y half plane
  }sidechain_donor_t;
  typedef std::map<atom_type, const sidechain_donor_t*> SC_tbl_atom_lvl;
  typedef std::map<residue_type, SC_tbl_atom_lvl> SC_tbl_res_lvl;

  //! Simple struct to hold protein sidechain acceptor atom information
  typedef struct{
    residue_type residue;
    atom_type acceptor;
    atom_type preacceptor;
  }sidechain_preacc_t;
  typedef std::map<atom_type, const sidechain_preacc_t*> preacc_tbl_atom_lvl;
  typedef std::map<residue_type, preacc_tbl_atom_lvl> preacc_tbl_res_lvl;

  //! do nothin constructor
  HbondGeometry()
  {
    A_prot = 0;
  }

  //! do nothin destructor
  ~HbondGeometry() {;}

  static inline bool
  salt_bridge(const atom_t &A, const atom_t &B, const my_float_t sq_dist)
  {
    // Due to the nature of His ND1, we do not know whether it is protonated.
    // However, it can only participate in a salt bridge if it has a proton,
    // and then as the positive component.
    interactionType A_act_type = A.act_type;
    if(A.res == HIS && A.name == ND1) A_act_type = DONOR;
    interactionType B_act_type = B.act_type;
    if(B.res == HIS && B.name == ND1) B_act_type = DONOR;
    // If the charges have opposite sign they will be negative when
    // multiplied.  Otherwise the product will be greater than or equal to 0.0
    if(MIN_SALT_BRIDGE_SQUARED_LENGTH <= sq_dist &&
       sq_dist <= MAX_SALT_BRIDGE_SQUARED_LENGTH &&
       is_hbond_interaction(A_act_type, B_act_type) && 
       A.charge * B.charge < 0.0)
      return true;
    return false;
  } 

  static inline bool
  is_hbond_interaction(interactionType A_act, interactionType B_act)
  {
    int act_sum = A_act + B_act;
    if(act_sum == ACCEPTOR_DONOR || act_sum == ACCEPTOR_DONEPTOR ||
       act_sum == DONOR_DONEPTOR || act_sum == DONEPTOR_DONEPTOR) return true;
    return false; 
  }

  //! For now the consensus is only acceptors should interact with metal ions
  static inline bool
  metal_hbond(const atom_t &polar_atom, const atom_t &metal, 
              const my_float_t sq_dist)
  {
    if((polar_atom.act_type == ACCEPTOR || polar_atom.act_type == DONEPTOR) &&
       ((metal.act_type == METAL_1 && 
         sq_dist <= PDB_metals::MAX_METAL_1_SQUARED_LENGTH) ||
        (metal.act_type == METAL_2 && 
         sq_dist <= PDB_metals::MAX_METAL_2_SQUARED_LENGTH)))
      return true;
    return false;
  }

  //! For now the consensus is only acceptors should interact with metal ions
  static inline bool
  metal_salt_bridge(const atom_t &polar_atom, const atom_t &metal, 
                    const my_float_t sq_dist)
  {
    // Due to the nature of His ND1, we do not know whether it is protonated.
    // However, it can only participate in a salt bridge if it has a proton,
    // and then as the positive component.
    if(polar_atom.res == HIS && polar_atom.name == ND1) return false;

    // the metal distance checks are not needed since we should always call
    // metal_hbond first and not call this function if metal_hbond is true,
    if(sq_dist <= MAX_SALT_BRIDGE_SQUARED_LENGTH && polar_atom.charge != 0 &&
       (polar_atom.act_type == ACCEPTOR || polar_atom.act_type == DONEPTOR) &&
       ((metal.act_type == METAL_1 && 
         sq_dist > PDB_metals::MAX_METAL_1_SQUARED_LENGTH) ||
        (metal.act_type == METAL_2 && 
         sq_dist > PDB_metals::MAX_METAL_2_SQUARED_LENGTH)))
      return true;
    return false;
  }

  bool
  intra_prot_hbond(atom_vci A, atom_vci D, my_float_t sq_dist,
                   PDBStructure &prot, my_float_t *cos_theta, 
                   my_float_t *cos_delta);

  bool 
  intra_lig_hbond(atom_vci A, atom_vci D, my_float_t sq_dist,
                  mol2File &ligand, my_float_t *cos_theta,
                  my_float_t *cos_delta);

  //! Atom type and distance chess as well as computing and checking the DHA 
  //! and H ... A-AA angles
  bool
  prot_lig_hbond(atom_vci prot_atom, atom_vci lig_atom, PDBStructure &prot, 
                 mol2File &ligand, my_float_t sq_dist, my_float_t *cos_theta, 
                 my_float_t *cos_delta);

  //! Compute the cosine of theta (DHA angle) when protein atom is the donor
  /*! Compute cosine of the donor-hydrogen-acceptor (DHA) angle.  This is
   * typically termed the theta angle.
   *  
   */
  bool
  protein_donor_cos_theta(atom_vci prot_atom, atom_vci acceptor, 
                          PDBStructure &prot, my_float_t *H_pos,
                          my_float_t *cos_theta);

  //! Compute the cosine of theta (DHA angle) when ligand atom is the donor
  /*! Compute cosine of the donor-hydrogen-acceptor (DHA) angle.  This is
   * typically termed the theta angle.
   *  
   */
  bool
  ligand_donor_cos_theta(atom_vci lig_atom, atom_vci acceptor, 
                         mol2File &ligand, my_float_t *H_pos, 
                         my_float_t *angle);

  //! Compute the cosine of delta (H ... A-AA angle) when the protein atom
  //! is the acceptor
  bool 
  prot_acceptor_cos_delta(atom_vci acceptor, PDBStructure &prot, 
                          my_float_t *H_pos, my_float_t *cos_delta);

  //! Compute the cosine of the delta angle (H ... A-AA angle) when the ligand
  //! atom is the acceptor
  bool
  lig_acceptor_cos_delta(atom_vci acceptor, mol2File &ligand,
                         my_float_t *H_pos, my_float_t *cos_delta);

  //! Compute the hydrogen position for a number of cases
  void
  compute_H_position(atom_vci hvy_nbr_A, atom_vci donor, atom_vci hvy_nbr_B, 
                     atom_vci acceptor, const my_float_t *local_H_pos, 
                     const local_coord_flip coord_flip, my_float_t *new_H_pos);

  //! Compute the local hydrogen position for a non sp3/sp4 ligand donor
  bool
  compute_local_H_pos(atom_vci hvy_nbr, atom_vci lig_atom, atom_vci H_atom, 
                      my_float_t *local_H_pos);

  //! Compute the cosine of the angle ABC
  my_float_t
  cosine_angle(my_float_t *A, my_float_t *B, my_float_t *C, 
               my_float_t *cos_angle);

  void
  set_protein(PDBStructure *prot) 
  {
    A_prot = prot; 
  }
 
  static const sidechain_donor_t*
  get_sidechain_donor(residue_type res, atom_type atom)
  {
    if(!A_SC_donors.size()) build_SC_donor_tables();
    SC_tbl_res_lvl::const_iterator res_iter = A_SC_donors.find(res);
    if(res_iter == A_SC_donors.end()) return 0;
    SC_tbl_atom_lvl::const_iterator atom_iter = res_iter->second.find(atom);
    if(atom_iter == res_iter->second.end()) return 0;
    return atom_iter->second;
  }

  static const sidechain_preacc_t*
  get_sidechain_preacc(residue_type res, atom_type atom)
  {
    if(!A_SC_preaccs.size()) build_SC_preacc_tables();
    preacc_tbl_res_lvl::const_iterator res_iter = A_SC_preaccs.find(res);
    if(res_iter == A_SC_preaccs.end()) return 0;
    preacc_tbl_atom_lvl::const_iterator atom_iter = res_iter->second.find(atom);
    if(atom_iter == res_iter->second.end()) return 0;
    return atom_iter->second;
  }
  
private:
  //! Initialize the A_SC_donors map
  static void
  build_SC_donor_tables();

  //! Initialize the A_SC_preaccs map
  static void
  build_SC_preacc_tables();

  PDBStructure *A_prot;

  static const my_float_t MIN_SALT_BRIDGE_SQUARED_LENGTH;
  static const my_float_t MAX_SALT_BRIDGE_SQUARED_LENGTH;
  static const my_float_t MIN_HBOND_SQUARED_LENGTH;
  static const my_float_t MAX_HBOND_SQUARED_LENGTH;
  static const my_float_t MIN_DHA_ANGLE;
  static const my_float_t MIN_H_A_AA_ANGLE;
  static const my_float_t COSINE_MIN_DHA_ANGLE;
  static const my_float_t COSINE_MIN_H_A_AA_ANGLE;

  static const my_float_t A_SC_H_positions[];  //!< Local hydrogen positions for protein side chain donor atoms
  static const my_float_t A_lig_H_positions[]; //!< Local hydrogen positions for ligand donor atoms
  static const my_float_t A_main_chain_NH_position[]; //!< local hydrogen position for main chain N-H
  static const size_t A_SC_donors_array_size;  //!< Number of elements in the A_SC_donors_array
  static const sidechain_donor_t A_SC_donors_array[];
  static SC_tbl_res_lvl A_SC_donors; //!< Map to look up protein sidechain donor info based on PDB residue and atom names
  static const size_t A_SC_preacc_array_size;  //!< Number of elements in the A_SC_preacc_array
  static const sidechain_preacc_t A_SC_preacc_array[];
  static preacc_tbl_res_lvl A_SC_preaccs;  //!< Map to help look up protein sidechain preacceptor info based on pdb residue and atom names
};
} // end namespace SimSite3D
#endif
