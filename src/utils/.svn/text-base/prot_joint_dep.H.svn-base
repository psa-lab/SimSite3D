#include <sstream>
#include <PDBStructure.H>

namespace SimSite3D{

typedef struct{
  const atom_type *atoms;
  int nrows;
  int num_atoms;
}joint_deps_t;

//! Data class to store the joints of one residue
class residue_joints{
public:
  //! Default constructor -- used by maps, etc before applying the assignment 
  //! operator
  residue_joints()
  {
    init();
  }

  //! Default constructor for a given residue
  residue_joints(const PDBStructure *prot, residue_vci res, 
                 const joint_deps_t &res_joints_info);

  //! Basic copy constructor
  residue_joints(const residue_joints &other)
  {
    do_copy(other);
  }

  //! Basic assignment operator
  const residue_joints&
  operator=(const residue_joints &other)
  {
    if(&other != this) do_copy(other);
    return *this;
  }

  //! Basic destructor
  ~residue_joints()
  {
    //if(axes_of_rotation) delete[] (axes_of_rotation);
//    if(A_init_end_effector_dir) delete[] (A_init_end_effector_dir);
    init();
  }

  //! Compute the so called phi, psi, and chi angles of the residue
  void
  compute_dihedral_angles(std::vector<my_float_t> *angles) const;

  //! Compute unit axes of rotation for side-chain joints using the current
  //! atomic positions
  /*! We need this function because we are relying on shallow copies (pointers
   * rather) to atomic positions.  Using this function makes it a bit less
   * likely to have stale axes of rotation, rather than requiring explicit
   * update of class variables after each moving of atoms (local and global).
   */
  bool
  compute_current_axes_of_rotation(my_float_t *rot_axes, int *num_axes,
                                   const int array_len) const;

  const int
  num_joints() const
  { return A_num_joints; }

  const std::vector<my_float_t>&
  orig_joint_angles() const
  { return A_orig_joint_angles; }

  const bool
  fail()
  { return A_fail; }

  // Variables to compute Jacobian of position with respect to change in joint
  // angles.  Although we don't necessarily need joint labels and joint
  // centers, it makes the processing loop much easier to follow.
  std::vector<atom_type> joint_labels;
  std::vector<const my_float_t *> joint_centers;  //!< Includes the CA position and one atom past the last joint
  //my_float_t *axes_of_rotation;  //!< An axis of rotation for each joint

private:

  bool
  check_atom(residue_vci res, atom_vci c_atom, atom_vci atoms_end, 
             atom_type expected_atom_name)
  {
    if(c_atom >= atoms_end  || c_atom->name != expected_atom_name){
      std::stringstream o_str;
      o_str << "The residue (" << res->chainID << ") " 
            << PDB_residues::residue_to_string(res->name) << res->number
            << " has missing atoms or does not have its atoms in canonical "
            << "order\nIt will not be modelled flexibly." << std::endl;
      warn("prot_joint_dep.H", "residue_joints::cstr", o_str.str());
      return false;
    }
    return true;
  }

  void
  do_copy(const residue_joints &other)
  {
    init();
    joint_centers.clear();
    joint_labels.clear();

    joint_centers.reserve(other.joint_centers.size());
    for(size_t i = 0; i < other.joint_centers.size(); ++i)
      joint_centers.push_back(other.joint_centers[i]);
    joint_labels.reserve(other.joint_labels.size());
    for(size_t i = 0; i < other.joint_labels.size(); ++i)
      joint_labels.push_back(other.joint_labels[i]);
/*
    if(other.A_num_joints > 0){
      axes_of_rotation = new my_float_t[3*other.A_num_joints];
      std::copy(other.axes_of_rotation, 
                other.axes_of_rotation + 3*other.A_num_joints, 
                axes_of_rotation); 
      A_init_end_effector_dir = new my_float_t[3*other.A_num_joints];
      std::copy(other.A_init_end_effector_dir, 
                other.A_init_end_effector_dir + 3*other.A_num_joints, 
                A_init_end_effector_dir); 
    }
*/
    A_num_joints = other.A_num_joints;
    A_resName = other.A_resName;
    A_atom_names = other.A_atom_names;
    A_atoms_pos_begin = other.A_atoms_pos_begin;
    A_atoms_pos_end = other.A_atoms_pos_end;
    A_prev_res_C_pos = other.A_prev_res_C_pos;
    A_next_res_N_pos = other.A_next_res_N_pos;
    A_orig_joint_angles = other.A_orig_joint_angles;
  }

  //! Initialize class variables
  void 
  init()
  {
//    axes_of_rotation = 0;
//    A_init_end_effector_dir = 0;
    A_fail = true;
    A_num_joints = 0;
    A_atom_names = 0;
    A_atoms_pos_begin = 0;
    A_atoms_pos_end = 0;
    A_prev_res_C_pos = 0;
    A_next_res_N_pos = 0;
  }

  //! Compute the dihedral angle for the directional bond chain A, B, C
  /*!
   * @return The dihedral angle in [-pi, pi]
   */
  my_float_t
  compute_dihedral_angle(const my_float_t *A, const my_float_t *B,
                         const my_float_t *C) const;

  //! A simple warning function to clean up the flow of the code
  void 
  warn_cannot_compute_dihedral(residue_vci res, std::string angle_name);

  // Variables to compute the value of the sidechain angles (chi angles)
  residue_type A_resName;
  int A_num_joints;
  bool A_fail;
//  my_float_t *A_init_end_effector_dir; //!< Unit vector from each joint to end effector that is not rotated by its corresponding joint or those after; used to measure amount of rotation about each joint relative to initial angles (positions)
  const atom_type *A_atom_names;
  const my_float_t *A_atoms_pos_begin;
  const my_float_t *A_atoms_pos_end;
  const my_float_t *A_prev_res_C_pos;
  const my_float_t *A_next_res_N_pos;
  std::vector<my_float_t> A_orig_joint_angles;

  static const std::string A_fname;
};

//! Class to handle mobile residues -- rotatable joints for each residue
/*! Build the data structure based using add_residue_info()
 *  Get the joints using get_joints() and check that the returned iterator
 *  is valid by comparing it against joints_end()
 */
class prot_joint_dep{
public:

  typedef std::map<atom_type, std::vector<atom_type> > joint_dep_tbl_atom_lvl;
  typedef std::map<residue_type, joint_dep_tbl_atom_lvl> joint_dep_tbl_res_lvl;
  typedef std::map<residue_vci, residue_joints>::const_iterator 
    res_to_joints_mci;
//  typedef std::map<int, residue_joints>::const_iterator 
//    res_to_joints_mci;

  prot_joint_dep()
  {
    init();
  }

  //! Add mobile atoms (based on bond network not other constraints)
  /*! Do not use this to add neighboring residues which we will hold fixed.
   *
   */
  bool
  add_residue_info(const PDBStructure *prot, residue_vci res, 
                   const residue_joints **j);

  res_to_joints_mci
  joints_begin() const
  {
    return A_joint_blocks.begin();
  }

  res_to_joints_mci
  joints_end() const
  {
    return A_joint_blocks.end();
  }

  const joint_deps_t&
  get_joints_info(residue_type res_name)
  {
    return A_res_map[res_name];
  }

  const bool
  fail()
  { return A_fail; }
  
private:
 
  void
  init();

  bool A_fail;
  joint_dep_tbl_res_lvl A_atoms_dep_on_joint;

  std::vector<atom_type> A_unknown_residue_vec;
  std::vector<atom_type> A_unknown_atom_vec;

  std::map<residue_type, joint_deps_t> A_res_map;

  //std::map<char, int> A_chainID_shift;
  std::map<residue_vci, residue_joints> A_joint_blocks;  //!< Groups of consecutive rows that correspond to joints from one residue -- ordered by residue number and chainID
//  std::map<int, residue_joints> A_joint_blocks;  //!< Groups of consecutive rows that correspond to joints from one residue -- ordered by residue number and chainID

  static const atom_type ALA_atoms[];
  static const atom_type ARG_atoms[];
  static const atom_type ASN_atoms[];
  static const atom_type ASP_atoms[];
  static const atom_type CYS_atoms[];
  static const atom_type GLN_atoms[];
  static const atom_type GLU_atoms[];
  static const atom_type GLY_atoms[];
  static const atom_type HIS_atoms[];
  static const atom_type ILE_atoms[];
  static const atom_type LEU_atoms[];
  static const atom_type LYS_atoms[];
  static const atom_type MET_atoms[];
  static const atom_type PHE_atoms[];
  static const atom_type PRO_atoms[];
  static const atom_type SER_atoms[];
  static const atom_type THR_atoms[];
  static const atom_type TRP_atoms[];
  static const atom_type TYR_atoms[];
  static const atom_type VAL_atoms[];

  static const joint_deps_t ALA_deps;
  static const joint_deps_t ARG_deps;
  static const joint_deps_t ASN_deps;
  static const joint_deps_t ASP_deps;
  static const joint_deps_t CYS_deps;
  static const joint_deps_t GLN_deps;
  static const joint_deps_t GLU_deps;
  static const joint_deps_t GLY_deps;
  static const joint_deps_t HIS_deps;
  static const joint_deps_t ILE_deps;
  static const joint_deps_t LEU_deps;
  static const joint_deps_t LYS_deps;
  static const joint_deps_t MET_deps;
  static const joint_deps_t PHE_deps;
  static const joint_deps_t PRO_deps;
  static const joint_deps_t SER_deps;
  static const joint_deps_t THR_deps;
  static const joint_deps_t TRP_deps;
  static const joint_deps_t TYR_deps;
  static const joint_deps_t VAL_deps;
};

}
