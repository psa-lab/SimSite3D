#include <prot_joint_dep.H>
#include <iostream>
#include <string>

using namespace ASCbase;

const std::string residue_joints::A_fname = "prot_joint_dep";

// Ignore GLY, ALA, and PRO at this time.  They will have to be revisited 
// if we get around to adding main chain flexiblity
const atom_type prot_joint_dep::ALA_atoms[] = {N, CA, C, O, CB};
const joint_deps_t prot_joint_dep::ALA_deps = { ALA_atoms, 0, 5 };

const atom_type prot_joint_dep::ARG_atoms[] = 
  {N, CA, C, O, CB, CG, CD, NE, CZ, NH1, NH2};
const joint_deps_t prot_joint_dep::ARG_deps = { ARG_atoms, 3, 11 };

const atom_type prot_joint_dep::ASN_atoms[] = {N, CA, C, O, CB, CG, OD1, ND2};
const joint_deps_t prot_joint_dep::ASN_deps = { ASN_atoms, 2, 8 };

const atom_type prot_joint_dep::ASP_atoms[] = {N, CA, C, O, CB, CG, OD1, OD2};
const joint_deps_t prot_joint_dep::ASP_deps = { ASP_atoms, 2, 8 };

const atom_type prot_joint_dep::CYS_atoms[] = {N, CA, C, O, CB, SG};
const joint_deps_t prot_joint_dep::CYS_deps = { CYS_atoms, 1, 6 };

const atom_type prot_joint_dep::GLN_atoms[] = 
  {N, CA, C, O, CB, CG, CD, OE1, NE2};
const joint_deps_t prot_joint_dep::GLN_deps = { GLN_atoms, 3, 9 };

const atom_type prot_joint_dep::GLU_atoms[] = 
  {N, CA, C, O, CB, CG, CD, OE1, OE2};
const joint_deps_t prot_joint_dep::GLU_deps = { GLU_atoms, 3, 9 };

const atom_type prot_joint_dep::GLY_atoms[] = {N, CA, C, O};
const joint_deps_t prot_joint_dep::GLY_deps = { GLY_atoms, 0, 4 };

const atom_type prot_joint_dep::HIS_atoms[] = 
  {N, CA, C, O, CB, CG, ND1, CD2, CE1, NE2};
const joint_deps_t prot_joint_dep::HIS_deps = { HIS_atoms, 2, 10 };

// Remember, ILE is different from the others in that there is a rotatable
// bond AFTER the first branching point
const atom_type prot_joint_dep::ILE_atoms[] = {N, CA, C, O, CB, CG1, CG2, CD1};
const joint_deps_t prot_joint_dep::ILE_deps = { ILE_atoms, 2, 8 };

const atom_type prot_joint_dep::LEU_atoms[] = {N, CA, C, O, CB, CG, CD1, CD2};
const joint_deps_t prot_joint_dep::LEU_deps = { LEU_atoms, 2, 8 };

const atom_type prot_joint_dep::LYS_atoms[] = {N, CA, C, O, CB, CG, CD, CE, NZ};
const joint_deps_t prot_joint_dep::LYS_deps = { LYS_atoms, 4, 9 };

const atom_type prot_joint_dep::MET_atoms[] = {N, CA, C, O, CB, CG, SD, CE};
const joint_deps_t prot_joint_dep::MET_deps = { MET_atoms, 3, 8 };

const atom_type prot_joint_dep::PHE_atoms[] = 
  {N, CA, C, O, CB, CG, CD1, CD2, CE1, CE2, CZ};
const joint_deps_t prot_joint_dep::PHE_deps = { PHE_atoms, 2, 11 };

const atom_type prot_joint_dep::PRO_atoms[] = {N, CA, C, O, CB, CG, CD};
const joint_deps_t prot_joint_dep::PRO_deps = { PRO_atoms, 0, 7 };

const atom_type prot_joint_dep::SER_atoms[] = {N, CA, C, O, CB, OG};
const joint_deps_t prot_joint_dep::SER_deps = { SER_atoms, 1, 6 };

const atom_type prot_joint_dep::THR_atoms[] = {N, CA, C, O, CB, OG1, CG2};
const joint_deps_t prot_joint_dep::THR_deps = { THR_atoms, 1, 7 };

const atom_type prot_joint_dep::TRP_atoms[] = 
  {N, CA, C, O, CB, CG, CD1, CD2, NE1, CE2, CE3, CZ2, CZ3, CH2};
const joint_deps_t prot_joint_dep::TRP_deps = { TRP_atoms, 2, 14 };

const atom_type prot_joint_dep::TYR_atoms[] = 
  {N, CA, C, O, CB, CG, CD1, CD2, CE1, CE2, CZ, OH};
const joint_deps_t prot_joint_dep::TYR_deps = { TYR_atoms, 2, 12 };

const atom_type prot_joint_dep::VAL_atoms[] = {N, CA, C, O, CB, CG1, CG2};
const joint_deps_t prot_joint_dep::VAL_deps = { VAL_atoms, 1, 7 };

residue_joints::residue_joints(const PDBStructure *prot, residue_vci res, 
                               const joint_deps_t &res_joints_info)
{
  init();
  if(res_joints_info.nrows == 0){
    //axes_of_rotation = 0;
    A_num_joints = 0;
    return;
  }

  A_num_joints = res_joints_info.nrows;
  std::cout << "Inside residue_joints, the number of joints is: "
            << A_num_joints << "\n";
  //axes_of_rotation = new my_float_t[3*A_num_joints];
  //A_init_end_effector_dir = new my_float_t[3*A_num_joints];
  joint_centers.reserve(A_num_joints + 1);
  joint_labels.reserve(A_num_joints);

  // Assume residue atoms are in canonical order and no atoms are missing
  atom_vci c_atom = res->atoms_begin;
  const atom_type* c_atom_name = res_joints_info.atoms;

  std::cout << "The residue (" << res->chainID << ") " 
            << PDB_residues::residue_to_string(res->name) << res->number
            << " will be allowed to move\n";
  
  // Main chain atoms
  for(int i = 0; i < 4; ++i, ++c_atom, ++c_atom_name){
    if(!check_atom(res, c_atom, res->atoms_end, *c_atom_name)) return;
    if(c_atom->name == CA) joint_centers.push_back(c_atom->pos);
  }

  // Side chain atoms
  for(int i = 0; i < A_num_joints; ++i, ++c_atom, ++c_atom_name){
    std::cout << "adding joints for SC atoms -- i: " <<  i << std::endl;
    if(!check_atom(res, c_atom, res->atoms_end, *c_atom_name)) return;
    std::cout << "Time to push back items\n";
    joint_labels.push_back(c_atom->name);
    joint_centers.push_back(c_atom->pos);
/*
    unit_vector(&axes_of_rotation[3*i], joint_centers[i+1], 
                joint_centers[i]);
    unit_vector(&A_init_end_effector_dir[3*i], , joint_centers[i]);
    std::cout << "axis of rot: " 
              << axes_of_rotation[3*i] << " "
              << axes_of_rotation[3*i + 1] << " "
              << axes_of_rotation[3*i + 2] << "\n";
*/
  }  

  // Keep track of initial directions (i.e. dihedral angles as measured by the
  // 3 vectors)
  joint_centers.push_back(c_atom->pos);
  
  // Finish checking the atom names (i.e. for those atoms which are not joints)
  for(int i = 0; i < res_joints_info.num_atoms - 4 - A_num_joints; ++i){
    if(!check_atom(res, c_atom, res->atoms_end, *c_atom_name)) return;
    ++c_atom;
    ++c_atom_name;
  }
  
  std::cout << "# joint_labels: " << joint_labels.size() << std::endl;
  std::cout << "# joint_centers: " << joint_centers.size() << std::endl;

  // Set values used to compute chi angles
  A_resName = res->name;
  A_atom_names = res_joints_info.atoms;
  A_atoms_pos_begin = res->atoms_begin->pos;
  A_atoms_pos_end = res->atoms_begin->pos + 3*res_joints_info.num_atoms;

  const my_float_t *C_pos = A_atoms_pos_begin + 6;
  residue_vci next_residue = res + 1;

  // check that next_residue is valid
  chain_const_iter chain = prot->chains_begin();
  for( ; chain < prot->chains_end(); ++chain)
    if(chain->chainID == res->chainID) break;

  if(next_residue >= chain->residues_end || chain->chainID != res->chainID)
    warn_cannot_compute_dihedral(res, "psi");
  // Check for chain breaks
  else{
    atom_vci N_atom = next_residue->get_atom(N);

    // N-C bond length is 1.32 (A) -- allow up to 1.5 (A) to be nice
    if(N_atom != atom_t::NULL_ATOM_VCI &&
       2.25 > dist_squared(C_pos, N_atom->pos)){
      A_next_res_N_pos = N_atom->pos;
    }else warn_cannot_compute_dihedral(res, "psi");
  }

  const my_float_t *N_pos = A_atoms_pos_begin;
  residue_vci prev_res = res - 1;

  // check that previous is valid
  for(chain = prot->chains_begin(); chain < prot->chains_end(); ++chain)
    if(chain->chainID == res->chainID) break;

  if(prev_res < chain->residues_begin || chain->chainID != res->chainID)
    warn_cannot_compute_dihedral(res, "phi");
  // Check for chain breaks
  else{
    atom_vci C_atom = prev_res->get_atom(C);

    // N-C bond length is 1.32 (A) -- allow up to 1.5 (A) to be nice
    if(C_atom != atom_t::NULL_ATOM_VCI &&
       2.25 > dist_squared(N_pos, C_atom->pos)){
      A_prev_res_C_pos = C_atom->pos;
    }else warn_cannot_compute_dihedral(res, "phi");
  }

  compute_dihedral_angles(&A_orig_joint_angles);
  A_fail = false;
}

bool
residue_joints::compute_current_axes_of_rotation(my_float_t *rot_axes,
                                                 int *num_axes,
                                                 const int array_len) const
{
  if(3*A_num_joints > array_len){
    std::cerr << "Cannot compute the current axes of rotation because too \n"
              << "little memory was given.\n";
    return false;
  }

  std::fill(rot_axes, rot_axes + 3*A_num_joints, 0.0);
  for(int i = 0; i < A_num_joints; ++i)
    unit_vector(&rot_axes[3*i], joint_centers[i+1], joint_centers[i]);
  *num_axes = A_num_joints;

  return true;
} 

void
residue_joints::compute_dihedral_angles(std::vector<my_float_t> *angles) const
{
  // In this function we explicitly use constant pointers in an effort
  // to make the logic a bit more clear -- the hope is the compiler will
  // optimize the extra pointers away

  // Compute the psi and phi bonds --
  // The bonds are C_{i-1} -> N -> CA -> C -> N_{i+1}
  const my_float_t *N_pos = A_atoms_pos_begin;
  const my_float_t *CA_pos = A_atoms_pos_begin + 3;
  const my_float_t *C_pos = A_atoms_pos_begin + 6;
  my_float_t bonds[18];
  vector(3, CA_pos, N_pos, bonds + 3);
  vector(3, C_pos, CA_pos, bonds + 6);

  // Compute phi as per Dunbrack orientation
  // Label the bonds:
  // A = C_{i-1} -> N = bonds
  // B = N -> CA = bonds + 3
  // C = CA -> C = bonds + 6
  if(A_prev_res_C_pos){
    vector(3, N_pos, A_prev_res_C_pos, bonds);
    my_float_t phi = compute_dihedral_angle(bonds, bonds + 3, bonds + 6);
    angles->push_back(phi);
  }else angles->push_back(my_float_max);
  
  // We need the reverse directions for psi as per Dunbrack
  // Label the bonds:
  // A = N_{i+1} -> C
  // B = C -> CA
  // C = CA -> N
  if(A_next_res_N_pos){
     vector(3, A_next_res_N_pos, C_pos, bonds + 9);
    for(int i = 0; i < 12; ++i) bonds[i] *= -1.0;
    my_float_t psi = compute_dihedral_angle(bonds + 9, bonds + 6, bonds + 3);
    angles->push_back(psi);

    // We need to "reuse" the N -> CA bond
    for(int i = 0; i < 3; ++i) bonds[i] = -1.0 * bonds[i+3];
  }else{
    std::copy(bonds + 3, bonds + 6, bonds);
    angles->push_back(my_float_max);
  }

  if(A_resName == ALA || A_resName == GLY || A_resName == PRO) return;

  // Compute the CA -> CB vector
  const my_float_t *CB_pos = A_atoms_pos_begin + 12;
  vector(3, CB_pos, CA_pos, bonds + 3); 

  // The order of the atoms was already checked so we know that the order 
  // of the atoms in the residue is the same as the order in A_atom_names
  if(A_resName == ILE){
    const my_float_t *CG1_pos = A_atoms_pos_begin + 15;
    vector(3, CG1_pos, CB_pos, bonds + 6);
    const my_float_t *CD1_pos = A_atoms_pos_begin + 21;
    vector(3, CD1_pos, CG1_pos, bonds + 9);
  }else{
    // Compute the rest of the vectors used to compute the side-chain dihedral 
    // angles.  We are assuming if there are two atoms at a certain level
    // (e.g. OD1 and OD2) that the atom with a '1' in its name is chosen 
    // (i.e. OD1, OG1, etc) in computation of the dihedral angles.
    const my_float_t *tail = A_atoms_pos_begin + 12;
    my_float_t *b = bonds + 6;
    for(int i = 0; i < A_num_joints; ++i){
      vector(3, tail + 3, tail, b);
      tail += 3;
      b += 3;
    }
  }

  // Compute the "chi" angles
  const my_float_t *first_bond = bonds;
  for(int i = 0; i < A_num_joints; ++i){
    my_float_t chi = 
      compute_dihedral_angle(first_bond, first_bond + 3, first_bond + 6);
    angles->push_back(chi);
    first_bond += 3;
  }

}

my_float_t
residue_joints::compute_dihedral_angle(const my_float_t *A, const my_float_t *B,
                                       const my_float_t *C) const
{
  // Note: since we are using atan2 we do not need unit vectors because
  // we are only interested in the ratio of the dot and cross products.
  // It is fairly straightforward to see that the dihedral angles specified
  // by the vectors A, B, C is given by
  // atan2(|B| A*(BxC), (AxB) * (BxC)) where * denotes inner product 
  // (dot product)  and x denotes the cross product.
  // Using identities for cross products (found in linear algebra books or 
  // Wolfram.com) we can reduce 
  // | (AxB) x (BxC) | to |B| * det(ABC) which is |B| A*(BxC)
  my_float_t A_x_B[3], B_x_C[3];
  cross(A, B, A_x_B);
  cross(B, C, B_x_C);
  my_float_t len_B = std::sqrt(B[3]*B[3] + B[4]*B[4] + B[5]*B[5]);
  return atan2(len_B * dot(A, B_x_C), dot(A_x_B, B_x_C));
}

void
residue_joints::warn_cannot_compute_dihedral(residue_vci res, 
                                             std::string angle_name)
{
  std::ostringstream ostr;
  ostr << "Break in the chain does not allow computation of the "
       << angle_name << " angle for (" << res->chainID << ") " 
       << PDB_residues::residue_to_string(res->name) << res->number << "\n";
  warn(A_fname, "residue_joints()", ostr.str());
}

void
prot_joint_dep::init()
{
  if(A_atoms_dep_on_joint.size()) return;

  A_res_map[ALA] = ALA_deps;
  A_res_map[ARG] = ARG_deps;
  A_res_map[ASN] = ASN_deps;
  A_res_map[ASP] = ASP_deps;
  A_res_map[CYS] = CYS_deps;
  A_res_map[GLN] = GLN_deps;
  A_res_map[GLU] = GLU_deps;
  A_res_map[GLY] = GLY_deps;
  A_res_map[HIS] = HIS_deps;
  A_res_map[ILE] = ILE_deps;
  A_res_map[LEU] = LEU_deps;
  A_res_map[LYS] = LYS_deps;
  A_res_map[MET] = MET_deps;
  A_res_map[PHE] = PHE_deps;
  A_res_map[PRO] = PRO_deps;
  A_res_map[SER] = SER_deps;
  A_res_map[THR] = THR_deps;
  A_res_map[TRP] = TRP_deps;
  A_res_map[TYR] = TYR_deps;
  A_res_map[VAL] = VAL_deps;
  A_fail = false;
}

bool
prot_joint_dep::add_residue_info(const PDBStructure *prot, residue_vci res, 
                                 const residue_joints **j)
{

  res_to_joints_mci res_map_iter = A_joint_blocks.find(res);
  if(res_map_iter == A_joint_blocks.end()){
    residue_joints my_joints(prot, res, A_res_map[res->name]);
    if(my_joints.fail()){
      *j = 0;
      return false;
    }
    A_joint_blocks[res] = my_joints;
    res_map_iter = A_joint_blocks.find(res);
  }
  *j = &(res_map_iter->second);
  return true;
}
