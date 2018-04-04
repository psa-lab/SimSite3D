#include <SurfDepsJoints.H>

using namespace SimSite3D;

SurfDepsJoints::SurfDepsJoints(prot_joint_dep &joints, PDBStructure *prot,
                               geometry::TransformableTrimesh *surf,
                               ModelHbondSurfaces* hbond_surfaces,
			       HbondPoints *hbond_pts) 
                                : A_joints(joints)
{
  init();
  if(surf) add_mol_surf_vertices(surf, prot);
  add_hbond_points(hbond_surfaces, hbond_pts, prot);
  if(!check_containers()) return;


  // compute number of mobile vertices and atoms
  int &num_verts = A_J_cols_in_blocks[0];
  int &num_cap_pts = A_J_cols_in_blocks[1];
  int &num_atoms = A_J_cols_in_blocks[2];
  int &num_joints = A_J_nrows;
  num_verts = num_cap_pts = num_atoms = num_joints = 0;
  prot_joint_dep::res_to_joints_mci res2j_iter = A_joints.joints_begin();
  res_to_verts_mi res2V_iter = A_verts_dep_on_joints.begin();
  res_to_verts_mi res2HB_iter = A_hbond_pts_dep_on_joints.begin();
  res_to_verts_mi res2A_iter = A_atoms_dep_on_joints.begin();
  for( ; res2A_iter != A_atoms_dep_on_joints.end(); ++res2A_iter){
    num_joints += res2j_iter->second.num_joints();
    num_atoms += res2A_iter->second.verts.size();
    res2A_iter->second.set_nrows(res2j_iter->second.num_joints());

    if(A_have_surf){
      res2V_iter->second.set_nrows(res2j_iter->second.num_joints());
      for(size_t i = 0; i < res2V_iter->second.verts.size(); ++i)
        A_vert_idz.push_back((res2V_iter->second.verts[i] -
                             surf->vertices_begin())/3);
      num_verts += res2V_iter->second.verts.size();
    } 

    if(A_have_hbond_surfs){
      res2HB_iter->second.set_nrows(res2j_iter->second.num_joints());
  std::cout << "HB cap for residue (" << res2HB_iter->first->chainID << ") " 
            << PDB_residues::residue_to_string(res2HB_iter->first->name) 
            << res2HB_iter->first->number
            << std::endl;
      std::cout << "number of hbond points: " 
                << res2HB_iter->second.verts.size() << std::endl;
      num_cap_pts += res2HB_iter->second.verts.size();
    }

#if 0 // check needs to be rewritten
    // This check needs to come last as we may have needed to adjust one
    // or more of the maps 
    if(res2V_iter->first != res2j_iter->first ||
       res2A_iter->first != res2j_iter->first ||
       res2HB_iter->first != res2j_iter->first){
      std::cerr << "Error in processing residue joints\n" << std::endl;
#if 1
    std::cout << "Vertex residue (" << res2V_iter->first->chainID << ") " 
              << PDB_residues::residue_to_string(res2V_iter->first->name) 
              << res2V_iter->first->number
              << std::endl;
    std::cout << "Joints residue (" << res2j_iter->first->chainID << ") " 
              << PDB_residues::residue_to_string(res2j_iter->first->name) 
              << res2j_iter->first->number
              << std::endl;
    std::cout << "Atoms residue (" << res2A_iter->first->chainID << ") " 
              << PDB_residues::residue_to_string(res2A_iter->first->name) 
              << res2A_iter->first->number
              << std::endl;
#endif
      return;
    }
#endif

    if(A_have_surf) ++res2V_iter;
    if(A_have_hbond_surfs) ++res2HB_iter;
    ++res2j_iter;
  }

  if(!A_have_surf && !A_have_hbond_surfs){
    std::cerr << "\nERROR!\nMust use either the molecular surfaces or hbond "
              << "caps (or both) to drive the morphing\n\n";
    return;
  }

  A_J_ncols = 3*(num_atoms + num_verts + num_cap_pts);
  A_J_nrows = num_joints;  // Each joint has only one degree of freedom
  A_JJtranspose = new my_float_t[A_J_nrows*A_J_nrows];
  A_fail = false;
}

void
SurfDepsJoints::compute_pinvJ_transpose()
{
  // compute JJ^t + lambda I
  compute_JJtranspose();
  const my_float_t lambda = 1E-04;
  for(int i = 0; i < A_J_nrows; ++i) A_JJtranspose[i*A_J_nrows + i] += lambda;

  std::cout << "num rows in J is: " << A_J_nrows << "\n\n" << std::endl;
  int *IPIV = new int[A_J_nrows];
  int rv = la_dgetrf(A_J_nrows, A_J_nrows, A_JJtranspose, A_J_nrows, IPIV);
//  std::cout << "return from dgetrf: " << rv << "\n";

  int lwork = 64*A_J_nrows;
  my_float_t *work = new my_float_t[lwork];
  rv = la_dgetri(A_J_nrows, A_JJtranspose, A_J_nrows, IPIV, work, lwork);
//  std::cout << "return from dgetri: " << rv << "\n";

#if 0
  for(int i = 0; i < A_J_nrows; ++i){
    std::cout << "\n";
    for(int j = 0; j < A_J_nrows; ++j) 
      std::cout << A_JJtranspose[i*A_J_nrows + j] << " ";
  }
  std::cout << "\n";
#endif

  compute_portion_of_pinvJtranspose(A_atoms_dep_on_joints.begin(),
                                    A_atoms_dep_on_joints.end(),
                                    A_JJtranspose);
  if(A_have_surf)
    compute_portion_of_pinvJtranspose(A_verts_dep_on_joints.begin(),
                                      A_verts_dep_on_joints.end(),
                                      A_JJtranspose);
  if(A_have_hbond_surfs)
    compute_portion_of_pinvJtranspose(A_hbond_pts_dep_on_joints.begin(),
                                      A_hbond_pts_dep_on_joints.end(),
                                      A_JJtranspose);
  if(IPIV) delete[] IPIV;
  if(work) delete[] work;
  IPIV = 0;
  work = 0;
}

// Rewrite this function to be more useful in debugging/testing
// Currently its original purpose is unknown
void
SurfDepsJoints::write_point_correspondences(std::ostream &out) const
{
  prot_joint_dep::res_to_joints_mci res2j_iter = A_joints.joints_begin();
  res_to_verts_mci res2pt_iter = A_verts_dep_on_joints.begin();
  for( ; res2pt_iter != A_verts_dep_on_joints.end(); ++res2pt_iter){
    residue_vci res = res2pt_iter->first;
    const pos_dep_on_joints &dep_points = res2pt_iter->second;

    for(size_t i = 0; i < dep_points.num_joints_used.size(); ++i){
      out << (res->atoms_begin + 4 + dep_points.num_joints_used[i])->atom_num
          << " " 
          << dep_points.verts[i][0] << " "
          << dep_points.verts[i][1] << " "
          << dep_points.verts[i][2] << "\n";
    }
  }
}

void 
SurfDepsJoints::compute_portion_of_JJtranspose(res_to_verts_mi pts_begin, 
                                               res_to_verts_mi pts_end,
                                               my_float_t *JJtrans)
{
  my_float_t *JJtrans_block_start = JJtrans;

  // Add to each block of JJtranspose the contribution of this portion of
  // the Jacobian
  res_to_verts_mi res2pt_iter;
  for(res2pt_iter = pts_begin; res2pt_iter != pts_end; ++res2pt_iter){
    pos_dep_on_joints *dep_points = &(res2pt_iter->second);
    const my_float_t *block_row_start = dep_points->J_block();
    const int nrows = dep_points->nrows();
    const int ncols = dep_points->ncols();
    AAtranspose(nrows, ncols, block_row_start, ncols, JJtrans_block_start,
                A_J_nrows, 1.0);
    JJtrans_block_start += nrows*A_J_nrows + nrows;
  }
}

void
SurfDepsJoints::compute_portion_of_pinvJtranspose(res_to_verts_mi pts_begin,
                                                  res_to_verts_mi pts_end,
                                                  my_float_t *JJtrans_inv)
{
  my_float_t *inv_JJtrans_block_start = JJtrans_inv;

  // Compute each block of this portion of pinv(J)^T
  res_to_verts_mi res2pt_iter;
  for(res2pt_iter = pts_begin; res2pt_iter != pts_end; ++res2pt_iter){
    pos_dep_on_joints *dep_points = &(res2pt_iter->second);
    const my_float_t *block_row_start = dep_points->J_block();
    const int nrows = dep_points->nrows();
    const int ncols = dep_points->ncols();

    my_float_t *pinv_Jtrans_block = dep_points->mutable_pinv_Jtrans_block();
    my_gemm(nrows, ncols, nrows, 1.0, inv_JJtrans_block_start, A_J_nrows,
            block_row_start, ncols, pinv_Jtrans_block, ncols, 0.0);

    inv_JJtrans_block_start += nrows*A_J_nrows + nrows;
  }
}

void
SurfDepsJoints::compute_portion_pinvJtrans_times_grad(res_to_verts_mi pts_begin,
                                                      res_to_verts_mi pts_end,
                                                      const my_float_t *grad_F,
                                                      my_float_t *Q)
{
  my_float_t *q = Q; 
  const my_float_t *grad_F_p = grad_F;
  res_to_verts_mi res2pt_iter;
  for(res2pt_iter = pts_begin; res2pt_iter != pts_end; ++res2pt_iter){
    pos_dep_on_joints *dep_points = &(res2pt_iter->second);
    const int nrows = dep_points->nrows();
    const int ncols = dep_points->ncols();
    const my_float_t *pinv_Jtrans_block = dep_points->pinv_Jtrans_block();

    // Need to add to q since we may have contributions from more
    // than one concept (surface pulling, repulsive atomic overlap, etc)
    my_gemm(nrows, 1, ncols, 1.0, pinv_Jtrans_block, ncols, grad_F_p, 1,
            q, 1, 1.0);

    grad_F_p += ncols;
    q += nrows;
  }
}

void
SurfDepsJoints::stop_rotation_given_atoms(const std::vector<bool> &overlap_vec,
                                          std::vector<bool> *fix_Q)
{
// In the future one may wish to check, for each overlapping atom 
// which could move, whether one could rotate in the positive or 
// negative direction. 


//Given the atoms that overlap -- determine which joints need to be fixed
  fix_Q->resize(A_J_nrows);
  std::fill(fix_Q->begin(), fix_Q->end(), false);
  std::vector<bool>::const_iterator overlap_flag = overlap_vec.begin();
  std::vector<bool>::iterator fix_flag = fix_Q->begin();
  res_to_verts_mi res2pt_iter = A_atoms_dep_on_joints.begin();
  for( ; res2pt_iter != A_atoms_dep_on_joints.end(); ++res2pt_iter){
    pos_dep_on_joints *dep_points = &(res2pt_iter->second);
    const int nrows = dep_points->nrows();

    int num_fixed_joints = 0;
    for(size_t i = 0; i < dep_points->num_joints_used.size(); ++i){
      if(*overlap_flag) num_fixed_joints = dep_points->num_joints_used[i];
      ++overlap_flag; 
    }

    for(int i = 0; i < num_fixed_joints; ++i) *(fix_flag + i) = true;
   
    fix_flag += nrows;
  }
}

void
SurfDepsJoints::apply_atoms_pinvJ_to_grad(const my_float_t *grad_F,
                                          my_float_t *Q, 
                                          const std::vector<bool> &overlap_vec)
{
  my_float_t *q = Q; 
  const my_float_t *grad_F_p = grad_F;
  std::vector<bool>::const_iterator overlap_flag = overlap_vec.begin();
  res_to_verts_mi res2pt_iter = A_atoms_dep_on_joints.begin();
  for( ; res2pt_iter != A_atoms_dep_on_joints.end(); ++res2pt_iter){
    pos_dep_on_joints *dep_points = &(res2pt_iter->second);
    const int nrows = dep_points->nrows();
    const int ncols = dep_points->ncols();
    const my_float_t *pinv_Jtrans_block = dep_points->pinv_Jtrans_block();

//    std::cout << "\n";
//    for(int i = 0; i < nrows; ++i){
//      for(int j = 0; j < ncols; ++j){
//        std::cout << pinv_Jtrans_block[i*ncols + j] << " ";
//      }
//      std::cout << "\n";
//    }
//    std::cout << "\n";

std::cout << "need to keep track of which directions to move to alleviate overlaps and if a joint needs to move in opposing directions just hold that joint fixed? \n\n\n";


    //////////////////////////////////////////////////////////////////////
    // If more than the tolerated overlap exists, we use a weighted
    // averaging scheme to attempt to alleviate the overlap
    // For each atom that has overlap, divide the rotations to 
    // remove the overlap amongst the joints that determine its position.
    // If more than one atom has overlap, average each joint rotation 
    // by the number of overlaps that joint is attempting to resolve
    //////////////////////////////////////////////////////////////////////

    residue_vci res = res2pt_iter->first;
    std::cout << "Inside matrix mult: (" << res->chainID << ") "
              << PDB_residues::residue_to_string(res->name) << res->number << "\n";
 //             << " " << PDB_residues::atom_to_string(res_atom->name) << " -- ";
    
    // Matrix - Vector multiplication 
    for(int i = 0; i < nrows; ++i){
      for(int j = 0; j < ncols/3; ++j){
        if(*(overlap_flag + j) == false) continue;
        if(dep_points->num_joints_used[j] < i+1) continue;

        std::cout << "(overlap flag + " << j << ") is true\n";

        // For each atom that has overlap, divide the rotations to 
        // remove the overlap amongst the joints that determine its position.
        for(int zz = 0; zz < 3; ++zz){
          q[i] += pinv_Jtrans_block[i*ncols + 3*j + zz] * grad_F_p[3*j + zz]
                    / dep_points->num_joints_used[j];
        }
      }
    }

    // If more than one atom has overlap, average each joint rotation 
    // by the number of overlaps that joint is attempting to resolve
    std::vector<int> num_overlaps(nrows);
    std::fill(num_overlaps.begin(), num_overlaps.end(), 0);  
    for(int i = 0; i < ncols/3; ++i){
      if(*(overlap_flag + i))
        for(int j = 0; j < dep_points->num_joints_used[i]; ++i)
          ++num_overlaps[j];
    }
    for(size_t i = 0; i < num_overlaps.size(); ++i)
      if(num_overlaps[i] > 1) q[i] /= num_overlaps[i];

    std::cout << "Q: ";
    for(size_t i = 0; i < num_overlaps.size(); ++i)
      std::cout << q[i] << " ";
    std::cout << "\n\n";
      

    grad_F_p += ncols;
    overlap_flag += ncols/3;
    q += nrows;
  }
}


void
SurfDepsJoints::update_positions(const my_float_t *delta_q)
{
  // Move each group by its change in angle
  const my_float_t *c_del_q = delta_q;
  prot_joint_dep::res_to_joints_mci res2j_iter = A_joints.joints_begin();
  res_to_verts_mi res2V_iter = A_verts_dep_on_joints.begin();
  res_to_verts_mi res2A_iter = A_atoms_dep_on_joints.begin();
  res_to_verts_mi res2HB_iter = A_hbond_pts_dep_on_joints.begin();
  for( ; res2A_iter != A_atoms_dep_on_joints.end(); ++res2A_iter, ++res2j_iter){
    const residue_joints &joints = res2j_iter->second;

    for(int j_num = 0; j_num < joints.num_joints(); ++j_num, ++c_del_q){
      // Get the rotation matrix for change in joint angle about the
      // current joint axis
      my_float_t axis[3];
      unit_vector(axis, joints.joint_centers[j_num + 1],
                  joints.joint_centers[j_num]);
      Quaternion Q(cos(*c_del_q), sin(*c_del_q), axis);
      my_float_t R[9];
      Q.get_ortho_rot_mat(R);

      // Move the atoms
      rotate_points(&(res2A_iter->second), R, joints.joint_centers[j_num], 
                    j_num);

      // Move the molecular surface patch vertices
      if(A_have_surf){
        rotate_points(&(res2V_iter->second), R, joints.joint_centers[j_num], 
                      j_num);
      }
      // Move the hbond points or hbond cap vertices
//      if(A_have_hbond_surfs){
      rotate_points(&(res2HB_iter->second), R, joints.joint_centers[j_num], 
                    j_num);
//      }
    }
    if(A_have_surf) ++res2V_iter;
    ++res2HB_iter;
    
  }
}

void
SurfDepsJoints::compute_dihedral_angles(std::vector< std::vector<my_float_t> > *angles) const
{
  // Move each group by its change in angle
  prot_joint_dep::res_to_joints_mci res2j_iter = A_joints.joints_begin();
  for( ; res2j_iter != A_joints.joints_end(); ++res2j_iter){
    std::vector<my_float_t> tmp;
    angles->push_back(tmp);
    res2j_iter->second.compute_dihedral_angles(&(*(angles->end() - 1)));
  }
}

void
SurfDepsJoints::add_dependent_atoms(const PDBStructure *prot, residue_vci res)
{
  A_site_residues[res] = true;
  if(res->name == ALA || res->name == GLY || res->name == PRO) return;

  // Check if we have already included the joints for this residue

  // testing stupid bug
//  prot->const_write("prot_in_add_dep_0.pdb");

  const residue_joints *c_joints;
  if(!A_joints.add_residue_info(prot, res, &c_joints)) return;
  res_to_verts_mci map_iter = A_atoms_dep_on_joints.find(res);
  if(map_iter != A_atoms_dep_on_joints.end()) return;

  // Store dependencies of atomic positions wrt joints
  const joint_deps_t joints_info = A_joints.get_joints_info(res->name);
  pos_dep_on_joints atoms_struct;
  
  // Block size is nrows by 3*(num_atoms - 5)
  if(res->name == ILE){
    atom_vci a = res->atoms_begin + 5;
    // CG1
    atoms_struct.verts.push_back(a->pos);
    atoms_struct.num_joints_used.push_back(1);
    ++a; 
    // CG2
    atoms_struct.verts.push_back(a->pos);
    atoms_struct.num_joints_used.push_back(1);
    ++a;
    // CD1
    atoms_struct.verts.push_back(a->pos);
    atoms_struct.num_joints_used.push_back(2);
  }else{
    int num_joints = 1;
    for(atom_vci a = res->atoms_begin + 5; a < res->atoms_end; ++a){
      atoms_struct.verts.push_back(a->pos);
      atoms_struct.num_joints_used.push_back(num_joints);
      if(num_joints < joints_info.nrows) ++num_joints;
    }
  }
  A_atoms_dep_on_joints[res] = atoms_struct;
} 

bool
SurfDepsJoints::determine_pt_dependence(my_float_t *pt, 
                                        const PDBStructure *prot,
                                        residue_vci res, atom_vci closest_atom)
{
  // At the present we are keeping these atoms fixed
  if(closest_atom->name == N || closest_atom->name == CA ||
     closest_atom->name == C || closest_atom->name == O ||
     (closest_atom->name == CB && res->name == ALA) ||
     res->name == PRO){
//     std::cout << "point is closest to a main chain atom\n";
    return false;
  }

  // testing stupid bug
//  prot->const_write("prot_in_add_pt_dep_0.pdb");

  const residue_joints *c_joints;
  if(!A_joints.add_residue_info(prot, res, &c_joints)) return false;

  // Currently the maximum number of joint is 4 -- lysine
  my_float_t axes_of_rotation[12];
  int num_joints;
  c_joints->compute_current_axes_of_rotation(axes_of_rotation, &num_joints,
                                             12);
#if 0
  std::cout << "Processing point for residue (" << res->chainID << ") " 
            << PDB_residues::residue_to_string(res->name) << res->number
            << std::endl;
  std::cout << "number of joints in residue: " << c_joints->num_joints() 
            << "\n";
  std::cout << "num labels: " << c_joints->joint_labels.size() 
            << std::endl;
  std::cout << "num joint centers: "<< c_joints->joint_centers.size()
            << std::endl;
#endif
  int num_joints_used = 0;
  for(int i = 0; i < c_joints->num_joints(); ++i, ++num_joints_used){
#if 0
    std::cout << "i: " <<  i << std::endl;
    
    std::cout << "Joint label: " 
              << PDB_residues::atom_to_string(c_joints->joint_labels[i]) << "\n";
#endif
    if(closest_atom->name == c_joints->joint_labels[i]){

      // determine plane at joint
      // Normal is the joint axis and the joint center is point on plane
      geometry::plane_t joint_plane(&axes_of_rotation[3*i],
                                    c_joints->joint_centers[i+1]);

      // check if pt is above or below the plane
      my_float_t d = joint_plane.signed_dist(pt);

      // If below, then it is not dependent on this joint, 
      // otherwise it is -- in either case we do not need to
      // consider any more joints
      if(d > 0) ++num_joints_used;
      break;
    }
  }
//  std::cout << "Num joints used: " << num_joints_used << "\n";

  // add info 
  if(num_joints_used){
    res_to_verts_mi pts2joints_iter = A_verts_dep_on_joints.find(res);
    if(pts2joints_iter == A_verts_dep_on_joints.end()){
      pos_dep_on_joints tmp;
      tmp.verts.push_back(pt);
      tmp.num_joints_used.push_back(num_joints_used);
      A_verts_dep_on_joints[res] = tmp;
    }else{
      pts2joints_iter->second.verts.push_back(pt);
      pts2joints_iter->second.num_joints_used.push_back(num_joints_used);
    }
    return true;
  }
  return false;
}

bool
SurfDepsJoints::determine_hbond_dep(HbondSurfaces<model_hbond_surf_t>::surfaces_vci cap, 
                                    const PDBStructure *prot, 
                                    residue_vci *res_out)
{
  residue_vci res = prot->get_residue(cap->atom);
  if(res == PDBStructure::NULL_RESIDUE_VCI) return false;

  if(cap->atom->name == O or cap->atom->name == N){
    std::cout << "Not adding HB cap for residue (" << res->chainID << ") " 
              << PDB_residues::residue_to_string(res->name) 
              << res->number
	      << " -- currently phi and psi angles are fixed"
              << std::endl;
    return false;
  }


  *res_out = res;
  const residue_joints *c_joints;
  std::cout << "Adding HB cap for residue (" << res->chainID << ") " 
            << PDB_residues::residue_to_string(res->name) 
            << res->number
            << std::endl;
  if(!A_joints.add_residue_info(prot, res, &c_joints)) return false;

  // Currently the maximum number of joint is 4 -- lysine
  my_float_t axes_of_rotation[12];
  int num_joints;
  c_joints->compute_current_axes_of_rotation(axes_of_rotation, &num_joints,
                                             12);

  // For each point in the cap add it to the "vector" 
  //std::map<residue_vci, pos_dep_on_joints> A_hbond_pts_dep_on_joints;
  res_to_verts_mi pts2joints_iter = A_hbond_pts_dep_on_joints.find(res);

  if(pts2joints_iter == A_hbond_pts_dep_on_joints.end()){
    pos_dep_on_joints tmp;
    A_hbond_pts_dep_on_joints.insert(std::pair<residue_vci, 
                                               pos_dep_on_joints>(res, tmp));
    pts2joints_iter = A_hbond_pts_dep_on_joints.find(res);
  }
  for(geometry::vertex_vci V = cap->verts_begin(); V < cap->verts_end(); ++V){
    pts2joints_iter->second.verts.push_back(V->pos);
    pts2joints_iter->second.normals.push_back(V->dir);
    pts2joints_iter->second.num_joints_used.push_back(num_joints); 
  }

  return true;
}

bool
SurfDepsJoints::determine_hbond_dep(hbond_ideal_pt_vci ideal_pt,
                                    const PDBStructure *prot, 
                                    residue_vci *res_out)
{
  residue_vci res = prot->get_residue(ideal_pt->atom);
  if(res == PDBStructure::NULL_RESIDUE_VCI){ 
    std::cout << "inside determine hbond dep (" << ideal_pt->atom->chainID 
              << ") " << PDB_residues::residue_to_string(ideal_pt->atom->res) 
              << ideal_pt->atom->res_num << " -- cannot find this?"
              << std::endl;
    return false;
  }

  if(ideal_pt->atom->name == O or ideal_pt->atom->name == N){
    std::cout << "Not adding HB points for residue (" << res->chainID << ") " 
              << PDB_residues::residue_to_string(res->name) 
              << res->number
	      << " -- currently phi and psi angles are fixed"
              << std::endl;
    return false;
  }

  // testing stupid bug
//  prot->const_write("prot_in_add_hbond_dep_0.pdb");


  *res_out = res;
  const residue_joints *c_joints;
  std::cout << "Adding HB cap for residue (" << res->chainID << ") " 
            << PDB_residues::residue_to_string(res->name) 
            << res->number
            << std::endl;
  if(!A_joints.add_residue_info(prot, res, &c_joints)) return false;

  // Currently the maximum number of joint is 4 -- lysine
  my_float_t axes_of_rotation[12];
  int num_joints;
  c_joints->compute_current_axes_of_rotation(axes_of_rotation, &num_joints,
                                             12);

  // For each point associated with the ideal point add it to the "vector" 
  res_to_verts_mi pts2joints_iter = A_hbond_pts_dep_on_joints.find(res);

  if(pts2joints_iter == A_hbond_pts_dep_on_joints.end()){
    pos_dep_on_joints tmp;
    A_hbond_pts_dep_on_joints.insert(std::pair<residue_vci, 
                                               pos_dep_on_joints>(res, tmp));
    pts2joints_iter = A_hbond_pts_dep_on_joints.find(res);
  }

  hbond_fit_pt_vci fit_pt = ideal_pt->fit_pts_beg;
  for( ; fit_pt < ideal_pt->fit_pts_end; ++fit_pt){
    pts2joints_iter->second.verts.push_back(fit_pt->pos);
    pts2joints_iter->second.normals.push_back(fit_pt->dir);
    pts2joints_iter->second.num_joints_used.push_back(num_joints); 
  }
  
  return true;
}

void
SurfDepsJoints::compute_portion_of_J(res_to_verts_mi pts_begin,
                                     res_to_verts_mi pts_end)
{
#ifdef __GLIBCXX__
  __glibcxx_requires_valid_range(pts_begin, pts_end);
#else
  __glibcpp_requires_valid_range(pts_begin, pts_end);
#endif

  // Compute each block of this portion of J
  prot_joint_dep::res_to_joints_mci res2j_iter = A_joints.joints_begin();
  res_to_verts_mi res2pt_iter;
  for(res2pt_iter = pts_begin; res2pt_iter != pts_end; ++res2pt_iter){
//    std::cout << "\n\n\n new residue \n\n\n";

    // Compute the current unit axes of rotation
    // Currently the maximum number of joint is 4 -- lysine
    my_float_t axes_of_rotation[12];
    int num_joints;
    const residue_joints &joints = res2j_iter->second;
    joints.compute_current_axes_of_rotation(axes_of_rotation, &num_joints,
                                            12);

    // For each joint, compute entries for Jacobian that depend on that
    // joint
    pos_dep_on_joints *dep_points = &(res2pt_iter->second);
    my_float_t *block_row_start = dep_points->mutable_J_block();
    size_t n = dep_points->ncols() * dep_points->nrows();
    std::fill(block_row_start, block_row_start + n, 0.0);
    for(int i = 0; i < joints.num_joints(); ++i){
//      std::cout << "\n new row \n";
      const my_float_t *joint_center = joints.joint_centers[i];
      const my_float_t *rot_axis = axes_of_rotation + 3*i;
//      std::cout << "block row start: " << block_row_start << "\n";
//      std::cout << "rot_axis: " << rot_axis[0] << " " << rot_axis[1] << " " << rot_axis[2] << "\n";
      for(size_t j = 0; j < dep_points->num_joints_used.size(); ++j){
//        std::cout << "[i,j]: " << i << ", " << j << "\n";
        // If the current point does not depend on the ith joint, skip it
        if(dep_points->num_joints_used[j] < i+1) continue;

        //const my_float_t *end_effector = dep_points->verts[j];
        my_float_t V[3];
        vector(3, dep_points->verts[j], joint_center, V);
/* NOTE if we have normals we will need to handle them here */

          // Assume rot_axis is a unit vector
        cross(V, rot_axis, block_row_start + 3*j);
#if 0
        std::cout << "V: " << V[0] << " " << V[1] << " " << V[2] << "\n";
        const my_float_t *b = block_row_start + 3*j;
        std::cout << "b addr: " << b << "\n";
        std::cout << "b: " << b[0] << " " << b[1] << " " << b[2] << "\n";
#endif
      }

      block_row_start += dep_points->ncols();
    }
    // reset block_row_start in case someone tries to use it in the 
    // future
    block_row_start = 0;

    ++res2j_iter;
  }
}


void
SurfDepsJoints::write_mat(std::ostream &out, const mat_type type) const
{
  if(A_have_surf)
    write_mat_component(A_verts_dep_on_joints.begin(), 
                        A_verts_dep_on_joints.end(), type, "V", out);
  write_mat_component(A_atoms_dep_on_joints.begin(), 
                      A_atoms_dep_on_joints.end(), type, "A", out);
  if(A_have_hbond_surfs)
    write_mat_component(A_hbond_pts_dep_on_joints.begin(), 
                        A_hbond_pts_dep_on_joints.end(), type, "HB", out);
}

void
SurfDepsJoints::write_mat_component(res_to_verts_mci pts_begin,
                                    res_to_verts_mci pts_end, 
                                    const mat_type type,
                                    const std::string name, std::ostream &out)
const
{
  res_to_verts_mci res2pt_iter;
  int cnt = 0;
  for(res2pt_iter = pts_begin; res2pt_iter != pts_end; ++res2pt_iter, ++cnt){
    out << name << cnt << " = [";
    const pos_dep_on_joints *dep_points = &(res2pt_iter->second);
    const my_float_t *block_row_start = 0;
    if(type == JACOBIAN) block_row_start = dep_points->J_block();
    else if(type == PINV_JACOBIAN) 
      block_row_start = dep_points->pinv_Jtrans_block();
    const size_t nrows = dep_points->nrows();
    const size_t ncols = dep_points->ncols();
    for(size_t i = 0; i < nrows; ++i){
      out << "[";
      for(size_t j = 0; j < ncols; ++j) 
        out << block_row_start[i*ncols + j] << ", ";
      out << "],\n";
    }
    out << "]\n";
  }
}

void
SurfDepsJoints::add_mol_surf_vertices(geometry::TransformableTrimesh *surf, 
                                      PDBStructure *prot)
{
  A_have_surf = true;
  const my_float_t *V_end =
    surf->vertices_begin() + 3*surf->number_of_vertices();
  for(my_float_t *V = surf->mutable_vertices_begin(); V < V_end; V += 3){
    // 1) find closest protein atom to the given point
    my_float_t d;
    atom_vci closest_atom = prot->closest_atom(V, &d);

    if(closest_atom->res == PDB_METAL){
      std::cout << "We Have a Metal!!!!!!!!!!!!!" << std::endl;
      continue;
    }

    if(closest_atom->res == HOH){
      std::cout << "Why are we mapping to WATER?!!!!!!!!!!!!!" << std::endl;
      continue;
    }

    if(closest_atom->res == 0 || closest_atom->res > 21){
      std::cout << "I don't know how to handle this residue type DAVE!!!!!!!!!!!!" << std::endl;
      continue;
    }

    // 2) look up the residue for the closest atom
    residue_vci close_res = prot->get_residue(closest_atom);

    if(closest_atom == atom_t::NULL_ATOM_VCI){
      std::cout << "Point doesn't have a close atom!!!!!!!!!!!!" << std::endl;;
      continue;
    }else if(close_res == PDBStructure::NULL_RESIDUE_VCI){
      std::cout << "atom doesn't have a residue!!!!!!!!!!!!" << std::endl;;
      continue;
    }

    // 3) get chain of joints for the point 
    // 4) Get the chain of joints for atoms in the residue
    if(determine_pt_dependence(V, prot, close_res, closest_atom))
      add_dependent_atoms(prot, close_res);
  }
}


void
SurfDepsJoints::add_hbond_points(ModelHbondSurfaces *hbond_surfaces,
                                 HbondPoints *hbond_pts,
                                 const PDBStructure *prot)
{
  if(hbond_surfaces){
    std::map<ModelHbondSurfaces::m_surf_vci, bool> hbond_caps_map;

    A_have_hbond_surfs = true;
    // For each cap, we know the atom already
    HbondSurfaces<model_hbond_surf_t>::surfaces_vi S;
    S = hbond_surfaces->surf_caps_begin();
    for( ; S != hbond_surfaces->surf_caps_end(); ++S){
      residue_vci res;

      // Determine the chain of joints for the cap points
      // If cap is added, then add the mobile side-chain atoms to the map
      if(determine_hbond_dep(S, prot, &res)){
        add_dependent_atoms(prot, res); 
        hbond_caps_map[S] = true;
      }
    }

    // move sorted caps to vector (they may already be sorted -- 
    // cannot remember at this point)
    std::map<ModelHbondSurfaces::m_surf_vci, bool>::iterator caps_map_iter;
    caps_map_iter = hbond_caps_map.begin();
    for( ; caps_map_iter != hbond_caps_map.end(); ++caps_map_iter){
      A_mobile_caps.push_back(caps_map_iter->first);
    }
  }else{
    A_have_hbond_surfs = false;    
  
    hbond_ideal_pt_vci HB_pt = hbond_pts->ideal_pts_beg();
    for( ; HB_pt != hbond_pts->ideal_pts_end(); ++HB_pt){
      residue_vci res;

      // Determine the chain of joints for the cap points
      // If cap is added, then add the mobile side-chain atoms to the map
      if(determine_hbond_dep(HB_pt, prot, &res)){
        add_dependent_atoms(prot, res); 
      }
    }
  }
}

bool
SurfDepsJoints::check_containers()
{
  // Try to make things easier by ensuring that each of the 3 maps has the 
  // same residues
  res_to_verts_mi res2V_iter = A_verts_dep_on_joints.begin();
  res_to_verts_mi res2HB_iter = A_hbond_pts_dep_on_joints.begin();
  res_to_verts_mi res2A_iter = A_atoms_dep_on_joints.begin();
  for( ; res2A_iter != A_atoms_dep_on_joints.end(); ++res2A_iter){
    if(A_have_surf){
      if(res2A_iter->first < res2V_iter->first ||
         res2V_iter == A_verts_dep_on_joints.end()){
        pos_dep_on_joints tmp;
        A_verts_dep_on_joints[res2A_iter->first] = tmp;
        res2V_iter = A_verts_dep_on_joints.find(res2A_iter->first);
      }else if(res2A_iter->first > res2V_iter->first){
        std::cerr << "Somehow have more residues that move based on vertices "
                  << "than atoms\n";
        return false;
      }
      ++res2V_iter;
    }
    
//    if(A_have_hbond_surfs){
    if(res2A_iter->first < res2HB_iter->first ||
       res2HB_iter == A_hbond_pts_dep_on_joints.end()){
      pos_dep_on_joints tmp;
      A_hbond_pts_dep_on_joints[res2A_iter->first] = tmp;
      res2HB_iter = A_hbond_pts_dep_on_joints.find(res2A_iter->first);
    }else if(res2A_iter->first > res2HB_iter->first){
      std::cerr << "Somehow have more residues that move based on hbond pts "
                << "than atoms\n";
      return false;
    }
    ++res2HB_iter;
//    }
  } 
  std::cout << "Size of surf vert map: " << A_verts_dep_on_joints.size();
  std::cout << "Size of atoms map: " << A_atoms_dep_on_joints.size();
  std::cout << "Size of hbond pts map: " << A_hbond_pts_dep_on_joints.size();

  return true;
}

