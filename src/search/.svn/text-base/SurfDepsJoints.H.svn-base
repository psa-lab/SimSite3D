#include <prot_joint_dep.H>
#include <TransformableTrimesh.H>
#include <sphere.H>
#include <ModelHbondSurfaces.H>
#include <HbondPoints.H>

namespace ASCbase{

//! A simple data class for a block of a block matrix
class pos_dep_on_joints{
public:

  pos_dep_on_joints()
  {
    init();
  }

  pos_dep_on_joints(const pos_dep_on_joints& other)
  { 
    init();
    do_copy(other); 
  }

  ~pos_dep_on_joints()
  {
    if(A_J_block) delete [](A_J_block);
    if(A_pinv_Jtrans_block) delete [](A_pinv_Jtrans_block);
    init();
  }

  const pos_dep_on_joints&
  operator=(const pos_dep_on_joints& other)
  {
    if(&other != this){
      do_copy(other);
    }
    return *this;
  }

  my_float_t*
  mutable_J_block()
  {
//    std::cout << "Need to design a method to handle only verts or verts "
//              << "and normals\n";
    if(A_J_block == 0){
      A_ncols = 3*(verts.size());
      A_J_block = new my_float_t[A_nrows * A_ncols];
    }
//    std::cout << "J block size: [" << A_nrows << " x " << A_ncols << "]\n";
    return A_J_block;
  }

  const my_float_t*
  J_block() const
  { return A_J_block; }
    

  my_float_t*
  mutable_pinv_Jtrans_block()
  {
    if(A_pinv_Jtrans_block == 0){
      A_pinv_Jtrans_block = new my_float_t[A_nrows * 3*verts.size()];
    }
    return A_pinv_Jtrans_block;
  }

  const my_float_t*
  pinv_Jtrans_block() const
  { return A_pinv_Jtrans_block; }
 
  void
  set_nrows(size_t nrows_in)
  { A_nrows = nrows_in; }

  const size_t
  nrows() const
  { return A_nrows; }

  const size_t
  ncols() const
  { return A_ncols; }

  //  int ncols; -- num cols is given by 3*verts.size()
  std::vector< my_float_t* > verts;
  std::vector< my_float_t* > normals;
  std::vector< int > num_joints_used;

private:

  void
  init()
  {
    A_nrows = 0;
    A_ncols = 0;
    A_J_block = 0;
    A_pinv_Jtrans_block = 0;
  }

  //! Do a copy of the other class instance
  /*! NOTE: block is a pointer to the first entry of a block with in a
   * large block matrix.  Thus this copy does not make a deep copy of the
   * block matrix nor the entire matrix
   */
  void
  do_copy(const pos_dep_on_joints& other)
  {
    A_nrows = other.A_nrows;
    A_ncols = other.A_ncols;
    verts.clear();
    verts = other.verts; 
    normals.clear();
    normals = other.normals;
    num_joints_used.clear();
    num_joints_used = other.num_joints_used;

    // Copy J block
    if(other.A_J_block){
      if(A_J_block) delete [](A_J_block);
      const int n = A_nrows * A_ncols;
      A_J_block = new my_float_t[n];
      std::copy(other.A_J_block, other.A_J_block + n, A_J_block);
    }else A_J_block = 0;

    // Copy pinv(J)^T block
    if(other.A_pinv_Jtrans_block){
      if(A_pinv_Jtrans_block) delete [](A_pinv_Jtrans_block);
      const int n = A_nrows * A_ncols;
      A_pinv_Jtrans_block = new my_float_t[n];
      std::copy(other.A_pinv_Jtrans_block, other.A_pinv_Jtrans_block + n,
                A_pinv_Jtrans_block);
    }else A_pinv_Jtrans_block = 0;
  }

  size_t A_nrows;
  size_t A_ncols;
  my_float_t *A_J_block;
  my_float_t *A_pinv_Jtrans_block;

};


//! Atoms (joints) and surface vertices that depend on rotatable joints
class SurfDepsJoints{
public:
  typedef std::map<residue_vci, pos_dep_on_joints> res_to_verts_map;
  typedef res_to_verts_map::iterator res_to_verts_mi;
  typedef res_to_verts_map::const_iterator res_to_verts_mci;

  typedef enum{
    JACOBIAN,
    PINV_JACOBIAN
  }mat_type;

  SurfDepsJoints(prot_joint_dep &joints, PDBStructure *prot,
                 geometry::TransformableTrimesh *surf, 
                 ModelHbondSurfaces* hbond_surfaces, HbondPoints *hbond_pts);

  ~SurfDepsJoints()
  {
    if(A_JJtranspose) delete [](A_JJtranspose);
    init();
  }

  //! Jacobian that give dependence of joint angles on small motions of vertices
  void
  compute_dihedrals_and_verts_J()
  {
    // compute the columns of J corresponding to surface vertex positions
    if(A_have_surf)
      compute_portion_of_J(A_verts_dep_on_joints.begin(), 
                           A_verts_dep_on_joints.end());
    if(A_have_hbond_surfs)
      compute_portion_of_J(A_hbond_pts_dep_on_joints.begin(), 
                           A_hbond_pts_dep_on_joints.end());
    // compute the columns of J corresponding to protein atom positions
    compute_portion_of_J(A_atoms_dep_on_joints.begin(), 
                         A_atoms_dep_on_joints.end());
  }

  void 
  compute_JJtranspose()
  {
    std::fill(A_JJtranspose, A_JJtranspose + A_J_nrows*A_J_nrows, 0.0);
    if(A_have_surf)
      compute_portion_of_JJtranspose(A_verts_dep_on_joints.begin(), 
                                     A_verts_dep_on_joints.end(), 
                                     A_JJtranspose);
    if(A_have_hbond_surfs)
      compute_portion_of_JJtranspose(A_hbond_pts_dep_on_joints.begin(), 
                                     A_hbond_pts_dep_on_joints.end(), 
                                     A_JJtranspose);
    compute_portion_of_JJtranspose(A_atoms_dep_on_joints.begin(), 
                                   A_atoms_dep_on_joints.end(), A_JJtranspose);
  }

  void
  update_positions(const my_float_t *delta_q);

  void
  compute_pinvJ_transpose();

  void
  write_point_correspondences(std::ostream &out) const;

  //! Compute the dihedral angles for those sidechains that are not fixed
  void
  compute_dihedral_angles(std::vector< std::vector<my_float_t> > *angles) const;

#if 0
  void 
  apply_pinvJ_to_grad(const my_float_t *grad_F, my_float_t *Q);
#endif

  //! Wrapper to compute_portion_pinvJtrans_times_grad for vertices
  /*!
   * Assumptions:
   *   len(grad_F) >= 3*(A_J_cols_in_blocks[0])
   *   len(Q) >= A_J_nrows
   */
  inline void
  apply_verts_pinvJ_to_grad(const my_float_t *grad_F, my_float_t *Q)
  {
    compute_portion_pinvJtrans_times_grad(A_verts_dep_on_joints.begin(),
                                          A_verts_dep_on_joints.end(),
                                          grad_F, Q);
  }

  //! Wrapper to compute_portion_pinvJtrans_times_grad for hbond cap pts
  /*!
   * Assumptions:
   *   len(grad_F) >= 3*(A_J_cols_in_blocks[1])
   *   len(Q) >= A_J_nrows
   */
  inline void
  apply_hbond_pts_pinvJ_to_grad(const my_float_t *grad_F, my_float_t *Q)
  {
    compute_portion_pinvJtrans_times_grad(A_hbond_pts_dep_on_joints.begin(),
                                          A_hbond_pts_dep_on_joints.end(),
                                          grad_F, Q);
  }

#if 0
  //! Wrapper to compute_portion_pinvJtrans_times_grad for atoms
  /*!
   * Assumptions:
   *   len(grad_F) >= 3*(A_J_cols_in_blocks[2])
   *   len(Q) >= A_J_nrows
   */
  inline void
  apply_atoms_pinvJ_to_grad(const my_float_t *grad_F, my_float_t *Q)
  {
    compute_portion_pinvJtrans_times_grad(A_atoms_dep_on_joints.begin(),
                                          A_atoms_dep_on_joints.end(),
                                          grad_F, Q);
  }
#endif

  //! Averaging method to handle overlapping atoms
  void
  apply_atoms_pinvJ_to_grad(const my_float_t *grad_F, my_float_t *Q, 
                            const std::vector<bool> &overlap_vec);

  //! Get the number of rows in the Jacobian
  const int
  nrows_J() const
  { return A_J_nrows; }

  //! Get the number of columns in the Jacobian
  const int
  ncols_J() const
  { return A_J_ncols; }

  //! Get the number of atom columns in the Jacobian
  const int
  num_atom_cols_J() const
  { return 3*A_J_cols_in_blocks[2]; }

  //! Get the number of vertex columns in the Jacobian
  const int
  num_vert_cols_J() const
  { return 3*A_J_cols_in_blocks[0]; }

  //! Get the number of hbond pts columns in the Jacobian
  const int
  num_hbond_pts_cols_J() const
  { return 3*A_J_cols_in_blocks[1]; }

  res_to_verts_mci
  atoms_dep_on_joints_begin() const
  { return  A_atoms_dep_on_joints.begin(); }

  res_to_verts_mci
  atoms_dep_on_joints_end() const
  { return  A_atoms_dep_on_joints.end(); }

  std::vector<int>::const_iterator
  vertex_idz_begin() const
  { return A_vert_idz.begin(); }

  std::vector<int>::const_iterator
  vertex_idz_end() const
  { return A_vert_idz.end(); }

  void
  write_mat(std::ostream &out, const mat_type type) const;

  bool
  fail() const
  { return A_fail; }

  std::vector<ModelHbondSurfaces::m_surf_vci>::const_iterator
  mobile_cap_iters_begin() const
  { return A_mobile_caps.begin(); }

  std::vector<ModelHbondSurfaces::m_surf_vci>::const_iterator
  mobile_cap_iters_end() const
  { return A_mobile_caps.end(); }

  void
  stop_rotation_given_atoms(const std::vector<bool> &overlap_vec,
                                          std::vector<bool> *fix_Q);

  //! Compute the atomic rmsd for the query and dataset residues 
  /*! WARNING! at the present this is only valid for proteins with the
   * same residues ...
   *
   * Note: At the present it is very inefficient, but I don't want
   * to miss anything and don't have the time to code it correctly.
   * The problem is PDBStructure does not place any guarentees on the
   * ordering of residues with respect to chainID, residue number, and iCode
   */
  my_float_t
  site_atomic_rmsd(const PDBStructure &dset_prot) const
  {
    my_float_t rmsd = 0;
    size_t num_atoms = 0;
    std::map<residue_vci, bool>::const_iterator q_res_iter;
    q_res_iter = A_site_residues.begin();
    for( ; q_res_iter != A_site_residues.end(); ++q_res_iter){
      const residue_vci q_res = q_res_iter->first;
      residue_vci dset_res = dset_prot.residues_begin();
      bool found_res = false;
      for( ; dset_res < dset_prot.residues_end() && !found_res; ++dset_res){
        if(dset_res->number == q_res->number &&
           dset_res->chainID == q_res->chainID &&
           dset_res->icode == q_res->icode){
          atom_vci Q = q_res->atoms_begin;
          atom_vci D = dset_res->atoms_begin;
          for( ; Q < q_res->atoms_end && D < dset_res->atoms_end; ++Q, ++D){
            rmsd += dist_squared(Q->pos, D->pos);
            ++num_atoms; 
          }
          found_res = true;
        }
      }
    }
    
    if(num_atoms == 0) return my_float_max;
    rmsd = std::sqrt(rmsd / num_atoms);
    return rmsd;
  }
  

private:
  //! Add the vertices for the given molecular surface to the points that can 
  //! move
  void
  add_mol_surf_vertices(geometry::TransformableTrimesh *surf, 
                        PDBStructure *prot);

  //! Add the points for the hbond groups in the binding site
  /*! If hbond_surfaces is not null, use the given hbond_surfaces, otherwise
   *  use the hbond points.  
   *
   *  NOTE: hbond_surfaces will be used to move side chains, hbond_pts are 
   *  NOT used to move sidechains but will "go along for the ride" if 
   *  hbond_surfaces are not used.
   */
  void
  add_hbond_points(ModelHbondSurfaces *hbond_surfaces, HbondPoints *hbond_pts,
                   const PDBStructure *prot);

  //! Check the 3 containers to make sure they have the same residues
  bool
  check_containers();

  //! Add the mobile atoms for those residues that can move
  void
  add_dependent_atoms(const PDBStructure *prot, residue_vci res);

  //! Add dependencies for each surface patch vertex
  bool
  determine_pt_dependence(my_float_t *pt, const PDBStructure *prot,
                          residue_vci res, atom_vci closest_atom);

  //! Add joint dependencies for the given hbond cap
  bool
  determine_hbond_dep(HbondSurfaces<model_hbond_surf_t>::surfaces_vci cap, 
                      const PDBStructure *prot, residue_vci *res_out);

  //! Add joint dependencies for the points associated with the hbond ideal 
  //! point
  bool
  determine_hbond_dep(hbond_ideal_pt_vci ideal_pt, const PDBStructure *prot, 
                      residue_vci *res_out);

  //! Compute the large portion of the Jacobian that pertains to the 
  //! given map
  void
  compute_portion_of_J(res_to_verts_mi pts_begin, res_to_verts_mi pts_end);

  //! Compute the Jacobian times its transpose for the given iterators
  void 
  compute_portion_of_JJtranspose(res_to_verts_mi pts_begin, 
                                 res_to_verts_mi pts_end,
                                 my_float_t *JJtrans);

  void
  compute_portion_of_pinvJtranspose(res_to_verts_mi pts_begin, 
                                    res_to_verts_mi pts_end, 
                                    my_float_t *JJtrans_inv);

  //! Apply the pinv of J to the gradient to determine change in Q due to 
  //! the type of actions
  /*!
   * Note: this function does not initialize Q as we may wish to accumulate
   * contributions to the "force" on Q
   */
  void 
  compute_portion_pinvJtrans_times_grad(res_to_verts_mi pts_begin,
                                        res_to_verts_mi pts_end,
                                        const my_float_t *grad_F, 
                                        my_float_t *Q);
  
  void
  write_mat_component(res_to_verts_mci pts_begin, res_to_verts_mci pts_end, 
                      const mat_type type, const std::string name, 
                      std::ostream &out) const;

  //! Rotate the dependent points using R about the axis of the joint 
  void
  rotate_points(pos_dep_on_joints *points, const my_float_t *R,
                const my_float_t *joint_center, const int joint_number)
  {
    for(size_t i = 0; i < points->num_joints_used.size(); ++i){
      if(points->num_joints_used[i] < joint_number + 1) continue;

      // Set joint center as the origin
      my_float_t tmp[3];
      std::copy(points->verts[i], points->verts[i] + 3, tmp);
      my_axpy(3, -1.0, joint_center, 1, tmp, 1);

      // rotate the point about the axis and add the joint offset
      std::copy(joint_center, joint_center + 3, points->verts[i]);
      my_gemv(3, 3, 1.0, R, 3, tmp, 1, 1.0, points->verts[i], 1);
    }
  }

  //! Initialize class variables
  void
  init()
  {
    A_have_surf = false;
    A_have_hbond_surfs = false;
    A_fail = true;
    A_J_nrows = 0;
    A_J_ncols = 0;
    A_JJtranspose = 0;
    std::fill(A_J_cols_in_blocks, A_J_cols_in_blocks + 3, 0);
  }

  bool A_have_surf;
  bool A_have_hbond_surfs;
  bool A_fail;

  prot_joint_dep &A_joints;
  std::vector<int> A_vert_idz;

  std::map<residue_vci, pos_dep_on_joints> A_verts_dep_on_joints; //!< Mobile molecular surface vertices
  std::map<residue_vci, pos_dep_on_joints> A_atoms_dep_on_joints; //!< Mobile protein atoms
  std::map<residue_vci, pos_dep_on_joints> A_hbond_pts_dep_on_joints; //!< Mobile hydrogen bond points (either from HbondPoints or HbondSurfaces)
  std::vector<ModelHbondSurfaces::m_surf_vci> A_mobile_caps;

  std::map<residue_vci, bool> A_site_residues;  //!< Used in testing ArtSurf binding site atomic RMSD 

  int A_J_nrows;  //!< Number of rows in the full Jacobian
  int A_J_ncols;
  int A_J_cols_in_blocks[3];
  my_float_t *A_JJtranspose;  //!< Work space variable -- do not use ouside of class 
};
}
