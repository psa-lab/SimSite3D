/******************************************************************************
 * Copyright (c) 2006,2007, Michigan State University (MSU) Board of Trustees.
 *   All rights reserved.
 *
 * This file is part of the ASCbase Software project.
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

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <ctime>
#include <basics.H>
#include <param_tools.H>
#include <mat_ops.H>
#include <IK_tests.H>
#include <NoTier1Score.H>
#include <SurfDepsJoints.H>
#include <HbondGeometry.H>
#include <ICPSites.H>

using namespace ASCbase;
const my_float_t IK_tests::MAX_ATOM_DIST_TO_CHECK_VDW_OVERLAP = 4.5;
const std::string IK_tests::A_fname = "IK_tests.C";

IK_tests::IK_tests(const SearchParameters& args_in)
{
  if(args_in.status() != BaseParameters::READY) return;

  A_site_score_p = 0;
  A_uses_hbond_surfaces = false;
  A_uses_site_surface = false;
  A_ICP_surf_pt_W = -1.0;
  A_ICP_hb_cap_pt_W = -1.0;
  A_hbond_caps_ICP = false;

  if(args_in.use_hbond_surfaces){
    std::cout << "using hbond caps score!!!!!\n";
    A_site_score_p = new HbondSurfacesScore();

    // do nothing -- use default values
    if(args_in.fine_tune_ratio == 0.0); 
    // use on hbond caps
    else if(args_in.fine_tune_ratio < 0.0){
      A_ICP_surf_pt_W = 0.0;  
      A_ICP_hb_cap_pt_W = 1.0;  
      A_hbond_caps_ICP = true;
    }else{
      if(args_in.fine_tune_ratio < 1.0){
        A_ICP_surf_pt_W = args_in.fine_tune_ratio;
        A_ICP_hb_cap_pt_W = 1.0;  
      }else if(args_in.fine_tune_ratio > 1.0){
        A_ICP_surf_pt_W = 0.0;  
        A_ICP_hb_cap_pt_W = 1.0 / (args_in.fine_tune_ratio);  
      }else{
        A_ICP_surf_pt_W = 1.0;  
        A_ICP_hb_cap_pt_W = 1.0;  
      }
      A_hbond_caps_ICP = true;
    }

    // do nothing -- use default values
    if(args_in.IK_type == SearchParameters::IK_SURFACES) 
      A_uses_site_surface = true;
    else if(args_in.IK_type == SearchParameters::IK_HB_CAPS) 
      A_uses_hbond_surfaces = true;
    else if(args_in.IK_type == SearchParameters::IK_BOTH_REPS){
      A_uses_hbond_surfaces = true;
      A_uses_site_surface = true;
    }

  }else{
    std::cout << "using point and surf score!!!!!\n";
    A_site_score_p = new point_and_surf_score();
  }


  std::cout << "A_uses_hbond_surfaces: " 
            << (A_uses_hbond_surfaces ? "true" : "false") << "\n";
  std::cout << "A_uses_site_surface: " 
            << (A_uses_site_surface ? "true" : "false") << "\n";

  A_site_score_p->set_scale_terms(args_in.scale_terms);
}

IK_tests::~IK_tests()
{
  if(A_site_score_p) delete A_site_score_p;
  A_site_score_p = 0;
}


void
IK_tests::test_moving_stuff(const SearchParameters &args_in)
{
#if 0
  if(fail()){
    warn(A_fname, "run", "Cannot run because an error has occured");
    return;
  }

  std::string dbfile_path, dbfile_struct_id;
  if(!args_in->dbase_file_name.length()){
    std::cerr << "Deformable surface tests require a database file name.\n"
              << "I don't have the time or will at the present to look into\n"
              << "allowing a database of site maps\n\n";
    return;
  }else{
    get_path_and_struct_id(args_in->dbase_file_name, &dbfile_path, 
                           &dbfile_struct_id);
  }
  std::string ofile_prefix = args_in->proj_output + "/" + dbfile_struct_id;
  ofile_prefix += "_surfs";
  if(!dir_exists(ofile_prefix, false) && !my_mkdir(ofile_prefix.c_str(), 0770)){
    std::cerr << "Could not create surf file directory\n";
    return;
  }  
  
  std::string model_path, model_struct_id;
  get_path_and_struct_id(args_in->model_file_name, &model_path, 
                         &model_struct_id);
  ofile_prefix += "/" + model_struct_id + "_surf_";

  // results file -- not the files used to write out the moved mesh
  std::ofstream out;
  if(!open_ofstream(out, args_in->ofname)) return;
  args_in->report_parameters(out);

  Timer my_timer;
  my_timer.start();

  DbaseSitemap dbase_site(dbfile_path, dbfile_struct_id, *args_in,
                          args_in->check_all_triangles,
                          args_in->use_hbond_surfaces);

  
  std::vector<align_w_rigid_align_t> test_align_vec;
  align_w_rigid_align_t test_align;
  test_align_vec.push_back(test_align);
  test_align_vec[0].match_print.resize(model_site->fit_points_size());
  align_w_rigid_align_t tmp_align = test_align_vec.begin();

  ICP_mobile_pts(*A_site_score_p, tmp_align, dbase_site, 100, 1E-06, true);
  run(dbase_site, tmp_align, dbfile_struct_id, out);

////////////////////////////////////////////////////////////////////////////////
  double real, virt, prof;
  my_timer.get(&real, &virt, &prof);
  std::cout << "\nreal time: " << real << "\n";
  std::cout << "virt time: " << virt << "\n";
  std::cout << "prof time: " << prof << "\n";
#endif
}

my_float_t
IK_tests::compute_initial_atomic_rmsd(ModelSitemap *model_site, 
                                      const DbaseSitemap &dbase_site, 
                                      rigid_align_t *align)
{
  // NOTE: it is computationally expensive to recompute all of the
  // correspondences, but easier to code.  At the present, the rmsd only
  // "makes sense" for sites with the same protein -- in the future we can
  // move this to a derived class or eliminate it.

  // Assign object (surface vertices,hbond vertices or points) to joint chains
  prot_joint_dep the_joints;
  SurfDepsJoints 
    query_surf_to_joints(the_joints, &(model_site->bind_site_atoms()), 
                         &(model_site->mutable_binding_site_mesh_handle()), 
			 model_site->mutable_model_hbond_surfaces(),
			 &(model_site->hbond_points()));
  if(query_surf_to_joints.fail()){
    std::cerr << "Failed to initialize the joints class" << std::endl;
    return my_float_max;
  }
  my_float_t atomic_rmsd = 
    query_surf_to_joints.site_atomic_rmsd(dbase_site.bind_site_atoms());
  align->set_site_atomic_rmsd(atomic_rmsd, 0);
  return atomic_rmsd;  
}

void
IK_tests::run(ModelSitemap *model_site, const DbaseSitemap &dbase_site, 
              rigid_align_t *align, const std::string db_struct_id, 
              std::ostream& out)
{
#if 0
  // for writing out surf files
  std::string ofile_prefix =  "./";
  ofile_prefix += dbase_site.get_struct_id() + "_surfs";
  if(!dir_exists(ofile_prefix, false) && !my_mkdir(ofile_prefix.c_str(), 0770)){
    std::cerr << "Could not create surf file directory\n";
    return;
  }  
  ofile_prefix += "/" + model_site->get_struct_id() + "_surf_";
#endif

  std::cout << "model hbond surfaces ptr: " 
            << model_site->mutable_model_hbond_surfaces() << std::endl;

  // Assign object (surface vertices,hbond vertices or points) to joint chains
  prot_joint_dep the_joints;
  SurfDepsJoints 
    query_surf_to_joints(the_joints, &(model_site->bind_site_atoms()), 
                         &(model_site->mutable_binding_site_mesh_handle()), 
			 model_site->mutable_model_hbond_surfaces(),
			 &(model_site->hbond_points()));
  if(query_surf_to_joints.fail()){
    std::cerr << "Failed to initialize the joints class" << std::endl;
    return;
  }

  my_float_t atomic_rmsd = 
    query_surf_to_joints.site_atomic_rmsd(dbase_site.bind_site_atoms());
  align->set_site_atomic_rmsd(atomic_rmsd, 1);
  align->save_preIK_vals();
                                
  //my_float_t prev_RMSE = my_float_max;

  // Compute the gradient of the objective function with respect to different
  // global positions moving
  const int nrows = query_surf_to_joints.nrows_J();
  my_float_t *grad_F_mol_surf_verts = 
    new my_float_t[query_surf_to_joints.num_vert_cols_J()];
  my_float_t *grad_F_atoms = 
    new my_float_t[query_surf_to_joints.num_atom_cols_J()];
  my_float_t *grad_F_hbond_surfs_verts = 
    new my_float_t[query_surf_to_joints.num_hbond_pts_cols_J()];

  my_float_t *q_vals = new my_float_t[nrows];
  my_float_t *surf_q_vals = new my_float_t[nrows];
  my_float_t *hbond_q_vals = new my_float_t[nrows];
  my_float_t *overlap_q_vals = new my_float_t[nrows];

for(int rounds = 0; rounds < 100; ++rounds){
  query_surf_to_joints.compute_dihedrals_and_verts_J();
  //const int ncols = query_surf_to_joints.ncols_J();
  //std::cout << "Number of columns: " << ncols << std::endl;

  //query_surf_to_joints.write_mat(std::cout, SurfDepsJoints::JACOBIAN);


  std::fill(grad_F_mol_surf_verts, 
            grad_F_mol_surf_verts + query_surf_to_joints.num_vert_cols_J(), 
            0.0);
  std::fill(grad_F_atoms, 
            grad_F_atoms + query_surf_to_joints.num_atom_cols_J(), 0.0);
  std::fill(grad_F_hbond_surfs_verts, 
            grad_F_hbond_surfs_verts + 
            query_surf_to_joints.num_hbond_pts_cols_J(), 0.0);

  // 1) Minimize the distance between corresponding surface points
  std::cout << "\nTime to score in IK + ICP\n\n";
  A_site_score_p->score(*model_site, dbase_site, align);
  //A_site_score_p->score(*model_site, dbase_site, test_align_vec.begin());
  // NOTE: surface area is junk since the area of triangles are not 
  // recomputed and triangles could intersect, etc.
  std::cout << "Surf terms: " << align->terms[10] << " "
            << align->terms[11] << " " << align->terms[12] << " "
            << align->terms[13] << "\n";
  if(A_uses_hbond_surfaces)
    std::cout << align->terms[16] << " " << align->terms[18] << "\n";
  else
    std::cout << "hbond terms: " << align->terms[2] << " "
              << align->terms[3] << "\n";

  // We have a problem of book keeping in that we need to have a mapping between
  // these surface points and those surface points that can move
  std::vector<const my_float_t *> q_pts;
  std::vector<const my_float_t *> corr_pts;
  A_site_score_p->correspondences(*model_site, &q_pts, &corr_pts);
  std::cout << "\n\nTime to check correspondences\n\n";
  std::cout << "Number of query surface points: " << q_pts.size() << " "
            << corr_pts.size() << "\n\n";

#if 0
  std::cout << "Corresponding surf pts -- for images\n";
  for(int i = 0 ; i < q_pts.size(); ++i){
    if(q_pts[i] == 0)
      std::cout << i << " NONE \n";
    else
      std::cout << i << " " << q_pts[i][0] << " " << q_pts[i][1] << " "
                << q_pts[i][2] << " " << corr_pts[i][0] << " "
                << corr_pts[i][1] << " " << corr_pts[i][2] << "\n";
  }
#endif

  std::vector<int>::const_iterator vert_idx;
  my_float_t *grad_F_p = grad_F_mol_surf_verts;
  for(vert_idx = query_surf_to_joints.vertex_idz_begin();
      vert_idx < query_surf_to_joints.vertex_idz_end(); 
      ++vert_idx, grad_F_p += 3){
    if(q_pts[*vert_idx] == 0) continue;

//    std::cout << *vert_idx << " " << q_pts[*vert_idx] << " " 
//              << corr_pts[*vert_idx]<< std::endl;

    vector(3, q_pts[*vert_idx], corr_pts[*vert_idx], grad_F_p);
    for(int i = 0; i < 3; ++i) grad_F_p[i] *= 2.0;
#if 0
    std::cout  << "pt idx: " << *vert_idx << "   " 
               << "pt dist: " << std::sqrt(dot(grad_F_p, grad_F_p)) << std::endl;    

#endif

  }

  const my_float_t min_atomic_percent_dist = 0.95;  

/******************************************************************************
/ We plan to try 3 different methods
/ 1) Allow a max of 5% overlap (or whatever low threshold we want) and stop
/    rotations of any joints for which overlapping atoms are end effectors
/ 2) Allow a max of 30% overlap (or whatever reasonable intermediate threshold
/    we weant.  
/    A) Use Sybzi to resolve the "bad" structures at the end
/    B) Reduce the overlap at the end by keeping those atoms "above" the
/       overlaping atoms as fixed.  IN this case we may need several iterations
/       because some atoms may need to move first inorder to remove overlap
/       if one or more atoms is overlapped by more than 1 atom
/ 
******************************************************************************/


  std::ostringstream less_ostr;
//  less_ostr << "Amount of lesser overlap in (A) ";
//  std::cout << "Amount of significant overlap in (A) " ;

  // 2) Minimize the overlap penalty for atoms
  std::vector<bool> overlap_vec(query_surf_to_joints.num_atom_cols_J() / 3);
  std::fill(overlap_vec.begin(), overlap_vec.end(), false);
  std::vector<bool>::iterator overlap_flag = overlap_vec.begin();

  SurfDepsJoints::res_to_verts_mci AJ_iter;
  AJ_iter = query_surf_to_joints.atoms_dep_on_joints_begin();
  grad_F_p = grad_F_atoms;
  for( ; AJ_iter != query_surf_to_joints.atoms_dep_on_joints_end(); ++AJ_iter){
    residue_vci res = AJ_iter->first;
    atom_vci res_atom;
    for(res_atom = res->atoms_begin + 5; res_atom < res->atoms_end; ++res_atom){
      point_storage<atom_t>::bins_vci_type nbrs, moved_nbrs;
      model_site->bind_site_atoms().get_atom_bin(res_atom, &nbrs, &moved_nbrs); 

   
      for(size_t i = 0; i < nbrs->positions.size(); ++i){
        // Ignore atoms from same side chain? -- good enough approx for now
 
        // Still too close (wrt vdw overlap to CA)
        atom_vci nbr = nbrs->assoc_data[i];
        if(nbr->name != O && nbr->name != N && nbr->name != C){
          if(res_atom->res_num == nbr->res_num &&
             res_atom->chainID == nbr->chainID &&
             res_atom->altLoc == nbr->altLoc) continue;
        }

        my_float_t diff[3];
        //vector(3, nbr->pos, res_atom->pos, diff);
        vector(3, res_atom->pos, nbr->pos, diff);
        my_float_t atomic_dist_sq = dot(diff, diff);
/*
      std::cout << "(" << nbr->chainID << ") "
                << PDB_residues::residue_to_string(nbr->res) 
                << nbr->res_num
                << " " << PDB_residues::atom_to_string(nbr->name) << " ";
*/

        // If the two atoms are farther apart than the sum of their vdw 
        // radii, they do not contribute to the energy function or gradient
       
        my_float_t vdw_rad_sum = res_atom->vdw_radius + nbr->vdw_radius;
        if(atomic_dist_sq >= vdw_rad_sum*vdw_rad_sum) continue;
 
        // Don't move atoms that overlap because of potential hbond 
        // interactions
        if(HbondGeometry::is_hbond_interaction(res_atom->act_type,
                                               nbr->act_type) &&
           atomic_dist_sq >= 2.5*2.5) continue;
        else if(HbondGeometry::is_hbond_interaction(res_atom->act_type,
                                                    nbr->act_type))
          std::cout << "hbond distance is: " << sqrt(atomic_dist_sq) << "\n";

        // best direction to move the current atom is along the direction
        // from neighbor to self
        my_float_t atomic_dist = normalize(diff); 
        my_float_t shortest_dist = min_atomic_percent_dist * vdw_rad_sum;

      std::cout << "(" << res->chainID << ") "
                << PDB_residues::residue_to_string(res->name) << res->number 
                << " " << PDB_residues::atom_to_string(res_atom->name) << " -- ";

      std::cout << "   (" << nbr->chainID << ") "
                << PDB_residues::residue_to_string(nbr->res) 
                << nbr->res_num
                << " " << PDB_residues::atom_to_string(nbr->name) << "  "
                << vdw_rad_sum - atomic_dist << " : "
                << round(100.0 - 100.0 * atomic_dist/vdw_rad_sum) << "\%\n";

      // AT the time it seemed like a good idea to move only small amount
      // if there are small overlaps, but we need to compete with the
      // forces moving in the opposite direction 

// CHECK this -- it was previously broken -- it may need to be turned off or
// possibly fixed in a more correct manner
#if 1
        // If the overlap is less than 20% move about 10% of the overlap
        if(atomic_dist >= shortest_dist){
          // atomic_dist is >= 4/5 of vdw sum
          my_axpy(3, atomic_dist - shortest_dist, diff, 1, grad_F_p, 1);
//          less_ostr << vdw_rad_sum - atomic_dist << " "; 

        // If the overlap is more than 20% move to at most 20% overlap
        }else{
          // atomic_dist is < 4/5 of vdw sum -- move to 4/5 of vdw sum
          my_axpy(3, shortest_dist - atomic_dist, diff, 1, grad_F_p, 1);
///          std::cout << vdw_rad_sum - atomic_dist << " "; 
          *overlap_flag = true;
        }
#endif

      }

      for(size_t i = 0; i < moved_nbrs->positions.size(); ++i){
        // Ignore atoms from same side chain? -- good enough approx for now
 
        // Still too close (wrt vdw overlap to CA)
        atom_vci nbr = moved_nbrs->assoc_data[i];
        if(nbr->name != O && nbr->name != N && nbr->name != C){
          if(res_atom->res_num == nbr->res_num &&
             res_atom->chainID == nbr->chainID &&
             res_atom->altLoc == nbr->altLoc) continue;
        }

        my_float_t diff[3];
        //vector(3, nbr->pos, res_atom->pos, diff);
        vector(3, res_atom->pos, nbr->pos, diff);
        my_float_t atomic_dist_sq = dot(diff, diff);
/*
      std::cout << "(" << nbr->chainID << ") "
                << PDB_residues::residue_to_string(nbr->res) 
                << nbr->res_num
                << " " << PDB_residues::atom_to_string(nbr->name) << " ";
*/

        // If the two atoms are farther apart than the sum of their vdw 
        // radii, they do not contribute to the energy function or gradient
        my_float_t vdw_rad_sum = res_atom->vdw_radius + nbr->vdw_radius;
        if(atomic_dist_sq >= vdw_rad_sum*vdw_rad_sum) continue;

        // Don't move atoms that overlap because of potential hbond 
        // interactions
        if(HbondGeometry::is_hbond_interaction(res_atom->act_type,
                                               nbr->act_type) &&
           atomic_dist_sq >= 2.5*2.5) continue;

        // best direction to move the current atom is along the direction
        // from neighbor to self
        my_float_t atomic_dist = normalize(diff); 
        my_float_t shortest_dist = min_atomic_percent_dist * vdw_rad_sum;

      std::cout << "(" << res->chainID << ") "
                << PDB_residues::residue_to_string(res->name) << res->number 
                << " " << PDB_residues::atom_to_string(res_atom->name) << " -- ";
      std::cout << "   (" << nbr->chainID << ") "
                << PDB_residues::residue_to_string(nbr->res) 
                << nbr->res_num
                << " " << PDB_residues::atom_to_string(nbr->name) << "  "
                << vdw_rad_sum - atomic_dist << " : "
                << round(100.0 - 100.0 * atomic_dist/vdw_rad_sum) << "\%\n";

        // If the overlap is less than 20% move about 10% of the overlap
        if(atomic_dist >= shortest_dist){
          // atomic_dist is >= 4/5 of vdw sum
          my_axpy(3, atomic_dist - shortest_dist, diff, 1, grad_F_p, 1);
//          less_ostr << vdw_rad_sum - atomic_dist << " "; 

        // If the overlap is more than 20% move to at most 20% overlap
        }else{
          // atomic_dist is < 4/5 of vdw sum -- move to 4/5 of vdw sum
          my_axpy(3, shortest_dist - atomic_dist, diff, 1, grad_F_p, 1);
//          std::cout << vdw_rad_sum - atomic_dist << " "; 
          *overlap_flag = true;
        }

      }

      ++overlap_flag;
      grad_F_p += 3;
    }
  }

//  std::cout << "\n";
//  std::cout << less_ostr.str() << "\n";


  // 3) Minimize angle deviations -- need to think about this one more
  // Yiying suggests enforcing constraint in joint space and ignoring it
  // in calculation of the gradient in cartesian space

#if 0
  // print out the gradient
  grad_F_p = grad_F;
  for(int i = 0; i < ncols/3; ++i){
    std::cout << "\n"; 
    for(int j = 0; j < 3; ++j, ++grad_F_p)
       std::cout << *grad_F_p << " ";
  }
  std::cout << "\n"; 
#endif

  // compute pseudo inverse using J^T(JJ^T + \lamdba I)
  query_surf_to_joints.compute_pinvJ_transpose();
  //query_surf_to_joints.write_mat(std::cout, SurfDepsJoints::PINV_JACOBIAN);
  //

  std::cout << "hbond pass: " << rounds << std::endl;

#if 1
  if(A_uses_hbond_surfaces){
    // First stab at polar term -- now compute gradient
    A_site_score_p->polar_correspondences(*model_site, 
                                     query_surf_to_joints.mobile_cap_iters_begin(),
                                     query_surf_to_joints.mobile_cap_iters_end(),
                                     &q_pts, &corr_pts);
//    std::cout << "HB caps pt check: " 
//              << query_surf_to_joints.num_hbond_pts_cols_J() << " "
//	      << q_pts.size() << " " << corr_pts.size() << "\n";
  
    for(int i = 0; i < query_surf_to_joints.num_hbond_pts_cols_J() / 3; ++i){
      if(q_pts[i] && corr_pts[i]){
        vector(3, q_pts[i], corr_pts[i], grad_F_hbond_surfs_verts + 3*i);
        std::cout << "hbond pts: " 
                  << q_pts[i][0] << " "
                  << q_pts[i][1] << " "
                  << q_pts[i][2] << " -- "
                  << corr_pts[i][0] << " "
                  << corr_pts[i][1] << " "
                  << corr_pts[i][2] << " -- "
                  << dist(q_pts[i], corr_pts[i]) << "\n";
      }else std::cout << "hbond pts: 0 0\n";
    } 
  }
#endif

  std::fill(q_vals, q_vals + nrows, 0.0);
  std::fill(surf_q_vals, surf_q_vals + nrows, 0.0);
  std::fill(hbond_q_vals, hbond_q_vals + nrows, 0.0);

//  A_uses_hbond_surfaces = false;
//  A_uses_site_surface = false;
  if(A_uses_hbond_surfaces){
    query_surf_to_joints.apply_hbond_pts_pinvJ_to_grad(grad_F_hbond_surfs_verts,
                                                       hbond_q_vals); 
    // multiply by alpha for hbond cap moves
    my_axpy(nrows, 0.4, hbond_q_vals, 1, q_vals, 1);
  }
  if(A_uses_site_surface){
    query_surf_to_joints.apply_verts_pinvJ_to_grad(grad_F_mol_surf_verts, 
                                                   surf_q_vals);
    // multiply by alpha for surf moves
    my_axpy(nrows, 0.1, surf_q_vals, 1, q_vals, 1);
  }

  // Apply "forces" to reduce overlap
  std::cout << "q vals before checking for overlap: ";
  for(int i = 0; i < nrows; ++i) std::cout << q_vals[i] << " ";
  std::cout << "\n";

  std::fill(overlap_q_vals, overlap_q_vals + nrows, 0.0);

#if 0
  // This method will be looked at last -- at the present it suffers from
  // the issue if there are multiple overlaps that a joint is trying to 
  // resolve it will resolve some at the expense of others
  query_surf_to_joints.apply_atoms_pinvJ_to_grad(grad_F_atoms, overlap_q_vals, 
                                                 overlap_vec);
#endif

  std::vector<bool> fixed_Q;
  query_surf_to_joints.stop_rotation_given_atoms(overlap_vec, &fixed_Q);

  std::cout << "Final Q vals:\n";
  for(int i = 0; i < nrows; ++i){
#if 0
    // We need to move this joint to resolve some overlap
    // move this joint only to resolve the overlap
    if(overlap_q_vals[i] < -1.0E-10 || 1.0E-10 < overlap_q_vals[i])
      q_vals[i] = overlap_q_vals[i];
//      q_vals[i] = 0.5 * overlap_q_vals[i];
    else q_vals[i] *= alpha;
#endif

    if(fixed_Q[i]) q_vals[i] = 0.0;

    // Cap rotation at max of +/- 5 degrees per timestep
    if(q_vals[i] > 0.087266462599716474) q_vals[i] = 0.087266462599716474;
    if(q_vals[i] < -0.087266462599716474) q_vals[i] = -0.087266462599716474;
    std::cout << q_vals[i] << " ";
  }
  std::cout << "\n";
 

  std::cout << "linear q vals: ";
  for(int i = 0; i < nrows; ++i)
    std::cout << q_vals[i] / M_PI * 180.0 << " ";
  std::cout << std::endl;  


#if 0
  // WARNING: this can quickly write out gigabytes of surface mesh files.
  // Write out files so we can visualize the contribution at each time step.
  std::ostringstream qmesh_fname;
  qmesh_fname << ofile_prefix << rounds;
  model_site->binding_site_mesh_handle().write(qmesh_fname.str());

  if(A_uses_hbond_surfaces){
    std::ostringstream caps_fname;
    caps_fname << ofile_prefix << "hbond_caps_" << rounds;
    std::string blah = caps_fname.str();
    model_site->write_msms_caps(blah);
  }
  qmesh_fname << ".pdb";
//  full_protein()->write(qmesh_fname.str());
  model_site->bind_site_atoms().write(qmesh_fname.str());
#endif


  query_surf_to_joints.update_positions(q_vals);

  // Just do 1 LSE (ICP) round
  if(A_hbond_caps_ICP == false)
    ICP::fine_tune_surfaces(model_site, &dbase_site, align, A_site_score_p,
                            1E-06, 1, true);
  else
    ICP::fine_tune_caps_and_surf(model_site, &dbase_site, align, A_site_score_p,
                                 A_ICP_surf_pt_W, A_ICP_hb_cap_pt_W, 
                                 1E-06, 1, true);

} // end of for( ; rounds ; )

  std::cout << "Final score call\n\n";
  align->score = A_site_score_p->score(*model_site, dbase_site, align);

  // NOTE: rotation and translation is currently computed based on sitemap
  // points -- as points may have moved relative to each other, they are no
  // longer valid.  One could use the main chain protein atoms at this 
  // point since they are currenlty fixed relative to the global orientation
  std::fill(align->R, align->R + 9, 9999999999);
  std::fill(align->T, align->T + 3, 9999999999);

  std::cout << "Note: IK score is not normalized\n";

  atomic_rmsd = 
    query_surf_to_joints.site_atomic_rmsd(dbase_site.bind_site_atoms());
  align->set_site_atomic_rmsd(atomic_rmsd, 2);


  if(grad_F_mol_surf_verts) delete[](grad_F_mol_surf_verts);
  if(grad_F_atoms) delete[](grad_F_atoms);
  if(grad_F_hbond_surfs_verts) delete[](grad_F_hbond_surfs_verts);
  if(q_vals) delete [](q_vals);
  if(overlap_q_vals) delete[](overlap_q_vals);
  if(surf_q_vals) delete[] surf_q_vals;
  if(hbond_q_vals) delete[] hbond_q_vals;
}
