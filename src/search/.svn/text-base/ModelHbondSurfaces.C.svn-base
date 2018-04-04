#include <ModelHbondSurfaces.H>
#include <HbondGeometry.H>

using namespace SimSite3D;

const uint model_hbond_surf_t::A_num_terms = 6;
  
model_hbond_surf_t::model_hbond_surf_t(const std::string& data_line, 
                                       PDBBase& prot_atoms, const uint level)
  : hbond_surface_t(data_line, prot_atoms)
{
  init();

  // 1) Get a tesselated cap -- precomputed so that we don't duplicate points
  // Get transformation to move from local coords to global
  my_float_t trans_A[3], trans_C[3];
  for(uint i = 0; i < 3; i++){
    trans_A[i] = carbon_nbr->pos[i] - atom->pos[i];
    trans_C[i] = second_nbr->pos[i] - atom->pos[i];
  }

  my_float_t R[9];
  get_local_orientation(trans_A, trans_C, R);

  // Quick & dirty method to get the model surface for the metal to be
  // transformation invariant (up to numerical precision of course).
  if(act_type == METAL_1){
    my_float_t my_dir[3];
    unit_vector(my_dir, carbon_nbr->pos, atom->pos);
    const pdb_metal_info_t* metal_info = PDB_metals::lookup(atom->name);
    A_mesh_surf = 
      geometry::TriMeshSphere(R, atom->pos, my_dir, 
                              geometry::TriMeshSphere::FULL_SPHERE, 
                              metal_info->interact_pt_rad, level);
  }else
    A_mesh_surf = 
      geometry::TriMeshSphere(R, atom->pos, ideal_dir(),
                              geometry::TriMeshSphere::SPHERICAL_CAP, 3.0, 
                              level);
  
  // 2) Remove all triangles that have 3 points outside of the surface --- we
  //    can use (hbond_surface_t::)closest_point since this is initialization
  A_mesh_surf.adjust_points(*this);

  // Allocate local memory
  A_closest_pts = new my_float_t[3*A_mesh_surf.num_verts()];
  A_dists = new my_float_t[A_mesh_surf.num_verts()];
}

//! Compute the best terms for each vertex in the model surface
/*!
 * For now we will throw the AA, DD, MA, MM into one basket and the
 * N* into another.
 *
 * This function updates the values of A_closest_pts and A_dists
 */
void
model_hbond_surf_t::compute_best_terms(std::vector<hbond_surface_t>::iterator 
                                       surfs_begin, 
                                       std::vector<hbond_surface_t>::iterator 
                                       surfs_end, my_float_t *AA_DD_terms, 
                                       my_float_t *doneptor_terms, 
                                       my_float_t *surface_area, 
                                       const my_float_t dist_tol) const
{
  std::fill(A_closest_pts, A_closest_pts + 3*A_mesh_surf.num_verts(),
            my_float_max);
  std::fill(A_dists, A_dists + A_mesh_surf.num_verts(), my_float_max);
  my_float_t AA_DD_cp[3], doneptor_cp[3]; 
  std::fill(AA_DD_cp, AA_DD_cp + 3, my_float_max);
  std::fill(doneptor_cp, doneptor_cp + 3, my_float_max);
  my_float_t AA_DD_cp_dist = my_float_max;
  my_float_t doneptor_cp_dist = my_float_max;

  *surface_area = 0.0;
  uint N = A_num_terms * A_mesh_surf.num_verts();
  std::fill(AA_DD_terms, AA_DD_terms + N, 0.0);
  std::fill(doneptor_terms, doneptor_terms + N, 0.0);
  my_float_t* AA_DD_row = AA_DD_terms;
  my_float_t* doneptor_row = doneptor_terms;

  std::map<geometry::vertex_vci, bool> hit_vertices;
  my_float_t *saved_cp = A_closest_pts;
  my_float_t *saved_dist = A_dists;
  geometry::vertex_vci q_vert = A_mesh_surf.verts_begin();
  for( ; q_vert < A_mesh_surf.verts_end(); ++q_vert, ++saved_dist, saved_cp+=3){
    std::fill(AA_DD_cp, AA_DD_cp + 3, my_float_max);
    std::fill(doneptor_cp, doneptor_cp + 3, my_float_max);
    AA_DD_cp_dist = my_float_max;
    doneptor_cp_dist = my_float_max;

/*
    std::cout << "query pt num: " << q_vert - A_mesh_surf.verts_begin() 
              << "\n";
*/

    hbond_surf_vi d_surf_iter;
    for(d_surf_iter = surfs_begin; d_surf_iter != surfs_end; ++d_surf_iter){
      hbond_surface_t& d_surf = *d_surf_iter;

      // Interaction filter
      if(act_type == ACCEPTOR && d_surf.act_type == DONOR) continue;
      else if(act_type == DONOR && (d_surf.act_type == ACCEPTOR ||
                                    d_surf.act_type == METAL_1)) continue;
      else if(act_type == METAL_1 && d_surf.act_type == DONOR) continue; 

/*
      std::cout << "DB HB cap for residue (" << d_surf.atom->chainID << ") "
	        << PDB_residues::residue_to_string(d_surf.atom->res)
		<< d_surf.atom->res_num << " "
		<< PDB_residues::atom_to_string(d_surf.atom->name) << std::endl;
*/

      //std::cout << "db surf num: " << d_surf_iter - surfs_begin << std::endl;

      // Distance tolerance filter
      my_float_t cp[3];
      if(!d_surf.closest_point(q_vert->pos, cp, dist_tol)){
//        std::cout << "pt too far away\n";
        continue;
      }

      // Direction filter
      my_float_t cp_dir[3];
      unit_vector(cp_dir, cp, d_surf.atom->pos);
      my_float_t dot_prod = dot(q_vert->dir, cp_dir);
      if(dot_prod <= 0.0){
//        std::cout << "dot product is negative\n";
        continue;
      }

      my_float_t tmp_dist = dist(q_vert->pos, cp);
//      std::cout << "distance is: " << tmp_dist << std::endl;
      if(act_type == DONEPTOR || d_surf.act_type == DONEPTOR){
        if(adjust_terms(dot_prod, tmp_dist, dist_tol, doneptor_row)){
          std::copy(cp, cp + 3, doneptor_cp);
          doneptor_cp_dist = tmp_dist;
        }
      }else{
        if(adjust_terms(dot_prod, tmp_dist, dist_tol, AA_DD_row)){
          std::copy(cp, cp + 3, AA_DD_cp);
          AA_DD_cp_dist = tmp_dist;
        }
      }
    }

    // For each point, we want only the best score and zero out the other 
    // point type score
    for(uint i = 1; i < A_num_terms; ++i){
      if(doneptor_row[i] > AA_DD_row[i]) AA_DD_row[i] = 0.0;
      else doneptor_row[i] = 0.0;
    }

    // Save correpsonding point and distance
    
    if(AA_DD_row[4] > doneptor_row[4] && AA_DD_row[4] != 0.0){
      std::copy(AA_DD_cp, AA_DD_cp + 3, saved_cp);
      *saved_dist = AA_DD_cp_dist;
    }else if(doneptor_row[4] != 0.0){
      // doneptors do not seem to be useful at this point to drive matches
      //std::copy(doneptor_cp, doneptor_cp + 3, saved_cp);
      //*saved_dist = doneptor_cp_dist;
    }

    if(AA_DD_row[0] || doneptor_row[0]){
      hit_vertices[q_vert] = true;
//      std::cout << "vertex was hit\n";
    }
/*
      std::cout << "saved dist: " << *saved_dist << " ";
      std::cout << "q vert: " << q_vert->pos[0] << " "
                << q_vert->pos[1] << " " << q_vert->pos[2] << " -- ";
      std::cout << "saved_cp: " << saved_cp[0] << " " << saved_cp[1] << " "
                << saved_cp[2] << "\n";
*/

    AA_DD_row += A_num_terms;
    doneptor_row += A_num_terms;
  }

#if 0
  std::cout << "Polar correspondances\n";
  q_vert = A_mesh_surf.verts_begin();
  saved_dist = A_dists;
  saved_cp = A_closest_pts;
  for( ; q_vert < A_mesh_surf.verts_end(); ++q_vert, ++saved_dist, saved_cp+=3){
    std::cout << *saved_dist << "|"
              << q_vert->pos[0] << " " << q_vert->pos[1] << " "
              << q_vert->pos[2] << "|"
              << saved_cp[0] << " " << saved_cp[1] << " "
              << saved_cp[2] << "|\n";
  }
#endif

  *surface_area = A_mesh_surf.complementary_surface_area(hit_vertices);
}

ModelHbondSurfaces::ModelHbondSurfaces(const std::string fname, 
                                       PDBBase& rad_atoms,
                                       const verbose_level_t verbosity)
  : HbondSurfaces<model_hbond_surf_t>(fname, rad_atoms, verbosity)
{
  init();
  if(verbosity != VERBOSE_SILENT){
    surfaces_vci zz;  
    for(zz = surf_caps_begin(); zz < surf_caps_end(); ++zz)
      zz->write(std::cout, zz->act_type);
  }

  // Get a count of the number of points
  for(surfaces_vci zz = surf_caps_begin(); zz < surf_caps_end(); ++zz){
    if(zz->act_type == ACCEPTOR || zz->act_type == DONOR ||
       zz->act_type == METAL_1 || zz->act_type == METAL_2) 
      A_num_A_D_points += zz->num_verts();
    A_num_polar_points += zz->num_verts();
  }
}
