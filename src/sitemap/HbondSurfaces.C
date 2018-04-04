#include <HbondSurfaces.H>
#include <sstream>

using namespace ASCbase;

const my_float_t hbond_surface_t::SURF_SPHERE_RAD = 3.0;

hbond_surface_t::hbond_surface_t(const std::string& data_line, 
                                 PDBBase& prot_atoms)
{
  // This is just a guess -- and keeps the geometry simple to compute
  // analytically
  my_float_t h = 1.5;

  std::vector<std::string> toks;
  string_tok(data_line, &toks, '|');

  // Load the point type
  if(toks[0] == "ACCEPTOR") act_type = ACCEPTOR;
  else if(toks[0] == "DONOR") act_type = DONOR;
  else if(toks[0] == "DONEPTOR") act_type = DONEPTOR;
  else if(toks[0] == "METAL") act_type = METAL_1;
  else{
    std::cerr << "Unknown point type in hbond_surface_t" << std::endl;
    return;
  }

  // Get the corresponding atom constant iterators
  uint atom_num = 0;
  my_strtoui(toks[1], &atom_num);
  atom = prot_atoms.get_atom(atom_num);
  if(atom == atom_t::NULL_ATOM_VCI){
    std::stringstream sstr;
    sstr << "Received atom number (" << atom_num 
         << ") for atom not in rad file";
    err_msg("", "hbond_surface_t::cstr()", sstr.str());
    return;
  }

  // Needed to get a unique orientation for the cap
  my_strtoui(toks[2], &atom_num);
  carbon_nbr = prot_atoms.get_atom(atom_num);
  my_strtoui(toks[3], &atom_num);
  second_nbr = prot_atoms.get_atom(atom_num);
  if(carbon_nbr == atom_t::NULL_ATOM_VCI || 
     second_nbr == atom_t::NULL_ATOM_VCI)
  {
    std::stringstream sstr;
    if(carbon_nbr == atom_t::NULL_ATOM_VCI)
      sstr << "Received atom number (" << toks[2] 
           << ") for carbon_nbr that is not in rad file";
    else
      sstr << "Received atom number (" << atom_num 
           << ") for atom not in rad file";
    err_msg("", "hbond_surface_t::cstr()", sstr.str());
    return;
  }
  
  // Point number for the hbond atom
  my_strtoui(toks[4], &pt_num);

  // Sphere is defined by center of hbond atom & some radius
  if(toks[0] != "METAL"){
    A_surf = geometry::sphere_t(atom->pos, SURF_SPHERE_RAD);

    // Load the plane for this cap
    my_float_t pvals[6];
    std::vector<std::string> stoks;
    string_tok(toks[5], &stoks, ' '); 
    for(int i = 0; i < 6; ++i) my_strtof(stoks[i], pvals + i);
    A_cut_plane = geometry::plane_t(pvals, pvals + 3);            
    A_cap_plane_rad = std::sqrt(A_surf.radius() * A_surf.radius() - h*h);
  }else{
    const pdb_metal_info_t* metal_info = PDB_metals::lookup(atom->name);
    A_surf = geometry::sphere_t(atom->pos, metal_info->interact_pt_rad);
  }

  
  // For the remaining toks
  // Parse each tok as comma delimited (denoted as ctok) where the first 
  // ctok is the sphere, second is the plane, following ctoks are the 
  // arcs.
  for(uint i = 6; i < toks.size(); ++i){
    if(toks[i].length() == 0) continue;

//    std::cout << "toks[" << i << "]: " << toks[i] << "\n";
    std::vector<std::string> ctoks;
    string_tok(toks[i], &ctoks, ',');
//    for(uint j = 0; j < ctoks.size(); ++j)
//      std::cout << i << ":" << j << " " << ctoks[j] << "\n";
    
    // Get sphere
    std::vector<std::string> stoks;
    string_tok(ctoks[0], &stoks, ' '); 
//    std::cout << "stoks " << stoks[0] << " " << stoks[1] << " "
//              << stoks[2] << " " << stoks[3] << "\n";
    my_float_t svals[4];
    for(int zz = 0; zz < 4; ++zz) my_strtof(stoks[zz], svals + zz);
    geometry::sphere_t circ_surf(svals, svals[3]);            

    // Get plane
    my_float_t pvals[6];
    stoks.clear();
    string_tok(ctoks[1], &stoks, ' '); 
    for(int zz = 0; zz < 6; ++zz) my_strtof(stoks[zz], pvals + zz);
    geometry::plane_t circ_plane(pvals, pvals + 3);            
    
    // Get arcs (if any)
    std::vector<geometry::arc_t> circ_arcs;
    for(uint j = 2; j < ctoks.size(); ++j){
//      std::cout << "J:" << j;
      stoks.clear();
      string_tok(ctoks[j], &stoks, ' '); 
      my_float_t arc_vals[9];
//      std::cout << "   num:" << stoks.size() << std::endl;
      for(int zz = 0; zz < 9; ++zz) my_strtof(stoks[zz], arc_vals + zz);
      circ_arcs.push_back(geometry::arc_t(circ_surf, arc_vals, arc_vals + 3, 
                                          arc_vals + 6));

    
    }
//    std::cout << "svals: " << svals[0] << " " << svals[1] << " "
//              << svals[2] << " " << svals[3] << "\n";
    geometry::iCircle my_circ(svals, svals[3], circ_plane.normal(), 
                              circ_arcs);
    A_circles.push_back(my_circ);
  }
}

bool
hbond_surface_t::closest_point(const my_float_t* pt, my_float_t* closest_pt, 
                               const my_float_t tol)
{
  my_float_t sq_tol = tol*tol;
  my_float_t proj_pt[3];
  std::fill(closest_pt, closest_pt + 3, my_float_max);
  
//  std::cout << "projecting point onto sphere\n";

  // Project point onto sphere
  A_surf.project_point(pt, proj_pt);
  if(dist_squared(proj_pt, pt) > sq_tol) return false;

  if(act_type == ACCEPTOR || act_type == DONOR || act_type == DONEPTOR){
//  std::cout << "Determine if sphere is above or below cap plane\n";
    // Determine if the sphere of intersection is above, below or is cut
    // by the plane used to define the spherical cap
    my_float_t dist_to_plane = A_cut_plane.signed_dist(proj_pt);
    if(dist_to_plane < -1.0 * tol) return false;

//  std::cout << "signed dist to plane: " << dist_to_plane << "\n";

    // If the projected point is below the plane, project it onto the
    // intersection of the cutting plane & the sphere (i.e. the circle of the
    // cap in the plane).
    if(dist_to_plane < 0.0){
//    std::cout << "Point is below plane ... Project to plane\n";
    // 1st, project pt to plane
//    std::cout << "plane normal: " << A_cut_plane.normal()[0] << " "
//              << A_cut_plane.normal()[1] << " " << A_cut_plane.normal()[2] << "\n";
      my_axpy(3, -1.0 * dist_to_plane, A_cut_plane.normal(), 1, proj_pt, 1);
//    std::cout << "proj_pt: " << proj_pt[0] << " "
//              << proj_pt[1] << " " << proj_pt[2] << "\n";
      if(dist_squared(proj_pt, pt) > sq_tol) return false;

      // 2nd, project the pt in plane to the circle (sphere restricted to the
      // plane
      my_float_t in_plane_dir[3];
      unit_vector(in_plane_dir, proj_pt, A_cut_plane.point());
  //    std::cout << "in plane dir: " << in_plane_dir[0] << " "
  //              << in_plane_dir[1] << " " << in_plane_dir[2] << "\n";
      std::copy(A_cut_plane.point(), A_cut_plane.point() + 3, proj_pt);
      my_axpy(3, A_cap_plane_rad, in_plane_dir, 1, proj_pt, 1);
      if(dist_squared(proj_pt, pt) > sq_tol) return false;
    }
  }

  // We now step through each of the kept circles and determine if the
  // projected point falls inside any of them

//  std::cout << "Now project the point to arcs\n";
//  std::cout << "Proj pt: " << proj_pt[0] << " " << proj_pt[1] << " "
//            << proj_pt[2] << "\n";

  std::vector<geometry::iCircle::point_type> arc_pts;
  std::vector<my_float_t> arc_pts_sq_dists;
  std::vector<geometry::iCircle>::iterator C;
  bool inside_a_circle = false;
  for(C = A_circles.begin(); C < A_circles.end(); ++C){
    if(!C->contains(proj_pt)) continue;
    inside_a_circle = true;
//    std::cout << "C " << *C << "\n";
    
    geometry::iCircle::point_type pt_on_C;
    my_float_t d2;
    if(!C->project_point(proj_pt, pt_on_C.pt, &d2, tol)) continue;
    if(d2 <= sq_tol){
      arc_pts.push_back(pt_on_C);
      arc_pts_sq_dists.push_back(d2);
    }

#if 0
    my_float_t d2 = dist_squared(pt_on_C.pt, pt);
    std::cout << "sq_dist: " << d2 << "\n";
    if(d2 > sq_tol) continue;
    else if(C->full_circle()){
      std::cout << "Full circle -- appending to vector\n";
      arc_pts.push_back(pt_on_C);
      arc_pts_sq_dists.push_back(d2);
    // Test each arc to determine the closest point
    }else if(C->num_final_arcs()){
      std::vector<geometry::arc_t>::const_iterator A;
      for(A = C->final_arcs_begin(); A < C->final_arcs_end(); ++A){
        std::cout << "  arc end points " << A->end_pts()[0] << " "
                  << A->end_pts()[1] << " " << A->end_pts()[2] << "   "
                  << A->end_pts()[3] << " " << A->end_pts()[4] << " "
                  << A->end_pts()[5] << "\n";
        std::cout << "  arc in dir: " << A->in_dir()[0] << " "
                  << A->in_dir()[1] << " " << A->in_dir()[2] << "\n";
        // If the projected point is in the arc, keep it
        if(A->contains(pt_on_C.pt)){
          std::cout << "  Arc contains pt\n";
          arc_pts.push_back(pt_on_C);
          arc_pts_sq_dists.push_back(d2);
        // Otherwise choose the closest end point of the arc
        }else{
          std::cout << "  Arc does not contain pt\n";
          my_float_t d2z[] = { dist_squared(A->end_pts(), pt), 
                               dist_squared(A->end_pts() + 3, pt) };
          my_float_t min_d2;
          int idx = argmin<my_float_t>(d2z, 2, &min_d2);
          if(min_d2 <= sq_tol){
            std::copy(A->end_pts() + 3*idx, A->end_pts() + 3*idx + 3, 
                      pt_on_C.pt);
            arc_pts.push_back(pt_on_C);
            arc_pts_sq_dists.push_back(min_d2);
          }
        }
      }
    }
#endif
  }
//  printf("arc_pts size: %d\n", arc_pts.size());

  // Find the point that is the closest to the original point and return that
  // one
  if(inside_a_circle){
    if(!arc_pts.size()) return false;

    my_float_t min_dist = arc_pts_sq_dists[0];
    uint idx = 0;
    for(uint i = 1; i < arc_pts_sq_dists.size(); ++i){
      if(arc_pts_sq_dists[i] < min_dist){
        min_dist = arc_pts_sq_dists[i];
        idx = i;
      }
    }
    std::copy(arc_pts[idx].pt, arc_pts[idx].pt + 3, closest_pt);

  }else std::copy(proj_pt, proj_pt + 3, closest_pt);

//  printf("closest pt: %f %f %f\n", closest_pt[0], closest_pt[1], closest_pt[2]);

  return true;
}
