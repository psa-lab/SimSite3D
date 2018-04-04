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

#include <errno.h>
#include <unistd.h>
#include <iomanip>
#include <sstream>
#include <param_tools.H>
#include <Timer.H>
#include <mol2File.H>
#include <Sitemap.H>
#include <interact_point_t.H>
#include <normalize_sitemap.H>

using namespace ASCbase;
const std::string Sitemap::A_fname = "Sitemap.C";
const my_float_t Sitemap::TOLERANCE = 2.0;
const int Sitemap::MIN_NUMBER_POINTS__POLAR_SITE = 20;
const int Sitemap::MIN_NUMBER_POINTS__HPHOB_SITE = 30;


Sitemap::Sitemap(const std::string path, const std::string struct_id,
                 const BaseParameters& args, const bool normalize,
                 const bool load_hbond_surfaces, const bool hydrophobic_query, 
                 const verbose_level_t verbosity)
{
  init();
  verbose = verbosity;
  A_site_path = path;
  A_dbase_ligs = args.dbase_ligs;
  A_struct_id = struct_id;

  // Verify that path to sitemap exists
  std::cout << "Loading site map for: " << A_struct_id << "\n";
  std::string my_struct = path + "/" + A_struct_id;
  std::string my_site = my_struct + "_s.csv";
  if(!normal_file_exists(my_site, false)){
    // If the sitemap as given does not exist, check in the dbase_sites dir
    std::cout << "Unable to find the file: " << my_site << "\n";
    A_site_path = args.dbase_sites;
    my_struct = args.dbase_sites + "/" + A_struct_id;
    my_site = my_struct + "_s.csv";
    std::cout << "Checking for " << A_struct_id 
              << "_s.csv in $ASCBASE_DBASE_SITES\n" 
              << "\t(" << args.dbase_sites << ")\n";
    if(!normal_file_exists(my_site, false)){
      std::cerr << "Unable to find the file: " << my_site << "\n"
                << "Skipping the binding site of " << A_struct_id << "\n";
      A_fail = true;
      return;
    }   
  }

  A_atoms_fname = my_struct + "_a.pdb";
  A_site_atoms = new SitemapAtomsFile(A_atoms_fname, verbosity);
  if(A_site_atoms == 0 || A_site_atoms->fail()){
    A_fail = true;
    return;
  } 

  // Changed so that DrugScore has enough atoms to work with
  A_atoms_fname = my_struct + "_rad.pdb";
  A_bind_site_res = new PDBStructure(A_atoms_fname, false, verbosity);
  if(A_bind_site_res->fail()){
    A_fail = true;
    return;
  } 

  points = new SitemapPointsFile(my_struct + "_s.pdb", verbosity);
  if(points == 0 || points->fail()){
    A_fail = true;
    return;
  } 

  // Move this to query sitemap since that is where it belongs
  if(normalize && 
     !get_norm_stats(my_struct + "_s.csv", args, hydrophobic_query))
  {
    A_fail = true;
    return;
  }else if(!normalize){
    // Init the mean and stdev to default values
    A_norm_moments[0] = 0.0;
    A_norm_moments[1] = -1.0;
  }

  std::ifstream labs_file;
  if(!open_ifstream(labs_file, my_struct + "_s.csv")){
    A_fail = true;
    return;
  }

#if 0
  // These items might not be that useful after all
  if(!load_volume(path, A_struct_id, args.dbase_sites, "<site_vol_est>", 
                  &A_site_vol_est, &A_site_vol_est_as_sphere)){
    std::string msg = "Unable to get the site map volume estimate\n";
    err_msg(A_fname, "cstr", msg);
    set_fail_flag(true);
    return;
  }

  // Find the mol_surf parameters section
  if(!find_section(labs_file, "<mol_surf>")){
    warn(A_fname, "Cstr()", std::string("Could not find a mol_surf section ")
         +  "in: " + A_struct_id);
  }else{
    // need to implement this
    std::cerr << "Need to implement the reading of msms_surf paramters "
              << "section\n\n";
  }
#endif
  
  // Find the ligand
  if(!find_section(labs_file, "<ligand_file>")){
    warn(A_fname, "Cstr()", std::string("Could not find a ligand file section ")
         +  "in: " + A_struct_id);
    A_fail = true;
    return;
  }
  labs_file >> A_lig_fname;

  // Look for point counts section
  if(!find_section(labs_file, "<point_counts>")){
    warn(A_fname, "Cstr()", std::string("Could not find the point count ") +
         "section in: " + A_struct_id);
    A_fail = true;
    return;
  }
  size_t num_hbond_pts, num_hphob_pts;
  labs_file >> num_hbond_pts >> num_hphob_pts;

  // Find the hbond section
  if(!find_section(labs_file, "<hydrogen_bond_points>")){
    warn(A_fname, "Cstr()", std::string("Could not find the hydrogen bond ") +
         "points section in: " + A_struct_id);
    A_fail = true;
    return;
  }
  hbond_pts = 
    new HbondPoints(labs_file, num_hbond_pts, *A_bind_site_res, *points);
  if(hbond_pts == 0 || hbond_pts->fail()){
    A_fail = true;
    return;
  } 
  std::cout << "\tloaded " 
            << hbond_pts->fit_pts_end() - hbond_pts->fit_pts_beg()
	    << " hbond points\n";

  // Find the hphob section
  if(!find_section(labs_file, "<hydrophobic_points>")){
    warn(A_fname, "Cstr()", std::string("Could not find the hydrophobic ") +
         "points section in: " + A_struct_id);
    A_fail = true;
    return;
  }
  hphob_pts = new HphobPoints(labs_file, num_hphob_pts, *A_site_atoms, *points);
  if(hphob_pts == 0){
    A_fail = true;
    return;
  } 
  std::cout << "\tloaded " << hphob_pts->size() << " hydrophobic points\n";

  if(args.require_min_npts && !has_enough_points()){
    A_fail = true;
    return;
  }
  
  // Generate the so called split points from the hbond fit points
  interact_point_t pt;
  hbond_fit_pt_vci hb = hbond_pts->fit_pts_beg();
  for( ; hb != hbond_pts->fit_pts_end(); ++hb){
    std::copy(hb->pos, hb->pos + 3, pt.pos);
    std::copy(hb->dir, hb->dir + 3, pt.dir);
    if(hb->act_type == DONEPTOR){
      pt.act_type = ACCEPTOR;
      interact_points.push_back(pt);
      pt.act_type = DONOR;
    }else pt.act_type = hb->act_type;
    interact_points.push_back(pt);
  }

  // hydrophobic interaction points
  pt.act_type = HYDROPHOB;
  std::fill(pt.dir, pt.dir + 3, 0);
  hphob_point_lci hp;
  for(hp = hphob_pts->begin(); hp != hphob_pts->end(); ++hp){
    std::copy(hp->pos, hp->pos + 3, pt.pos);
    interact_points.push_back(pt);    
  }

  if(load_hbond_surfaces){
    // Verify that path to sitemap exists
    std::cout << "Loading hbond surfaces for: " << A_struct_id << "\n";
    my_site = my_struct + "_surf_caps.csv";
    if(!normal_file_exists(my_site, false)){
      std::cerr << "Unable to find the file: " << my_site << "\n"
                << "Skipping the binding site of " << A_struct_id << "\n";
      A_fail = true;
      return;
    }   
    A_hbond_surfaces = 
      new HbondSurfaces<hbond_surface_t>(my_site, *A_bind_site_res);
  }


  A_fail = false;
}

Sitemap::Sitemap(const GenPointsParameters& args)
{
  init();
  A_include_metals = args.include_metals;
  A_min_res_per_chain = args.min_res_per_chain;
  verbose = args.verbose_level;
  // 0) Init the mean and stdev to default values
  A_norm_moments[0] = 0.0;
  A_norm_moments[1] = -1.0;

  // 1) Load protein structure -- ignoring altLocs
  prot_atoms = new PDBStructure(args.prot_fname, true, verbose);
  std::cout << "number of protein atoms: " << prot_atoms->num_atoms() << std::endl;
  if(prot_atoms->fail()){
    A_fail = true;
    std::cerr << "Failed to read the protein file\n";
    return;
  }

  // 2) Define site volume, discretize volume and remove the grid points which
  //    are too far from the protein or are too close to the protein
  A_site_vol = setup_volume(args.lig_fname, args.grid_spacing, args.sphere_str,
                            args.prune_to_lig);
  if(!A_site_vol){
    A_fail = true;
    std::cerr << "Failed to get the bounding volume\n";
    return;
  }

  // 2.5) Determine the "binding site" interacting atoms
  std::vector<atom_vci> bind_atoms;
  std::vector<atom_vci> rad_atoms; 
  interact_atoms_vec hbond_atoms;
  interact_atoms_vec hphob_atoms;
  hphob_triad_vec hphob_triads;
  get_prot_interact_atoms(prot_atoms, A_site_vol, &bind_atoms, &rad_atoms,
                          &hbond_atoms, &hphob_atoms, &hphob_triads,
                          args.min_res_per_chain);
  if(rad_atoms.empty()){
    std::cerr << "Failed to get the binding site hbonds\n";
    A_fail = true;
    return;
  }

  // 3) Compute the interaction points and keep track of their associated atoms
  points = new SitemapPointsFile;
  build_points_header(args);
  find_binding_site_hbonds(args.hbond_method, 
                           args.install_dir + "/ASCbaseSoftParams", args.waters,
                           A_site_vol, hbond_atoms, &rad_atoms, &bind_atoms);
  find_binding_site_hphobs(args.hphob_method, args.cluster_diameter,
                           hphob_atoms, hphob_triads, bind_atoms, A_site_vol, 
                           args.sphere_sample_level); 

  // We are only interested in the binding site atoms and can ignore those
  // atoms which are farther away
  if(args.lig_fname.length())
    A_site_atoms = 
      new SitemapAtomsFile(prot_atoms, interact_atoms, A_lig_fname,
                           args.min_res_per_chain);
  else{
    my_float_t center[4];
    std::istringstream sphere_args(args.sphere_str);
    if(sphere_args.fail()){
      warn(A_fname, "get_user_surf_patches", 
           "Error parsing center and radius for sphere", std::cerr);
      A_fail = true;
      return;
    }
    for(uint i = 0; i < 4; i++) sphere_args >> center[i];
    Ball tmp_vol(center, center[3] + 12.0);
    A_site_atoms = new SitemapAtomsFile(prot_atoms, &tmp_vol, 
                                        args.min_res_per_chain, interact_atoms);
  }

  // 4) Compute the MSMS surface for the binding site
  if(args.call_msms){
    A_probe_radius = args.probe_radius;
    A_num_pts_per_area = args.num_pts_per_area;
    if(get_MSMS_surf_patches(args.lig_fname, args.scratch_dir, args.sphere_str,
                             args.install_dir, args.include_metals, 
                             args.msms_binary, args.waters))
      A_MSMS_surf = true;
    else return;
  }else if(args.user_provided_surf.length()){
    if(get_user_surf_patches(args.lig_fname, args.scratch_dir,
                          args.sphere_str, args.user_provided_surf))
      A_user_surf = true;
    else return;
  }

  A_fail = false;
}

bool
Sitemap::get_user_surf_patches(const std::string &lig_fname,
                               const std::string &scratch_dir,
                               const std::string &sphere_str,
                               const std::string &user_vert_fname)
{
  std::cout << "Pruning the given molecular surface to be near the binding site\n";

  // Create the volume representation of the site
  // Test different volumetric distances here
  BoundingVolume* mesh_bv = 0;
  if(lig_fname.size() < 1){
    my_float_t center[4];
    std::istringstream sphere_args(sphere_str);
    for(uint i = 0; i < 4; i++) sphere_args >> center[i];
    if(sphere_args.fail()){
      warn(A_fname, "get_user_surf_patches", 
           "Error parsing center and radius for sphere", std::cerr);
      return false;
    }
    mesh_bv = new Ball(center, center[3] + 1.0);

  }else mesh_bv = new UnionOfBalls(mol2File(lig_fname), 4.0);

  // Load the given mesh -- remove the ".vert" from the filename (if it is there)
  SimpleTrimesh full_mesh(user_vert_fname, true);
  if(full_mesh.fail()) return false;

  // Prune it back to the ligand/sphere
  SimpleTrimesh *curr_comp_mesh = new SimpleTrimesh(full_mesh, *mesh_bv);
  if(curr_comp_mesh == 0 || curr_comp_mesh->fail()){
    if(curr_comp_mesh) delete curr_comp_mesh;
    return false;
  }
  old_binding_site_mesh = curr_comp_mesh;

  return true;
}

bool 
Sitemap::get_MSMS_surf_patches(const std::string &lig_fname,
                               const std::string &scratch_dir,
                               const std::string &sphere_str,
                               const std::string &install_dir,
                               const bool include_metals,
                               const std::string &msms_binary,
                               const std::vector<std::string> &waters)
{
  std::cout << "Generating the MSMS surface for the binding site" << std::endl;

  ////////////////////////////////////////////////////////////////////////
  // It may be that we have a mostly closed pocket, we don't know.  To get
  // around this limitation, we have msms generate all of the surface
  // components, cull each component to the faces near the binding site,
  // keep the component with the highest remaining surface area.
  ////////////////////////////////////////////////////////////////////////

  // Test different volumetric distances here
  BoundingVolume* mesh_bv = 0;
  if(lig_fname.size() < 1){
    my_float_t center[4];
    std::istringstream sphere_args(sphere_str);
    for(uint i = 0; i < 4; i++) sphere_args >> center[i];
    if(sphere_args.fail()){
      warn(A_fname, "sitemap cstr", 
           "Error parsing center and radius for sphere", std::cerr);
      return false;
    }

    // We may need to increase the tolerance
    mesh_bv = new Ball(center, center[3] + 1.0);
  }else mesh_bv = new UnionOfBalls(mol2File(lig_fname), 4.0);

  // one thing at a time -- after testing using the "rad" file,
  // we can use a atom the starting probe must touch if we are still
  // getting too many "RS components"
  
  // Get a temporary file name to help convert PDB prot to XYZR
  std::string xyzr_str = std::string(scratch_dir) + "/my_xyzr_XXXXXX";
  char *xyzr_fname = new char[xyzr_str.length() + 1];
  strncpy(xyzr_fname, xyzr_str.c_str(), xyzr_str.length() + 1);
  if(mkstemp(xyzr_fname) == -1){
    int errsv = errno;
    std::ostringstream msg;
    msg << "Unable to create the temporary xyzr file for MSMS: " 
        << strerror(errsv);
    err_msg(A_fname, "get_MSMS_surf_patches()", msg.str());
    return false;
  }

  // Convert PDB prot to XYZR
  std::ofstream ofile;
  uint dummy_var;
  if(!open_ofstream(ofile, xyzr_fname)) return false;
  A_site_atoms->write_xyzr(ofile, atom_t::NULL_ATOM_VCI, &dummy_var);
  ofile.close();

  // Get a temporary file name to help convert PDB prot to XYZR
  std::string msms_log_str = std::string(scratch_dir) + "/msms_stdout_XXXXXX";
  char *msms_log_fname = new char[msms_log_str.length() + 1];
  strncpy(msms_log_fname, msms_log_str.c_str(), msms_log_str.length() + 1);
  if(mkstemp(msms_log_fname) == -1){
    int errsv = errno;
    std::ostringstream msg;
    msg << "Unable to create the temporary log file for MSMS: " 
        << strerror(errsv);
    err_msg(A_fname, "get_MSMS_surf_patches()", msg.str());
    return false;
  }

  // Create "full" MSMS surface -- default probe radius of 1.5 is a bit too
  // big, use 1.4 instead (1.4 is closer to the 1.36 we use for water anyhow)
  std::stringstream cmd_line;
  if(msms_binary == "default") cmd_line << install_dir <<  "/bin/linux_msms";
  else cmd_line << msms_binary;
  cmd_line << " -if " << xyzr_fname << " -of " << xyzr_fname
           << " -probe_radius " << A_probe_radius << " "
           << " -density " << A_num_pts_per_area << " "
           << " -all_components >> " << msms_log_fname;

  if(system(cmd_line.str().c_str()) == -1){
    err_msg(A_fname, "get_MSMS_surf_patches()", 
            "Could not create protein MSMS surface");
    return false;
  }

  // parse the log file and determine the number of component files
  std::ifstream msms_log_file;
  if(!open_ifstream(msms_log_file, msms_log_fname)){
    err_msg(A_fname, "get_MSMS_surf_patches()", 
            "Unable to read the MSMS log file");
    return false;
  }
  std::string head_line_start = 
    "    RS component  #faces  #edges #free_edges  #vertices   genus";
  bool found_line = false;
  for(std::string line; std::getline(msms_log_file, line); ){
    size_t pos = line.find(head_line_start);
    if(pos != std::string::npos){
      found_line = true;
      break;
    }
  }
  if(!found_line){
    err_msg(A_fname, "get_MSMS_surf_patches()", 
            "The MSMS log file does not have an \"RS component\" table");
    msms_log_file.close();
    return false;
  }
  
  std::string end_line_start = "    Time Reduced Surface real:";
  int num_components = 0;
  found_line = false;
  for(std::string line; std::getline(msms_log_file, line); ++num_components){
    size_t pos = line.find(end_line_start);
    if(pos != std::string::npos){
      found_line = true;
      break;
    }
  }
  if(!found_line){
    err_msg(A_fname, "get_MSMS_surf_patches()", 
            "The MSMS log file does not have a complete \"RS component\" table");
    msms_log_file.close();
    return false;
  }
  std::cout << "Found " << num_components << " MSMS surface components" 
            << std::endl;
  msms_log_file.close();
  unlink(msms_log_fname);

  // Keep the component with the largest surface area in the binding site
  my_float_t max_SA = 0.0;
  for(int i = 0; i < num_components; ++i){
    std::ostringstream xyzr_stream;
    xyzr_stream << xyzr_fname;
    if(i) xyzr_stream << "_" << i;
    SimpleTrimesh full_mesh(xyzr_stream.str(), true);
    if(full_mesh.fail()) continue;
//    std::cout << "the full surface area for component " << i << " is: " 
//              << full_mesh.get_total_SA() << "\n";

    // determine the verts and faces near the binding site
    SimpleTrimesh *curr_comp_mesh = new SimpleTrimesh(full_mesh, *mesh_bv);
    if(curr_comp_mesh == 0 || curr_comp_mesh->fail()){
      if(curr_comp_mesh) delete curr_comp_mesh;
      continue;
    }

    // get surface area left for each component
    my_float_t curr_SA = curr_comp_mesh->get_total_SA();
//    std::cout << "the nearby surface area for component " << i << " is: " 
//              << curr_SA << "\n";
    if(curr_SA > max_SA){
      max_SA = curr_SA;
      if(old_binding_site_mesh) delete old_binding_site_mesh;
      old_binding_site_mesh = curr_comp_mesh;
    }else delete curr_comp_mesh;
    curr_comp_mesh = 0;
 
    // delete the component files
    unlink((xyzr_stream.str() + ".vert").c_str());
    unlink((xyzr_stream.str() + ".face").c_str());

  }
  // unlink the xyzr file  
//  unlink(xyzr_fname);

  delete mesh_bv;
  delete [] xyzr_fname;
  delete [] msms_log_fname;
  xyzr_fname = 0;
  msms_log_fname = 0;

  if(old_binding_site_mesh == 0) return false;
  return true;
}

Sitemap::~Sitemap()
{
  if(A_site_atoms) delete(A_site_atoms);
  if(A_bind_site_res) delete(A_bind_site_res);
  if(points) delete(points);
  if(prot_atoms) delete(prot_atoms);
  if(hbond_pts) delete(hbond_pts);
  if(hphob_pts) delete(hphob_pts);
  if(A_site_vol) delete(A_site_vol);
  if(A_site_vol_est) delete(A_site_vol_est);
  if(old_binding_site_mesh) delete(old_binding_site_mesh);
//  if(A_binding_site_mesh) delete(A_binding_site_mesh);
  if(A_hbond_surfaces) delete(A_hbond_surfaces);

  init();
}

void
Sitemap::init()
{
  A_site_atoms = 0;
  A_bind_site_res = 0;
  points = 0;
  prot_atoms = 0;
  hbond_pts = 0;
  hphob_pts = 0;
  A_site_vol = 0;
  A_site_vol_est = 0;
//  A_binding_site_mesh = 0;
  old_binding_site_mesh = 0;
  A_MSMS_surf = false;
  A_user_surf = false;
  A_hbond_surfaces = 0;
  A_num_pts_per_area = 1;
  A_min_res_per_chain = 1;
}

bool
Sitemap::get_norm_stats(const std::string csv_fname, const BaseParameters& args,
                        bool hydrophobic_query)
{
  std::ifstream labs_file;
  if(!open_ifstream(labs_file, csv_fname)) return false;

  // Get the normalization stats
  if(hydrophobic_query){
    for(std::string line; std::getline(labs_file, line); )
      if(line[0] == '<' && line.substr(0,18) == "<hphob_norm_stats>")
        break;
  }else{
    for(std::string line; std::getline(labs_file, line); )
      if(line[0] == '<' && line.substr(0,21) == "<normalization_stats>")
        break;
  }

  labs_file >> A_norm_moments[0] >> A_norm_moments[1];
  labs_file.close();
  if(A_norm_moments[1] > 0) return true;

  if(!normalize_sitemap(csv_fname, args, hydrophobic_query)) return false;
  return get_norm_stats(csv_fname, args, hydrophobic_query); 
}

BoundingVolume*
Sitemap::setup_volume(const std::string ligand_fname,
                      const my_float_t grid_spacing, 
                      const std::string sphere_str, 
                      const bool prune_to_lig)
{
  BoundingVolume *rv = NULL;  
  if(ligand_fname.length()){
    A_lig_fname = ligand_fname;
    my_float_t min[6];
    my_float_t *max = min + 3;
    if(prune_to_lig)
      rv = new UnionOfBalls(mol2File(ligand_fname));
    else if(gen_lig_box(ligand_fname, min))
      rv = new RectangularSolid(min, max, grid_spacing);
  }else if(sphere_str.length()){
    if(prune_to_lig)
      warn(A_fname, "get_bounding_volume", 
           "Cannot prune to ligand voluem if a spherical site volume is given", 
           std::cerr);
    A_lig_fname = "NONE";
    my_float_t tmp[4];
    std::istringstream sphere_args(sphere_str);
    for(uint i = 0; i < 4; i++) sphere_args >> tmp[i];
    if(sphere_args.fail())
      warn(A_fname, "get_bounding_volume", 
           "Error parsing center and radius for sphere", std::cerr);
    else rv = new Ball(tmp, tmp[3], grid_spacing);
  }

  return rv;
}

void
Sitemap::get_prot_interact_atoms(const PDBStructure* prot,
                                 BoundingVolume* site_vol,
                                 std::vector<atom_vci>* bind_atoms,
                                 std::vector<atom_vci>* rad_atoms, 
                                 interact_atoms_vec* hbond_atoms,
                                 interact_atoms_vec* hphob_atoms,
                                 hphob_triad_vec* hphob_triads,
                                 const int min_chain_size)
{
  // Find the protein atoms which are within binding and/or rad distance of a 
  // grid point.  If an atom can make hbonds and is within binding distance, 
  // add that atom to the list of possible hbonding atoms.
  chain_const_iter chain;
  for(chain = prot->chains_begin(); chain < prot->chains_end(); chain++){
    if(chain->residues_end - chain->residues_begin < min_chain_size){
      std::cout << "Chain " << chain->chainID << " has fewer than "
                << min_chain_size << " residues.  It is not considered to be "
                << "part of the structure when generating sitemaps\n";
      continue;
    }

    residue_vci residue = chain->residues_begin;
    for( ; residue < chain->residues_end; residue++){
      atom_vci atom = residue->atoms_begin;
      for( ; atom < residue->atoms_end; atom++){
        // %$%@%&*($&*^(&# hydrogen atoms
        if(atom->name == H || atom->name == D) continue;

        if(site_vol->BIND_vol_contains(atom->pos)){
          rad_atoms->push_back(atom);
          bind_atoms->push_back(atom);
          interact_atom_t atom_info;
          atom_info.chain = chain;
          atom_info.residue = residue;
          atom_info.atom = atom;
          if(atom->is_hbonder()) hbond_atoms->push_back(atom_info);
          else{
            const hphob_triad_t* triad = 0;
            if(HphobPoints::atom_is_hphob(atom->res, atom->name, &triad)){
              hphob_atoms->push_back(atom_info);
              hphob_triads->push_back(triad);
            }
          }
        }else if(site_vol->RAD_vol_contains(atom->pos)) 
          rad_atoms->push_back(atom);
      }
    }
  }
}

void
Sitemap::find_binding_site_hbonds(const GenPointsParameters::hbond_method_t 
                                  hbond_density,
                                  const std::string params_dir,
                                  const std::vector<std::string>& waters,
                                  BoundingVolume* site_vol,
                                  const interact_atoms_vec& hbond_atoms,
                                  std::vector<atom_vci>* bind_atoms,
                                  std::vector<atom_vci>* rad_atoms)
{
  hbond_pts = new HbondPoints(hbond_density, params_dir, site_vol);
  hbond_pts->find_binding_site_hbonds(prot_atoms, waters, A_include_metals, 
                                      hbond_atoms, bind_atoms, rad_atoms);
  points->add_points(hbond_pts->fit_pts_beg(), hbond_pts->fit_pts_end());
  hbond_ideal_pt_vci p = hbond_pts->ideal_pts_beg();
  for( ; p < hbond_pts->ideal_pts_end(); ++p){
    interact_atoms[p->atom->atom_num] = p->atom;
    if(p->carbon_nbr != HbondPoints::NO_CARBON_NEIGHBOR)
      interact_atoms[p->carbon_nbr->atom_num] = p->carbon_nbr;
  }
}

bool
Sitemap::find_binding_site_hphobs(const GenPointsParameters::hphob_method_t 
                                  hphob_method,
                                  const my_float_t cluster_diameter,
                                  const interact_atoms_vec& hphob_atoms,
                                  const hphob_triad_vec& hphob_triads,
                                  const std::vector<atom_vci>& bind_atoms, 
                                  BoundingVolume* site_vol,
                                  const sphere_sample_level_t sample_level)
{
  hphob_pts = new HphobPoints(hphob_method, cluster_diameter, sample_level);
  if(!hphob_pts->gen_points(hphob_atoms, hphob_triads, bind_atoms, site_vol))
    return false;
  if(verbose){
    std::cout << "  " << hphob_pts->size()
              << " hydrophobic points before culling those too close to polar "
              << "points\n\n";
  }
  hphob_pts->cull_too_close_to_polar(hbond_pts->fit_pts_beg(), 
                                     hbond_pts->fit_pts_end());
  points->add_points(hphob_pts->begin(), hphob_pts->end());
  for(hphob_point_lci p = hphob_pts->begin(); p != hphob_pts->end(); ++p){
    std::list<atom_vci>::const_iterator li; 
    for(li = p->atoms.begin(); li != p->atoms.end(); ++li)
      interact_atoms[(*li)->atom_num] = *li;
  }
  return true;
}

bool 
Sitemap::has_enough_points()
{
  int num_hphob_pts = hphob_pts->size();
  int num_hphill_pts = hbond_pts->fit_pts_end() - hbond_pts->fit_pts_beg();
  int npts = num_hphob_pts + num_hphill_pts;
  if(npts < MIN_NUMBER_POINTS__POLAR_SITE ||
     (num_hphill_pts < 6 && npts < MIN_NUMBER_POINTS__HPHOB_SITE)){
    std::string msg = "Site map has too few points.  This is likely to ";
    msg += "cause inconsistent results";
    err_msg(A_fname, "has_enough_points", msg);
    return false;
  }
  return true;
}

bool
Sitemap::write_files(const std::string struct_id, const std::string path) 
{
  if(A_fail){
    std::cerr << "Could not write the sitemap files -- an error has occured "
              << "while generating the sitemap\n";
    return false;
  }

  // Write the binding site mesh files 
  if((A_MSMS_surf || A_user_surf) && 
     !old_binding_site_mesh->write(path + "/" + struct_id + "_surf")){
    std::cerr << "could not write binding site mesh\n";
    return false;
  }

  // Write the sitemap points file (the sitemap file)
  if(!points->write(path + "/" + struct_id + "_s.pdb")) return false;

  // Open file in append mode and add "hbvecs"
  std::ofstream out;
  std::string sitemap_fname = path + "/" + struct_id + "_s.pdb"; 
  out.open(sitemap_fname.c_str(), std::ios_base::app|std::ios_base::ate); 
  if(out.fail()){
    std::cerr << "Unable to open the output file name: " << sitemap_fname 
              << "\n";
    return false;
  }
  out.setf(std::ios_base::fixed, std::ios_base::floatfield);
  out.precision(3);

  uint cnt = 4000; 
  my_float_t dir[] = {0, 0, 0};
  hbond_fit_pt_vci fit_pt = hbond_pts->fit_pts_beg();
  for( ; fit_pt < hbond_pts->fit_pts_end(); ++fit_pt){
    ++cnt;
    unit_vector(dir, fit_pt->pos, fit_pt->atom->pos);
    out << "#HBVEC" << std::setw(5) << cnt << "  O   HOH  "
        << std::setw(4) << cnt << "    "
        << std::setw(8) << dir[0] << std::setw(8) << dir[1]
        << std::setw(8) << dir[2] << "  1.00200.00\n";
  } 
  out.close();

  // Write the sitemap "rad" file
  if(A_site_atoms == 0){
    std::cerr << "Site atoms object is NULL -- RAD file may contain short peptides if they are near the binding sites\n";
    if(A_lig_fname != "NONE") 
      A_site_atoms = 
        new SitemapAtomsFile(prot_atoms, interact_atoms, A_lig_fname, 1);
    else 
      A_site_atoms = 
        new SitemapAtomsFile(prot_atoms, A_site_vol, 1, interact_atoms);
  }
  if(!A_site_atoms->write(path + "/" + struct_id + "_rad.pdb")) return false;
//  delete A_site_atoms;

  // Write the sitemap atoms file 
  SitemapAtomsFile act_atoms(interact_atoms);
  if(!act_atoms.write(path + "/" + struct_id + "_a.pdb")) return false;

  // Write the glue file
  if(!open_ofstream(out, path + "/" + struct_id + "_s.csv")) return false;
  out.setf(std::ios_base::fixed, std::ios_base::floatfield);
  out.precision(3);

  // Header
  out << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n" 
      << "<!-- Note: no xml namespace for now -->\n"
      << "<ASCbase_sitemap>\n"
      << "<!--\n"
      << "# gen_points (" << PACKAGE_NAME << ") " << PACKAGE_VERSION
      << " sitemap labels file\n";
  char time_buf[80];
  out << "# Date: " << Timer::get_local_time(time_buf, 80) << "\n"
      << "-->\n";

  out << "<protein_structure_file>\n  ";
  std::string::size_type pos = prot_atoms->name().find_last_of("/");
  if(pos == std::string::npos) out << prot_atoms->name() << "\n";
  else out << prot_atoms->name().substr(pos + 1) << "\n";
  out << "</protein_structure_file>\n";

  out << "<ligand_file>\n  ";
  pos = A_lig_fname.find_last_of("/");
  if(pos == std::string::npos) out << A_lig_fname << "\n";
  else out << A_lig_fname.substr(pos + 1) << "\n";
  out << "</ligand_file>\n";

  out << "<mol_surf>\n";
  if(A_MSMS_surf)
    out << "  <probe_radius>" << A_probe_radius << "</probe_radius>\n"    
        << "  <density>" << A_num_pts_per_area << "</density>\n";
  out << "</mol_surf>\n";

  out << "<site_vol_est>\n";
  size_t npts = hbond_pts->fit_pts_end() - hbond_pts->fit_pts_beg();
  npts += hphob_pts->size();
  get_site_volume_estimate(hbond_pts->fit_pts_beg(), 
                           hbond_pts->fit_pts_end(), hphob_pts->begin(), 
                           hphob_pts->end(), npts, 1.5);
  out << A_site_vol_est->xml_str() << "\n"
      << "</site_vol_est>\n";

  out << "<min_res_per_chain>\n" << A_min_res_per_chain 
      << "\n</min_res_per_chain>\n";

  // change this
  out << "<volume_type>\n" << A_site_vol->xml_str() << "\n</volume_type>\n";

  out << "<!--\n"
      << "# Normalization stats:\n"
      << "# 1) Mean score of sitemap against normalization database\n"
      << "# 2) Stdev of scores of sitemap against normalization database\n#\n"
      << "# Note:  If the mean is 0 and the stdev is -1, the sitemap "
      << "normalization\n#        statistics were not computed.\n"
      << "-->\n";
  out << "<normalization_stats>\n  " 
      << A_norm_moments[0] << " " << A_norm_moments[1] << "\n"
      << "</normalization_stats>\n";

  // Hmm, looks like some or all XML parsers don't like "--" inside comment
  // sections
  out << "<!--\n"
      << "# Hydrophobic pocket normalization stats:\n"
      << "# Normalization stats using hydrophobic_query flag\n"
      << "# 1) Mean score of sitemap against normalization database\n"
      << "# 2) Stdev of scores of sitemap against normalization database\n#\n"
      << "# Note:  If the mean is 0 and the stdev is -1, the sitemap "
      << "normalization\n#        statistics were not computed using the "
      << "hydrophobic_query flag.\n"
      << "-->\n";
  out << "<hphob_norm_stats>\n" 
      << "  0.0 -1.0\n"
      << "</hphob_norm_stats>\n";

  out << "<!--\n"
      << "# Point counts:\n# 1) Number of polar points\n"
      << "# 2) Number of hydrophobic points\n"
      << "-->\n";
  out << "<point_counts>\n  " 
      << points->num_polar_points()
      << " " << points->num_hphobic_points() << "\n"
      << "</point_counts>\n";

  // Hbond "fit" points
  out << "<!--\n"
      << "# Fields in hydrogen bond section\n"
      << "# 1) Sitemap point number\n"
      << "# 2) Associated atom string (res, res num, iCode, chain id, atom "
      << "name and\n#    sitemap point descriptor)\n"
      << "# 3) Chain ID\n"
      << "# 4) Residue name\n"
      << "# 5) Residue number + insertion code (if it is not a space char)\n"
      << "# 6) Atom name\n"
      << "# 7) Atom number (serial)\n"
      << "# 8) Hbond point number (for the given atom) - ASCbase internal "
      << "use\n"
      << "# 9) Point position  (ASCbase internal use)\n"
      << "# 10) Carbon neighbor atom number  (ASCbase internal use)\n"
      << "-->\n";
  out  << "<hydrogen_bond_points>\n";
  uint counts[] = {0, 0, 0};
  uint serial_adds[] = {1000, 2000, 3000};

  atom_vci prev_atom = atom_t::NULL_ATOM_VCI;
  uint pts_num[] = {1, 1, 1}; 
  hbond_ideal_pt_vci ideal_pt = hbond_pts->ideal_pts_beg();
  for( ; ideal_pt < hbond_pts->ideal_pts_end(); ++ideal_pt){
    const atom_t& a = *(ideal_pt->atom);
    if(prev_atom != ideal_pt->atom){
      std::fill(pts_num, pts_num + 3, 1);
      prev_atom = ideal_pt->atom;
    }

    std::string atom_name;
    std::string res_name;
    uint carbon_nbr_atom_num = 0;
    if(a.is_hetero){
      atom_name = a.name_str;
      res_name = a.res_str;
    }else{
      atom_name = PDB_residues::atom_to_string(a.name);
      res_name = PDB_residues::residue_to_string(a.res);
      carbon_nbr_atom_num = ideal_pt->carbon_nbr->atom_num; 
    }
    clean_up_string(atom_name);

    uint* serial = 0; 
    uint serial_add = 0;
    char lab = ' ';
    uint* fit_pt_num = 0;
    if(ideal_pt->act_type == ACCEPTOR){
      serial = &(counts[0]);
      serial_add = serial_adds[0];
      lab = 'A';
      fit_pt_num = pts_num;
    }else if(ideal_pt->act_type == DONOR){
      serial = &(counts[1]);
      serial_add = serial_adds[1];
      lab = 'D';
      fit_pt_num = pts_num + 1;
    }else if(ideal_pt->act_type == DONEPTOR){
      serial = &(counts[2]);
      serial_add = serial_adds[2];
      lab = 'N';
      fit_pt_num = pts_num + 2;
    }else{
      warn(A_fname, "write_files()", "Unknown hydrogen bond point type");
      return false;
    }

    fit_pt = ideal_pt->fit_pts_beg;
    for( ; fit_pt < ideal_pt->fit_pts_end; ++fit_pt){
      ++(*serial);
      uint pt_num = *serial + serial_add;
      std::string iCode, chainID;
      if(a.iCode != ' ') iCode += a.iCode;
      if(a.chainID != ' ') chainID += a.chainID;

      out << pt_num << "|" 
          << res_name << a.res_num << iCode << "_" << chainID << "_" 
          << atom_name << "/" << *fit_pt_num << lab << "|" 
          << a.chainID << "|" << res_name << "|" << a.res_num << iCode << "|" 
          << atom_name << "|" << a.atom_num << "|" << ideal_pt->pt_num << "|";
      for(uint j = 0; j < 3; ++j) out << ideal_pt->pos[j] << " "; 
      out << "|" << carbon_nbr_atom_num << "|\n";
      ++(*fit_pt_num);
    }
  }
  out << "</hydrogen_bond_points>\n";

  // Hydrophobic points
  out << "<!--\n"
      << "# Fields in hydrophobic points section\n"
      << "# 1) Sitemap point number\n"
      << "# 2) Associated atoms' strings\n"
      << "-->\n";
  out << "<hydrophobic_points>\n";
  cnt = 0;
  for(hphob_point_lci p = hphob_pts->begin(); p != hphob_pts->end(); ++p){
    std::list<atom_vci>::const_iterator li;
    ++cnt;
    bool first = true;

    out << cnt << "|";
    for(li = p->atoms.begin(); li != p->atoms.end(); ++li){
      const atom_t& a = **li;
      std::string iCode, chainID;
      if(a.iCode != ' ') iCode += a.iCode;
      if(a.chainID != ' ') chainID += a.chainID;
      std::string atom_name;
      std::string res_name;
      if(a.is_hetero){
        atom_name = a.name_str;
        res_name = a.res_str;
      }else{
        atom_name = PDB_residues::atom_to_string(a.name);
        res_name = PDB_residues::residue_to_string(a.res);
      }
      clean_up_string(atom_name);

      if(first) first = false;
      else out << "_";
      out << res_name << a.res_num << iCode << "_" << chainID << "_" 
          << atom_name;
    }
    out << "|\n";
  }
  out << "</hydrophobic_points>\n";
  write_feature_scales(out);
  out << "</ASCbase_sitemap>\n";
  return true;
}

void
Sitemap::clean_up_string(std::string& S)
{
  for(std::string::iterator s = S.begin(); s < S.end(); ++s)
    if(isspace(*s)){
      s = S.erase(s);
      --s;
    }
}

void
Sitemap::build_points_header(const GenPointsParameters& args)
{
  std::ostringstream header;
  std::string REMARK_10 = "REMARK  10 ";
  // Print a short header 
  header << REMARK_10 << "Sitemap generated by ASCbase\n"
           << REMARK_10 << "Protein: " << args.prot_fname << "\n";
  if(args.lig_fname.length())
    header << REMARK_10 << "Sitemap Generation method: ligand bounding box\n"
             << REMARK_10 << "  Ligand: " << args.lig_fname << "\n";
  else if(args.sphere_str.length())
    header << REMARK_10 
             << "Sitemap Generation method: ball (point and radius)\n"
             << REMARK_10 << "  Params: " << args.sphere_str << "\n";

  header << REMARK_10 << "Hbond points density: ";
  switch(args.hbond_method){
  case GenPointsParameters::OPTIMUM_HBONDS:
    header << "optimum (one point per H or lone pair)";
    break; 
  case GenPointsParameters::MIN_HBONDS:
    header << "minimal";
    break; 
  case GenPointsParameters::SPARSE_HBONDS:
    header << "sparse";
    break; 
  case GenPointsParameters::DENSE_HBONDS:
    header <<  "dense";
    break; 
  default:
    // should never get here
    break;
  } 
  header << "\n";
  
  if(args.hphob_method == GenPointsParameters::THREE_PROTEIN_ATOMS)
    header << REMARK_10 << "Hydrophobic method: 3 nearby protein atoms\n"
             << REMARK_10 << "Grid spacing: " << args.grid_spacing 
             << " (A)\n"
             << REMARK_10 << "Cluster threshold: "
             << args.cluster_diameter << " (A)\n";
  else if(args.hphob_method == GenPointsParameters::THREE_MORE_HPHOB)
    header << REMARK_10 << "Hydrophobic method: 3 more hphob neighbors than "
             << "hphill\n" 
             << REMARK_10 << "Grid spacing: " << args.grid_spacing
             << " (A)\n"
             << REMARK_10 << "Cluster threshold: "
             << args.cluster_diameter << " (A)\n";
  else if(args.hphob_method == GenPointsParameters::ATOM_CENTERS)
    header << REMARK_10 << "Hydrophobic method: Hydrophobic atom centers\n";
  else if(args.hphob_method == GenPointsParameters::PSEUDO_SURFACE)
    header << REMARK_10 
           << "Hydrohobic method: hydrophobic spheres (sample level " 
           << args.sphere_sample_level << ")\n"; 

  points->set_header(header.str());
}

void 
Sitemap::write_feature_scales(std::ostream &out)
{
  out << "<!--\n"
      << "# Note: These are the largest values for the features of the sitemap"
      << " when\n # compared with itself\n"
      << "# Features\n"
      << "# 0) Closest polar sum\n"
      << "# 1) polar mismatch sum\n"
      << "# 2) AA DD sum\n"
      << "# 3) NA AN ND DN NN sum\n"
      << "# 4) hphob count\n"
      << "# 5) Unsat query polar \n"
      << "# 6) best polar sum\n"
      << "# 7) best AA DD sum\n"
      << "# 8) best NA AN ND DN NN sum\n"
      << "# 9) best unsat polar\n";
  if(A_MSMS_surf || A_user_surf)
    out << "# 10) number of query surface points within 1.5 (A) of db surface\n"
        << "# 11) number of query faces with all 3 vertices within 1.5 (A) of "
        << "db surface\n"
        << "# 12) RMSE of query points to db surf (each point must be within "
        << "1.5 (A) )\n"
        << "# 13) area of faces in 11\n"
        << "# 14) RMSE of query normals to db surf normals\n";
  out << "-->\n"
      << "<features_scales>\n";

  my_float_t num_hbonds = hbond_pts->fit_pts_end() - hbond_pts->fit_pts_beg();
  if(num_hbonds < 1.0) num_hbonds = 1.0;
  out << num_hbonds << " " << num_hbonds << " ";

  my_float_t num_AA_DD = 0.0;
  hbond_fit_pt_vci fit_pt = hbond_pts->fit_pts_beg();
  for( ; fit_pt < hbond_pts->fit_pts_end(); ++fit_pt)
    if(fit_pt->act_type == ACCEPTOR || fit_pt->act_type == DONOR) 
      num_AA_DD += 1.0;
  if(num_AA_DD < 1.0) num_AA_DD = 1.0;
  out << num_AA_DD << " " << num_hbonds << " ";
  if(hphob_pts->size() < 1) out << "1.0 ";
  else out << hphob_pts->size() << " " ;
  out << num_hbonds << " " << num_hbonds << " " << num_AA_DD << " " 
      << num_hbonds << " " << num_hbonds << "\n";
  if(A_MSMS_surf || A_user_surf){
    // The feature 3rd feature used to have a max val of 1.0, currently
    // it is 1.5 and it probably should go up
    out << old_binding_site_mesh->number_of_vertices() << " "
        << old_binding_site_mesh->number_of_faces() << " 1.5 "
        << old_binding_site_mesh->get_total_SA();
    // Currently max value for norm penalty is 4.0
    out << " 4.0\n";

  }
  out << "</features_scales>\n";
}

bool
Sitemap::gen_lig_box(const std::string fname, my_float_t* corners)
{
  corners[0] = corners[1] = corners[2] = my_float_max;
  corners[3] = corners[4] = corners[5] = -corners[0];

  CoordFile* lig = 0;
  if(!load_ligand(fname, &lig)) return false;

  // Get the ligand bounding box -- skip the hydrogen atoms
  for(atom_vci a = lig->atoms_begin(); a < lig->atoms_end(); ++a){
    // Skip hydrogen atoms
    if(a->name == H || a->name == D) continue;

    my_float_t* min = corners;
    my_float_t* max = corners + 3;
    const my_float_t* p = a->pos;
    for(uint i = 0; i < 3; ++i, ++p, ++min, ++max){ 
      if(*min > *p) *min = *p;
      if(*max < *p) *max = *p;
    }
  }
  delete lig;

  // Add the tolerance to the bounding box
  for(uint i = 0; i < 3; i++){
    corners[i] -= TOLERANCE;
    corners[3 + i] += TOLERANCE;
  }
  return true;
}

bool
Sitemap::load_volume(std::string path, std::string struct_id, 
                     std::string dbase_sites, std::string vol_tag,
                     BoundingVolume **vol, geometry::sphere_t *S)
{
  if(*vol) delete *vol;
  *vol = 0;

  std::string my_site = path + "/" + struct_id + "_s.csv"; 
  if(!normal_file_exists(my_site, false))
    my_site = dbase_sites + "/" + struct_id + "_s.csv";

  std::ifstream labs_file;
  if(!open_ifstream(labs_file, my_site)){
    A_fail = true;
    return false;
  }
  
  for(std::string line; std::getline(labs_file, line); )
    if(line[0] == '<' && line.substr(0,vol_tag.length()) == vol_tag){
      if(!std::getline(labs_file, line)) break;

      if("<sphere>" == line.substr(0, 8)){
        // Get the center
        size_t pos = line.find("<center>");
        if(pos == std::string::npos) break;
        pos += 8;
        std::string s = line.substr(pos);
        eat_leading_whitespace(&s);
        std::vector<std::string> toks; 
        string_tok(s, &toks, ' ');
        my_float_t C[3];
        for(size_t i = 0; i < 3; ++i) my_strtof(toks[i], &(C[i])); 
        
        // Get the radius
        pos = line.find("<radius>");
        if(pos == std::string::npos) break;
        pos += 8;
        s = line.substr(pos);
        eat_leading_whitespace(&s);
        my_float_t r;      
        my_strtof(s, &r);
  
        // If a sitemap point happens to be on the surface of the ball,
        // it can match stuff upto 1.5 (A) outside the ball -- but this
        // tolerance is already added when the volume was estimated
        *vol = new Ball(C, r);
        if(S) *S = geometry::sphere_t(C,r);
        break;
      }else if("<hyperrectangle>" == line.substr(0, 16)){
        // Our initial estimate is that 3.0 is a good balance between
        // the query ligand volume and ligand fragments
        *vol = new UnionOfBalls(mol2File(ligand_file_name()), 3.0);
        break;
      }

    }

  if(*vol) return true;
  err_msg(A_fname, "load_volume", "Unable to load the site volume");
  return false;
}

void
Sitemap::get_site_volume_estimate(hbond_fit_pt_vci hbond_pts_beg,
                                  hbond_fit_pt_vci hbond_pts_end, 
                                  hphob_point_lci hphob_points_beg,
                                  hphob_point_lci hphob_points_end, size_t npts,
                                  my_float_t vol_est_tol)
{
  my_float_t *positions = new my_float_t[3*npts];
  my_float_t *pos = positions;

  hbond_fit_pt_vci fit_pt = hbond_pts_beg;
  for( ; fit_pt < hbond_pts_end; ++fit_pt, pos += 3)
    std::copy(fit_pt->pos, fit_pt->pos + 3, pos);

  hphob_point_lci p;
  for(p = hphob_points_beg; p != hphob_points_end; ++p, pos += 3)
    std::copy(p->pos, p->pos + 3, pos);

  size_t i, j;
  my_float_t dist;
  max_pair_dist(positions, npts, &i, &j, &dist);

  my_float_t center[3];
  for(size_t z = 0; z < 3; ++z)
    center[z] = (positions[3*i + z] + positions[3*j + z]) / 2.0;
  my_float_t radius = dist / 2.0 + vol_est_tol;
  if(A_site_vol_est){
    std::cerr << "The estimate of the site volume has been compute "
              << "previously\n";
    delete A_site_vol_est;
    A_site_vol_est = 0;
  }
  A_site_vol_est = new Ball(center, radius);

  delete [] positions;
  pos = positions = 0;
}

bool
Sitemap::load_ligand(const std::string &fname, CoordFile** lig)
{
  // Check if this is pdb or mol2 file
  std::string::size_type pos = fname.rfind(".");
  if(pos == std::string::npos){
    std::cerr << "The ligand file name must include either the .pdb or .mol2 "
              << "file extension as part of its name\n";
    return false;
  }

  // Load the ligand file
  *lig = 0;
  std::string ext = fname.substr(pos + 1);
  if(ext == "pdb") *lig = new PDBBase(fname);
  else if(ext == "mol2") *lig = new mol2File(fname);
  else{
    std::cerr << "The ligand file name must include either the .pdb or .mol2 "
              << "file extension as part of its name\n";
    return false;
  }if(*lig == 0){
    std::cerr << "Unable to allocate memory to open the ligand file\n";
    return false;
  }

  return true;
}

bool
Sitemap::get_prot_atom_close_to_ligand(const std::string ligand_fname,
                                       std::vector<atom_vci>& bind_atoms,
                                       atom_vci* close_atom)
{
  CoordFile* lig = 0;
  if(!load_ligand(ligand_fname, &lig)){
    err_msg(A_fname, "get_prot_atom_close_to_ligand",
            "Unable to create surface, cannot load the ligand: " + ligand_fname);
    return false;
  }

  std::vector<atom_vci>::const_iterator b_iter;
  *close_atom = *(bind_atoms.begin());
  my_float_t best_sq_dist = my_float_max;
  const my_float_t sq_tol = 2.8*2.8;
  for(b_iter = bind_atoms.begin(); b_iter < bind_atoms.end(); ++b_iter){
    for(atom_vci a = lig->atoms_begin(); a < lig->atoms_end(); ++a){
      // Skip hydrogen atoms
      if(a->name == H || a->name == D) continue;

      my_float_t sq_dist = dist_squared((*b_iter)->pos, a->pos);
      if(sq_dist < best_sq_dist){
        best_sq_dist = sq_dist;
        *close_atom = *b_iter; 
        if(sq_dist <= sq_tol) return true;
      }
    }
  }

  delete lig;
  return true;
}

void
Sitemap::compute_max_feature_values(std::vector<my_float_t> *max_feat_vals) 
const
{
  max_feat_vals->reserve(14);
  my_float_t num_hbonds = hbond_pts->fit_pts_end() - hbond_pts->fit_pts_beg();
  if(num_hbonds < 1.0) num_hbonds = 1.0;
  max_feat_vals->push_back(num_hbonds);
  max_feat_vals->push_back(num_hbonds);

  my_float_t num_AA_DD = 0.0;
  hbond_fit_pt_vci fit_pt = hbond_pts->fit_pts_beg();
  for( ; fit_pt < hbond_pts->fit_pts_end(); ++fit_pt)
    if(fit_pt->act_type == ACCEPTOR || fit_pt->act_type == DONOR) 
      num_AA_DD += 1.0;
  if(num_AA_DD < 1.0) num_AA_DD = 1.0;
  max_feat_vals->push_back(num_AA_DD);
  max_feat_vals->push_back(num_hbonds);

  if(hphob_pts->size() < 1) max_feat_vals->push_back(1.0);
  else max_feat_vals->push_back(hphob_pts->size());
  max_feat_vals->push_back(num_hbonds);
  max_feat_vals->push_back(num_hbonds);
  max_feat_vals->push_back(num_AA_DD);
  max_feat_vals->push_back(num_hbonds);
  max_feat_vals->push_back(num_hbonds);
}

my_float_t
Sitemap::compute_site_vol_intersect(const Sitemap& other) const
{ 
  if(A_site_vol_est_as_sphere.radius() <= 0.0 && 
     other.A_site_vol_est_as_sphere.radius() <= 0.0) return -10.0;
  
  my_float_t vol;
  geometry::sphere_t::intersectionType rv = 
    A_site_vol_est_as_sphere.intersects(other.A_site_vol_est_as_sphere, &vol);
  if(rv == geometry::sphere_t::NONE) return -10.0;
  const my_float_t &r = A_site_vol_est_as_sphere.radius();
  return vol / (4.0/3.0 * M_PI * r*r*r);
}
