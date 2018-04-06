/******************************************************************************
 * Copyright (c) 2006,2007, Michigan State University (MSU) Board of Trustees.
 *   All rights reserved.
 *
 * This file is part of the SimSite3D Software project.
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

#include <sstream>
#include <cstring>
#include <popt.h>
#include <GenPointsParameters.H>
#include <param_tools.H>

using namespace SimSite3D;
const std::string GenPointsParameters::A_fname = "GenPointsParameters.C";

GenPointsParameters::GenPointsParameters() : BaseParameters() 
{
  init();
}

GenPointsParameters::GenPointsParameters(const int argc, const char** argv)
 : BaseParameters()
{
  init();
  // Check if anything has failed in BaseParameters
  if(A_status == FATAL_ERROR) return;

  A_status = get_opts(argc, argv);
  if(A_status != INITIALIZING) return;
  A_status = verify_params();
}

GenPointsParameters::~GenPointsParameters()
{
  free_cstrings();
}

void
GenPointsParameters::init()
{
  verbose_level = VERBOSE_ONE;
  prots_dir = NULL;
  ligs_dir = NULL;
  A_scratch_dir = NULL;
  A_proj_output = NULL;
  A_waters_str = NULL;
  hbond_dens_str = NULL;
  cluster_diameter = 3.5;
  grid_spacing = 0.5;
  hphob_method = PSEUDO_SURFACE;
  // Need more testing to determin the default level
  sphere_sample_level = DISCRETE_SPHERE_LEVEL_TWO;
  //sphere_sample_level_int = 0;
  hbond_method = SPARSE_HBONDS;
  normalize = true;
  include_metals = false;
  call_msms = false;
  prune_to_lig = true;
  probe_radius = 1.4;
  num_pts_per_area = 1;
  min_res_per_chain = 1;
  msms_binary = "default";
}

void 
GenPointsParameters::free_cstrings()
{
  if(prots_dir) free(prots_dir);
  if(ligs_dir) free(ligs_dir);
  if(A_scratch_dir) free(A_scratch_dir);
  if(A_proj_output) free(A_proj_output);
  if(A_waters_str) free(A_waters_str);
  if(hbond_dens_str) free(hbond_dens_str);
  init();
}

BaseParameters::status_t
GenPointsParameters::get_opts(const int argc, const char** argv)
{
  char* pts_fname_Cstr = 0;
  char* prot_fname_Cstr = 0;
  char* lig_fname_Cstr = 0;
  char* sphere_Cstr = 0;
  char* user_provided_surf_Cstr = 0;
  char* msms_bin_Cstr = 0;

  // Setup the popt options table
  // At present, the flags with no arg are only those used to determine
  // the hydrophobic method; we start counting from 1 since 0 implies
  // no return.
  struct poptOption supportedOptionsTable[] = {
    { "version", 'v', POPT_ARG_NONE, 0, DISPLAY_VERSION, 
      "Prints version of software and exits"},
    { "protein", 'p', POPT_ARG_STRING, &(prot_fname_Cstr), 0,
      "Path to protein PDB file used to generate the sitemap", "XXXXXX_p.pdb"},
    { "ligand", 'l', POPT_ARG_STRING, &(lig_fname_Cstr), 0,
      "mol2 file of ligand to use to define the sitemap volume", 
      "XXXXXX_l.mol2"},
    { "sphere", '\0', POPT_ARG_STRING, &(sphere_Cstr), 0,
      "Generate sitemap volume as intersection of sphere and solvent volume",
      "\"X Y Z radius\""},
    { "use_lig_bound_box", '\0', POPT_ARG_NONE, 0, 
      USE_BOUNDING_BOX_FOR_SITE_VOLUME,
      "Prune site map points, so that the remaining are within 3.0 (A) of the minimal axis aligned bounding box containing all ligand atoms"},
    { "surf_vert_file", '\0', POPT_ARG_STRING, &user_provided_surf_Cstr, 0,
      "path to the vertex file for the site's molecular surface"},
    { "proj_output", '\0', POPT_ARG_STRING, &A_proj_output, 0,
      "Set directory to hold the generated sitemap", "path/to/sitemap"},
    { "msms_surf", '\0', POPT_ARG_NONE, 0, MSMS_SURFACE,
      "Use Michel Sanner's msms to create pocket surface ($SIMSITE3D_INSTALL_DIR/bin/linux_msms) -- noncommerical use ONLY"},
    { "probe_radius", '\0', POPT_ARG_DOUBLE, &(probe_radius), 0, 
      "The probe radius to use for the MSMS surface" },
    { "num_pts_per_area", '\0', POPT_ARG_DOUBLE, &(num_pts_per_area), 0, 
      "The number of vertices per (A)^2"},
    POPT_TABLEEND
  };

  struct poptOption otherOptionsTable[] = {
    { "msms_binary", '\0', POPT_ARG_STRING, &(msms_bin_Cstr), 0, 
      "Path to your msms binary", "/path/to/your/msms"},
    { "waters", '\0', POPT_ARG_STRING, &(A_waters_str), 0,
      "Consider the listed water residues in the protein structure as \"part\" of the protein",
      "\"H2O residues\""},
    { "include_metals", '\0', POPT_ARG_NONE, 0, INCLUDE_METALS,
      "Consider ALL metal atoms in the protein structure as \"part\" of the protein"},
    { "allow_small_site_maps", '\0', POPT_ARG_NONE, 0, ALLOW_SMALL_SITEMAPS,
      "Allow site maps to have fewer than the recommended minimum number of points"},
    {"min_res_per_chain", '\0', POPT_ARG_INT, &(min_res_per_chain), 0,
     "Minimum number of residues per protein chain -- default is 1"},
    { "no_normalization", '\0', POPT_ARG_NONE, 0, NO_NORMALIZATION,
      "Do not normalize the generated sitemap as a query with respect to the diverse sitemaps"},
    { "scratch_dir", '\0', POPT_ARG_STRING, &A_scratch_dir, 0,
      "Set name of scratch directory", "path/to/tmp"},
    POPT_TABLEEND
  };

#if 0
  struct poptOption testingOptionsTable[] = {
    { "hbond_density", '\0', POPT_ARG_STRING, &(hbond_dens_str), 0,
      "Hueristic for number of Hbond points using the SLIDE method",
      "density"},
    { "grid_spacing", '\0', POPT_ARG_DOUBLE, &(grid_spacing), 0,
      "The grid spacing (default is 0.5(A))"},
    { "cluster_diameter", '\0', POPT_ARG_DOUBLE, &(cluster_diameter),
      0, "The maximal cluster diameter for complete linkage of hphobs"},
    { "hphob_baseline", '\0', POPT_ARG_NONE, 0, BASELINE_HPHOB_REP,
      "Use the absense of polar atoms as hydrophobic locations"},
    { "hphob_pseudo_surf", '\0', POPT_ARG_INT | POPT_ARGFLAG_OPTIONAL,
       &(sphere_sample_level), PSEUDO_SURF_HPHOB_REP,
      "Use sphere around each hydrophobic atom"},
    POPT_TABLEEND
  };
#endif

  struct poptOption mainOptsTable[] = {
    {NULL, '\0', POPT_ARG_INCLUDE_TABLE, &supportedOptionsTable, 0,
     "Supported options", NULL},
    {NULL, '\0', POPT_ARG_INCLUDE_TABLE, &otherOptionsTable, 0,
     "Expert user options", NULL},
// All of these options no longer serve a purpose unless one is seeking to 
// change the point density (hbond_density or hphob_pseudo_surf), but
// our results have shown that increasing the point density has little payoff
// and decreasing the point density has very drastic consequences.
//    {NULL, '\0', POPT_ARG_INCLUDE_TABLE, &testingOptionsTable, 0,
//     "Software testing and development options\n\t** USE AT YOUR OWN RISK **", NULL},
    POPT_AUTOHELP
    POPT_TABLEEND
  };

  std::ostringstream ostr;
  ostr << "\nThat is " << *argv << " [OPTIONS]* [XXXXXX_s.pdb]\n"
       << "Typically\n  gen_points -p <XXXXXX_p>.pdb -l <XXXXXX_l>.mol2 "
       << "[path_to_store_sitemap]\n\t\t- OR -\n"
       << "  gen_points -p <XXXXXX_p>.pdb --sphere \"Cx Cy Cz r\" "
       << "[path_to_store_sitemap]\n"
       << "    where C = (Cx,Cy,Cz) is the center of the sphere and "
       << "r is the radius"
       << "When a ligand is given, each sitemap point is required to be "
       << "within\n3.0 (A) of at least one ligand heavy atom.\n"
       << "If you wish to override this functionality, use the\n"
       << "--use_lig_bound_box flag or provide a sphere instead\n";
  poptContext optCon = poptGetContext(argv[0], argc, argv, mainOptsTable, 0);
  poptSetOtherOptionHelp(optCon, ostr.str().c_str());

  if(argc < 2) {
    poptPrintUsage(optCon, stderr, 0);
    return DISPLAY_HELP_ONLY;
  }

  int rc;
  // Process the options (hydrophobic point method)
  for(rc = 0; (rc = poptGetNextOpt(optCon)) >= 0; ){
    switch(rc){
    case POPT_NULL_FLAG:
      // Should not get zero returned
      fprintf(stderr, "Error in processing command line arguments\n");
      return FATAL_ERROR;
    case DISPLAY_VERSION:
      std::cout << "\n\nCopyright (C) 2006-2008, Michigan State University "
                << "(MSU) Board of Trustees.\nAll rights reserved.\n"
                << "\nWritten by Jeffrey R. Van Voorst and Leslie A. Kuhn\n";
      return DISPLAY_HELP_ONLY;
      break;
    case NO_NORMALIZATION:
      normalize = false;
      break;
    case BASELINE_HPHOB_REP:
      hphob_method = THREE_PROTEIN_ATOMS;
      break;
    case PSEUDO_SURF_HPHOB_REP:
      hphob_method = PSEUDO_SURFACE;
      break;
    case INCLUDE_METALS:
      include_metals = true;
      break;
    case MSMS_SURFACE:
      call_msms = true;
      break;
    case ALLOW_SMALL_SITEMAPS:
      require_min_npts = false;
      break;
    case USE_BOUNDING_BOX_FOR_SITE_VOLUME:
      prune_to_lig = false;
      break;
    default:
      fprintf(stderr, "Error in processing command line arguments\n");
      return FATAL_ERROR;
    }
  }

  // An error occurred during option processing
  if (rc < -1) {
    fprintf(stderr, "%s: %s\n", poptBadOption(optCon, POPT_BADOPTION_NOALIAS),
            poptStrerror(rc));
    return FATAL_ERROR;
  }

  // We must have the protein pdb file
  if(!prot_fname_Cstr){
    err_msg(A_fname, "get_opts", "A protein PDB structure file is required\n");
    poptPrintUsage(optCon, stderr, 0);
    return INVALID_PARAMETER;
  }else{
    prot_fname = prot_fname_Cstr; 
    free(prot_fname_Cstr);
    prot_fname_Cstr = 0;
  }

  // Either an explicit volume or a ligand file is required
  if(!lig_fname_Cstr && !sphere_Cstr){
    err_msg(A_fname, "get_opts", 
            "Either a ligand file or an explicit volume must be provided\n");
    poptPrintUsage(optCon, stderr, 0);
    return INVALID_PARAMETER;
  }if(lig_fname_Cstr){
    lig_fname = lig_fname_Cstr;
    free(lig_fname_Cstr);
    lig_fname_Cstr = 0;
  }if(sphere_Cstr){
    sphere_str = sphere_Cstr;
    free(sphere_Cstr);
    sphere_Cstr = 0;
  }

  // Get the sitemap PDB file and check that we have at most 1 sitemap output 
  // file name 
  get_popt_arg(optCon, &(pts_fname_Cstr));
  if(poptPeekArg(optCon) != NULL){
    err_msg(A_fname, "get_opts", 
            "Please specify at most 1 output sitemap file\n");
    poptPrintUsage(optCon, stderr, 0);
    return INVALID_PARAMETER;
  }if(pts_fname_Cstr){
    pts_fname = pts_fname_Cstr;
    free(pts_fname_Cstr);
    pts_fname_Cstr = 0;
  }

  if(user_provided_surf_Cstr){
    if(call_msms){
      std::string msg = 
"Please provide either an input molecular surface (in MSMS format) for the \n \
site -OR- use the --msms_surf flag to have SimSite3D call MSMS to generate\n \
the site's molecular surface, but not both\n\n";
      err_msg(A_fname, "get_opts", msg);
      return INVALID_PARAMETER; 
    }

    user_provided_surf = user_provided_surf_Cstr;
    free(user_provided_surf_Cstr);
    user_provided_surf_Cstr = 0;
  }

  if(msms_bin_Cstr){
    msms_binary = msms_bin_Cstr;
    free(msms_bin_Cstr);
    msms_bin_Cstr = 0;
  }

  poptFreeContext(optCon);
  return INITIALIZING; 
}

BaseParameters::status_t
GenPointsParameters::verify_params()
{
  if(!normal_file_exists(prot_fname)) return INVALID_PARAMETER;
  if(lig_fname.length() && !normal_file_exists(lig_fname)) 
    return INVALID_PARAMETER;

  if(A_scratch_dir) scratch_dir = A_scratch_dir;
  if(A_proj_output) proj_output = A_proj_output;
  if(scratch_dir == ".")
    std::cerr << "Using the current directory as scratch space\n";

  // If an output name is not given, use the prefix of the protein PDB file
  // and put the file in the working directory
  if(pts_fname.length() == 0){
    std::string tmp = prot_fname;
    tmp = tmp.substr(tmp.rfind('/') + 1);
    size_t len = tmp.length() + 1;
    tmp.at(len - 6) = 's';
    pts_fname = proj_output +  "/" + tmp;
  // If an output name is given, we should write it to the proj_output dir or 
  // current dir? -- at present it is written to the work dir
  }else if(pts_fname[0] != '/'){
    std::string tmp = proj_output + "/" + pts_fname;
    pts_fname = tmp;
  }

  if(hbond_dens_str){ 
    std::string tmp = hbond_dens_str;
    if(tmp == "optimum_only") hbond_method = OPTIMUM_HBONDS;
    else if(tmp == "minimal") hbond_method = MIN_HBONDS;
    else if(tmp == "sparse") hbond_method = SPARSE_HBONDS;
    else if(tmp == "dense") hbond_method = DENSE_HBONDS;
    else{
      std::string msg = tmp;
      msg += " is an unknown hydrogen bond points density";
      err_msg(A_fname, "verify_params()", msg);
      return INVALID_PARAMETER;
    }
  }

  if(A_waters_str) string_tok(A_waters_str, &waters, ' ');

  if(probe_radius < 1.0 || 5.0 < probe_radius){
    err_msg(A_fname, "verify_params()", 
            "probe radius must be between 1.0 and 5.0 (A)");
    return INVALID_PARAMETER;
  }if(num_pts_per_area < 0.1 || 10 < num_pts_per_area){
    err_msg(A_fname, "verify_params()",
            "The number of points per (A)^2 must be between 0.1 and 10");
    return INVALID_PARAMETER;
  }

  if(min_res_per_chain < 1){
    err_msg(A_fname, "verify_params()",
            "Minimum number of residues per protein chain must be positive");
    return INVALID_PARAMETER;
  }

  return READY;
}
