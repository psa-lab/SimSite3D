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

#include <SearchParameters.H>
#include <param_tools.H>
#include <basics.H>
#include <popt.h>
#include <iomanip>
#include <Timer.H>

using namespace SimSite3D;
const std::string SearchParameters::A_fname = "SearchParameters.C";

SearchParameters::SearchParameters() : BaseParameters() 
{
  init_vars();
}

SearchParameters::SearchParameters(const int argc, const char** argv)
 : BaseParameters() 
{
  init_vars();
  // Check if anything has failed in BaseParameters
  if(A_status == FATAL_ERROR) return;

  A_status = get_opts(argc, argv);
  if(A_status != INITIALIZING) return;
  A_status = verify_parameters();
}

void 
SearchParameters::init_vars()
{
  db_start = 0;
  db_stop = 0;
  normalize = true;
  num_scores_to_keep = 1;
  max_tier1_aligns = 10;
  score_cutoff = -1.5;
  min_num_atoms = 5;
  dmetol = 0.3;
  lsetol = 0.3; 
  ligand_rmsd = false;
  sitemap_rmsd = false;
  write_ligands = true;
  align_to_query = true;
  num_rand_aligns = 0;
  time_process = true;
  allow_hphob_triangles = false;
  add_struct_id_field = false;
  do_internal_prot_lig_score = false;
  use_hbond_surfaces = false;
  A_proj_output = 0;
  A_scratch_dir = 0;
  A_dbase_dir = 0;
  A_ligs_dir = 0;
  load_surf_files = true;
  fine_tune_tier2_alignments = false;
  fine_tune_best_tier2_alignment = true;
  check_all_triangles = false;
  save_rigid_scores = false;
  scale_terms = false;
  do_IK = false;
  max_corr_surf_pt_dist = 1.5;
  fine_tune_ratio = 0.0;
  IK_type = IK_BOTH_REPS;
}

void
SearchParameters::free_cstrings()
{
  // According to the manpages for linux we do not need to check if these
  // are null since it is supposed to be handled by free.  However, if I recall
  // correctly some old unix verions of free do not check if the pointer is 
  // NULL.
  if(A_proj_output) free(A_proj_output);
  if(A_scratch_dir) free(A_scratch_dir);
  if(A_dbase_dir) free(A_dbase_dir);
  if(A_ligs_dir) free(A_ligs_dir);
  init_vars();
}

BaseParameters::status_t
SearchParameters::get_opts(const int argc, const char** argv)
{   
  char *ext_score_method_Cstr = 0;
  char *ofname_Cstr = 0;
  char *score_Cstr = 0;
  char *db_idx_fname_Cstr = 0;
  char *query_prot_Cstr = 0;
  char *alignments_fname_Cstr = 0;

  // Setup the supported popt options table
  struct poptOption supportedOptionsTable[] = {
    { "version", '\0', POPT_ARG_NONE, 0, PRINT_VERSION,
      "Print version of software and exit"},
    { '\0', 'o', POPT_ARG_STRING, &ofname_Cstr, 0,
      "Set name of results (output) file", "my_results.txt"},
    { "prot_lig_score", '\0', POPT_ARG_STRING, &ext_score_method_Cstr, 0,
      "Set external prot-lig scoring function to assess prot-lig interactions",
      "SF name (ID)"},
    { "lig_frag_size", '\0', POPT_ARG_INT, &min_num_atoms, 0,
      "Disregard ligand fragments with a size less than the given size", 
      "integer"},
    { "keep_n_scores", '\0', POPT_ARG_INT, &num_scores_to_keep,
      0, "Max # of aligns with scores better than the threshold to keep per dbase site",
      "integer"},
    { "score_threshold", '\0', POPT_ARG_DOUBLE, &score_cutoff, 0,
      "Set the score threshold", "real number"},
    { "no_normalization", '\0', POPT_ARG_NONE, 0, DO_NOT_NORMALIZE,
      "Do not nornalize the sitemap scores with respect the query's scores against the diverse sitemaps"},
    { "SCORE_ONLY", '\0', POPT_ARG_NONE, 0, SCORE_ONLY,
      "Score sitemaps only--assumes sitemaps are already aligned"}, 
    { "hydrophobic_query", '\0', POPT_ARG_NONE, 0, ALLOW_HPHOB_ONLY_TRIANGLES,
      "Use if pocket is highly hydrophobic", 0},
    { "db_index_file", '\0', POPT_ARG_STRING, &db_idx_fname_Cstr, 0,
      "Set the database index file", 0},
    { "dbstart", '\0', POPT_ARG_INT, &db_start, 0,
      "Line number of starting database site (in db index file)", "NUMBER"},
    { "dbstop", '\0', POPT_ARG_INT, &db_stop, 0,
      "Line number of final database site (in db index file)", "NUMBER"},
    { "allow_small_site_maps", '\0', POPT_ARG_NONE, 0, ALLOW_SMALL_SITEMAPS,
      "Allow site maps to have fewer than the recommended minimum number of points"},
    { "no_fine_tuning", '\0', POPT_ARG_NONE, 0, NO_FINE_TUNING,
      "Do not adjust the candidate alignments"},
    { "omit_surfaces", '\0', POPT_ARG_NONE, 0, OMIT_SURFACES,
      "Score aligned sites using only site map points" },
    { "add_struct_id", '\0', POPT_ARG_NONE, 0, ADD_STRUCT_ID_FIELD,
      "Add struct id as the first field in a results line -- useful for empty pocket searches"}, 
    { "prot_lig_scoring", '\0', POPT_ARG_NONE, 0, DO_PROT_LIG_SCORING,
      "Compute the affiscore and orientscore of query protein and database ligand atoms (internal to SimSite3D)"},
    POPT_TABLEEND
  };

  struct poptOption directoryOptionsTable[] = {
    { "proj_output", '\0', POPT_ARG_STRING, &A_proj_output, 0,
      "Set directory to hold the search results", "path/to/store/results"},
    { "dbase_sites", '\0', POPT_ARG_STRING, &A_dbase_dir, 0,
      "Set directory of templates to search", "path/to/dbase/sites"},
    { "dbase_ligs", '\0', POPT_ARG_STRING, &A_ligs_dir, 0,
      "Set directory of ligands corresponding to the dbase to search",
      "path/to/dbase/ligands"},
    { "scratch_dir", '\0', POPT_ARG_STRING, &A_scratch_dir, 0,
      "Set directory for temporary files", "path/to/scratch/dir"},
    POPT_TABLEEND
  };

  // Setup the usupported popt options table -- mainly for testing different
  // things
  struct poptOption UNsupportedOptionsTable[] = {
    { "IK_surf", '\0', POPT_ARG_NONE, 0, IK_SURFACES,
      "test IK surf implementation", 0},
    { "IK_hb_caps", '\0', POPT_ARG_NONE, 0, IK_HB_CAPS,
      "test IK hbond caps implementation", 0},
    { "IK_both", '\0', POPT_ARG_NONE, 0, IK_BOTH_REPS,
      "test IK surf AND hbond caps implementation", 0},
    { "score_method", '\0', POPT_ARG_STRING, &(score_Cstr), 0,
      "Scoring functions", "DO NOT USE"},
    { "DO_NOT_TIME", 't', POPT_ARG_NONE, 0, DO_NOT_TIME_PROCESS,
      "Do not use process timers to time the process", 0},
    { "ligand_rmsd", '\0', POPT_ARG_NONE, 0, LIGAND_RMSD,
      "RMSD of moved ligand with respect to initial alignment"},
    { "sitemap_rmsd", '\0', POPT_ARG_NONE, 0, SITEMAP_RMSD,
      "RMSD of moved sitemap with respect to initial alignment"},
    { "DO_NOT_WRITE_LIGANDS", '\0', POPT_ARG_NONE, 0, DO_NOT_WRITE_LIGANDS,
      "Do not write out any ligand files (full or fragments)"},
    { "RAND_ALIGNS", '\0', POPT_ARG_INT, &num_rand_aligns, SCORE_ONLY,
      "generate N random alignments -- implies SCORE_ONLY"},
    { "hbond_surfaces", '\0', POPT_ARG_NONE, 0, USE_HBOND_SURFACES_MODEL,
      "Use the hbond surfaces modeling to compare binding sites"},
    { "query_prot", '\0', POPT_ARG_STRING, &(query_prot_Cstr), 0,
      "Query protein file -- all residues near those that move"},
    { "fine_tune_tier2", '\0', POPT_ARG_NONE, 0, FINE_TUNE_TIER2_ALIGNMENTS,
      "Do a rigid body surface refinement of second tier alignments before scoring them" },
//    { "fine_tune_best_tier2", '\0', POPT_ARG_NONE, 0, 
//      FINE_TUNE_BEST_TIER2_ALIGNMENT,
//      "Score tier2 alignments first and then use ICP to fine tune the best "
//      "scoring tier2 alignment" },
    { "max_tier1_aligns", '\0', POPT_ARG_INT, &(max_tier1_aligns), 0,
      "max # of scores that may pass tier1 sieve (continue to tier2)", "10" },
    { "scale_terms", '\0', POPT_ARG_NONE, 0, SCALE_SF_TERMS,
      "Scale SF terms to be between [0.0, 1.0] and use weights for scaled SF"},
    { "check_all_triangles", '\0', POPT_ARG_NONE, 0, CHECK_ALL_TRIANGLES,
      "Test surface method by checking every dbase triangle for closest point to query each query vertex -- extremely slow (100x or more)" },
    { "save_rigid_scores", '\0', POPT_ARG_NONE, 0, SAVE_RIGID_SCORES,
      "Save the scores before and after fine tuning of alignments"},
    { "max_corr_surf_pt_dist", '\0', POPT_ARG_DOUBLE, &(max_corr_surf_pt_dist),
      0, "The max distance for corresponding points between 2 binding site surfaces (should be between 1.5 and 4.0)"},
    { "fine_tune_ratio", '\0', POPT_ARG_DOUBLE, &(fine_tune_ratio), 0,
      "The ratio of weights for shape versus chemistry points [default=0.0 --"
      "implies do not use chemistry points]"},
    { "alignments_fname", '\0', POPT_ARG_STRING, &alignments_fname_Cstr, 0,
      "Path to the alignments file -- basically a SimSite3D results file" 
      "this exists primarily for testing purposes"},
    POPT_TABLEEND
  };

  struct poptOption mainTable[] = {
    {NULL, '\0', POPT_ARG_INCLUDE_TABLE, &supportedOptionsTable, 0, 
     "Supported options", NULL},
    {NULL, '\0', POPT_ARG_INCLUDE_TABLE, &directoryOptionsTable, 0, 
     "Directory options", NULL},
    POPT_AUTOHELP
    {NULL, '\0', POPT_ARG_INCLUDE_TABLE, &UNsupportedOptionsTable, 0, 
     "Unsupported options -- developer use only -- results may be undefined",
     NULL},
    POPT_TABLEEND
  };

  poptContext optCon = poptGetContext(argv[0], argc, argv, mainTable, 0);
  poptSetOtherOptionHelp(optCon, "[OPTIONS]* <query sitemap> [<sitemap 2>]");
  if(argc < 2) {
    poptPrintUsage(optCon, stderr, 0);
    return DISPLAY_HELP_ONLY;
  }

  int rc;
  // Process the options 
  for(rc = 0; (rc = poptGetNextOpt(optCon)) >= 0; ){
    switch(rc){
    case PRINT_VERSION:
      print_version(std::cout, argv[0]);
      return DISPLAY_HELP_ONLY;
      break;
    case PRINT_HELP:
      print_version(std::cout, argv[0]);
      return DISPLAY_HELP_ONLY;
      break;
    case DO_NOT_NORMALIZE:
      normalize = false;
      break;
    case LIGAND_RMSD:
      ligand_rmsd = true;
      break;
    case SITEMAP_RMSD:
      sitemap_rmsd = true;
      break;
    case DO_NOT_WRITE_LIGANDS: 
      write_ligands = false;
      fprintf(stderr, "\n -- WARNING --\nNo ligands will be written during "
              "this search.\n");
      break;
    case SCORE_ONLY:
      align_to_query = false; 
      if(num_rand_aligns)
        fprintf(stderr, "\n -- WARNING --\nRandom alignments will be generated.\n\tThis method requires the sitemaps to be already aligned if it is to make sense.\n");
      else
        fprintf(stderr, "\n -- WARNING --\nNo alignments will be generated.\n\tSitemaps will be scored as given.\n");
      break;
    case DO_NOT_TIME_PROCESS:
      time_process = false;
      break;
    case ALLOW_HPHOB_ONLY_TRIANGLES:
      allow_hphob_triangles = true;
      break;
    case ADD_STRUCT_ID_FIELD:
      add_struct_id_field = true;
      break;
    case DO_PROT_LIG_SCORING:
      do_internal_prot_lig_score = true;
      break;
    case USE_HBOND_SURFACES_MODEL:
      use_hbond_surfaces = true;
      score_str = "HbondSurfacesScore";
      load_surf_files = true;
      break;
    case OMIT_SURFACES:
      use_hbond_surfaces = false;
      load_surf_files = false;
      fine_tune_best_tier2_alignment = false;
      score_str = "WeightedSumsScore";
      break;
    case ALLOW_SMALL_SITEMAPS:
      require_min_npts = false;
      break;
    case FINE_TUNE_TIER2_ALIGNMENTS:
      fine_tune_tier2_alignments = true;
      fine_tune_best_tier2_alignment = false;
      break;
    case NO_FINE_TUNING:
      fine_tune_tier2_alignments = false;
      fine_tune_best_tier2_alignment = false;
      break;
    case SAVE_RIGID_SCORES:
      save_rigid_scores = true;
      break;
    case CHECK_ALL_TRIANGLES:
      check_all_triangles = true;
      break;
    case SCALE_SF_TERMS:
      scale_terms = true;
      break;
    case  IK_SURFACES:
      do_IK = true;
      //use_hbond_surfaces = true;
      IK_type = IK_SURFACES;
      break;
    case  IK_HB_CAPS:
      do_IK = true;
      use_hbond_surfaces = true;
      IK_type = IK_HB_CAPS;
      score_str = "HbondSurfacesScore";
      break;
    case  IK_BOTH_REPS:
      do_IK = true;
      use_hbond_surfaces = true;
      IK_type = IK_BOTH_REPS;
      score_str = "HbondSurfacesScore";
      break;
    default:
      fprintf(stderr, "Error in processing command line arguments\n");
      return FATAL_ERROR;
      break;
    }
  }

  // An error occurred during option processing
  if (rc < -1) {
    fprintf(stderr, "%s: %s\n", poptBadOption(optCon, POPT_BADOPTION_NOALIAS),
            poptStrerror(rc));
    return FATAL_ERROR;
  }

  // Get the query sitemap -- we require exactly one query.
  char *model_fname_Cstr = 0;
  get_popt_arg(optCon, &(model_fname_Cstr));
  // Check if we are given a second sitemap
  char *dbase_fname_Cstr = 0;
  get_popt_arg(optCon, &(dbase_fname_Cstr));
  // Make sure we have exactly 1 or 2 sitemaps (no more, no less)
  if(!model_fname_Cstr || poptPeekArg(optCon) != NULL){
    err_msg(A_fname, "get_opts", "Please specify 1 or 2 sitemap files\n");
    poptPrintUsage(optCon, stderr, 0);
    return INVALID_PARAMETER;
  }

  // Copy the input Cstrings to std::strings and free the memory allocated by
  // popt.  std::strings are preferred since the are easier to convert 
  // using Boost.Python than Cstrings.
  model_file_name = model_fname_Cstr;
  free(model_fname_Cstr);
  if(dbase_fname_Cstr){
    dbase_file_name = dbase_fname_Cstr;
    free(dbase_fname_Cstr);
  }
  // Added for testing of surfaces
  if(score_Cstr){ 
    score_str = score_Cstr;
    free(score_Cstr);
  }

  if(score_str == "ModelSiteRMSD") load_surf_files = false;
  if(score_str == "WeightedSumsScore") load_surf_files = false;
  if(score_str == "point_and_surf_score") load_surf_files = true;
  if(ofname_Cstr){
    ofname = ofname_Cstr;
    free(ofname_Cstr);
  }
  if(db_idx_fname_Cstr){
    db_index_fname = db_idx_fname_Cstr; 
    free(db_idx_fname_Cstr);
  }
  if(query_prot_Cstr){
    query_prot = query_prot_Cstr;
    free(query_prot_Cstr);
  }
  if(alignments_fname_Cstr){
    alignments_fname = alignments_fname_Cstr;
    free(alignments_fname_Cstr);
  }  

  poptFreeContext(optCon);
  return INITIALIZING;
}

BaseParameters::status_t 
SearchParameters::verify_parameters()
{   
  if(A_scratch_dir) scratch_dir = A_scratch_dir;

  std::string path;
  std::string struct_id;
  get_path_and_struct_id(model_file_name, &path, &struct_id);
  if(!A_proj_output) proj_output += "/" + struct_id + "_results";
  else proj_output = A_proj_output; 

  if(A_dbase_dir) dbase_sites = A_dbase_dir;
  if(A_ligs_dir) dbase_ligs = A_ligs_dir;

  if(ofname.length() == 0){
    ofname = proj_output + "/" + struct_id;
    if(dbase_file_name.length()){
      std::string dbfile_path, dbfile_struct_id;
      get_path_and_struct_id(dbase_file_name, &dbfile_path, &dbfile_struct_id);
      ofname += "_to_" + dbfile_struct_id + ".out";
    }else ofname += "_blast.out";
  }

  // Check
  if(dbase_sites.length() == 0){
    std::string msg = "Directory for searchable sitemaps was not given";
    msg += " (i.e. $SIMSITE3D_DBASE_SITES)";
    err_msg(A_fname, "verify_parameters", msg);
    return INVALID_PARAMETER;
  }if(dbase_ligs.length() == 0){
    std::string msg = "Directory for ligands corresponding to the searchable ";
    msg += "sitemaps was not given\n (i.e. $SIMSITE3D_DBASE_LIGS)";
    err_msg(A_fname, "verify_parameters", msg);
    return INVALID_PARAMETER;
  }if(num_scores_to_keep < 1){
    std::string msg = "Number of scores to keep must be positive, ";
    msg += "keeping only the best.";
    warn(A_fname, "verify_parameters", msg); 
    num_scores_to_keep = 1;
  }if(min_num_atoms < 0){
    std::string msg = "The minimum number of heavy, ligand atoms in a fragment";
    msg += " must be nonnegative";
    err_msg(A_fname, "verify_parameters", msg); 
    return INVALID_PARAMETER;
  }if(normalize && diverse_sites.length() == 0){
    std::string msg = "Directory for diverse sitemaps was not given";
    msg += " (i.e. $SIMSITE3D_DIVERSE_SITES)";
    err_msg(A_fname, "verify_parameters", msg);
    return INVALID_PARAMETER;
  }if(normalize && diverse_ligs.length() == 0){
    std::string msg = "Directory for ligands corresponding to the diverse ";
    msg += "sitemaps was not given\n (i.e. $SIMSITE3D_DIVERSE_LIGS)";
    err_msg(A_fname, "verify_parameters", msg);
    return INVALID_PARAMETER;
  }if(db_start > 0 && db_stop > 0 && !db_index_fname.length()){
    err_msg(A_fname, "verify_parameters", std::string("If you wish to use ")
            + "db_start/db_stop, please specify the db_index_file");
    return INVALID_PARAMETER;
  }if(db_index_fname.length() && db_start > db_stop){
    err_msg(A_fname, "verify_parameters",
            "The db_start must be less than or equal to db_stop");
    return INVALID_PARAMETER;
// omitting this check since we now are testing fine tuning using site map
// points
//  }if((fine_tune_tier2_alignments || fine_tune_best_tier2_alignment) &&
//      !load_surf_files){
//    err_msg(A_fname, "verify_parameters",
//            "Fine tuning of alignments requires molecular surface files");
//    return INVALID_PARAMETER;
  }if(fine_tune_tier2_alignments && fine_tune_best_tier2_alignment){
    err_msg(A_fname, "verify_parameters",
            "Please choose at most 1 of fine tuning all tier2 alignemnts or "
            "scoring first and then fine tuning only the best tier2 alignment");
    return INVALID_PARAMETER;
  }
  
  if(max_corr_surf_pt_dist < 1.5 || 4.0 < max_corr_surf_pt_dist){
    err_msg(A_fname, "verify_parameters",
            "Please choose a maximum distance for corresponding points to be "
            "in [1.5, 4.0] (A)");
    return INVALID_PARAMETER;
  }
    
// Allow Sitemap() to handle the checking -- why do it twice?
#if 0
  // We require a model file name
  std::string tmp_name = path + "/" + struct_id + "_s.csv";
  if(!normal_file_exists(tmp_name, false)){
    std::cout << "Unable to find: " << tmp_name << "\nlooking for " 
              << struct_id << " in $SIMSITE3D_DBASE_SITES\n"
              << "\t(" <<  dbase_sites <<  ")\n";
    tmp_name = dbase_sites + "/" + struct_id + "_s.csv";
    if(!normal_file_exists(tmp_name)) return INVALID_PARAMETER;
  }
  model_file_name = tmp_name;
#endif
  return READY;
}   

void
SearchParameters::report_parameters(std::ostream& out) const
{
  char buf[257];
  Timer::get_local_time(buf, 80);
  out << "# Parameters for search_sitemaps (" << PACKAGE_NAME << ") "
      << PACKAGE_VERSION << " run\n";
  out << std::left;
  out << std::setw(50) << "# Local start time:" << buf << "\n";
 
  buf[256] = 0; 
  out << std::setw(50) << "# Working directory (via getcwd()):";
  if(getcwd(buf, 256) == 0){
    warn(A_fname, "report_parameters", 
         "Could not get working directory using getcwd()");
    out << "ERROR\n";
  }else out << buf << "\n";

  out << std::setw(50) << "# Searchable sitemaps directory:"
      << dbase_sites << "\n";
  out << std::setw(50) << "# Corresponding ligands directory:"
      << dbase_ligs  << "\n";
  if(db_index_fname.length()){
    out << std::setw(50) << "# Database index file:"
        << db_index_fname << "\n";
    if(db_start > 0 && db_stop > 0)
      out << std::setw(50) << "#   Database start line number:"
          << db_start << "\n"
          << std::setw(50) << "#   Database stop line number:"
          << db_stop << "\n";
  }
  out << std::setw(50) << "# Query (Model) Sitemap Name:"
      << model_file_name << "\n";
  if(dbase_file_name.length())
    out << std::setw(50) << "# Second Sitemap Name:" << dbase_file_name << "\n";
  out << std::setw(50) << "# Max number of scores to keep for each sitemap:"
      << num_scores_to_keep << "\n";
  out << std::setw(50) << "# Score threshold:"
      << score_cutoff << "\n";
  out << std::setw(50) << "# Minimum number of atoms required in a fragment:"
      << min_num_atoms << "\n";
  out << std::setw(50) << "# Average distance metric error tolerance:"
      << dmetol << "\n";
  out << std::setw(50) << "# Average least squares error tolerance:"
      << lsetol << "\n";
  out << std::setw(50) << "# Highly hydrophobic query pocket:";
  if(allow_hphob_triangles) out << "Yes\n";
  else out << "No\n";
  out << std::setw(50) << "# Scoring function terms were scaled:"
      << (scale_terms ? "Yes" : "No") << "\n";
  out << std::setw(50) << "# Max dist between corresponding surface points:"
      << max_corr_surf_pt_dist << "\n";
}
