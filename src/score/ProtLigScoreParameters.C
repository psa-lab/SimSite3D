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
 * Authors: Jeffrey Van Voorst, jeff.vanvoorst@gmail.com
 *          Leslie Kuhn, Ph.D., KuhnL@msu.edu 
 *****************************************************************************/
  
#include <popt.h>
#include <ProtLigScoreParameters.H>

using namespace SimSite3D;

const std::string ProtLigScoreParameters::A_fname = "ProtLigScoreParameters.C";

ProtLigScoreParameters::ProtLigScoreParameters()
{
  init();
  // Check if anything has failed in BaseParameters
  if(A_status == FATAL_ERROR) return;
}

ProtLigScoreParameters::ProtLigScoreParameters(const int argc, 
                                               const char** argv)
{
  init();
  // Check if anything has failed in BaseParameters
  if(A_status == FATAL_ERROR) return;

  A_status = get_opts(argc, argv);
  if(A_status != INITIALIZING) return;
  A_status = verify_parameters();
}

ProtLigScoreParameters::~ProtLigScoreParameters()
{

}

void
ProtLigScoreParameters::init()
{
  build_interact_tbl = false;
  print_interactions = false;
}

BaseParameters::status_t 
ProtLigScoreParameters::get_opts(const int argc, const char** argv)
{
  char *prot_fname_Cstr = 0;
  char *lig_fname_Cstr = 0;
  char *lig_list_fname_Cstr = 0;

  struct poptOption singleOptionsTable[] = {
    { "version", 'v', POPT_ARG_NONE, 0, PRINT_VERSION,
      "Print version and exit", 0},
//    { "side_chains", '\0', POPT_ARG_STRING, &opts->sc_fname, 0,
//      "Path to PDB file holding moved side chains", 0},
    { "ligand", 'l', POPT_ARG_STRING, &lig_fname_Cstr, 0,
      "Path to ligand mol2 file (in docked conformation)", 0},
//    { "waters", 'w', POPT_ARG_STRING, &opts->waters_fname, 0,
//      "Path to conserved waters pdb file", 0},
    { "interactions", '\0', POPT_ARG_NONE, 0, PRINT_INTERACTIONS,
      "Print a list of the protein interactions satisifed by the ligand", 0},
    POPT_TABLEEND
  };

  struct poptOption multiOptionsTable[] = {
//    { "sc_list", '\0', POPT_ARG_STRING, &opts->sc_list_fname, 0,
//      "File listing the PDB files for moved side chains", 0},
    { "lig_list", '\0', POPT_ARG_STRING, &lig_list_fname_Cstr, 0,
      "File listing the ligand mol2 files (in docked conformation)", 0},
//    { "water_list", '\0', POPT_ARG_STRING, &opts->water_list_fname, 0,
//      "File listing the conserved waters pdb files", 0},
    { "interactions_table", '\0', POPT_ARG_NONE, 0, BUILD_INTERACTIONS_TABLE,
      "Build a table of protein-ligand interactions", 0},
    POPT_TABLEEND
  };

  struct poptOption mainOptionsTable[] = {
    { "protein", 'p', POPT_ARG_STRING, &prot_fname_Cstr, 0,
      "Path to protein PDB file or protein RAD file", 0},
    {NULL, '\0', POPT_ARG_INCLUDE_TABLE, &singleOptionsTable, 0,
      "Options for single target and ligand", 0},
    {NULL, '\0', POPT_ARG_INCLUDE_TABLE, &multiOptionsTable, 0,
      "Options for a single target and multiple ligand and side chain files",
      0},
    POPT_AUTOHELP
    POPT_TABLEEND
  };
  poptContext optCon = poptGetContext(argv[0], argc, argv, mainOptionsTable, 0);

  if(argc < 2) {
    poptPrintUsage(optCon, stderr, 0);
    return DISPLAY_HELP_ONLY;
  }

  int rc;
  /* Process the options */
  for(rc = 0; (rc = poptGetNextOpt(optCon)) >= 0; ){
    switch(rc){
    case PRINT_VERSION:
      print_version(std::cout, argv[0]);
      return DISPLAY_HELP_ONLY;
      break;
    case BUILD_INTERACTIONS_TABLE:
      build_interact_tbl = true;
      break;
    case PRINT_INTERACTIONS:
      print_interactions = true;
      break;
    default:
      fprintf(stderr, "Error in processing command line arguments\n");
      return FATAL_ERROR;
    }
  }

  /* An error occurred during option processing */
  if (rc < -1) {
    fprintf(stderr, "%s: %s\n", poptBadOption(optCon, POPT_BADOPTION_NOALIAS),
            poptStrerror(rc));
    return FATAL_ERROR;
  }

  if(!prot_fname_Cstr){
    std::cerr << "A protein filename is required" << std::endl;
    return FATAL_ERROR;
  }else{
    prot_fname = prot_fname_Cstr;
    free(prot_fname_Cstr);
    prot_fname_Cstr = 0;
  }

  if(lig_fname_Cstr && lig_list_fname_Cstr){
    std::cerr << "Please choose either a single docked ligand file or a text "
              << "file listing\ndocked ligand files and not both" << std::endl;
    return FATAL_ERROR;
  }else if(lig_fname_Cstr){
    lig_fname = lig_fname_Cstr;
    free(lig_fname_Cstr);
    lig_fname_Cstr = 0;
  }else if(lig_list_fname_Cstr){
    lig_list_fname = lig_list_fname_Cstr;
    free(lig_list_fname_Cstr);
    lig_list_fname_Cstr = 0;
  }else{
    std::cerr << "Please choose either a single docked ligand file or a text "
              << "file listing\ndocked ligand files and not both" << std::endl;
    return FATAL_ERROR;
  }

  poptFreeContext(optCon);
  return INITIALIZING;
}
  
BaseParameters::status_t 
ProtLigScoreParameters::verify_parameters()
{
  return READY;
}    
