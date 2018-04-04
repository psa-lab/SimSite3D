#include <sstream>
#include <iomanip>
#include <popt.h>
#include <ImmovableTrimesh.H>
#include <TransformableTrimesh.H>
#include <stream_basics.H>
#include <PDBBase.H>
#include <mol2File.H>

using namespace SimSite3D;
using namespace SimSite3D::geometry;

class comp_surf_parameters{
public:
  std::string query_in_fname; //!< Path to the input query MSMS vert (surf) file
  std::string dbase_in_fname; //!< Path to the input dbase MSMS vert (surf) file
  std::string dbase_out_fname;//!< Path to the output dbase MSMS vert file
  bool fine_tune;             //!< If true, use ICP to fine tune surfaces
  std::string prot_in_fname;  //!< Path to the input protein dbase file
  std::string prot_out_fname; //!< Path to the output protein dbase file
  std::string lig_in_fname;   //!< Path to the input ligand dbase file
  std::string lig_out_fname;  //!< Path to the output ligand dbase file

  typedef enum{
    FINE_TUNE = 1
  }popt_args_t;

  comp_surf_parameters(const int argc, const char **argv)
  {
    init_vars();
    get_opts(argc, argv);
    if(verify() == false) A_fail = true;
  }

  const bool
  fail() const
  { return A_fail; }

private:

  void 
  init_vars()
  { 
    fine_tune = false; 
    
    A_fail = false;
  }

  bool
  get_opts(const int argc, const char **argv)
  {
    char *query_in_Cstr = 0;
    char *dbase_in_Cstr = 0;
    char *dbase_out_Cstr = 0;
    char *prot_in_Cstr = 0;
    char *prot_out_Cstr = 0;
    char *lig_in_Cstr = 0;
    char *lig_out_Cstr = 0;

    struct poptOption OptionsTable[] = {
      { "query", 'q', POPT_ARG_STRING, &query_in_Cstr, 0,
        "The name of the query MSMS surface file", "<QUERY_PATH>.vert"},
      { "dbase_in", 'd', POPT_ARG_STRING, &dbase_in_Cstr, 0,
        "The name of the dbase MSMS surface file", "<DBASE_PATH>.vert"},
      { "fine_tune", '\0', POPT_ARG_NONE, 0, FINE_TUNE,
        "Fine tune the surfaces using ICP"},
      { "dbase_out", 'o', POPT_ARG_STRING, &dbase_out_Cstr, 0,
        "The name of the output dbase MSMS surface file (after fine tuning)"},
      { "prot_in", '\0', POPT_ARG_STRING, &prot_in_Cstr, 0,
        "Path to the dbase protein to transform (using final fine tuning transform)", "<prot>.pdb"},
      { "prot_out", '\0', POPT_ARG_STRING, &prot_out_Cstr, 0,
        "Path to store the transformed dbase protein", "<prot>.pdb"},
      { "ligand_in", '\0', POPT_ARG_STRING, &lig_in_Cstr, 0,
        "Path to the dbase ligand to transform (using final fine tuning transform)", "<lig>.mol2"},
      { "ligand_out", '\0', POPT_ARG_STRING, &lig_out_Cstr, 0,
        "Path to store the transformed dbase ligand", "<lig>.mol2"},
      POPT_AUTOHELP
      POPT_TABLEEND
    };

    std::ostringstream ostr;
    ostr << "\nThe purpose of this program is to compare the query surface to "
         << "the database\n surface.  In particular query points are matched "
         << "to the closest points on the\n database surface.\n";
    poptContext optCon = poptGetContext(argv[0], argc, argv, OptionsTable, 0);
    poptSetOtherOptionHelp(optCon, ostr.str().c_str());
    if(argc < 2) {
      poptPrintUsage(optCon, stderr, 0);
      return false;
    }

    int rc;
    // Process the options 
    for(rc = 0; (rc = poptGetNextOpt(optCon)) >= 0; ){
      switch(rc){
      case FINE_TUNE:
        fine_tune = true;
        break;
      default:
        std::cerr <<  "Error in processing command line arguments\n";
        return false;
        break;
      }
    }

    Cstring_to_string(query_in_Cstr, &query_in_fname);
    Cstring_to_string(dbase_in_Cstr, &dbase_in_fname);
    Cstring_to_string(dbase_out_Cstr, &dbase_in_fname);
    Cstring_to_string(prot_in_Cstr, &prot_in_fname); 
    Cstring_to_string(prot_out_Cstr, &prot_out_fname); 
    Cstring_to_string(lig_in_Cstr, &lig_in_fname); 
    Cstring_to_string(lig_out_Cstr, &lig_out_fname); 

  return true; 
  }

  bool 
  Cstring_to_string(char *Cstr, std::string *str)
  {
    if(Cstr){
      *str = Cstr;
      free(Cstr);
      Cstr = 0;
      return true;
     }
    return false;
  }

  bool 
  verify()
  {
    std::vector<std::string> msms_surf_extensions;
    msms_surf_extensions.push_back(".vert");
    msms_surf_extensions.push_back(".surf");

    if(query_in_fname.length() == 0){
      std::cerr << "Error:  A query surface file is required as input"
                << std::endl;
      return false;
    }else{
      // Check that the file exists
      if(!normal_file_exists(query_in_fname, true)) return false;

      // Check for file extension
      if(!eat_file_extension(&query_in_fname, msms_surf_extensions)){
        std::cerr << "Warning: Expected a query surface file name that ends "
                  << "with \".vert\"\n Instead received: " 
                  << query_in_fname << "\n";
      }
    }

    if(dbase_in_fname.length() == 0){
      std::cerr << "Error:  A dbase surface file is required as input"
                << std::endl;
      return false;
    }else{
      // Check that the file exists
      if(!normal_file_exists(dbase_in_fname, true)) return false;

      // Check for file extension
      if(!eat_file_extension(&dbase_in_fname, msms_surf_extensions)){
        std::cerr << "Warning: Expected a dbase surface file name that ends "
                  << "with \".vert\"\n Instead received: " 
                  << dbase_in_fname << "\n";
      }
    }

    if(prot_out_fname.length() && prot_in_fname.length() == 0){
      std::cerr << "Error: a transformed protein file requires an input "
                << "protein file to transform\n";
      return false;
    }
    if(lig_out_fname.length() && lig_in_fname.length() == 0){
      std::cerr << "Error: a transformed ligand file requires an input "
                << "ligand file to transform\n";
      return false;
    }

    if(!fine_tune){
      if(dbase_out_fname.length() || prot_out_fname.length() ||
         lig_out_fname.length())
        std::cerr << "Warning:  Since fine tuning was not requested, "
                  << "output files will contain the\nsame positions as their "
                  << "corresponding input files\n";
    }
    return true;
  }

  bool
  eat_file_extension(std::string *path, 
                     const std::vector<std::string> &extensions)
  {
    bool found = false;
    std::vector<std::string>::const_iterator ext;
    for(ext = extensions.begin(); ext < extensions.end(); ++ext){
      size_t pos = path->rfind(*ext);
      if(pos != std::string::npos && 
         pos + ext->length() == path->length()){
        std::string tmp = path->substr(0, pos);
        *path = tmp;
        found = true;
        break;
      }
    }
    return found;
  }

  bool A_fail; //!< True implies something went too wrong to continue

};

/*
get_opts(const int argc, const char **argv, std::string *query_str,
         std::string *dbase_str, std::string *dbase_ofname, bool *fine_tune)

bool
get_opts(const int argc, const char **argv, std::string *query_str,
         std::string *dbase_str, std::string *dbase_ofname, bool *fine_tune);
*/

my_float_t
correspondences(const TransformableTrimesh &mesh, const my_float_t *dists,
                const my_float_t *surf_closest_pts, 
                const my_float_t max_surf_pt_dist, my_float_t **vertices_ptr, 
                my_float_t **closest_pts_ptr, size_t *npts);

bool
fine_tune_surfaces(TransformableTrimesh *q_mesh, ImmovableTrimesh *db_mesh,
                   my_float_t *closest_pts, my_float_t *dists,
                   const my_float_t max_surf_dist);


int main(const int argc, const char **argv)
{
  const my_float_t max_surf_dist = 1.5;

  std::cout << "#\n# " << argv[0] << " (" << PACKAGE_NAME << ") "
            << PACKAGE_VERSION << std::endl;

  comp_surf_parameters surf_params(argc, argv);
  if(surf_params.fail()){
    std::cout << "Errors in command-line arguments; cannot continue\n";
    return -1;
  }

  TransformableTrimesh q_mesh(surf_params.query_in_fname);
  ImmovableTrimesh db_mesh(surf_params.dbase_in_fname);
/*
  std::string query_fname, dbase_fname, dbase_ofname;
  bool fine_tune = false;
  if(!get_opts(argc, argv, &query_fname, &dbase_fname, &dbase_ofname, 
               &fine_tune)){
    std::cerr << "Error parsing command line\n";
    exit(-1);
  }

  // Check that the files exist
  if(!normal_file_exists(query_fname, true)) exit(-1);
  if(!normal_file_exists(dbase_fname, true)) exit(-1);

  // Read in the surfaces
  size_t pos = query_fname.rfind(".vert");
  if(pos == std::string::npos) pos = query_fname.rfind(".surf");
  if(pos == std::string::npos){
    std::cerr << "Expected a query surface file name that ends with \".vert\"\n"
              << "Instead received: " << query_fname << "\n";
    exit(-1);
  }

  pos = dbase_fname.rfind(".vert");
  if(pos == std::string::npos) pos = dbase_fname.rfind(".surf");
  if(pos == std::string::npos){
    std::cerr << "Expected a dbase surface file name that ends with \".vert\"\n"
              << "Instead received: " << dbase_fname << "\n";
    exit(-1);
  }
  std::string orig_dbase_fname = dbase_fname.substr(0, pos);
  ImmovableTrimesh db_mesh(orig_dbase_fname);
*/

  if(q_mesh.fail() || db_mesh.fail()){
    std::cerr << "Failed to load surfaces" << std::endl;
    exit(-1);
  }

  my_float_t *closest_pts = new my_float_t[4 * q_mesh.number_of_vertices()];
  my_float_t *dists = closest_pts + 3 * q_mesh.number_of_vertices();

  // Compute 
  if(surf_params.fine_tune)
    fine_tune_surfaces(&q_mesh, &db_mesh, closest_pts, dists, max_surf_dist);
  else
    db_mesh.compare(q_mesh.vertices_begin(), q_mesh.number_of_vertices(),
                    closest_pts, dists, 0, 0, 0, max_surf_dist);

  my_float_t num_surf_pts, num_faces, RMSE, area, RMS_norm_err, max_q_pt_dist;
  q_mesh.compute_features(dists, max_surf_dist, 0, &num_surf_pts, &num_faces, 
                          &RMSE, &RMS_norm_err, &area, &max_q_pt_dist);

  std::cout << "# Maximum distance between close points: " << max_surf_dist 
            << "\n";
  std::cout << "# num close surface points|"
            << "num close faces|"
            << "RMSE of surfaces|"
            << "area of close faces|\n";
  std::cout << std::fixed << std::setprecision(0)
            << num_surf_pts << "|" << num_faces << "|" 
            << std::fixed << std::setprecision(3) << RMSE << "|"
            << std::fixed << std::setprecision(1) << area << "|\n";

  // Get the overall transformation to go from original to current positions
  my_float_t T[3], R[9];
  T[0] = T[1] = T[2] = 0.0;
  R[0] = R[4] = R[8] = 1.0;
  R[1] = R[2] = R[3] = R[5] = R[6] = R[7] = 0.0;
  
  if(surf_params.fine_tune){
    Quaternion Q;
    q_mesh.get_current_inverse_3D_transform(&Q, T);
    Q.get_ortho_rot_mat(R);
    std::cout << "The transformation to move the dbase surface to the query "
              << "surface (y = Rx + T)\n\tis:\n";
    std::cout << std::fixed << std::setprecision(6) 
              << "\tR = [ " << R[0] << ", " << R[3] << ", " << R[6] << ",\n"
              << "\t      " << R[1] << ", " << R[4] << ", " << R[7] << ",\n"
              << "\t      " << R[2] << ", " << R[5] << ", " << R[8] << "]\n";
    std::cout << "\tT = [ " << T[0] << ", " << T[1] << ", " << T[2] << "]\n";

    if(surf_params.dbase_out_fname.length()){
      std::string dbase_ofname = surf_params.dbase_out_fname;
      size_t pos = dbase_ofname.rfind(".vert");
      if(pos == std::string::npos) pos = dbase_ofname.rfind(".surf");
      // Assume user didnt' supply .vert or .surf in the file name
      if(pos == std::string::npos) pos = dbase_ofname.length();

      TransformableTrimesh temp_mesh(surf_params.dbase_in_fname);
      temp_mesh.transform(R, T);
      temp_mesh.write(dbase_ofname.substr(0, pos));      
    }
  }

  // Write out any given output files
  if(surf_params.prot_out_fname.length()){
    PDBBase prot(surf_params.prot_in_fname);
    prot.transform(R, T);
    prot.write(surf_params.prot_out_fname);
  }

  if(surf_params.lig_out_fname.length()){
    mol2File lig(surf_params.lig_in_fname);
    lig.transform(R, T);
    lig.write(surf_params.lig_out_fname);
  }

  delete [] closest_pts;
}

bool
fine_tune_surfaces(TransformableTrimesh *q_mesh, ImmovableTrimesh *db_mesh,
                   my_float_t *closest_pts, my_float_t *dists,
                   const my_float_t max_surf_dist)
{
  db_mesh->compare(q_mesh->vertices_begin(), q_mesh->number_of_vertices(),
                   closest_pts, dists, 0, 0, 0, max_surf_dist);

  // Correspondences needed for ICP
  my_float_t *q_verts = 0, *corr_pts = 0;
  size_t num_corr = 0;
  correspondences(*q_mesh, dists, closest_pts, max_surf_dist, &q_verts,
                  &corr_pts, &num_corr);
  if(num_corr == 0){
    std::cerr << "Warning: No correspondences between the two surfaces\n";
    return false;
  }
  std::cout << "Fine tuning surfaces " << std::flush;

  const my_float_t RMSE_TOL = 0.05;
  const size_t N = 100;

  for(size_t n = 0; n < N; ++n){
    if(n % 10 == 0) std::cout << "." << std::flush;

    // Get transformation 
    Quaternion Q;
    my_float_t T[3];
    lse_3D_fit(q_verts, corr_pts, num_corr, &Q, T);

    // Clean up mem
    if(q_verts) delete [] q_verts;
    if(corr_pts) delete [] corr_pts;
    q_verts = 0;
    corr_pts = 0;

    // Apply transformation
    my_float_t R[9];
    Q.get_ortho_rot_mat(R);
    q_mesh->inverse_transform(R, T);
/*
    std::cout << std::scientific << std::setprecision(6) 
              << "R = [ " << R[0] << ", " << R[3] << ", " << R[6] << ",\n"
              << "      " << R[1] << ", " << R[4] << ", " << R[7] << ",\n"
              << "      "  << R[2] << ", " << R[5] << ", " << R[8] << "]\n";
    std::cout << "T = [ " << T[0] << ", " << T[1] << ", " << T[2] << "]\n";
*/

    // Get current correspondences
    db_mesh->compare(q_mesh->vertices_begin(), q_mesh->number_of_vertices(),
                     closest_pts, dists, 0, 0, 0, max_surf_dist);
    my_float_t RMSE = correspondences(*q_mesh, dists, closest_pts, 
                                      max_surf_dist, &q_verts, &corr_pts, 
                                      &num_corr);
    // RMSE tol has been reached
    if(RMSE <= RMSE_TOL) break;

    // Within numerical accuracy
    my_float_t eps = 1e-8;
    if(0 - eps < T[0] && T[0] < 0 + eps &&
       0 - eps < T[1] && T[1] < 0 + eps &&
       0 - eps < T[2] && T[2] < 0 + eps &&
       1.0 - eps < R[0] && R[0] < 1.0 + eps &&
       1.0 - eps < R[4] && R[4] < 1.0 + eps &&
       1.0 - eps < R[8] && R[8] < 1.0 + eps &&
       0.0 - eps < R[1] && R[1] < 0.0 + eps &&
       0.0 - eps < R[2] && R[2] < 0.0 + eps &&
       0.0 - eps < R[3] && R[3] < 0.0 + eps &&
       0.0 - eps < R[5] && R[5] < 0.0 + eps &&
       0.0 - eps < R[6] && R[6] < 0.0 + eps &&
       0.0 - eps < R[7] && R[7] < 0.0 + eps) break;
  }
  std::cout << std::endl;

  if(q_verts) delete [] q_verts;
  if(corr_pts) delete [] corr_pts;
  q_verts = 0;
  corr_pts = 0;
  return true;
}
#if 0
bool
get_opts(const int argc, const char **argv, std::string *query_str,
         std::string *dbase_str, std::string *dbase_ofname, bool *fine_tune)
{
  *fine_tune = false;
  char *query_Cstr = 0;
  char *dbase_Cstr = 0;
  char *dbase_out_Cstr = 0;
  char *prot_Cstr = 0;
  char *lig_Cstr = 0;

  struct poptOption OptionsTable[] = {
    { "query", 'q', POPT_ARG_STRING, &query_Cstr, 0,
      "The name of the query MSMS surface file", "<QUERY_PATH>.vert"},
    { "dbase", 'd', POPT_ARG_STRING, &dbase_Cstr, 0,
      "The name of the dbase MSMS surface file", "<DBASE_PATH>.vert"},
    { "fine_tune", '\0', POPT_ARG_NONE, 0, FINE_TUNE,
      "Fine tune the surfaces using ICP"},
    { "dbase_out", 'o', POPT_ARG_STRING, &dbase_out_Cstr, 0,
      "The name of the output dbase MSMS surface file (after fine tuning)"},
    { "prot_in", 'p', POPT_ARG_STRING, &prot_Cstr, 0,
      "Path to the dbase protein to transform (using final fine tuning transform)", "<prot>.pdb"},
    { "ligand", 'p', POPT_ARG_STRING, &lig_Cstr, 0,
      "Path to the dbase ligand to transform (using final fine tuning transform)", "<lig>.mol2"},
    POPT_AUTOHELP
    POPT_TABLEEND
  };

  std::ostringstream ostr;
  ostr << "\nThe purpose of this program is to compare the query surface to "
       << "the database\n surface.  In particular query points are matched to "
       << "the closest points on the\n database surface.\n";
  poptContext optCon = poptGetContext(argv[0], argc, argv, OptionsTable, 0);
  poptSetOtherOptionHelp(optCon, ostr.str().c_str());
  if(argc < 2) {
    poptPrintUsage(optCon, stderr, 0);
    return false;
  }

  int rc;
  // Process the options 
  for(rc = 0; (rc = poptGetNextOpt(optCon)) >= 0; ){
    switch(rc){
    case FINE_TUNE:
      *fine_tune = true;
      break;
    default:
      std::cerr <<  "Error in processing command line arguments\n";
      return false;
      break;
    }
  }
 

  if(query_Cstr){
    *query_str = query_Cstr;
    free(query_Cstr);
    query_Cstr = 0;
  }else{
    std::cerr << "Error:  Did not receive a query surface file\n";
    return false;
  }
  if(dbase_Cstr){
    *dbase_str = dbase_Cstr;
    free(dbase_Cstr);
    dbase_Cstr = 0;
  }else{
    std::cerr << "Error:  Did not receive a dbase surface file\n";
    return false;
  }
  if(dbase_out_Cstr){
    *dbase_ofname = dbase_out_Cstr;
    free(dbase_out_Cstr);
    dbase_out_Cstr = 0;
  }
//  char *prot_Cstr = 0;
//  char *lig_Cstr = 0;

  return true; 
}
#endif

// This is a copy of the correspondences from src/search/point_and_surf_score.C
my_float_t
correspondences(const TransformableTrimesh &mesh, const my_float_t *dists,
                const my_float_t *surf_closest_pts, 
                const my_float_t max_surf_pt_dist, my_float_t **vertices_ptr, 
                my_float_t **closest_pts_ptr, size_t *npts)
{
  const int N = mesh.number_of_vertices();
  *vertices_ptr = new my_float_t[3*N];
  *closest_pts_ptr = new my_float_t[3*N];
  *npts = 0;

  my_float_t SE = 0.0;
  my_float_t *vert_out = *vertices_ptr;
  my_float_t *cp_out = *closest_pts_ptr;
  const my_float_t *dist = dists;
  const my_float_t *saved_cp = surf_closest_pts;
  const my_float_t *vert = mesh.vertices_begin();
  const my_float_t *verts_end = mesh.vertices_begin() + 3*N;
  for( ; vert < verts_end; vert += 3, saved_cp += 3, ++dist){
    if(*dist <= max_surf_pt_dist){
      std::copy(vert, vert + 3, vert_out);
      std::copy(saved_cp, saved_cp + 3, cp_out);
      vert_out += 3;
      cp_out += 3;
      SE += (*dist)*(*dist);
      *npts += 1;
    }else{
      SE += max_surf_pt_dist*max_surf_pt_dist;
    }
  }

  return std::sqrt(SE/N);
}






