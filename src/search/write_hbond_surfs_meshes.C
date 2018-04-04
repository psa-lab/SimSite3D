#include <Search.H>

using namespace SimSite3D;


int main(const int argc, const char **argv)
{
  std::string A_fname = "write_hbond_surfs_meshes.C";

  std::cout << "\n" << argv[0] << " (" << PACKAGE_NAME << ") " 
            << PACKAGE_VERSION << "\n\n";

  // Do not return -1 here since system, fork, etc return -1 on failure, and we
  // wish to distinguish between system and program failure
  SearchParameters args(argc, argv);
  BaseParameters::status_t status = args.status();
  if(status == BaseParameters::DISPLAY_HELP_ONLY) return 0;
  else if(status != BaseParameters::READY){
    std::cerr << "\n" << argv[0]
              << " *FAILED* \n\tCould not initialize parameters\n";
    return 1;
  }

  std::cout << "This program writes out the caps to a modified MSMS surf\n";

  std::string path;
  std::string struct_id;
  get_path_and_struct_id(args.model_file_name, &path, &struct_id);

  // for testing using articulated side chains
  PDBStructure *full_prot = 0;
  if(args.query_prot.size() == 0)
  {
    err_msg(A_fname, "Cstr()",
            "Currently hbond surfaces requires a query prot file as input");
    return -1;
  }else{
    full_prot = new PDBStructure(args.query_prot);
    if(full_prot->fail()){
      err_msg(A_fname, "Cstr()",
              "Currently hbond surfaces requires a query prot file as input");
      return -1;
    }
  }

  ModelSitemap my_site(path, struct_id, args, false,
                       args.use_hbond_surfaces,
                       args.allow_hphob_triangles);
  if(my_site.fail()){
    err_msg(A_fname, "Cstr()", "Failed to intialize the model sitemap");
    return -1;
  }

  my_site.write_msms_caps(args.ofname); 


}
