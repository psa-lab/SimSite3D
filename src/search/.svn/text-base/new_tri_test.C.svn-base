#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <ctime>
#include <basics.H>
#include <param_tools.H>
#include <mat_ops.H>
#include <Search.H>
#include <HbondSurfacesScore.H>
#include <NoTier1Score.H>
#include <ModelSiteRMSD.H>
#include <AllPairsSitemapTest.H>


using namespace SimSite3D;

class test_Search : public Search{
public:
  test_Search(const SearchParameters* args_in) : Search(args_in)
  { ; }

  ~test_Search() 
  { ; }

  bool
  run()
  {
    if(fail()){
      warn("new_tri_test", "run", "Cannot run because an error has occured");
      return false;
    }

    if(args()->time_process) start_timer();
    align_method_t align_method;
    if(!set_alignment_method(&align_method)) return false;
    std::vector<align_w_tri_params_t> alignments;


    std::cout << "This is a RESEARCH METHOD: it computes the RMSD for each all reasonable triangle matchs relative to the query and database starting positions -- it only makes sense in when the input sites are already aligned\n\n\n";
    score(&alignments, align_method,
          new ScoreRigidAlignments<NoTier1Score, ModelSiteRMSD,
                                   align_w_tri_params_t>(model(), *args()));

    get_timer_and_write_to_file();
    return true;
  }
private:

};


int main(const int argc, const char **argv)
{
  std::cout << "\n" << argv[0] << " (" << PACKAGE_NAME << ") " 
            << PACKAGE_VERSION << "\n\n";

  // Do not return -1 here since system, fork, etc return -1 on failure, and we
  // wish to distinguish between system and program failure
  SearchParameters my_params(argc, argv);
  BaseParameters::status_t status = my_params.status();
  if(status == BaseParameters::DISPLAY_HELP_ONLY) return 0;
  else if(status != BaseParameters::READY){
    std::cerr << "\n" << argv[0]
              << " *FAILED* \n\tCould not initialize parameters\n";
    return 1;
  }
  
  //my_params.use_hbond_surfaces = true;
  test_Search my_search(&my_params);
  if(my_search.fail()){
    std::cerr << "\n" << argv[0]
              << " *FAILED* \n\tCould not initialize search\n";
    return 1;
  }

  if(!my_search.run()) return 1;
  return 0; 
}

