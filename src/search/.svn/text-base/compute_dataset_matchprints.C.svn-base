/******************************************************************************
 * Copyright (c) 2011, Jeffrey Van Voorst
 * This file is part of the SimSite3D software project.
 *
 * Authors: Jeffrey Van Voorst, vanvoor4@msu.edu
 *          Leslie Kuhn, Ph.D., KuhnL@msu.edu 
 *
 *  SimSite3D is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  SimSite3D is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc.,
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *  You may also visit http://www.gnu.org/licenses/gpl-2.0.html.
 * 
 * $Date: $: Date of last commit
 * $Author: $: Author of last commit
 * $Rev: $: svn revision of last commit
 *
 * svn file: $Id: $
 * file location: $URL: $
 *****************************************************************************/
#include <RigidAlignWithDsetMP.H>
#include <WSS_and_dataset_MP.H>
#include <NoTier1Score.H>
#include <Search.H>

using namespace SimSite3D;

class compute_dataset_matchprints : public Search{
public:

  //! I am doing this quickly, so most parameters will be ignored
  compute_dataset_matchprints(const SearchParameters* args_in): Search(args_in)
  { ; }

  virtual ~compute_dataset_matchprints()
  { ; }

  bool
  run()
  {
    if(fail()){      
      warn("compute_dataset_matchprints", "run",
           "Cannot run because an error has occured");
      return false;
    }

    if(args()->time_process) start_timer();
    align_method_t align_method = UNKNOWN_ALIGNMENT_METHOD;
    if(args()->alignments_fname.length()){
      align_method = FROM_RESULTS_FILE;
      std::cout << "NOTE: using given results file ("
                << args()->alignments_fname << ") for alignments\n";
    }else{
      warn("compute_dataset_matchprints", "run",
           "Cannot run because this method only supports alignments from a results file\nPlease use the --alignments_fname flag with the results file for which you seek the dataset matchprints.\n");
      return false;
    }

    std::vector<RigidAlignWithDsetMP> alignments;
    Search::score(&alignments, align_method,
                  new ScoreRigidAlignments< NoTier1Score, WSS_and_dataset_MP,
                                            RigidAlignWithDsetMP >(model(), *(args())));
    if(!Search::run(&alignments, align_method)) return false;

    get_timer_and_write_to_file();
    return true;
  } 

};

int main(const int argc, const char** argv)
{
  std::cout << "\n" << argv[0] << " (" << PACKAGE_NAME << ") "
            << PACKAGE_VERSION << "\n"
            << "Note: this program was written quickly.  The command line flags for search options will not be used (this includes protein-ligand scoring, fine_tuning etc.).  The score and terms correspond to WeightedSumsScore (aka SF8)\n\n";

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

  compute_dataset_matchprints my_search(&my_params);
  if(my_search.fail()){
    std::cerr << "\n" << argv[0]
              << " *FAILED* \n\tCould not initialize search\n";
    return 1;
  }

  if(!my_search.run()) return 1;
  return 0;

}
