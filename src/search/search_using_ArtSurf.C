/******************************************************************************
 * Copyright (c) 2006-2011, Michigan State University (MSU) Board of Trustees.
 * This file is part of the SimSite3D software project.
 *
 * Authors: Jeffrey Van Voorst, jeff.vanvoorst@gmail.com
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
 * $Date: 2011-02-25 10:12:03 -0500 (Fri, 25 Feb 2011) $: Date of last commit
 * $Author: vanvoor4 $: Author of last commit
 * $Rev: 1620 $: svn revision of last commit
 *
 * svn file: $Id: search_using_ArtSurf.C 1620 2011-02-25 15:12:03Z vanvoor4 $
 * file location: $URL: file:///psa/share/repository/SimSite3D/branches/surfaces-branch/src/search/search_using_ArtSurf.C $
 *****************************************************************************/
#include <Search.H>
#include <NoTier1Score.H>
#include <HbondSurfacesScore.H>
#include <ICPSites.H>

using namespace SimSite3D;

template < class tier1_SF, class tier2_SF, typename align_T >
class ScoreWithWeightedICP : 
  public ScoreRigidAlignments<tier1_SF, tier2_SF, align_T> {
public:

  typedef ScoreRigidAlignments<tier1_SF, tier2_SF, align_T>
    ScoreBase;
  typedef typename ScoreBase::rigid_align_obj            rigid_align_obj;


  ScoreWithWeightedICP(ModelSitemap* model_in, const SearchParameters& params,
                       const my_float_t surf_pt_w, const my_float_t hb_cap_w)
    : ScoreRigidAlignments<tier1_SF, tier2_SF, align_T>(model_in, params)
  {
    A_model = ScoreBase::model();
    A_fine_tune_tier2_alignments = params.fine_tune_tier2_alignments;
    A_mu = ScoreBase::mu();
    A_sigma = ScoreBase::sigma();
    A_save_rigid_results = params.save_rigid_scores;
    A_surf_pt_w = surf_pt_w;
    A_hb_cap_w = hb_cap_w;
  }

  // First a straight copy to compare results with previous method -- should
  // be identical.  Next update to include cap surf points
  void 
  fine_tune_align(DbaseSitemap* search, rigid_align_obj *align)
  {
    std::cout << "Inside test_caps_N_surf::fine_tune_align()" << std::endl;
    align->save_rigid_alignment_vals();
    if(ScoreBase::do_IK())
      ScoreBase::IK_tests_handle()
        .compute_initial_atomic_rmsd(A_model, *search, &(*align));

    tier2_SF &my_SF = ScoreBase::A_tier2_score_class;
    ICP::fine_tune_caps_and_surf(A_model, search, &(*align), &my_SF,
                                 A_surf_pt_w, A_hb_cap_w);

    // Get a single transformation that moves the dbase site to the query
    A_model->get_current_inverse_3D_transform(&(align->Q), align->T);
    align->Q.get_ortho_rot_mat(align->R);
    align->score = (align->score - A_mu) / A_sigma;
  }

  void
  write_score_header(std::ostream& out)
  {
    ScoreBase::write_score_header(out); 
    out << std::setw(50) << "# Mol surf point weight:" << A_surf_pt_w << "\n"
        << std::setw(50) << "# Hbond cap surf point weight:" << A_hb_cap_w 
        << "\n";
  }


private:

  ModelSitemap *A_model;
  bool A_fine_tune_tier2_alignments;
  bool A_save_rigid_results;
  my_float_t A_surf_pt_w;
  my_float_t A_hb_cap_w;
  my_float_t A_mu;
  my_float_t A_sigma;

};


class test_caps_N_surf_ICP : public Search{
public:
  test_caps_N_surf_ICP(const SearchParameters* args_in) : Search(args_in)
  { ; }

  virtual ~test_caps_N_surf_ICP()
  { ; }

  bool
  run()
  {
    if(fail()){
      warn("test_caps_N_surf_ICP", "run", 
           "Cannot run because an error has occured");
      return false;
    }

    if(args()->do_internal_prot_lig_score){
      warn("search_using_ArtSurf.C", "run", "Currently the saving of rigid scores or ArtSurf is mutually exclusive with protein-ligand scoring.");
      return false;
    }

    if(args()->time_process) start_timer();
    align_method_t align_method;
    if(!set_alignment_method(&align_method)) return false;

    my_float_t surf_pt_w = 1.0, hbond_pt_w = 0.0;
    // Given our that our implementation of ICP uses a closed form to solve
    // for the optimum transformation, we require that the max weight for
    // each be <= 1.0 and one should be 1.0
    if(args()->fine_tune_ratio > 0.0){
      if(args()->fine_tune_ratio > 1.0){
        surf_pt_w = 1.0;
        hbond_pt_w = 1.0 / (args()->fine_tune_ratio);
      }else{
        surf_pt_w = args()->fine_tune_ratio;
        hbond_pt_w = 1.0;
      }
    }if(args()->fine_tune_ratio < 0.0){
      surf_pt_w = 0.0;
      hbond_pt_w = 1.0;
    }
    std::cout << "Surf, hbond weights (resp.): " << surf_pt_w << " " 
              << hbond_pt_w << "\n";
  
    std::cout << "Using hbond caps+ scoring method -- training method\n";
    if(args()->do_IK){
      std::vector<align_w_preIK_data_t> alignments;
      score(&alignments, align_method, 
            new ScoreWithWeightedICP< NoTier1Score, HbondSurfacesScore, align_w_preIK_data_t >
            (model(), *(args()), surf_pt_w, hbond_pt_w));
    }else{
      std::vector<align_w_rigid_data_t> alignments;
      score(&alignments, align_method, 
            new ScoreWithWeightedICP< NoTier1Score, HbondSurfacesScore, align_w_rigid_data_t >
            (model(), *(args()), surf_pt_w, hbond_pt_w));

   }


    // no if statement here -- we always run the test case
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
  
  my_params.use_hbond_surfaces = true;
  test_caps_N_surf_ICP my_search(&my_params);
  if(my_search.fail()){
    std::cerr << "\n" << argv[0]
              << " *FAILED* \n\tCould not initialize search\n";
    return 1;
  }

  if(!my_search.run()) return 1;
  return 0; 
}
