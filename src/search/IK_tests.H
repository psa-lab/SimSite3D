
#ifndef INVERSE_KINEMATICS_TEST_HEADER_FILE_INCLUDED
#define INVERSE_KINEMATICS_TEST_HEADER_FILE_INCLUDED

#include <HbondSurfacesScore.H>
#include <ScoreMapBase.H>


namespace SimSite3D{

//! Testing/Validation class for surface deformations and hbonds moving
class IK_tests{
public:
  //! Only allowed Cstor for search_sitemaps
  IK_tests(const SearchParameters &args);

  //! Free up the user_args structure
  ~IK_tests();

  //! Big messy test function 
  /*! Assumes that args() contains both query and db sites
   *  (i.e. no dataset searching)
   */
  void
  test_moving_stuff(const SearchParameters &args_in);

  static const my_float_t
  max_atom_dist_for_checking_vdw_overlap()
  { return MAX_ATOM_DIST_TO_CHECK_VDW_OVERLAP; }

  //! Run flexibility using the current orientation of the model site
  /*! Assumption: ICP steps have been done to get binding site into a 
   *  good correspondence
   */
  void
  run(ModelSitemap *model_site, const DbaseSitemap &dbase_site, 
      rigid_align_t *align, const std::string db_struct_id, 
      std::ostream& out);

  my_float_t
  compute_initial_atomic_rmsd(ModelSitemap *model_site, 
                              const DbaseSitemap &dbase_site, 
                              rigid_align_t *align);

private:

/*
  my_float_t
  ICP_mobile_pts(ModelSitemap *model_site, const DbaseSitemap &dbase_site, 
                 rigid_align_t* align, const size_t number_of_iterations,
                 const my_float_t min_delta, const bool verbose);
*/

  SitesScoreBase *A_site_score_p;
  bool A_uses_site_surface;        //! ArtSurf includes site surface
  bool A_uses_hbond_surfaces;      //! ArtSurt includes hbond caps 

  bool A_hbond_caps_ICP;
  my_float_t A_ICP_surf_pt_W;
  my_float_t A_ICP_hb_cap_pt_W;

  static const my_float_t MAX_ATOM_DIST_TO_CHECK_VDW_OVERLAP; //!< Maximum distance between any 2 to require checking of intersection of Van der Waals spheres.
  static const std::string A_fname;

};
} // end namespace SimSite3D
#endif
