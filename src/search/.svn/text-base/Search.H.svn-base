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

#ifndef ASCBASE_MAIN_HEADER_FILE
#define ASCBASE_MAIN_HEADER_FILE

#include <errno.h>
#include <dirent.h>
#include <functional>
#include <ScoreRigidAlignments.H>
#include <ExternalScoringFunction.H>
#include <SearchParameters.H>
#include <MatchTriangles.H>
#include <IdentityAlignment.H>
#include <AlignmentsFromResultsFile.H>
#include <Timer.H>
#include <HbondSurfacesScore.H>
#include <NoTier1Score.H>
#include <ModelSiteRMSD.H>
#include <AllPairsSitemapTest.H>

namespace SimSite3D{
 
#if 0 
  template < typename align_T >
  struct align_impl;
#endif

  //! Main (Driver) class for searching sitemaps
  class Search{
  public:

    typedef enum align_method_enum{
      UNKNOWN_ALIGNMENT_METHOD = 0,
      NO_ALIGNMENT,
      FROM_RESULTS_FILE,
      TRIANGLE_MATCHES,
      RANDOM_SAMPLING,
      GRID_SAMPLING
    }align_method_t;

    //! Only allowed Cstor for search_sitemaps
    Search(const SearchParameters* args_in);

    //! Free up the user_args structure
    virtual ~Search();

    //! Do the search
    /*! Given the template gunk as well as gcc issues, we need a function
     * to call score with the desired scoring method (template parameter).
     *
     * @return True if no errors have occured, else false
     */
    bool 
    run();

    //! If true, running will yield undefined results or core dump
    /*!
     * @return True if an error has occured, otherwise false
     */
    bool 
    fail() const 
    { return A_fail; }

    //! Setup the directories for output
    static bool 
    setup_directories(std::string base_dir, mode_t mode);

    const SearchParameters*
    args() const 
    { return A_args; } 

#if 0
    //! Don't modify this -- instead specialize the class below
    template < typename align_T > bool
    align(ModelSitemap *query_site, Sitemap &dset_site, 
          std::vector<align_T> *alignments)
    {
      align_impl<align_T>::align(query_site, dset_site, alignments);
    }
#endif

  protected:
    ModelSitemap* 
    model()
    { return A_model; }

    void
    start_timer()
    { A_timer.start(); }

    bool
    get_timer_and_write_to_file();

    template < typename align_T> bool
    run(std::vector<align_T> *alignments, const align_method_t align_method)
    {
      // Call the scoring method with the desired scoring class
      if(A_args->score_str.length() == 0 ||
         A_args->score_str == "default"){
        std::cout << "using default scoring method\n"; 
        score(alignments, align_method, 
              new ScoreRigidAlignments< WeightedSumsScore, point_and_surf_score,
                                        align_T >(A_model, *A_args));
      }else if(A_args->score_str == "WeightedSumsScore"){
        std::cout << "using site map points only scoring method\n";
        score(alignments, align_method,
              new ScoreRigidAlignments< NoTier1Score, WeightedSumsScore,
                                        align_T >(A_model, *A_args));
      }else if(A_args->score_str == "point_and_surf_score"){
        std::cout << "using site points and surface scoring method and fine tuning of the best alignment\n";
        score(alignments, align_method,
              new ScoreRigidAlignments< NoTier1Score, point_and_surf_score,
                                        align_T >(A_model, *A_args));
      }else if(A_args->score_str == "HbondSurfacesScore" ||
               A_args->use_hbond_surfaces){
        std::cout << "Using hbond caps+ scoring method -- training method\n";
        score(alignments, align_method,
              new ScoreRigidAlignments< NoTier1Score, HbondSurfacesScore,
                                        align_T >(A_model, *A_args));
      }else if(A_args->score_str == "ModelSiteRMSD"){
        std::cout << "Only computing the RMSD of the alignments\n"
                  << "NOTE: this requires that the sites are already aligned\n\n";
        score(alignments, align_method,
              new ScoreRigidAlignments< NoTier1Score, ModelSiteRMSD,
                                        align_T >(A_model, *A_args));
      }else if(A_args->score_str == "sample_all_triangles"){
        std::cout << "This is a RESEARCH METHOD: it computes the RMSD for each all reasonable triangle matchs relative to the query and database starting positions -- it only makes sense in when the input sites are already aligned\n\n\n";
        score(alignments, align_method,
              new ScoreRigidAlignments< NoTier1Score, ModelSiteRMSD,
                                        align_T >(A_model, *A_args));
      }else if(A_args->score_str == "SiteIntersectionVolume"){
        std::cout << "THIS WAS REMOVED -- volume estimate must be much faster than site point score and spherical volumes are too much of an over estimate\n";
    //    std::cout << "Only computes the percentage of the query site volume estimate in the intersection of query site and dbase site\n\n";
    //    score(new ScoreRigidAlignments<SiteIntersectionVolume,
    //                                   CopyTier1Score< SiteIntersectionVolume > >(A_model, *A_args));
      }else{
        std::ostringstream ostr;
        ostr << "Unknown scoring method -- terminating\n"
             << "  Supported methods are:\n"
             << "    point_and_surf_score (Default)\n"
             << "    WeightedSumsScore (Faster than default, but less accurate)\n";
        err_msg(A_fname, "run()", ostr.str());
        return false;
      }
      return true;
    }

    //! Wrapper method to allow pairwise comparisons and database 
    //! (directory) blasts using the given scoring method.
    template <typename scoring_T, typename align_T> bool
    score(std::vector<align_T> *alignments, const align_method_t align_method, 
          scoring_T* score_method)
    {
      std::ofstream results_out;
      if(!open_ofstream(results_out, A_args->ofname)) return false;
//      if(A_args->save_rigid_scores && 
//         !open_ofstream(rigid_results_out, A_args->ofname + ".rigid")) 
//        return false;
      A_args->report_parameters(results_out);
//      A_args->report_parameters(rigid_results_out);
    
      std::string sf_file = A_args->install_dir;
      sf_file += "/SimSite3DSoftParams/ext_prot_lig_score_fcns.conf";
      if(A_args->ext_score_method.length() == 0){
        if(!score_method->set_external_scoring_method(sf_file, "default")) 
          return false;
      }else if(!score_method->set_external_scoring_method(sf_file, A_args->ext_score_method)) return false;
      score_method->write_score_header(results_out);
//      score_method->write_score_header(rigid_results_out);
   
      // Compare a pair of sitemaps
      if(A_args->dbase_file_name.length()){
        std::string path;
        std::string struct_id;
        get_path_and_struct_id(A_args->dbase_file_name, &path, &struct_id);
        alignments->clear();
        align_and_score_dset_site(alignments, score_method, align_method, path, 
                                  struct_id, results_out);
 
      // Blast model against structures in the database index file
      }else if(A_args->db_index_fname.length()){
        std::ifstream db_idx_file;
        if(!open_ifstream(db_idx_file, A_args->db_index_fname)) return false;
        std::cout << "Searching the database: (" << A_args->dbase_sites 
                  << ")\nUsing the database index file: (" 
                  << A_args->db_index_fname << ")\n";
       
        // Note: if this ever gets used in a nonsecure environment, we will
        // want to limit the number of lines read 
        std::vector<std::string> struct_ids;
        // Use the entire file
        if(A_args->db_start == 0 && A_args->db_stop == 0){
          std::cout << "Using all of the structures in the index file\n";
          for(std::string line; std::getline(db_idx_file, line); )
            struct_ids.push_back(line);
        // Use dbstart && dbstop
        }else{
          const int dbstart = A_args->db_start;
          const int dbstop =  A_args->db_stop;
          std::string line;
          int lineno = 1;
          for( ; std::getline(db_idx_file, line) && lineno <= dbstop; ++lineno)
            if(dbstart <= lineno) struct_ids.push_back(line);
          std::cout << "Using structures " << dbstart << " thru " << dbstop
                    << " in the index file\n";
        }

        std::vector<std::string>::const_iterator st_iter = struct_ids.begin();
        for( ; st_iter < struct_ids.end(); ++st_iter){
          alignments->clear();
          align_and_score_dset_site(alignments, score_method, align_method, 
                                    A_args->dbase_sites, *st_iter,
                                    results_out);
        }

      // Blast model against the database (actually a directory)
      }else{
        // Get a listing of all the files in dbase_dir
        struct dirent **dbfilelist;
        int n = scandir(A_args->dbase_sites.c_str(), &dbfilelist, 0, 0);
        if(n < 0){
          int errsv = errno;
          std::cerr << "Unable to get a listing of the files to search.\n"
                    << "Scandir failed: " << strerror(errsv);
          if(errsv == ENOENT)
            std::cerr << " (" << A_args->dbase_sites.c_str() << ")";
          std::cerr << "\n";
          return false;
        }

        std::cout << "Searching the database: (" << A_args->dbase_sites 
                  << ")\n";
    
        // The tag denoting a sitemap labels file is the suffix "_s.csv".
        while(n--){ 
          // Check if file in dbase_dir is a sitemap and remove the last 6 chars
          // to arrive at the structure id
          std::string struct_id = dbfilelist[n]->d_name;
          if(struct_id.length() < 6 || 
             struct_id.substr(struct_id.length() - 6) != "_s.csv"){
            free(dbfilelist[n]);
            continue;
          }
          struct_id.erase(struct_id.end() - 6, struct_id.end());
    
          alignments->clear();
          align_and_score_dset_site(alignments, score_method, align_method, 
                                    A_args->dbase_sites, struct_id,
                                    results_out);
          free(dbfilelist[n]);
        }
         free(dbfilelist);
      }
      results_out.close();
    
      // Need to check if this is a "normalization" run or "actual" run
      std::vector<std::string> tokens;
      string_tok(A_args->proj_output, &tokens, '/');
      if(tokens.back().substr(0,13) != ".SimSite3D_tmp_")
        std::cout << "You may find the results of your search in\n  "
                  << A_args->proj_output << "\n";
     
      return true;
    }

    //! Set the alignment method enum based on the arguments
    bool
    set_alignment_method(align_method_t *align_method)
    {
      *align_method = UNKNOWN_ALIGNMENT_METHOD;
      if(A_args->alignments_fname.length()){
        *align_method = FROM_RESULTS_FILE;
        std::cout << "NOTE: using given results file ("
                  << A_args->alignments_fname << ") for alignments\n";
      }else if(A_args->align_to_query){
        *align_method = TRIANGLE_MATCHES;
      }else if( A_args->num_rand_aligns){
        *align_method = RANDOM_SAMPLING;
        return false;
      }else{
        *align_method = NO_ALIGNMENT;
      } 
      return true;
    }

  private:
    //! Initialize class variables
    void init_vars();

    //! Given the target sitemap and a scoring method, generate alignments
    //! using the desired generation method and score the transformed target
    //! with respect to the query sitemap.
    /*! 
     * Note that this function no longer clears the alignments container
     */
    template < typename scoring_T, typename align_T > bool
    align_and_score_dset_site(std::vector<align_T> *alignments, 
                              scoring_T *score_method,
                              const align_method_t align_method,
                              const std::string target_sitemap_path,
                              const std::string struct_id, 
                              std::ofstream& results_out)
    {
      // 1 Load target sitemap
      DbaseSitemap dbase_site(target_sitemap_path, struct_id, *A_args,
                              A_args->max_corr_surf_pt_dist,
                              A_args->check_all_triangles, 
                              A_args->use_hbond_surfaces);
      if(dbase_site.fail()){
        std::string msg = "Failed to initialize the sitemap:\n\t";
        err_msg(A_fname, "align_and_score_target", 
                msg + target_sitemap_path + "/" + struct_id);
        return false;
      }

      // 2 Generate alignments to move target to query
      A_model->revert(true, true);
      align(A_model, dbase_site, alignments, align_method);

      // 3 Score all alignments
      score_method->score_alignments(*alignments, &dbase_site, results_out);
      return true;
    }

    template < typename align_T > bool
    align(ModelSitemap *query_site, Sitemap &dset_site, 
          std::vector<align_T> *alignments, const align_method_t align_method)
    {
      if(align_method == NO_ALIGNMENT){
        IdentityAlignment my_align;
        my_align.align(query_site, dset_site, alignments);
      }else if(align_method == TRIANGLE_MATCHES){
        MatchTriangles my_align(A_args->dmetol, A_args->lsetol, 
                                A_args->allow_hphob_triangles);
        my_align.align(query_site, dset_site, alignments);
      }else if(align_method == FROM_RESULTS_FILE){
        AlignmentsFromResultsFile my_align(A_args->alignments_fname);
        my_align.align(query_site, dset_site, alignments);
      }else if(align_method == RANDOM_SAMPLING){
        std::cout << "Unimplemented sampling method -- cannot continue\n";
        return false;
      }else if(align_method == GRID_SAMPLING){
        std::cout << "Unimplemented sampling method -- cannot continue\n";
        return false;
      }else{
        std::cout << "Unknown sampling method -- cannot continue\n";
        return false;
      }
      if(alignments->size() > 0) return true;
      return false;
    }


    bool A_fail;              //!< If true, running will yield undefined results
    const SearchParameters* A_args;  //! Cmdline and config file args
    ModelSitemap* A_model;   //!< Model to which we wish to align dbase sitemaps
    Timer A_timer;

    static const std::string A_fname;  //!< Source file name
  };

#if 0
  //! specialize this class template to get your alignment method
  template < typename align_T >
  struct align_impl
  {
    // Overload this function
    static bool align(Sitemap& dset_site, std::vector<align_T> *alignments);
  };
#endif

}
#endif
