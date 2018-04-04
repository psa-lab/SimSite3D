/******************************************************************************
 * Copyright (c) 2006,2007, Michigan State University (MSU) Board of Trustees.
 *   All rights reserved.
 *
 * This file is part of the ASCbase Software project.
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
  
/*
 * $Source: /psa/share/repository/pfizer_proj/src/search/ASCbaseSearchParameters.H,v $
 * $Revision: 1.5 $
 * $Author: vanvoor4 $
 * $Date: 2009-01-21 01:52:37 $
 * 
 * $Log: not supported by cvs2svn $
 * Revision 1.4  2009/01/12 21:29:23  vanvoor4
 * Changed to make it more easier to build python extensions as well
 * as changed numeric values from popt to enums.
 *
 * Revision 1.3  2008/07/29 16:29:17  vanvoor4
 * Updated the command line status flag to support more values than 2.
 *
 * Revision 1.2  2008/05/13 16:32:10  vanvoor4
 * Added command line flags for score only (sitemaps aligned by an
 * outside method) and to generate random alignments
 *
 * Revision 1.1  2007/11/01 16:31:23  vanvoor4
 * Split from ASCbaseMain.H
 * 
 *
 */

#ifndef ASCBASE_SEARCH_PARAMETERS_HEADER_INCLUDED
#define ASCBASE_SEARCH_PARAMETERS_HEADER_INCLUDED

#include <BaseParameters.H>
#include <types.H>

namespace ASCbase{

  //! Use popt to parse a search's commandline parameters 
  class SearchParameters : public BaseParameters{
  public:
    typedef enum{
      PRINT_VERSION = 1,
      PRINT_HELP,
      DO_NOT_NORMALIZE,
      LIGAND_RMSD,
      SITEMAP_RMSD,
      DO_NOT_WRITE_LIGANDS,
      SCORE_ONLY,
      DO_NOT_TIME_PROCESS,
      ALLOW_HPHOB_ONLY_TRIANGLES,
      ADD_STRUCT_ID_FIELD,
      DO_PROT_LIG_SCORING,
      USE_HBOND_SURFACES_MODEL,
      ALLOW_SMALL_SITEMAPS,
      OMIT_SURFACES,
      USE_ONLY_SURFACE_TO_RANK,
      NO_FINE_TUNING,
      FINE_TUNE_TIER2_ALIGNMENTS,
      FINE_TUNE_BEST_TIER2_ALIGNMENT,
      SAVE_RIGID_SCORES,
      CHECK_ALL_TRIANGLES,
      SCALE_SF_TERMS,
    }popt_args_t;

    typedef enum{
      IK_SURFACES = 1000,
      IK_HB_CAPS,
      IK_BOTH_REPS
    }IK_type_t;

    //! To be used by the ASCbasePy search.parameters module
    SearchParameters();

    //! Set the search parameters using the values in argv
    SearchParameters(const int argc, const char** argv);

    ~SearchParameters()
    { free_cstrings(); }

    void report_parameters(std::ostream& out) const;

    std::string ext_score_method; //!< Name of external scoring method to use
    std::string model_file_name;  //!< Name of the model points file
    std::string dbase_file_name;  //!< Name of the current dbase points file
    std::string ofname;           //!< Name of the output file
    std::string score_str;        //!< Scoring method 
    std::string db_index_fname;   //!< Name of database index file
    int db_start;                 //!< Line number of initial db site
    int db_stop;                  //!< Line number of final db site
    bool normalize;         //!< false implies no normalization
    int num_scores_to_keep; //!< max # of final alignments to keep for each search sitemap
    int max_tier1_aligns;   //!< max # of tier1 alignments to pass to tier2 scoring
    my_float_t score_cutoff; //!< Keep alignments with score better than cutoff
    int min_num_atoms;      //!< Number of ligand atoms required in query pocket
    my_float_t dmetol;      //!< Tolerance for the average dist. metric error
    my_float_t lsetol;      //!< Tolerance for the average least squares error
    bool ligand_rmsd;
    bool sitemap_rmsd;
    bool write_ligands;     //!< True implies write ligands, else do not write
    bool align_to_query;    //!< True implies align dbase site to query, else score as given
    size_t num_rand_aligns; //!< If not align_to_query, use this many random alignments
    bool time_process;      //!< If true, use itimers to time program
    bool allow_hphob_triangles; //!< If false, triangle need at least 1 hbond
    bool add_struct_id_field;   //!< If false & empty pockets searches, write struct id into ligand field for empty pocket hits
    bool do_internal_prot_lig_score;  //!< If true compute affi and orient score for query prot and db lig interactions
    bool use_hbond_surfaces;  //!< If true, use the hbond surfaces modeling of hydrogen bonding regions (volumes)
    bool fine_tune_tier2_alignments; //!< If true, use an ICP like method to fine tune all tier2 alignments
    bool fine_tune_best_tier2_alignment;  //!< If true, use an ICP-like method to fine tune the best scoring tier2 alignment
    bool check_all_triangles;  //!< If true, check all triangles in dbase mesh for closest point to each query vertex -- useful for testing
    std::string query_prot;    //!< Currently required for hbond_surfaces
    bool save_rigid_scores;    //!< When using ICP, save the rigid scores to another output file
    bool scale_terms;         //!< If true, scale each term of the scoring function(s) linearly from [0.0, max] to [0.0, 1.0] where max is the max value the query site can have for that term
    bool do_IK;
    my_float_t max_corr_surf_pt_dist; //!< Maximum distance to consider corresponding surface points
    my_float_t fine_tune_ratio;  //!< Ratio of shape point weight (msms surf) to hbond cap point weight (1.0 imples equal, 2.0 implies surface points get weighted twice that of hbond caps, etc.)  Default is Zero, which implies use only the MSMS surface for ICP
    std::string alignments_fname; //!< Path to alignments file (in SimSite3D results file format with at least 4 fields per line -- struct_id, score, R, T; and in that order and pipe delimited fields with numbers single space delimited)
    IK_type_t IK_type; //!< Type of IK method to use
  private:

    void init_vars();

    void free_cstrings();

    status_t get_opts(const int argc, const char** argv);
    
    status_t verify_parameters();

    char* A_conf_fname;     //!< Location of environment configure file
    char* A_proj_output;    //!< Directory to store results
    char* A_scratch_dir;    //!< Directory to use as a scratch space
    char* A_dbase_dir;      //!< Location of the database to query
    char* A_ligs_dir;       //!< Location of ligands associated with dbase_dir
    static const std::string A_fname;  //!< Name of the source file
  };
}
#endif
