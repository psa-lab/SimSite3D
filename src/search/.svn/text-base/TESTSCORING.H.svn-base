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

/*
 * $Source: /psa/share/repository/pfizer_proj/src/search/TESTSCORING.H,v $
 * $Revision: 1.3 $
 * $Author: vanvoor4 $
 * $Date: 2008-01-04 15:41:39 $
 * 
 * $Log: not supported by cvs2svn $
 * Revision 1.2  2007/11/01 16:35:11  vanvoor4
 * Now pass in the parameters obj instead of those needed
 *
 * Revision 1.1  2007/10/11 16:21:08  vanvoor4
 * Initial checkin
 *
 * 
 */

#ifndef TEST_SCORING_SUMS_SCORE_HEADER_FILE
#define TEST_SCORING_SUMS_SCORE_HEADER_FILE

#include <ScoreRigidAlignments.H>

namespace SimSite3D{

  //! Base class for scoring rigid alignments of dbase sitemaps to a given model
  class TestScoring : public ScoreRigidAlignments<std::less<my_float_t> >{
  public:
    //! Default constructor for scoring methods of rigidly aligned sitemaps 
    /*!
     * @param model_in Pointer to the model points (sitemap class)
     * @param params Reference to the search's runtime parameters
     */
    TestScoring(ModelSitemap* model_in, 
                const SearchParameters& params);

    //! basic destruction
    ~TestScoring(); 

  protected:
    //! Given an alignment of search to the query, score said alignment
    /*!
     * @param search Pointer to the sitemap aligned to the model sitemap
     * @return The score of the alignment
     */
    my_float_t score(const DbaseSitemap& search, rigid_align_vi scores);

    my_float_t
    correspondences(my_float_t **query_pts_ptr, my_float_t **db_pts_ptr,
                    size_t *npts)
    {
      *query_pts_ptr = 0;
      *db_pts_ptr = 0;
      *npts = 0;
      return my_float_max;
    }

  private:
    static const std::string _fname;          //!< "WeightedSumsScore.C"
    static const my_float_t POLAR_SUM_W;      //!< -10370
    static const my_float_t HPHOBIC_COUNT_W;  //!< -4130
    static const my_float_t POLAR_MISMATCH_W; //!< 11428
    static const my_float_t UNSAT_POLAR_W;    //!<   963
    static const my_float_t CONSTANT_TERM;    //!< -58070
    static const my_float_t max_dist;         //!< 1.5 angstroms
  };
}
#endif
