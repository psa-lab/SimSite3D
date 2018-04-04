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
 * $Source: /psa/share/repository/pfizer_proj/src/search/SimplifiedMayoHbondsScore.H,v $
 * $Revision: 1.1 $
 * $Author: vanvoor4 $
 * $Date: 2008-05-13 16:48:25 $
 * 
 * $Log: not supported by cvs2svn $
 * 
 */

#ifndef SIMPLIFIED_MAYO_HBONDS_SCORE_HEADER_FILE
#define SIMPLIFIED_MAYO_HBONDS_SCORE_HEADER_FILE

#include <ScoreRigidAlignments.H>

namespace ASCbase{

  //! Base class for scoring rigid alignments of dbase sitemaps to a given model
  class SimplifiedMayoHbondsScore  
   : public ScoreRigidAlignments<std::less<my_float_t> >{

    typedef struct{
      hbond_ideal_pt_vci S;
      hbond_ideal_pt_vci M;
    }S_and_M;

  public:
    //! Default constructor for scoring methods of rigidly aligned sitemaps 
    /*!
     * @param model_in Pointer to the model points (sitemap class)
     * @param params Reference to the search's runtime parameters
     */
    SimplifiedMayoHbondsScore(ModelSitemap* model_in, 
                              const SearchParameters& params);

    //! basic destruction
    ~SimplifiedMayoHbondsScore(); 

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
    my_float_t simp_mayo_fun(hbond_ideal_pt_vci M, hbond_ideal_pt_vci S);

    static const my_float_t V_zero;
    static const my_float_t d_zero;
    static const my_float_t d_zero_squared;
    static const std::string A_fname;          //!< Source file name
  };
}
#endif
