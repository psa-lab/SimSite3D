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
 * $Source: /psa/share/repository/pfizer_proj/src/search/EllipsoidalHbondsScore.H,v $
 * $Revision: 1.4 $
 * $Author: vanvoor4 $
 * $Date: 2008-05-13 16:01:40 $
 * 
 * $Log: not supported by cvs2svn $
 * Revision 1.3  2008/01/04 15:50:09  vanvoor4
 * Added template gunk
 *
 * Revision 1.2  2007/11/01 16:36:24  vanvoor4
 * Now pass in the parameters obj instead of those needed
 *
 * Revision 1.1  2007/10/11 16:21:28  vanvoor4
 * Initial checkin
 *
 * 
 */

#ifndef ELLISPSOIDAL_HBONDS_SCORE_HEADER_FILE
#define ELLISPSOIDAL_HBONDS_SCORE_HEADER_FILE

#include <ScoreRigidAlignments.H>

namespace ASCbase{

  //! Base class for scoring rigid alignments of dbase sitemaps to a given model
  class EllipsoidalHbondsScore  
   : public ScoreRigidAlignments<std::greater<my_float_t> >{

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
    EllipsoidalHbondsScore(ModelSitemap* model_in, 
                           const SearchParameters& params);

    //! basic destruction
    ~EllipsoidalHbondsScore(); 

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
    my_float_t ellipsoidal_hbonds(hbond_ideal_pt_vci M, hbond_ideal_pt_vci S);

    static const std::string _fname;          //!< "EllipsoidalHbondsScore.C"
  };
}
#endif
