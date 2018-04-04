/******************************************************************************
 * Copyright (c) 2010, Michigan State University (MSU) Board of Trustees.
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

#ifndef MODEL_SITE_RMSD_HEADER_FILE
#define MODEL_SITE_RMSD_HEADER_FILE

#include <SitesScoreBase.H>

namespace ASCbase{

//! Weight interaction NN interaction pair sums
class ModelSiteRMSD : public SitesScoreBase{
public:
  typedef std::less<my_float_t>		score_cmp; //! A lower score is more favorable

  //! Default constructor for site map point scoring of (aligned) sitemaps 
  ModelSiteRMSD() { ; }

  //! basic destruction
  virtual
  ~ModelSiteRMSD() { ; }

  //! Given an alignment of dbase site to the query, score the two sites 
  //! to the query relative using the site map points
  /*!
   * @param model Const ref to the model site
   * @param dbase Const ref to the sitemap aligned to the model sitemap
   * @param scores iterator to the score data class
   * @return The score of the alignment
   */
  //template< typename align_T_iterator > my_float_t 
  my_float_t
  score(const ModelSitemap &model, const DbaseSitemap &dbase,
        rigid_align_t *scores)
        //align_T_iterator scores)
  {
    std::fill(scores->match_print.begin(), scores->match_print.end(), false);
    return model.compute_RMSD();
  }
                   
  inline bool
  uses_surface_mesh() 
  { return false; }

  inline bool
  uses_hbond_surfaces() 
  { return false; }

private:
  static const std::string A_fname;             //!< Source file name
};
} // end namespace ASCbase
#endif
