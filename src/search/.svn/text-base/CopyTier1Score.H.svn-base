/******************************************************************************
 * Copyright (c) 2010, Michigan State University (MSU) Board of Trustees.
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

#ifndef COPY_TIER1_SCORE_HEADER_FILE
#define COPY_TIER1_SCORE_HEADER_FILE

#include <SitesScoreBase.H>

namespace SimSite3D{

//! Returns the tier1 score as the score for the current tier
/*! Used to fill in C++ template parameters if we have only 1 tiered scoring
 *  
 *  NOTE: the template parameter MUST be the scoring function class used
 *  for tier1
 */
template< class tier1_SF >
class CopyTier1Score : public SitesScoreBase{
public:
    typedef typename tier1_SF::score_cmp	score_cmp; //! Use same compare as tier1 SF

  //! Default constructor for site map point scoring of (aligned) sitemaps 
  CopyTier1Score() { ; }

  //! basic destruction
  ~CopyTier1Score() { ; }

  //! Returns score from tier1
  /*!
   * @param model Const ref to the model site
   * @param dbase Const ref to the sitemap aligned to the model sitemap
   * @param scores iterator to the score data class
   * @return score from tier1 (scores->tier1_score)
   */
  my_float_t 
  score(const ModelSitemap &model, const DbaseSitemap &dbase,
        rigid_align_vi scores)
  { return scores->tier1_score; }

  inline bool
  uses_surface_mesh()
  { return false; }

  inline bool
  uses_hbond_surfaces()
  { return false; }
                   
};
} // end namespace SimSite3D
#endif
