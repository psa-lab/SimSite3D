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

#ifndef NO_TIER1_SCORE_HEADER_FILE
#define NO_TIER1_SCORE_HEADER_FILE

#include <SitesScoreBase.H>

namespace ASCbase{

//! A do nothing score -- used to check the value of the parameter
/*! Used to fill in the first C++ template parameter if we have only 1 tiered 
 * scoring.  This is bogus, but it is currently the best I could think of 
 * in this short amoutn of time and still work well with the logic of 
 * fine tuning and not copy the vector of alignments
 *  
 */
class NoTier1Score : public SitesScoreBase{
public:
  typedef std::less<my_float_t>              score_cmp; //! Use std::less just because

  //! Default constructor for site map point scoring of (aligned) sitemaps 
  NoTier1Score() { ; }

  //! basic destruction
  ~NoTier1Score() { ; }

  //! Returns my_float_max
  /*!
   * @param model Const ref to the model site
   * @param dbase Const ref to the sitemap aligned to the model sitemap
   * @param scores iterator to the score data class
   * @return my_float_max
   */
  //template< class align_T_iterator > my_float_t 
  my_float_t
  score(const ModelSitemap &model, const DbaseSitemap &dbase,
        rigid_align_t *scores)
        //align_T_iterator scores)
  { return my_float_max; }

  inline bool
  uses_surface_mesh()
  { return false; }

  inline bool
  uses_hbond_surfaces()
  { return false; }

  //! This class's score() function score is a no op!!
  const bool
  score_is_noop() const
  { return true; }
                   
};
} // end namespace ASCbase
#endif
