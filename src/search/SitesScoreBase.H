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
 * Authors: Jeffrey Van Voorst, jeff.vanvoorst@gmail.com
 *          Leslie Kuhn, Ph.D., KuhnL@msu.edu 
 *****************************************************************************/

#ifndef SITES_SCORE_BASE_HEADER_FILE
#define SITES_SCORE_BASE_HEADER_FILE

#include <ModelSitemap.H>
#include <DbaseSitemap.H>
#include <ScoreMapBase.H>

namespace SimSite3D{

//! Base class for scoring two (aligned) sites 
/*! This class defines the interface required for a class to be able to be
 * used to score two sites.  A derived class can be used at either level of
 * the scoring (template parameter 1 or 2 of ScoreRigidAlignments)
 *
 * notes on template parameters:
 * align_T: must be rigid_align_t or derived class
 */
class SitesScoreBase{
public:

  typedef std::less<my_float_t>		score_cmp; //! A lower score is more favorable -- change this in derived class to std::greater if a higher score is more favorable

  typedef std::vector<ModelHbondSurfaces::m_surf_vci>::const_iterator
    m_cap_iters_vci;


  //! Default constructor for site map point scoring of (aligned) sitemaps 
  SitesScoreBase() 
  { A_scale_terms = false; }

  //! basic destruction
  virtual 
  ~SitesScoreBase() { ; }

  //! Given an alignment of dbase site to the query, score the two sites 
  //! to the query relative
  /*!
   * @param model Const ref to the model site
   * @param dbase Const ref to the sitemap aligned to the model sitemap
   * @param scores iterator to the score data class
   * @param save_correspondences should the score method save the points
   *        correspondences?  Saving is required for methods like ICP 
   * @return The score of the alignment
   */
  //template< class align_T_iterator > my_float_t 
  virtual my_float_t
  score(const ModelSitemap &model, const DbaseSitemap &dbase, 
        rigid_align_t *scores) = 0;

  //! Get the latest correspondences of the db to query sites
  virtual my_float_t
  correspondences(const ModelSitemap &model, my_float_t **query_pts_ptr, 
                  my_float_t **db_pts_ptr, size_t *npts)
  {
    std::cerr << "Correspondences are NULL\n";
    *query_pts_ptr = 0;
    *db_pts_ptr = 0;
    *npts = 0;
    return my_float_max;
  }

  //! Used for ICP where we want the correspondences in an array
  virtual my_float_t
  polar_correspondences(const ModelSitemap &model, my_float_t **vertices_ptr,
                        my_float_t **closest_pts_ptr, size_t *npts) const
  {
    std::cerr << "ICP polar correspondences are NULL\n";
    *vertices_ptr = 0;
    *closest_pts_ptr = 0;
    *npts = 0;
    return my_float_max;
  }

  //! Used for IK where not all of the caps may be mobile
  virtual void 
  polar_correspondences(const ModelSitemap &model,
                        m_cap_iters_vci cap_iters_beg,
                        m_cap_iters_vci cap_iters_end,
                        std::vector<const my_float_t *> *q_pts,
                        std::vector<const my_float_t *> *corr_pts) const
  {
    std::cerr << "IK polar correspondences are NULL\n";
  }

  virtual void
  correspondences(const ModelSitemap &model, 
                  std::vector<const my_float_t *> *q_pts,
                  std::vector<const my_float_t *> *corr_pts)
  {
    std::cerr << "this function is supposed to be overriden by a derived class"
              << std::endl;
  }

  //! Derived class should return false, unless, the site's surface needs to be
  //! to be rotated and translated
  inline virtual bool
  uses_surface_mesh()
  { return false; }
  
  //! Derived class should return false, unless the derived class is using the
  //! hydrogen bond surface caps
  inline virtual bool
  uses_hbond_surfaces()
  { return false; }

  //! Set the class variable scale_terms to the given value
  void
  set_scale_terms(bool scale_in)
  { A_scale_terms = scale_in; }

  //! Should each term be linearly scaled to from [0.0, max] to [0.0, 1.0] 
  //! where max is the maximum value the query site can have for that term?
  inline const bool
  scale_terms() const
  { return A_scale_terms; }

  //! Is the function score a no op?
  /*! It is safe to return false here since, this class is a pure virtual class
   */
  virtual const bool
  score_is_noop() const
  { return false; }

private:
  bool A_scale_terms;
};
} // end namespace SimSite3D
#endif
