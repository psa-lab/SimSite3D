/******************************************************************************
 * Copyright (c) 2011, Jeff Van Voorst
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
#ifndef WEIGHTED_SUMS_SCORE_AND_DATASET_MP_HEADER_FILE
#define WEIGHTED_SUMS_SCORE_AND_DATASET_MP_HEADER_FILE

#include <WeightedSumsScore.H>
#include <RigidAlignWithDsetMP.H>
#include <cassert>

namespace ASCbase{

//! Weight interaction NN interaction pair sums
class WSS_and_dataset_MP : public WeightedSumsScore{
public:
  typedef std::less<my_float_t>		score_cmp; //! A lower score is more favorable

  //! Default constructor for site map point scoring of (aligned) sitemaps 
  WSS_and_dataset_MP() 
  { 
  }

  //! basic destruction
  virtual ~WSS_and_dataset_MP() 
  { 
  }

  //! Given an alignment of dbase site to the query, score the two sites 
  //! to the query relative using the site map points
  /*!
   * @param model Const ref to the model site
   * @param dbase Const ref to the sitemap aligned to the model sitemap
   * @param scores iterator to the score data class
   * @param save_correspondences should the score method save the points
   *        correspondences?  Saving is required for methods like ICP 
   * @return The score of the alignment
   */
  my_float_t
  score(const ModelSitemap &model, const DbaseSitemap &dbase,
        RigidAlignWithDsetMP *scores)
        
//        rigid_align_t *scores)
  {
    my_float_t site_score = WeightedSumsScore::score(model, dbase, scores);

    scores->dataset_mp.resize(dbase.fit_points_size());
    std::fill(scores->dataset_mp.begin(), scores->dataset_mp.end(), false);
  
    //std::fill(A_closest_pts, A_closest_pts + 4*A_max_num_pts, my_float_max);
    // get A_closest_pts from base class

    std::vector<bool>::iterator D_MP_start = scores->dataset_mp.begin();

    const HbondPoints& S_hbonds = dbase.hbond_points();
    const my_float_t* S_hbonds_coords_start = S_hbonds.fit_pts_beg()->pos;
    const HbondPoints& M_hbonds = model.hbond_points();
    std::vector<bool>::const_iterator M_hit = scores->match_print.begin();
    const my_float_t* cp = WeightedSumsScore::closest_points();
    for(hbond_fit_pt_vci M = M_hbonds.fit_pts_beg(); M < M_hbonds.fit_pts_end();
        ++M, ++M_hit, cp += 3){
      if(*M_hit && *cp < my_float_max){
        int idx = (cp - S_hbonds_coords_start ) / 3;    
        assert(idx >= 0);
        *(D_MP_start + idx) = true;
      }
    }
  
    const HphobPoints& M_hphobs = model.hphob_points();
    const HphobPoints& S_hphobs = dbase.hphob_points();
    const my_float_t* S_hphob_coords_start = S_hphobs.begin()->pos;
    D_MP_start += 3 * S_hbonds.fit_pts_size(); 
    for(hphob_point_lci H = M_hphobs.begin(); H != M_hphobs.end();
        ++H, ++M_hit, cp += 3){
      if(*M_hit && *cp < my_float_max){
        int idx = (cp - S_hphob_coords_start ) / 3;    
        assert(idx >= 0);
        *(D_MP_start + idx) = true;
      }
    }

    return site_score;
  }

private:

};
} // end namespace ASCbase
#endif
