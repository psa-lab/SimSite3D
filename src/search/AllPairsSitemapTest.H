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
 * Authors: Jeffrey Van Voorst, jeff.vanvoorst@gmail.com
 *          Leslie Kuhn, Ph.D., KuhnL@msu.edu 
 *****************************************************************************/
#ifndef ALL_PAIRS_SITEMAP_TEST_HEADER_FILE_INCLUDED
#define ALL_PAIRS_SITEMAP_TEST_HEADER_FILE_INCLUDED

#include <ModelSitemap.H>
#include <hash_class.H>

namespace SimSite3D{

//! A test model sitemap used to determine how the "best" we can do
//! given certain tolerances on triangle matches (not triangle sizes)
/*!
 * This class is implemented in a kludge manner so that the interface for
 * get_bucket_iters is the same as that for ModelSitemaps.  We don't really
 * care about speed penalty here as it is not production code
 */
class AllPairsSitemapTest : public ModelSitemap{
public:
  AllPairsSitemapTest(const std::string path, const std::string struct_id,
               const BaseParameters& args, const bool normalize = false,
               const bool load_hbond_surfaces = false,
               const bool hydrophobic_query = false) :
    ModelSitemap(path, struct_id, args, false, false, false)
  {
    A_color_bins.resize(ModelSitemap::NUMBER_OF_TYPE_HASH_CLASSES);
    create_color_bins();
  }

  ~AllPairsSitemapTest() { ; }

  //! Given c,p,l,s get the iterators for the corresponding color bucket
  /*! Given a "color", get all the triangles with that color via iterators to
   * the begin and end of the bucket.
   *
   * @param c hash_class index 
   * @param p perimeter level index  -- not used in this class
   * @param l longest edge length index -- not used in this class
   * @param s shortest edge length index -- not used in this class
   * @param bucket_begin Pointer to the beginning of the bucket
   * @param bucket_end Pointer to the end of the bucket
   * @return true if bucket is nonempty
   */
  bool get_bucket_iters(const uint c, const my_float_t perimeter,
                        const my_float_t max_len, const my_float_t min_len,
                        bucket_vci* bucket_begin, bucket_vci* bucket_end)
  {
    if(c > ModelSitemap::NUMBER_OF_TYPE_HASH_CLASSES) return false;
    *bucket_begin = A_color_bins[c].begin(); 
    *bucket_end = A_color_bins[c].end(); 
    return true;
  }

private:
  void 
  create_color_bins()
  {
    interact_pts_vci A, B, C;
    triangle_t triangle;
    for(A = interact_pts_beg(); A < interact_pts_end(); ++A){
      triangle.vertices[0] = A;
      for(B = A + 1; B < interact_pts_end(); ++B){
        triangle.vertices[1] = B;
        for(C = B + 1; C < interact_pts_end(); ++C){
          triangle.vertices[2] = C;
          triangle.edge_len[0] = dist(A->pos, B->pos);
          triangle.edge_len[1] = dist(B->pos, C->pos);
          triangle.edge_len[2] = dist(C->pos, A->pos);
          
          // Need a minimum value -- cannot have points that are
          // almost colinear or duplicate points
          int long_idx = 0;
          int short_idx = 0;
          my_float_t perimeter = triangle.edge_len[0];
          for(int i = 1; i < 3; ++i){
            long_idx = (triangle.edge_len[i] > triangle.edge_len[long_idx] ? 
                        i : long_idx);
            short_idx = (triangle.edge_len[i] < triangle.edge_len[short_idx] ? 
                         i : short_idx);
            perimeter += triangle.edge_len[i];
          }
          
          if(triangle.edge_len[short_idx] > 0.1 &&
             (perimeter - 2*triangle.edge_len[long_idx]) > 0.1){
            int color_bin = hash_class[A->act_type + B->act_type + C->act_type];
            // Probably used a list in ModelSitemap to avoid invalid iterators
            A_why_use_a_list.push_front(triangle);
            A_color_bins[color_bin].push_back(A_why_use_a_list.begin());
          }
        }
      }
    }
    std::cout << "bin sizes:\n";
    for(int i = 0; i < 10; ++i)
      std::cout << "bin " << i << " size: " << A_color_bins[i].size() << "\n";
    std::cout << std::endl;
  }



//  std::vector< std::vector<triangle_t> > A_color_bins;
  idx_tbl_lvl_long A_color_bins;  //!< Abuse the name ...
  std::list<triangle_t> A_why_use_a_list; //!< A double pass is likely to be more efficient -- first populate the vector, then hold it constant while determining the iterators for each bin
  static const std::string A_fname;   //!< Source file name
};

} // end namespace SimSite3D

#endif
