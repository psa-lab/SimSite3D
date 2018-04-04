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

#ifndef _MATCH_TRIANGLES_HEADER_FILE
#define _MATCH_TRIANGLES_HEADER_FILE

#include <list>
#include <ModelSitemap.H>
#include <ScoreMapBase.H>
#include <hash_class.H>

namespace SimSite3D{

//! Registration of sitemaps based on close 3 point alignments
//class MatchTriangles : public InitialAlignmentsBase{
class MatchTriangles{
public:
  //! Cstor 
  /*!
   * @param dmetol_in Max tolerated average distance metric error
   * @param lsetol_in Max tolerated average least squares error
   * @param allow_hphob_triangles If false, triangles must contain at least 1
   *        hydrogen bond
   */
  MatchTriangles(const my_float_t dmetol_in, const my_float_t lsetol_in, 
                 const bool allow_hphob_triangles);

  //! basic destruction
  virtual
  ~MatchTriangles(); 
#if 0
  template < typename align_T >
  struct test_align : public align_impl< align_T>{

    bool
    align(Sitemap& dset_site, std::vector<align_T> *alignments)
    {
      if(model->fail() || dset_site.fail() || !alignments) return false;
    
      errorType error = NONE;
      targ_pts_beg = dset_site.interact_pts_beg();
      std::vector<uint> tmp_array(dset_site.interact_pts_end() - targ_pts_beg);
      for(uint i = 0; i < tmp_array.size(); ++i) tmp_array[i] = i;
      for(uint i = 2; !error && i < tmp_array.size(); i++){
        tmp_array[2] = tmp_array[i];
        compute_triangles(&tmp_array, 2, 3, alignments);
      }
  
      return true;
    }
    
  };
#endif

  //! Run to generate the 3pt rigid alignments.  
  /*!
   * @param dset_site Pointer to the "second" sitemap
   * @param alignments Pointer to the vector holding found correspondences
   */
  template< typename align_T > bool
  align(ModelSitemap *query_site, Sitemap& dset_site, 
        std::vector<align_T> *alignments)
  {
    if(query_site->fail() || dset_site.fail() || !alignments) return false;
    A_model = query_site;
  
    errorType error = NONE;
    A_targ_pts_beg = dset_site.interact_pts_beg();
    std::vector<uint> tmp_array(dset_site.interact_pts_end() - A_targ_pts_beg);
    for(uint i = 0; i < tmp_array.size(); ++i) tmp_array[i] = i;
    for(uint i = 2; !error && i < tmp_array.size(); i++){
      tmp_array[2] = tmp_array[i];
      compute_triangles(&tmp_array, 2, 3, alignments);
    }

    return true;
  }
  

private:
  ModelSitemap* A_model;
  my_float_t target_centroid[3];
  interact_pts_vci A_targ_pts_beg;
  my_float_t dmetol_2;  //!< Square of avg distance metric error tolerance
  my_float_t lsetol_2;  //!< Square of avg least squares pt correspondence tol
  //std::vector<rigid_align_t> *A_alignments; //!< candidate alignments
  bool A_allow_hphob_triangles; 

  // Constants for ease of use/readability
  static const std::string _fname;           //!< Name of the source file
  static const int permutations[6][3];       //!< Permutations of edges
  static const int point_permutations[6][3]; //!< Permutations of vertices
  static const int ACCEPTOR_ACCEPTOR;   //!< 2 * ACCEPTOR
  static const int DONOR_DONOR;         //!< 2 * DONOR
  static const int HYDROPHOB_HYDROPHOB; //!< 2 * HYDROPHOB
  static const int ACCEPTOR_DONOR;      //!< ACCEPTOR + DONOR
  static const int ACCEPTOR_DONEPTOR;   //!< ACCEPTOR + DONEPTOR
  static const int DONOR_DONEPTOR;      //!< DONOR + DONEPTOR
  static const int DONEPTOR_DONEPTOR;   //!< DONEPTOR + DONEPTOR


  template< typename align_T > void
  compute_triangles(std::vector<uint>* arr_p, int index, int length,
                    std::vector<align_T> *alignments)
  {
    std::vector<uint>& arr = *arr_p;
    triangle_t target_tri;

    // Get edge lengths and perimeter for the given vertices
    for(uint i = 0; i < 3; ++i) 
      target_tri.vertices[i] = A_targ_pts_beg + arr[i];
    my_float_t perimeter = 0;
    for(uint i = 0; i < 3; ++i){
      target_tri.edge_len[i] = 
        dist(target_tri.vertices[i]->pos, target_tri.vertices[(i+1) % 3]->pos);
      perimeter += target_tri.edge_len[i];
    }
    
    // Find min and max edge lengths -- assumes no eqilateral triangle
    my_float_t min_len = target_tri.edge_len[0];
    my_float_t max_len = target_tri.edge_len[0];
    for(uint i = 1; i < 3; ++i){
      if(target_tri.edge_len[i] > max_len) max_len = target_tri.edge_len[i]; 
      else if(target_tri.edge_len[i] < min_len) 
        min_len = target_tri.edge_len[i];
    }
  
    // Get the interaction types of the vertices
    bool have_Hbond = false;
    for(int i = 0; i < 3; ++i)
      if(target_tri.vertices[i]->act_type == ACCEPTOR 
         || target_tri.vertices[i]->act_type == DONOR) 
        have_Hbond = true;
    
    // Get iterators to the bucket holding the similar model triangles
    int c = 0;
    for(uint i = 0; i < 3; ++i) c += target_tri.vertices[i]->act_type;
    c = hash_class[c];
    if(c >= 10) std::cerr << "need to handle doneptors\n";
    bucket_vci bucket_begin, bucket_end;
    bool have_model_triangles = 
      A_model->get_bucket_iters(c, perimeter, max_len, min_len, 
                                &bucket_begin, &bucket_end);
  
    // Need a minimum value -- cannot have points that are
    // almost colinear or duplicate points
    if(min_len > 0.1 && (perimeter - 2*max_len) > 0.1 && 
       (have_Hbond || A_allow_hphob_triangles) && have_model_triangles){
      for(bucket_vci model_tri = bucket_begin; model_tri < bucket_end;
          ++model_tri){
        int idx = get_best_correspondence(**model_tri, target_tri);
        if(idx == -1) continue;
  
        // Get the permuted search triangle corresponding to idx and find a 
        // rigid rotation and translation (by a least squares fit) to move the 
        // permuted search triangle into alignment with the current model 
        // triangle.
        triangle_t c_tri;
        for(uint i = 0; i < 3; ++i)
          c_tri.vertices[i] = target_tri.vertices[point_permutations[idx][i]];
        align_T my_align;
        my_align.frag_file = NULL;
        if(ThreePtLeastSquares(c_tri, **model_tri, &(my_align.Q), my_align.R,
                               my_align.T)){
          // Potentially inefficient depending on vector design
          my_align.set_triangle_params(perimeter, max_len, min_len);
          alignments->push_back(my_align);
  
          //Kludge for testing sampling since some sites have an exorbant
          // number of alignments
          if(alignments->size() > 4000000){
            std::cout << "number of alignments: " << alignments->size()
                      << std::endl;
           return;
          }
        }
      }
    }

    if(index > 0){
      int orig = arr[index-1]; /* store original entry at position index-1 */
      while(arr[index-1] < arr[index] - 1) {
        /* increase value at position index - 1 as long as it is
           still smaller than the value at position index */
        arr[index-1]++;
        /* recursively increase all elements in the array at
           positions [0..index-2] */
        compute_triangles(arr_p, index - 1, length, alignments);
      }
      /* restore original value at position index-1 */
      arr[index-1] = orig;
    }
  }

  //! Find the best correspondence between a given model and search triangle
  /*! For each (6) permuation of the vertices of the search triangle see if 
   * the corresponding vertices (of the model triangle) are compatible.  If
   * they are, determine the average difference of the corresponding 
   * edge lengths.  The permutation with the lowest average difference is 
   * taken to be the one of choice.  If there is no permutation with 
   * corresponding vertices with an average difference of edge lengths under
   * the given tolerance, return -1 (no correspondence).
   *
   * @param model_tri Reference to a model triangle
   * @param search_tri Reference to a search triangle
   * @return 0--5 depending which permutation is closest match, -1 if none
   */
  int get_best_correspondence(const triangle_t& model_tri, 
                              const triangle_t& search_tri);

  //! Wrapper to get the least squares fit for two triangles
  /*! Get a weighted least squares fit to align A to B.  That is, A is moved
   * to the reference frame of B.
   *
   * @param A Triangle to move to the reference point of B
   * @param B Triangle to align A to
   * @param res Pointer to a results structure (holds the transformation)
   * @return True implies least squares had error below tolerance
   */
  bool 
  ThreePtLeastSquares(const triangle_t& A, const triangle_t& B, Quaternion *Q,
                      my_float_t *R, my_float_t *T);


  bool
  get_triangles(triangle_t *target_tri, const std::vector<uint>& arr, 
                bucket_vci *bucket_begin, bucket_vci *bucket_end);

};
} // end namespace SimSite3D
#endif
