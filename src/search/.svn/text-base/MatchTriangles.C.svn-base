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

//#define TRACE
#include <cmath>
#include <iomanip>
#include <iostream>
#include <basics.H>
#include <MatchTriangles.H>
//#include <InitialAlignmentsBase.H>

using namespace SimSite3D;

const std::string MatchTriangles::_fname = "MatchTriangles.C";
const int MatchTriangles::permutations[6][3] = 
{ { 0, 1, 2 },
  { 0, 2, 1 },
  { 1, 0, 2 },
  { 1, 2, 0 },
  { 2, 0, 1 },
  { 2, 1, 0 } };
const int MatchTriangles::point_permutations[6][3] = 
{ { 0, 1, 2 },
  { 1, 0, 2 },
  { 2, 1, 0 },
  { 1, 2, 0 },
  { 2, 0, 1 },
  { 0, 2, 1 } };

const int MatchTriangles::ACCEPTOR_ACCEPTOR   =   2;
const int MatchTriangles::DONOR_DONOR         =   8;
const int MatchTriangles::HYDROPHOB_HYDROPHOB =  32;
const int MatchTriangles::ACCEPTOR_DONOR      =   5;
const int MatchTriangles::ACCEPTOR_DONEPTOR   =  65;
const int MatchTriangles::DONOR_DONEPTOR      =  68;
const int MatchTriangles::DONEPTOR_DONEPTOR   = 128;


MatchTriangles::MatchTriangles(const my_float_t dmetol_in, 
                               const my_float_t lsetol_in, 
                               const bool allow_hphob_triangles)
{
  dmetol_2 = dmetol_in * dmetol_in;
  lsetol_2 = lsetol_in * lsetol_in;
  A_allow_hphob_triangles = allow_hphob_triangles;
}

MatchTriangles::~MatchTriangles()
{
  dmetol_2 = -1;
  lsetol_2 = -1;
}


#if 0
bool
MatchTriangles::align(Sitemap& target, align_vec *alignments)
{
  if(model->fail() || target.fail() || !alignments) return false;

  A_alignments = alignments;
  errorType error = NONE;
  targ_pts_beg = target.interact_pts_beg();
  std::vector<uint> tmp_array(target.interact_pts_end() - targ_pts_beg);
  for(uint i = 0; i < tmp_array.size(); ++i) tmp_array[i] = i;
  for(uint i = 2; !error && i < tmp_array.size(); i++){
    tmp_array[2] = tmp_array[i];
    compute_triangles(&tmp_array, 2, 3);
  }

  return true;
}
#endif

#if 0
void
MatchTriangles::compute_triangles(std::vector<uint>* arr_p, int index, 
                                  int length)
{
  std::vector<uint>& arr = *arr_p;
#if 0
  
  // Get edge lengths and perimeter for the given vertices
  triangle_t target_tri;
  for(uint i = 0; i < 3; ++i) target_tri.vertices[i] = targ_pts_beg + arr[i];
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
    else if(target_tri.edge_len[i] < min_len) min_len = target_tri.edge_len[i];
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
    model->get_bucket_iters(c, perimeter, max_len, min_len, 
                            &bucket_begin, &bucket_end);
#endif

  bucket_vci bucket_begin, bucket_end;
  // Find the model triangles which are close enough
  if(get_triangles(arr, &bucket_begin, &bucket_end)){
    for(bucket_vci model_tri = bucket_begin; model_tri < bucket_end; 
        ++model_tri){
      int idx = get_best_correspondence(**model_tri, target_tri);
      if(idx == -1) continue;

      // Get the permuted search triangle corresponding to idx and find a rigid
      // rotation and translation (by a least squares fit) to move the permuted
      // search triangle into alignment with the current model triangle.
      triangle_t c_tri;
      for(uint i = 0; i < 3; ++i)
        c_tri.vertices[i] = target_tri.vertices[point_permutations[idx][i]];
      rigid_align_t my_align;
      my_align.frag_file = NULL;
      if(ThreePtLeastSquares(c_tri, **model_tri, &my_align)){
#if 0
        // Added for testing of triangle sizes
        for(size_t iii = 0; iii < 3; ++iii) my_align.pts_idx[iii] = arr[iii];
        my_align.tri_params[0] = perimeter;
        my_align.tri_params[1] = max_len;
        my_align.tri_params[2] = min_len;
#endif
        // Potentially inefficient depending on vector design
        A_alignments->push_back(my_align);

        //Kludge for testing sampling since some sites have an exorbant
        // number of alignments
        
        if(A_alignments->size() > 2000000){
          std::cout << "number of alignments: " << A_alignments->size() 
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
      compute_triangles(arr_p, index - 1, length);
    }
    /* restore original value at position index-1 */
    arr[index-1] = orig;
  }
}
#endif

int 
MatchTriangles::get_best_correspondence(const triangle_t& model_tri, 
                                        const triangle_t& search_tri)
{
  int min = -1;
  my_float_t min_err_2 = my_float_max;
  // Keep the points of the model triangle fixed and find the best (if any)
  // correspondence of the 6 permutations of the (3) search triangle points 
  // to the (3) model triangle points.
  for(uint i = 0; i < 6; ++i){

    bool match = true;
    for(uint j = 0; j < 3 && match; ++j){
      int sum = model_tri.vertices[j]->act_type + 
        search_tri.vertices[point_permutations[i][j]]->act_type;
      if(!(sum == ACCEPTOR_ACCEPTOR || sum == DONOR_DONOR
         || sum == HYDROPHOB_HYDROPHOB || sum == DONOR_DONEPTOR 
	 || sum == ACCEPTOR_DONEPTOR)) match = false;
    }
  
    // Here we check only if the distances are relatively compatible.
    // If they are, we will do a weighted least squares fit to check if
    // it is possible to align the triangles well.
    if(match){
      std::vector<my_float_t> w(3, 1.0);
      // edge[0] is the edge between vertices 0 and 1, etc.  Thus, if point 0
      // is a hydrophobic point the matching of the edges[0] and [2] are given
      // more tolerance, etc.
      if(model_tri.vertices[0]->act_type == HYDROPHOB){
        w[2] = HPHOB_MATCH_WEIGHT; 
        w[0] = HPHOB_MATCH_WEIGHT; 
      }if(model_tri.vertices[1]->act_type == HYDROPHOB){
        w[0] = HPHOB_MATCH_WEIGHT; 
        w[1] = HPHOB_MATCH_WEIGHT; 
      }if(model_tri.vertices[2]->act_type == HYDROPHOB){
        w[1] = HPHOB_MATCH_WEIGHT; 
        w[2] = HPHOB_MATCH_WEIGHT; 
      }

      my_float_t err_2 = 0.0;
      for(uint j = 0; j < 3; ++j){
        my_float_t tmp = w[j] * (model_tri.edge_len[j] 
	                         - search_tri.edge_len[permutations[i][j]]);
        err_2 += tmp*tmp;
      }
      err_2 /= 3.0;
      if(err_2 < min_err_2){
        min_err_2 = err_2;
	min = i;
      }
    }
  }

  if(min_err_2 > dmetol_2) min = -1;
  return min;
}

bool
MatchTriangles::ThreePtLeastSquares(const triangle_t& A, const triangle_t& B, 
                                    Quaternion *Q, my_float_t *R, my_float_t *T)
{
  // Use a weighted superposition, where an hydrophobic point has a weight that
  // is HPHOB_MATCH_WEIGHT and hydrophillic points have a weight of 1.0.
  my_float_t w[3];
  for(uint i = 0; i < 3; ++i)
    if(A.vertices[i]->act_type == HYDROPHOB) w[i] = HPHOB_MATCH_WEIGHT;
    else w[i] = 1.0;

  my_float_t a_verts[12], b_verts[9];
  for(uint i = 0; i < 3; ++i){
    std::copy(A.vertices[i]->pos, A.vertices[i]->pos + 3, a_verts + 3*i);
    std::copy(B.vertices[i]->pos, B.vertices[i]->pos + 3, b_verts + 3*i);
  }
  std::copy(target_centroid, target_centroid + 3, a_verts + 9);

  my_float_t new_a_verts[12];
  three_pt_align(a_verts, b_verts, w, Q, T, new_a_verts);
  Q->get_ortho_rot_mat(R);

  // Calculate the average difference for positions of the corresponding 
  // vertices.
  my_float_t err_2 = 0;
  for(uint i = 0; i < 3; ++i){
    my_float_t tmp = w[i] * dist_squared(b_verts + 3*i, new_a_verts + 3*i);
    // NOTE: dissertation chap 4 & 5 results were computed without this error 
    // check
    err_2 += tmp*tmp;
  }
  err_2 /= (w[0] + w[1] + w[2]);
#if 0
  my_float_t err_2 = 0;
  for(uint i = 0; i < 3; ++i){
    my_float_t tmp = w[i] * dist(b_verts + 3*i, new_a_verts + 3*i);
    err_2 += tmp*tmp;
  }
  err_2 /= 3;
#endif
  if(err_2 > lsetol_2) return false;

  // Determine if the corresponding vectors (directions) are compatible
  my_float_t a_dirs[9], newa_dirs[9], no_T[3] = {0, 0, 0};
  for(uint i = 0; i < 3; ++i) 
    std::copy(A.vertices[i]->dir, A.vertices[i]->dir + 3, a_dirs + 3*i);
  move_positions(3, newa_dirs, a_dirs, R, no_T);

  // less than is more correct than less than or equal to since the current 
  // hphob method has zero vectors.
  for(size_t i = 0; i < 3; ++i)
    if(dot(newa_dirs + 3*i, B.vertices[i]->dir) < 0) return false;
 
  return true;
}

bool
MatchTriangles::get_triangles(triangle_t *target_tri, 
                              const std::vector<uint>& arr, 
                              bucket_vci *bucket_begin, bucket_vci *bucket_end)
{
  // Get edge lengths and perimeter for the given vertices
  for(uint i = 0; i < 3; ++i) target_tri->vertices[i] = A_targ_pts_beg + arr[i];
  my_float_t perimeter = 0;
  for(uint i = 0; i < 3; ++i){
    target_tri->edge_len[i] = 
      dist(target_tri->vertices[i]->pos, target_tri->vertices[(i+1) % 3]->pos);
    perimeter += target_tri->edge_len[i];
  }
  
  // Find min and max edge lengths -- assumes no eqilateral triangle
  my_float_t min_len = target_tri->edge_len[0];
  my_float_t max_len = target_tri->edge_len[0];
  for(uint i = 1; i < 3; ++i){
    if(target_tri->edge_len[i] > max_len) max_len = target_tri->edge_len[i]; 
    else if(target_tri->edge_len[i] < min_len) 
      min_len = target_tri->edge_len[i];
  }

  // Get the interaction types of the vertices
  bool have_Hbond = false;
  for(int i = 0; i < 3; ++i)
    if(target_tri->vertices[i]->act_type == ACCEPTOR 
       || target_tri->vertices[i]->act_type == DONOR) 
      have_Hbond = true;
  
  // Get iterators to the bucket holding the similar model triangles
  int c = 0;
  for(uint i = 0; i < 3; ++i) c += target_tri->vertices[i]->act_type;
  c = hash_class[c];
  if(c >= 10) std::cerr << "need to handle doneptors\n";
  //bucket_vci bucket_begin, bucket_end;
  bool have_model_triangles = 
    A_model->get_bucket_iters(c, perimeter, max_len, min_len, 
                              bucket_begin, bucket_end);

  return((have_Hbond || A_allow_hphob_triangles) && have_model_triangles);
}
