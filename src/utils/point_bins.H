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

#ifndef POINT_BINS_HEADER_INCLUDED
#define POINT_BINS_HEADER_INCLUDED

#include <map>
#include <string_basics.H>
#include <point.H>
#include <math_basics.H>
#include <my_float_array.H>

namespace ASCbase{

//! Divide up the space based on a coordinate aligned grid
/*! Assumptions:
 *    The positions may be transformed via a rigid transformation.
 *    The type _Tp has a field "pos" that holds the 3D position of the type.
 *    It is expected that this will be used for the query sites; for that
 *    reason it is not necessarily optimized for speed of construction.
 */
template<typename _Tp>
class point_bins{

  // Concept requirements
#ifdef __GLIBCXX__
  __glibcxx_class_requires(_Tp, _ForwardIteratorConcept);
#else
  __glibcpp_class_requires(_Tp, _ForwardIteratorConcept);
#endif

public:
  struct three_idz{
    int idz[3];
  };

  template<typename T>
  struct point_bin_t{
    std::vector<const my_float_t *> positions; //!< Positions stored in the bin
    std::vector<T> assoc_data;   //!< Data associated with the positions
  };

  typedef typename std::vector< point_bin_t< _Tp > > bin_vec;
  typedef typename bin_vec::iterator bin_vi;
  typedef typename bin_vec::const_iterator bin_vci;
  typedef typename std::multimap<bin_vci, const my_float_t*> bin2vert_mmap;

  //! Default constructor -- sometimes the array cannot be constructed at class
  //! construction
  point_bins()
  { init(); }

  point_bins(_Tp begin, _Tp end, my_float_t bin_width = 4.5)
  {
    init(); 
    setup_grid(begin, end, bin_width);
  }

  point_bins(const point_bins<_Tp> &other) 
  {
    do_copy(other);
  }

  const point_bins&
  operator=(const point_bins& other)
  {
    if(&other != this){
      do_copy(other);
    }
    return *this;
  }

  void
  setup_grid(_Tp begin, _Tp end, const my_float_t bin_width)
  {
  // Have STL (or is it g++?) check if the given range is valid
#ifdef __GLIBCXX__
    __glibcxx_requires_valid_range(begin, end);
#else
    __glibcpp_requires_valid_range(begin, end);
#endif

    if(A_bins_are_setup){
      std::cerr << "Bins are already setup - cannot set them up again\n";
      return;
    }
    A_begin = begin;
    A_end = end;
    A_bin_width = bin_width + 0.01;
    A_bin_width_inv = 1.0 / bin_width;
    if(begin == end) return;

    my_float_t min_pos[3], max_pos[3];
    get_min_max_coords(begin, end, min_pos, max_pos);
    determine_space(min_pos, max_pos, A_lower_corner, A_upper_corner,
                    A_num_bins);
    populate_grid(begin, end, A_bin_width, &A_bins);
//    A_bins_for_moved_pts->resize(A_num_bins[0]*A_num_bins[1]*A_num_bins[2]);

    A_plane_normals[0] = 1.0;
    A_plane_normals[4] = 1.0;
    A_plane_normals[8] = 1.0;
    std::copy(A_lower_corner, A_lower_corner + 3, A_orig_lower);
    std::copy(A_upper_corner, A_upper_corner + 3, A_orig_upper);

    A_bins_are_setup = true;
  }

  //! Transform the bins
  /*! NOTE: this does not transform the positions, and must be used in 
   * conjunction with moving of the positions 
   */
  void
  transform_bins(const my_float_t* R, const my_float_t* T)
  {
    if(!A_bins_are_setup) return;

    // Rotate the normals by R
    my_float_t scratch[9];
    std::copy(A_plane_normals, A_plane_normals + 9, scratch);
    my_gemm(3, 3, 3, 1.0, scratch, 3, R, 3, A_plane_normals, 3, 0.0);

    // Transform upper and lower corner
    std::copy(A_upper_corner, A_upper_corner + 3, scratch);
    std::copy(T, T + 3, A_upper_corner);
    my_gemm(1, 3, 3, 1.0, scratch, 3, R, 3, A_upper_corner, 3, 1.0);
    std::copy(A_lower_corner, A_lower_corner + 3, scratch);
    std::copy(T, T + 3, A_lower_corner);
    my_gemm(1, 3, 3, 1.0, scratch, 3, R, 3, A_lower_corner, 3, 1.0);

    A_bins_are_transformed = true;
  }

  void
  inverse_transform_bins(const my_float_t* R, const my_float_t* T)
  {
    if(!A_bins_are_setup) return;

    my_float_t R_transpose[9];
    for(int i = 0; i < 3; ++i)
      for(int j = 0; j < 3; ++j) R_transpose[3*j + i] = R[3*i + j];
    my_float_t T_inverse[3];
    for(int i = 0; i < 3; ++i) T_inverse[i] = -1.0 * T[i];

    // Rotate the normals by the transpose of R
    my_float_t scratch[9];
    std::copy(A_plane_normals, A_plane_normals + 9, scratch);
    my_gemm(3, 3, 3, 1.0, scratch, 3, R_transpose, 3, A_plane_normals, 3, 0.0);

    // Transform upper and lower corner (C_prime = (C - T)R^t)
    std::copy(A_upper_corner, A_upper_corner + 3, scratch);
    my_axpy(3, 1.0, T_inverse, 1, scratch, 1);
    my_gemm(1, 3, 3, 1.0, scratch, 3, R_transpose, 3, A_upper_corner, 3, 0.0);
    std::copy(A_lower_corner, A_lower_corner + 3, scratch);
    my_axpy(3, 1.0, T_inverse, 1, scratch, 1);
    my_gemm(1, 3, 3, 1.0, scratch, 3, R_transpose, 3, A_lower_corner, 3, 0.0);

    A_bins_are_transformed = true;
  }

  void
  revert()
  {
    if(!A_bins_are_setup) return;

    std::copy(A_orig_lower, A_orig_lower + 3, A_lower_corner);
    std::copy(A_orig_upper, A_orig_upper + 3, A_upper_corner);
    std::fill(A_plane_normals, A_plane_normals + 9, 0.0);
    A_plane_normals[0] = 1.0;
    A_plane_normals[4] = 1.0;
    A_plane_normals[8] = 1.0;
    A_bins_are_transformed = false;
  }

  void
  get_min_max_coords(_Tp begin, _Tp end, my_float_t *min_pos, 
                     my_float_t *max_pos)
  {
    std::copy(begin->pos, begin->pos + 3, min_pos);
    std::copy(begin->pos, begin->pos + 3, max_pos);
    for(_Tp my_iter = begin; my_iter != end; ++my_iter)
      for(int i = 0; i < 3; ++i){
        if(max_pos[i] < my_iter->pos[i]) max_pos[i] = my_iter->pos[i];
        if(my_iter->pos[i] < min_pos[i]) min_pos[i] = my_iter->pos[i];
      }
  }

  void
  determine_space(const my_float_t *min_pos, const my_float_t *max_pos,
                  my_float_t *min_corner, my_float_t *max_corner, int *num_bins)
  {
    /// We want a large enough buffer so that moved points do not move into
    // any of the "edge" bins
    for(int i = 0; i < 3; ++i){
      min_corner[i] =
        A_bin_width * std::floor(min_pos[i] / A_bin_width) - 3*A_bin_width;
      max_corner[i] =
        A_bin_width * std::ceil(max_pos[i] / A_bin_width) + 3*A_bin_width;
      num_bins[i] =
        static_cast<int>(std::ceil((max_corner[i] - min_corner[i])
                                   / A_bin_width));
    }
  }

  void
  populate_grid(_Tp begin, _Tp end, const my_float_t bin_width,
                bin_vec *bins)
  {
    bins->resize(A_num_bins[0]*A_num_bins[1]*A_num_bins[2]);
    A_bins_of_moved_pts.resize(A_num_bins[0]*A_num_bins[1]*A_num_bins[2]);
//    A_current_bin_idz.resize(A_num_bins[0]*A_num_bins[1]*A_num_bins[2]);

    for(_Tp curr = begin; curr != end; ++curr){
      three_idz my_idz;
      int *idx = my_idz.idz;
      for(int i = 0; i < 3; ++i)
        idx[i] = 
          static_cast<int>(std::floor((curr->pos[i] - A_lower_corner[i]) *
                                      A_bin_width_inv));
      
      A_orig_bin_idz.push_back(my_idz);

      for(int i = -1; i < 2; ++i){
        bin_vi cbin = A_bins.begin();
        cbin += (idx[0] + i) * A_num_bins[1] * A_num_bins[2];
        cbin += (idx[1] - 1) * A_num_bins[2];
        for(int j = -1; j < 2; ++j){
          for(int k = -1; k < 2; ++k){
            (cbin + idx[2] + k)->positions.push_back(curr->pos);
            (cbin + idx[2] + k)->assoc_data.push_back(curr);
          }
          cbin += A_num_bins[2];
        }
      }
    }
    A_current_bin_idz.resize(A_orig_bin_idz.size());
    std::copy(A_orig_bin_idz.begin(), A_orig_bin_idz.end(),
              A_current_bin_idz.begin());
  }

  ~point_bins() 
  { 
    init();
  }

  inline void
  get_bin(const my_float_t *pt, bin_vci *nbrs, bin_vci *moved_nbrs) const
  {
    if(!A_bins_are_setup){
      *nbrs = A_bins.end();
      *moved_nbrs = A_bins.end();
      warn("point_bins", "get_bin", 
           "Cannot get neighborhood before setting up grid");
      return;
    }

    int idx[3];
    if(get_bin_idz(pt, idx)){
      *nbrs = A_bins.begin() + (idx[0]*A_num_bins[1]*A_num_bins[2] + 
                                idx[1]*A_num_bins[2] + idx[2]);
      *moved_nbrs = A_bins_of_moved_pts.begin() + 
                    (idx[0]*A_num_bins[1]*A_num_bins[2] + 
                     idx[1]*A_num_bins[2] + idx[2]);
    }else{
     *nbrs = A_bins.end();
     *moved_nbrs = A_bins.end();
    }
  } 

/*
  void
  get_bins(const my_float_t *pts, const uint npts,
           bin2vert_mmap *pts_in_bins) const
  {
    const my_float_t *pts_end = pts + 3*npts;
    for(const my_float_t *pt = pts; pt < pts_end; pt += 3)
      pts_in_bins->insert(std::pair<bin_vci, const my_float_t*>(get_bin(pt), pt));
  }
*/

  //! Constant iterator to one past the last bin in the vector
  inline bin_vci
  bins_end() const
  { return A_bins.end(); }

  //! Constant iterator to one past the last bin in the vector
  inline bin_vci
  moved_bins_end() const
  { return A_bins_of_moved_pts.end(); }

  const bool
  bins_are_setup() const
  { return A_bins_are_setup; }

  //! Check if the given data point needs to have its bins updated
  void
  update_point_bin(_Tp &data)
  {
    // At the present we only want to support updating if the template 
    // parameter is a random access iterator
#ifdef __GLIBCXX__
    __glibcxx_function_requires(_RandomAccessIteratorConcept<_Tp>);
#else
    __glibcpp_function_requires(_RandomAccessIteratorConcept<_Tp>);
#endif

    if(data < A_begin || data >= A_end){
      std::cerr << "Data point is not in the bins, cannot update bin "
                << "locations\n";
      return;
    }

    // Compute the bin for the pt, given the current position
    size_t idx = data - A_begin;
    int idz[3];
    if(!get_bin_idz(data->pos, idz)){
      std::cerr << "A point has moved out side of the bins, please increase "
                << "the buffer around the\ninitial positions and rerun the "
                << "computations\n";
      return;
    }

    // First check if we have not moved from the current bin
    const int *curr_bin_idz = A_current_bin_idz[idx].idz;
    if(idz[0] == curr_bin_idz[0] && idz[1] == curr_bin_idz[1] &&
       idz[2] == curr_bin_idz[2]) return;

    // Next check if we have moved to the point's original bin
    const int *orig_bin_idz = A_orig_bin_idz[idx].idz;
    if(idz[0] == orig_bin_idz[0] && idz[1] == orig_bin_idz[1] &&
       idz[2] == orig_bin_idz[2]){
      std::copy(orig_bin_idz, orig_bin_idz + 3, curr_bin_idz);
      return;
    }

    for(int i = 0; i < 3; ++i){
      if(idz[i] <= 0 or (idz[i] + 1) >= A_num_bins[i]){
        std::cerr << "A point has moved almost out side of the bins, please increase "
                  << "the buffer around the\ninitial positions and rerun the "
                  << "computations\n";
        return;
      }
    }

    for(int i = -1; i < 2; ++i){
      bin_vi cbin = A_bins_of_moved_pts.begin();
      cbin += (idz[0] + i) * A_num_bins[1] * A_num_bins[2];
      cbin += (idz[1] - 1) * A_num_bins[2];
      for(int j = -1; j < 2; ++j){
        for(int k = -1; k < 2; ++k){
          (cbin + idz[2] + k)->positions.push_back(data->pos);
          (cbin + idz[2] + k)->assoc_data.push_back(data);
        }
        cbin += A_num_bins[2];
      }
    }

#if 0
    // this isn't difficult, but will take some time to write or think
    // about it -- at the present I need the prototype of the method more
    // than testing this piece of code

    int idz_diff[3];
    for(int i = 0; i < 3; ++i) idz_diff[i] = idz[i] - curr_bin_idz[i];

    if((idz[0] == -1 || idz[0] == 1) && idz[1] == 0 && idz[2] == 0){
    
    }else if(idz[0] == 0 && (idz[1] == -1 || idz[1] == 1) && idz[2] == 0){

    }else if(idz[0] == 0 && idz[1] == 0 && (idz[2] == -1 || idz[2] == 1)){

    }else{
      // We aren't saving that much so just add it to everything?
      
    }
#endif

    std::copy(idz, idz + 3, A_current_bin_idz[idx].idz);
  }

  void
  remove_moved_points_from_bins()
  {
    bin_vi b;
    for(b = A_bins_of_moved_pts.begin(); b < A_bins_of_moved_pts.end(); ++b)
      b->clear();

    std::copy(A_orig_bin_idz.begin(), A_orig_bin_idz.end(), 
              A_current_bin_idz.begin());
  }

  void
  write_bounds(std::ostream &out) const
  {
    std::cout << "A_upper_corner: " << A_upper_corner[0] << " "
              << A_upper_corner[1] << " " << A_upper_corner[2] << "\n";
    std::cout << "A_lower_corner: " << A_lower_corner[0] << " "
              << A_lower_corner[1] << " " << A_lower_corner[2] << "\n";
    std::cout << "Normals:" ;
    for(int i = 0; i < 9; ++i) std::cout << " " << A_plane_normals[i];
  }

private:

  void
  init()
  {
    A_bin_width = my_float_max;
    std::fill(A_upper_corner, A_upper_corner + 3, 0.0);
    std::fill(A_lower_corner, A_lower_corner + 3, 0.0);
    std::fill(A_num_bins, A_num_bins + 3, 0);
    A_bins_are_setup = false;
    A_bins_are_transformed = false;
    std::fill(A_plane_normals, A_plane_normals + 9, 0.0);
  } 

  //! Copies could be expensive, but depend greatly on the template parameter
  /*! We are assuming that most uses will have template parameters that are
   * pointers or iterators (hence we assume _Tp->pos exists and is a valid
   * 3D point).
   */
  void
  do_copy(const point_bins<_Tp> &src)
  {
    init();
    A_bins_are_setup = src.A_bins_are_setup;
    A_bins_are_transformed = src.A_bins_are_transformed;
    A_bin_width = src.A_bin_width;
    A_bin_width_inv = src.A_bin_width_inv;
    std::copy(src.A_upper_corner, src.A_upper_corner + 3, A_upper_corner);
    std::copy(src.A_lower_corner, src.A_lower_corner + 3, A_lower_corner);
    std::copy(src.A_num_bins, src.A_num_bins + 3, A_num_bins);
    std::copy(src.A_plane_normals, src.A_plane_normals + 9, A_plane_normals);

    A_bins.reserve(src.A_bins.size());
    for(bin_vci b = src.A_bins.begin(); b < src.A_bins.end(); ++b)
      A_bins.push_back(*b);
  
    A_current_bin_idz.resize(src.A_current_bin_idz.size());
    std::copy(src.A_current_bin_idz.begin(), src.A_current_bin_idz.end(),
              A_current_bin_idz.begin());
    A_orig_bin_idz.resize(src.A_orig_bin_idz.size());
    std::copy(src.A_orig_bin_idz.begin(), src.A_orig_bin_idz.end(),
              A_orig_bin_idz.begin());
  }

  bool
  get_bin_idz(const my_float_t *pt, int *idx) const 
  {
    if(A_bins_are_transformed){
      // use the normals to each plane and the lower corner as point on 
      // the plane

      // Get signed distance of point to each "lower" plane
      my_float_t signed_dist[3];
      for(int i = 0; i < 3; ++i){
        my_float_t V[3];
        vector(3, pt, A_lower_corner, V);
        signed_dist[i] = dot(V, A_plane_normals + 3*i);
      }

      // Get the bin indices
      for(int i = 0; i < 3; ++i){
        if(signed_dist[i] < 0.0) return false;
        idx[i] = static_cast<int>(std::floor(signed_dist[i] / A_bin_width));
        if(idx[i] >= A_num_bins[i]) return false;
      }
    }else{
      // Check if point falls outside the range of the bins
      for(int i = 0; i < 3; ++i)
        if(pt[i] < A_lower_corner[i] || A_upper_corner[i] < pt[i])
          return false;
  
      // Compute the index for the bin containing the point
      for(int i = 0; i < 3; ++i){
        idx[i] = static_cast<int>(std::floor((pt[i] - A_lower_corner[i]) 
                                             / A_bin_width));
        if(idx[i] >= A_num_bins[i]) return false;
      }
    }
    return true;
  }

  _Tp A_begin;
  _Tp A_end;
  bool A_bins_are_setup;
  my_float_t A_bin_width; //!< Side length of a bin
  my_float_t A_bin_width_inv; //!< multiplicative inverse (1 over) side length
 
  my_float_t A_upper_corner[3];
  my_float_t A_lower_corner[3];
  my_float_t A_orig_lower[3];   //!< Save the original lower corner
  my_float_t A_orig_upper[3];   //!< Save the original upper corner
  my_float_t A_plane_normals[9];
  bool A_bins_are_transformed;

  int A_num_bins[3];
  bin_vec A_bins;
  bin_vec A_bins_of_moved_pts;
  std::vector<three_idz> A_current_bin_idz;
  std::vector<three_idz> A_orig_bin_idz;

};
} // end namespace ASCbase
#endif
