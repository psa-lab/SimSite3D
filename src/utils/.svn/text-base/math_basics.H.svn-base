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

/*
 * $Source: /psa/share/repository/pfizer_proj/src/basics/math_basics.H,v $
 * $Revision: 1.11 $
 * $Author: vanvoor4 $
 * $Date: 2009/01/12 20:58:39 $
 * 
 * $Log: math_basics.H,v $
 * Revision 1.11  2009/01/12 20:58:39  vanvoor4
 * Additional math functions to support ICP, etc.
 *
 * Revision 1.10  2008/07/29 14:28:09  vanvoor4
 * Updated the LSE fitting function and wrote some short comments
 *
 * Revision 1.9  2008/07/28 15:16:00  vanvoor4
 * Added a number of functions to help with ICP as well as projecting
 * points to closest points on moth eaten spheres.
 *
 * Revision 1.8  2008/05/13 15:40:15  vanvoor4
 * Removed compare_double() -- no longer used
 *
 * Revision 1.7  2008/04/28 17:48:37  vanvoor4
 * Allow corresponding point for point to line treat the line as either
 * a line segment or line (infinite).
 *
 * Revision 1.6  2008/01/08 19:06:00  vanvoor4
 * Updated the interfaces to corresponding points
 *
 * Revision 1.5  2007/12/17 21:05:16  vanvoor4
 * Added a few functions for pt to surface / curve
 *
 * Revision 1.4  2007/09/24 15:33:48  vanvoor4
 * Added support for quaternions rather than ortho matrices.
 * The advantage of quaternions is when multiple rotations are used.
 * The method of BKP Horn for 3pt align is a bit "faster" than
 * using homebrew code for eigen values.
 *
 * Revision 1.3  2007/08/29 20:16:14  vanvoor4
 * Added the dist squared function
 *
 * Revision 1.2  2007/08/21 15:35:41  vanvoor4
 * Added license header
 *
 * Revision 1.1  2007/02/07 15:59:32  vanvoor4
 * The math functions previously in utils
 *
 *
 * 
 * 
 */

#ifndef MATH_BASICS_HEADER_INCLUDED
#define MATH_BASICS_HEADER_INCLUDED

#include <cmath>
#include <vector>
#include <list>
#include <types.H>
#include <Quaternion.H>
#include <mat_ops.H>

static inline my_float_t
deg2rad(my_float_t deg)
{ return deg/180.0 * M_PI; }

//! Distance squared -- try to reduce number of calls to sqrt
inline my_float_t 
dist_squared(const my_float_t a[3], const my_float_t b[3])
{
  my_float_t sum = 0, tmp = 0;
  for(uint i = 0; i < 3; ++i){
    tmp = b[i] - a[i];
    sum += tmp*tmp;
  }
  return sum;
}

//! distance is sqrt(dist_squared())
inline my_float_t 
dist(const my_float_t a[3], const my_float_t b[3])
{ return std::sqrt(dist_squared(a,b)); }

//! Set X to a unit vector in the same direction and return the magitude of X
//! (before normalization)
inline my_float_t 
normalize(my_float_t* X, const uint ndim=3)
{
  my_float_t* x = X;
  my_float_t mag = 0;
  for(uint i = 0; i < ndim; ++i, ++x) mag += (*x) * (*x);
  mag = std::sqrt(mag);
  x = X;
  for(uint i = 0; i < ndim; ++i, ++x) *x /= mag;
  return mag;
}

//! Get the unit vector in the direction of the vector from tail to head
inline my_float_t
unit_vector(my_float_t* n, const my_float_t* head, const my_float_t* tail)
{
  for(uint i = 0; i < 3; ++i, ++head, ++tail, ++n) *n = *head - *tail;
  return normalize(n - 3);
}

//my_float_t dot(const my_float_t* U, const my_float_t* V);

//! Compute U o V for 3D vectors
inline my_float_t 
dot(const my_float_t* U, const my_float_t* V)
{
  my_float_t rv = 0;
//  const my_float_t *v = V;
//  const my_float_t *u = U;
  for(uint i = 0; i < 3; ++i, ++U, ++V) rv += *U * (*V);
  return rv;
}

//! Compute U x V
void cross(const my_float_t* U, const my_float_t* V, my_float_t* rv);


//! Get the mid point of 2 3D points
/*! Assumptions
 *    A, B, M are all of length at least N
 */ 
inline void
mid_point(const my_float_t* A, const my_float_t* B, my_float_t* M,
          const int N = 3)
{
  std::copy(A, A + N, M);
  my_axpy(N, 1.0, B, 1, M, 1);
  for(int i = 0; i < N; ++i) M[i] *= 0.5;
}

//! Argmin that also keeps track of the minimum value
template<class T> inline const int 
argmin(T* vals, int len, T* min_val)
{
  if(len < 1) return -1;
  int idx = 0;
  *min_val = *vals;
  for(int i = 0; i < len; ++i)
    if(*min_val > vals[i]){
      *min_val = vals[i];
      idx = i;
    } 
  return idx;
}

//! Argmax that also keeps track of the maximum value
template<class T> inline const int 
argmax(T* vals, int len, T* max_val)
{
  if(len < 1) return -1;
  int idx = 0;
  *max_val = *vals;
  for(int i = 0; i < len; ++i)
    if(*max_val < vals[i]){
      *max_val = vals[i];
      idx = i;
    } 
  return idx;
}

//! Argmax of the absolute value that also keeps track of the maximum value
template<class T> inline const int 
argmaxabs(T* vals, int len, T* max_val)
{
  if(len < 1) return -1;
  int idx = 0;
  my_float_t tmp = (*vals < 0 ? -1.0 * (*vals) : *vals);
  for(int i = 0; i < len; ++i){
    my_float_t m = (vals[i] < 0 ? -1.0 * vals[i] : vals[i]);
    if(tmp < m){
      tmp = m;
      idx = i;
    } 
  }
  *max_val = *(vals + idx);
  return idx;
}

//!< should be replaced by the matrix vector mult blas routine
/*!
 * Assumptions:<OL>
 * <LI> stride is 3
 * <LI> The rotmat is 3x3 (9 el)
 * <LI> tvec is 3x1 (3 el)
 * <LI> n is the length of pos/3 == length of new_pos/3
 * <LI> Memory for pos and new_pos must be grabbed prior to call.
 * </OL>
 *
 * Replace by BLAS or similar routine in future
 */
void move_positions(const uint n, my_float_t* new_pos, const my_float_t* pos, 
                    const my_float_t* rotmat, const my_float_t* tvec);

//! Get the rotation matrix to align translated global coordinates to the 
//! local coordinate system based on vectors U and V
/*! 
 * First subtract translation component, and then premultiply R^t by the
 * translated global coord to get to local coordinate system
 * Premulitply R by local coord and add the translation component to get to 
 * global coordinates.  Note: this function does not compute the translation
 * component since U & V are assumed to be unit vectors.
 * <OL>
 * <LI> place vectors U and V in the XY plane to get U' and V'
 * <LI> rotate U' and V' about Z axis so that U'' lies on the X axis
 * <LI> return the rotation matrix to undo the alignment
 * </OL>
 */
void get_local_orientation(const my_float_t *U, const my_float_t *V, 
                           my_float_t* R);

//! Align the vector V to the X-axis -- assumption is V lies in the XY-plane
void align_to_X_axis(const my_float_t *V, Quaternion* Q);

//! Find a quaternion Q which rotates X_norm to be the same as Y_norm
void align_planes(const my_float_t* X_norm, const my_float_t* Y_norm, 
                  Quaternion* Q);

//! get normal to plane defined by three points in X
void plane_normal(const my_float_t* X, my_float_t* N);

//! Given XA (X[0:3[), XB (X[3:6[) and XC (X[6:9[), set U as the vector from 
//! XB to XA and V as the vector from XB to XC.
void get_vectors(const my_float_t* X, my_float_t* U, my_float_t* V);

//! Compute the weighted three point alignment for X to Y using the method
//! of Berthold K. P. Horn for quaternions.
void three_pt_align(const my_float_t* X, const my_float_t *Y,
                    const my_float_t *W, Quaternion* Q, my_float_t* T, 
                    my_float_t* newX);

//! If the discriminant is positive return both real roots, otherwise return
//! false
bool quadratic_roots(const my_float_t *P, my_float_t *roots);

//! Get one of the three roots from a full cubic polynomial
my_float_t get_root_from_cubic(const my_float_t *P);

//! Assumes the reduced quartic (x^4 + px^2 + qx + r = 0) has 4 real roots
//! Returns all four roots
void quartic_roots(const my_float_t p, const my_float_t q, const my_float_t r,
                   my_float_t *roots);

//! Compute the rigid transformation moving the target points to the model
//! points by minimizing the least squared error between corresponding points
/*!
 * @param model_pts const pointer to an array of length 3*npts -- point 
 *        representation of the model/query object
 * @param target_pts const pointer to an array of length 3*npts -- point
 *        representation of the target/dataset object 
 * @param npts number of points to fit
 * @param Q quaternion to store the computed rotation
 * @param T array of length 3 to store the compute translation
 * @param W Null or array of length npts for a weighted fit
 * @return False if npts < 3
 */
bool lse_3D_fit(const my_float_t* model_pts, const my_float_t* target_pts, 
                const size_t npts, Quaternion *Q, my_float_t *T,
                const my_float_t* W=0);

//! Assume a stride of 3 
/*! Should this centroid be updated as with Laurentio's method for 
 * updating mean and variance after each observation or do we get enough
 * precision here?
 */
void centroid_3D(const my_float_t* pts, const size_t npts, my_float_t* C,
                 const my_float_t* W=0);

//! Compute the upper triangular pair wise squared distance of 3D points
void pair_dist_squared(const my_float_t* pts, size_t npts, 
                       const my_float_t *d_sq);

void max_pair_dist(const my_float_t* pts, const size_t npts, size_t *i,
                   size_t *j, my_float_t *d);

//! Find the best corresponding point on the triangle V
/*! NOTE: this function is very bad in that it returns true if
 *  abs(prev_best_d) < abs(d)
 */
bool corresponding_point(const my_float_t* pt, const my_float_t *Vi, 
                         const my_float_t *Vj, const my_float_t *Vk, 
		         const my_float_t prev_best_d,
                         my_float_t* d, my_float_t* pt_on_plane);

//! Find the best corresponding poin on the line V
bool corresponding_point(const my_float_t* pt, const my_float_t* Vi, 
                         const my_float_t* Vj, my_float_t* d, 
                         my_float_t* pt_on_line,
                         const bool LINE_SEGMENT = true);

//! 3D version to check if the point is in the triangle 
bool pt_inside_triangle(const my_float_t *N, const my_float_t *Vi,
                        const my_float_t *Vj, const my_float_t *Vk,
                        const my_float_t *pt);

//! 2D version to check if the point is in the triangle 
bool pt_inside_triangle(const my_float_t* x, const my_float_t* y);

//! Compute distance from point to plane
my_float_t pt2plane_dist(const my_float_t *N, const my_float_t *p0, 
                         const my_float_t *pt);

//! Compute projection of point on plane
my_float_t pt2plane_proj(const my_float_t *N, const my_float_t *p0, 
                         const my_float_t *pt, my_float_t* pt_on_plane);

//! Compute projection of point on plane
void pt2plane_proj(const my_float_t *N, const my_float_t d, 
                   const my_float_t *pt, my_float_t* pt_on_plane);

//! Compute projection of point on sphere
void pt2sphere_proj(const my_float_t *c, const my_float_t r, 
                    const my_float_t *pt, my_float_t *pt_on_sphere);

//! Determine the intersection of the two planes
void intersect_planes(const my_float_t *N1, const my_float_t d1,
                      const my_float_t *N2, const my_float_t d2,
                      my_float_t *m, my_float_t *b);

class sphere{
public:
  sphere(const my_float_t *center, const my_float_t radius)
  {
    std::copy(center, center + 3, A_c);
    A_r = (radius < 0 ? -1.0 * radius : radius);
    A_rsquared = radius*radius; 
    A_vol = 0.0;
  }

  bool
  operator==(const sphere &S) const 
  {
    if(S.A_r == A_r && S.A_c[0] == A_c[0] && S.A_c[1] == A_c[1] && 
       S.A_c[2] == A_c[2]) return true;
    return false;
  } 
  
  ~sphere() {;}

  my_float_t volume()
  {
    if(A_vol == 0.0) A_vol = 4.0/3.0 * M_PI * A_rsquared * A_r; 
    return A_vol;
  }

  void
  proj2surf(const my_float_t *pt, my_float_t *pt_on_sphere) const
  {
    my_float_t dir[3];
    unit_vector(dir, pt, A_c);
    std::copy(A_c, A_c + 3, pt_on_sphere);
    my_axpy(3, A_r, dir, 1, pt_on_sphere, 1);
  }

  my_float_t A_c[3];
  my_float_t A_r;
  my_float_t A_rsquared;
  my_float_t A_vol;
};

// Dump in sphere class at some point
bool
line_thru_sphere(const my_float_t *m, const my_float_t *b,
                 const sphere &S, my_float_t *t1, my_float_t *t2);

bool
intersect(const sphere &A, const sphere &B, my_float_t *N, my_float_t *p0,
          my_float_t *radius);

class circular_segment{
public:
  //friend bool merge(circular_segment *A, circular_segment *B);

  typedef enum{
    NO_POINT = 0,
    FIRST_POINT,
    SECOND_POINT,
    BOTH_POINTS
  }intersect_type; 

  circular_segment(const my_float_t *p0, const my_float_t *p1, 
                   const my_float_t *in_dir, const my_float_t *center)
  {
    std::copy(p0, p0 + 3, A_p0);
    std::copy(p1, p1 + 3, A_p1);
    std::copy(in_dir, in_dir + 3, A_in_dir);
    std::copy(center, center + 3, A_circle_center);
    for(size_t i = 0; i < 3; ++i) A_mid_pt[i] = 0.5 * (A_p0[i] + A_p1[i]);
  }

  circular_segment(const circular_segment &src)
  {
    if(this == &src) return;
    std::copy(src.A_p0, src.A_p0 + 3, A_p0);
    std::copy(src.A_p1, src.A_p1 + 3, A_p1);
    std::copy(src.A_circle_center, src.A_circle_center + 3, A_circle_center);
    std::copy(src.A_in_dir, src.A_in_dir + 3, A_in_dir);
    std::copy(src.A_mid_pt, src.A_mid_pt + 3, A_mid_pt);
  }

  // Idea -- if the vector from this chord's end point to the end point from 
  // the other chord has a positive projection onto the in direction, then
  // that end point from the other chord is inside this circular segment
  intersect_type
  intersects(const circular_segment &other_seg)
  {
#if 0
    my_float_t p0_to_x[3];
    std::copy(other_seg.A_p0, other_seg.A_p0 + 3, p0_to_x);
    my_axpy(3, -1.0, A_p0, 1, p0_to_x, 1);
    bool p0_inside = false;
    if(dot(A_in_dir, p0_to_x) > 0.0) p0_inside = true;

    std::copy(other_seg.A_p1, other_seg.A_p1 + 3, p0_to_x);
    my_axpy(3, -1.0, A_p0, 1, p0_to_x, 1);
    bool p1_inside = false;
    if(dot(A_in_dir, p0_to_x) > 0.0) p1_inside = true;
#endif

    // Still have a problem if both end points are the same

    bool p0_inside = false;
    bool p1_inside = false;
    if(contains(other_seg.A_p0)) p0_inside = true;
    if(contains(other_seg.A_p1)) p1_inside = true;

    intersect_type rv = NO_POINT;
    if(false == p0_inside && false == p1_inside) rv = NO_POINT;
    else if(true == p0_inside && true == p1_inside) rv = BOTH_POINTS;
    else if(p0_inside) rv = FIRST_POINT;
    else if(p1_inside) rv = SECOND_POINT;

    return rv;
  }

  bool merge(circular_segment *B)
  {
    circular_segment::intersect_type A_pts_in_B = intersects(*B);
    circular_segment::intersect_type B_pts_in_A = B->intersects(*this);
    if(NO_POINT == A_pts_in_B && NO_POINT == B_pts_in_A) return false;
    if(BOTH_POINTS == B_pts_in_A) return true;
    if(BOTH_POINTS == A_pts_in_B){
      std::cout << "Attempting to assign to *this\n";
      *this = *B;
      return true;
    }
  
    // else update the end points of the chord (and arc)
    if(A_pts_in_B == FIRST_POINT){
      if(B_pts_in_A == FIRST_POINT) std::copy(B->A_p1, B->A_p1 + 3, A_p0);
      else std::copy(B->A_p0, B->A_p0 + 3, A_p0);
    }else{ 
      if(B_pts_in_A == FIRST_POINT) std::copy(B->A_p1, B->A_p1 + 3, A_p1);
      else std::copy(B->A_p0, B->A_p0 + 3, A_p1);
    }
  
    // MUST UPDATE THE IN DIRECTION !!!!!!!!!!!!!!!!!!
    for(size_t i = 0; i < 3; ++i) A_mid_pt[i] = (A_p0[i] + A_p1[i]) / 2;
    unit_vector(A_in_dir, A_mid_pt, A_circle_center);
  
    my_float_t interior_vector[3];
    if(B_pts_in_A == FIRST_POINT)
      unit_vector(interior_vector, B->A_p0, A_circle_center);
    else unit_vector(interior_vector, B->A_p1, A_circle_center);
  
    if(dot(A_in_dir, interior_vector) < 0.0)
      for(size_t i = 0; i < 3; ++i) A_in_dir[i] *= -1.0;
  
    return true; 
  }

  // Ignore the radius for now
  bool
  contains(const my_float_t *pt)
  {
    std::cout << "\nCheck for Chord containment\n";
    std::cout << "pt: " << pt[0] << " " << pt[1] << " " << pt[2] << "\n";
    std::cout << "mid_pt: " << A_mid_pt[0] << " " << A_mid_pt[1] << " "
              << A_mid_pt[2] << "\n";

    // Doubt we need unit vector here
    my_float_t pt_dir[3];
    unit_vector(pt_dir, pt, A_mid_pt);
    std::cout << "pt_dir: " << pt_dir[0] << " " << pt_dir[1] << " "
              << pt_dir[2] << "\n";
    std::cout << "dot prod: " << dot(A_in_dir, pt_dir) << "\n\n";
    if(dot(A_in_dir, pt_dir) < 0.0) return false;
    return true;
  }    
  
  const my_float_t* 
  chord_endpt_one() const 
  { return A_p0; }

  const my_float_t* 
  chord_endpt_two() const 
  { return A_p1; }

private:
  my_float_t A_p0[3];
  my_float_t A_p1[3];
  my_float_t A_mid_pt[3];
  my_float_t A_in_dir[3];  //!< Vector lying in the plane and _|_ to the chord
  my_float_t A_circle_center[3];
};

bool
my_close_point(const sphere &S, const std::vector<sphere> &nbrs, 
               const my_float_t *pt, const my_float_t tolerance,
               my_float_t *close_pt);

// Must update one of the end points as well as the direction of the line
// and the in direction
/*
bool merge(circular_segment *A, circular_segment *B)
{
  circular_segment::intersect_type A_pts_in_B = A->intersects(*B);
  circular_segment::intersect_type B_pts_in_A = B->intersects(*A);
  if(circular_segment::NO_POINT == A_pts_in_B == B_pts_in_A) return false;
  if(A_pts_in_B == circular_segment::BOTH_POINTS) return true;
  if(B_pts_in_A == circular_segment::BOTH_POINTS){
    *B = *A;
    return true;
  }

  // else update the end points of the chord (and arc)
  if(A_pts_in_B == circular_segment::FIRST_POINT){
    if(B_pts_in_A == circular_segment::FIRST_POINT){
      std::copy(A->A_p1, A->A_p1 + 3, B->A_p0);
    }else{ 
      std::copy(A->A_p1, A->A_p1 + 3, B->A_p1);
    }
  }else{ 
    if(B_pts_in_A == circular_segment::FIRST_POINT){
      std::copy(A->A_p0, A->A_p0 + 3, B->A_p0);
    }else{ 
      std::copy(A->A_p0, A->A_p0 + 3, B->A_p1);
    }
  }

  // Update vector parallel to the chord
  unit_vector(B->A_p0p1, B->A_p1, B->A_p0);
  
  // MUST UPDATE THE DIRECTION FOR B!!!!!!!!!!!!!!!!!!
  my_float_t mid_pt[3];
  for(size_t i = 0; i < 3; ++i) mid_pt[i] = (B->A_p0[i] + B->A_p1[i]) / 2;
  unit_vector(B->A_in_dir, mid_pt, B->A_circle_center);

  my_float_t interior_vector[3];
  if(A_pts_in_B == circular_segment::FIRST_POINT)
    unit_vector(interior_vector, A->A_p0, B->A_circle_center);
  else
    unit_vector(interior_vector, A->A_p1, B->A_circle_center);

  if(dot(B->A_in_dir, interior_vector) < 0.0)
    for(size_t i = 0; i < 3; ++i) B->A_in_dir[i] *= -1.0;

  return true; 
}
*/

//! Find the min and max corners of the axis aligned box holding the points
/*! Assumptions:
 *    The variables begin and end are iterators (should be constant)
 *    The objects iterated over have a pos variable
 *    min and max_pt have enough memory allocated
 *    The caller has checked to ensure begin and end are a valid range
 */
template<typename Const_Iter> void
get_bounds(Const_Iter begin, Const_Iter end, my_float_t* min_pt, 
           my_float_t* max_pt, int ndim = 3)
{
  // Spin through the positions to find the min and max values
  std::copy(begin->pos, begin->pos + ndim, max_pt);
  std::copy(begin->pos, begin->pos + ndim, min_pt);
  for(Const_Iter p = begin; p < end; ++p){
    for(int i = 0; i < ndim; ++i){
      if(max_pt[i] < p->pos[i]) max_pt[i] = p->pos[i];
      if(min_pt[i] > p->pos[i]) min_pt[i] = p->pos[i];
    }
  }
}


#endif
