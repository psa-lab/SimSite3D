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
 * $Source: /psa/share/repository/pfizer_proj/src/basics/math_basics.C,v $
 * $Revision: 1.16 $
 * $Author: vanvoor4 $
 * $Date: 2009-01-12 20:58:40 $
 * 
 * $Log: not supported by cvs2svn $
 * Revision 1.15  2008/07/29 14:28:09  vanvoor4
 * Updated the LSE fitting function and wrote some short comments
 *
 * Revision 1.14  2008/07/28 15:16:00  vanvoor4
 * Added a number of functions to help with ICP as well as projecting
 * points to closest points on moth eaten spheres.
 *
 * Revision 1.13  2008/05/13 15:40:15  vanvoor4
 * Removed compare_double() -- no longer used
 *
 * Revision 1.12  2008/04/30 13:14:14  vanvoor4
 * Projection of point to line was previously performed in a round about
 * manner.
 *
 * Revision 1.11  2008/04/28 17:48:37  vanvoor4
 * Allow corresponding point for point to line treat the line as either
 * a line segment or line (infinite).
 *
 * Revision 1.10  2008/02/26 19:10:22  vanvoor4
 * Added a check to make sure we don't have divide by zero
 *
 * Revision 1.9  2008/01/08 19:06:35  vanvoor4
 * Numerous errors corrected in corresponding point funcitons
 *
 * Revision 1.8  2007/12/17 21:05:27  vanvoor4
 * Added functions for pt to surf /curve
 *
 * Revision 1.7  2007/10/05 17:55:14  vanvoor4
 * Checks for numerically identity alignments are more robust
 *
 * Revision 1.6  2007/10/04 15:05:18  vanvoor4
 * Added checks for cos/sin values to be near 1 or 0.
 *
 * Revision 1.5  2007/10/02 17:26:36  vanvoor4
 * Need to multiply 4 not 3 points
 *
 * Revision 1.4  2007/09/24 15:34:32  vanvoor4
 * Added support for quaternions rather than ortho matrices.
 * The advantage of quaternions is when multiple rotations are used.
 * The method of BKP Horn for 3pt align is a bit "faster" than
 * using homebrew code for eigen values.
 *
 * Revision 1.3  2007/08/29 20:16:04  vanvoor4
 * Added the dist_squared function
 *
 * Revision 1.2  2007/08/21 15:35:09  vanvoor4
 * Added license header
 *
 * Revision 1.1  2007/02/07 15:59:20  vanvoor4
 * The math functions previously in utils
 *
 *
 */

#include <algorithm>
#include <iostream>
#include <vector>
#include <list>

#include <math_basics.H>
#include <mat_ops.H>
#include <iomanip>

static const my_float_t min_significant_value = 1E-07;

void 
move_positions(const uint n, my_float_t* new_pos, const my_float_t* pos, 
               const my_float_t* rotmat, const my_float_t* tvec)
{
  for(uint i = 0; i < n; ++i, pos += 3, new_pos += 3){
    std::copy(tvec, tvec + 3, new_pos);
    for(uint j = 0; j < 3; ++j)
      for(uint k = 0; k < 3; ++k) new_pos[j] += rotmat[3*j + k] * pos[k];
  }
}

void
get_local_orientation(const my_float_t *U, const my_float_t *V, my_float_t* R)
{
  // Find a unit normal to the plane containing U and V
  my_float_t U_cross_V[3];
  cross(U, V, U_cross_V);
  normalize(U_cross_V);

  // Get the rotation that places U and V in the XY plane
  my_float_t unit_Z[] = {0, 0, 1};
  Quaternion Q1;
  align_planes(U_cross_V, unit_Z, &Q1);
  my_float_t R_Q[9];
  Q1.get_ortho_rot_mat(R_Q);
  my_float_t tmp[6];
  std::copy(U, U + 3, tmp);
  std::copy(V, V + 3, tmp + 3);
  // U' = vecs_prime[0:3] ;  V' = vecs_prime[3:6]
  my_float_t vecs_prime[6];
  my_gemm(2, 3, 3, 1.0, tmp, 3, R_Q, 3, vecs_prime, 3, 0.0);

  // Rotate U and V about the Z axis to align U with the X axis 
  Quaternion Q2;
  align_to_X_axis(vecs_prime, &Q2);
  Q2.get_ortho_rot_mat(R_Q);
  my_float_t V_prime_prime[3];
  my_gemm(1, 3, 3, 1.0, vecs_prime + 3, 3, R_Q, 3, V_prime_prime, 3, 0.0);

  // Require the Y component of V (in the local coordinates) to be negative
  if(V_prime_prime[1] > 0.0){
    my_float_t tmp_Q[] = {0.0, 1.0, 0.0, 0.0};
    Quaternion Q3(tmp_Q, 4);
    Quaternion Q_tmp = Q1 * Q2;
    (Q_tmp * Q3).conjugate().get_ortho_rot_mat(R);
  }else (Q1 * Q2).conjugate().get_ortho_rot_mat(R);
}

// Add output for the rotated and translated X -- rotate X_cent and X at
// the same time
void
three_pt_align(const my_float_t* X, const my_float_t *Y, const my_float_t *W,
               Quaternion* Q, my_float_t* T, my_float_t* newX)
{
  my_float_t X_cent[3], Y_cent[3];
  my_float_t W_sum = W[0] + W[1] + W[2];
  for(uint i = 0; i < 3; ++i){
    X_cent[i] = (W[0]*X[i] + W[1]*X[3 + i] + W[2]*X[6 + i]) / W_sum;
    Y_cent[i] = (W[0]*Y[i] + W[1]*Y[3 + i] + W[2]*Y[6 + i]) / W_sum;
  }

  my_float_t X_V[9], Y_V[9];
  for(uint i = 0; i < 3; ++i)
    for(uint j = 0; j < 3; ++j){
      X_V[3*i+j] = X[3*i+j] - X_cent[j];
      Y_V[3*i+j] = Y[3*i+j] - Y_cent[j];
    }
 
  my_float_t X_norm[3], Y_norm[3];
  cross(X_V, X_V + 3, X_norm);
  cross(Y_V, Y_V + 3, Y_norm);
  normalize(X_norm);
  normalize(Y_norm);

  // Rotate the X plane so that the X_norm and Y_norm are the same
  Quaternion Q1;
  align_planes(X_norm, Y_norm, &Q1);
  my_float_t R_Q[9];
  Q1.get_ortho_rot_mat(R_Q);
  my_float_t X_V_R1[9];
  my_gemm(3, 3, 3, 1.0, X_V, 3, R_Q, 3, X_V_R1, 3, 0.0); 

  // compute the weighted C and S
  my_float_t C = 0;
  for(uint i = 0; i < 3; ++i) C += W[i] * dot(Y_V + 3*i, X_V_R1 + 3*i);
  my_float_t my_sum[] = {0, 0, 0};
  for(uint i = 0; i < 3; ++i){
    my_float_t tmp[3];
    cross(Y_V + 3*i, X_V_R1 + 3*i, tmp); 
    for(uint j = 0; j < 3; ++j) my_sum[j] += W[i] * tmp[j];
  }
  my_float_t S = dot(my_sum, Y_norm);
  my_float_t denom = std::sqrt(S*S + C*C);  

  // Note, we are post multiplying by the rotation matrices from the 
  // quaternions.  This means, to chain rotations, we need to multiply the
  // previous quaternion on the right by the latest rotation (quaternion)
  Quaternion Q2(C / denom, S / denom, Y_norm);
  *Q = Q1 * Q2;
  Q->get_ortho_rot_mat(R_Q);

  // Translation component
  my_float_t X_cent_R[3];
  my_gemm(1, 3, 3, 1.0, X_cent, 3, R_Q, 3, X_cent_R, 3, 0.0);
  for(uint i = 0; i < 3; ++i) T[i] = Y_cent[i] - X_cent_R[i];

  // Apply transformation to X and store in newX
  for(uint i = 0; i < 4; ++i) std::copy(T, T+3, newX + 3*i); 
  my_gemm(4, 3, 3, 1, X, 3, R_Q, 3, newX, 3, 1);
}

// Get the vector denoting the direction of the line of intersection 
// between the two planes.  We also get the quaternion representig
// the rotation of the plane containing the points from X to the plane 
// containing the points from Y about the line of intersection.
void
align_planes(const my_float_t* X_norm, const my_float_t* Y_norm, Quaternion* Q)
{
  my_float_t cos_phi = dot(Y_norm, X_norm);
  my_float_t sin_phi = 0.0;
  my_float_t C[] = {1.0, 0.0, 0.0};

  // We don't need exact machine epsilon here, since we have at most three 
  // decimal digits of precision for atom coordinates.
  if(cos_phi > 1 - min_significant_value) cos_phi = 1.0;
  else if(cos_phi < -1 + min_significant_value)cos_phi = -1.0;
  else{
    cross(Y_norm, X_norm, C);
    sin_phi = normalize(C, 3);
  }
  *Q = Quaternion(cos_phi, sin_phi, C);
}

void
align_to_X_axis(const my_float_t *V, Quaternion* Q)
{
  my_float_t len_V = std::sqrt(dot(V,V));
  my_float_t cos_theta = *V / len_V;
  my_float_t sin_theta = 0.0;
  my_float_t Z_axis[] = { 0, 0, 1 };
  
  if(cos_theta > 1 - min_significant_value) cos_theta = 1.0;
  else if(cos_theta < -1 + min_significant_value) cos_theta = -1.0;
  else sin_theta = V[1] / len_V;
  *Q = Quaternion(cos_theta, sin_theta, Z_axis);
}

void
plane_normal(const my_float_t* X, my_float_t* N)
{
  my_float_t U[6];
  for(uint i = 0; i < 3; ++i){
    U[i] = X[i] - X[3 + i];
    U[3 + i] = X[6 + i] - X[3 + i];
  }
  cross(U, U+3, N);
  normalize(N);
}

/*
my_float_t 
dot(const my_float_t* U, const my_float_t* V)
{
  my_float_t rv = 0;
  const my_float_t *v = V;
  const my_float_t *u = U;
  for(uint i = 0; i < 3; ++i, ++u, ++v) rv += *u * (*v);
  return rv;
}
*/

void 
cross(const my_float_t* U, const my_float_t* V, my_float_t* rv)
{
  const my_float_t* u = U + 1;
  const my_float_t* v = V + 2;
  my_float_t* uxv = rv;

  // UxV[0] = U[1]*V[2] - U[2]*V[1]
  *uxv = (*u) * (*v);
  --v;
  ++u;
  *uxv -= (*u) * (*v);
  ++uxv;

  // UxV[1] = U[2]*V[0] - U[0]*V[2]
  --v;
  *uxv = (*u) * (*v);
  v += 2;
  u -= 2;
  *uxv -= (*u) * (*v);
  ++uxv;

  // UxV[2] = U[0]*V[1] - U[1]*V[0]
  --v;
  *uxv = (*u) * (*v);
  ++u;
  --v;
  *uxv -= (*u) * (*v);
}

#if 0
my_float_t 
normalize(my_float_t* X, const uint ndim)
{
  my_float_t* x = X;
  my_float_t mag = 0;
  for(uint i = 0; i < ndim; ++i, ++x) mag += (*x) * (*x);
  mag = std::sqrt(mag);
  x = X;
  for(uint i = 0; i < ndim; ++i, ++x) *x /= mag;
  return mag;
}
#endif

my_float_t
get_root_from_cubic(const my_float_t *P)
{
  // Rewrite as x^3 + bx^2 + cx + d = 0
  my_float_t b = P[1]/P[0];
  my_float_t c = P[2]/P[0];
  my_float_t d = P[3]/P[0];
//  std::cout << "1.0 " << b << " " << c << " " << d << "\n";

  // Write in the form x^3 + mx = n
  my_float_t m = c - b*b/3.0;
  my_float_t n = -1.0*d + b*c/3.0 -2.0*b*b*b/27;
//  std::cout << "m:  " << m << "\n";
//  std::cout << "n:  " << n << "\n";

  my_float_t M = m/3.0;
  my_float_t N = 0.5 * n;
//  std::cout << "M:  " << M << "\n";
//  std::cout << "N:  " << N << "\n";

  // Compute the discriminant D
  my_float_t D = N*N + M*M*M;

  my_float_t reduced_root = 0.0;
  // Compute the root given by S + T where:
  // S = cube root(N + sqrt(D)) & T = cube root(N - sqrt(D)) 
  if(D < 0){
    // There may be a better way to compute this, but for now we write
    // S = cube root(N + sqrt(D)) using complex numbers + De Moivre's formula
    // and taking the ?? branch

    // Note in this case D < 0 ==> -D > 0 thus sqrt(-D) = +/- sqrt(D).
    // Also the tangent of the slope in the complex plane is 
    // +/- sqrt(D) / N.  However, we need to be careful with atan() since
    // it returns a value in [-M_PI_2, M_PI_2].
  
    // Positive sqrt case:
    // theta =? atan( sqrt(D) / N)
    //   If N > 0, theta lies in the first quadrant
    //   If N < 0, theta actually lies in the second quadrant, but atan() 
    //     will return a value in the fourth quadrant (need to add M_PI to the
    //     value returned by atan )
   
    // Negative sqrt case:
    // We can ignore this case since the vector is reflected about the 
    // imaginary axis (or X axis if you prefer) and cos of the angle
    // is the same 
  
    my_float_t theta = std::atan(std::sqrt(-1.0 * D) / N);
    theta = ( theta < 0.0 ? theta + M_PI : theta);

    //std::cout << "Need to jump to complex world\n";
    //return 2.0 * std::pow(-1.0 * M*M*M, 1.0/6.0) * std::cos(theta / 3.0);
    reduced_root = 2.0 * std::pow(N*N - D, 1.0/6.0) * std::cos(theta / 3.0);
  }else{
    // Kludge to use pow even if the arguments would be negative since we
    // want the root with the imaginary component of zero
    my_float_t arg1 = N + std::sqrt(D);
    my_float_t arg2 = N - std::sqrt(D);
    if(arg1 < 0) reduced_root = -1.0 * std::pow(-1.0*arg1, 1.0/3.0);
    else reduced_root = std::pow(arg1, 1.0/3.0);
    if(arg2 < 0) reduced_root += -1.0 * std::pow(-1.0*arg2, 1.0/3.0);
    else reduced_root += std::pow(arg2, 1.0/3.0);
    //reduced_root = std::pow(N + std::sqrt(D), 1.0/3.0) + 
                   //std::pow(N - std::sqrt(D), 1.0/3.0);
  }
  return reduced_root - b/3.0;
}

bool
quadratic_roots(const my_float_t *P, my_float_t *roots)
{
  my_float_t radical = P[1]*P[1] - 4*P[0]*P[2];
  if(radical < -1.0E-14) return false;
  else if(radical < 1.0E-14) radical = 0.0;
  else radical = std::sqrt(radical);

  roots[0] = (-1.0 * P[1] + radical) / (2.0 * P[0]);
  if(radical != 0.0)
    roots[1] = (-1.0 * P[1] - radical) / (2.0 * P[0]);
  else 
    roots[1] = roots[0];
  return true;
}

// Solve for x in x^4 + p*x^2 + q*x + r = 0
void
quartic_roots(const my_float_t p, const my_float_t q, const my_float_t r,
              my_float_t *roots)
{
  std::fill(roots, roots + 4, 0.0);
  my_float_t cubic_P[] = {1.0, -1.0*p, -4.0*r, 4.0*r*p - q*q};
  my_float_t y = get_root_from_cubic(cubic_P);
 
  my_float_t R = std::sqrt(y - p);
  if(-1.0*my_float_eps <= R && R <= my_float_eps){
    std::cerr << "R is very small (" << R << ")\n";
    std::cout << "before return, eigenvals: " << roots[0] << " "
              << roots[1] << " " << roots[2] << " " << roots[3] << "\n";
    return;
  }

  my_float_t D = std::sqrt(-y - p - 2*q/R);
  roots[0] = 0.5 * (R+D);
  roots[1] = 0.5 * (R-D);

  my_float_t E = std::sqrt(-y - p + 2*q/R);
  roots[2] = 0.5 * (-R+E);
  roots[3] = 0.5 * (-R-E);
}

// Assume a bijection between model_pts and target_pts -- this is reasonable
// in the case of pt 2 surf since for each model_pt we will compute the surf
// pt (before this function).  the idea si to just get the model to target
// transformation and just use the inverse (which is easy to find for 3D 
// affine transformations.  
bool
lse_3D_fit(const my_float_t* model_pts, const my_float_t* target_pts, 
           const size_t npts, Quaternion *Q, my_float_t *T, const my_float_t *W)
{
  if(npts < 3){
    std::cerr << "A 3D alignment was requested with fewer than 3 points "
              << "(lse_3D_fit)\n";
    return false;
  }
  if(npts == 3){
    const my_float_t W_tmp[] = {1.0, 1.0, 1.0};
    my_float_t dummy[9];
    if(W) three_pt_align(target_pts, model_pts, W, Q, T, dummy);
    else three_pt_align(target_pts, model_pts, W_tmp, Q, T, dummy);
    return true;
  }
  //std::cout << "number of points in LSE fit: " << npts << "\n";

#if 0
  my_float_t rad = 0.25;

  std::cout << "COLOR, 0.9, 0.2, 0.2,\n";
  for(size_t i = 0; i < npts; ++i)
    if(W[i]) 
      std::cout << "SPHERE, " << model_pts[3*i] << ", "
                << model_pts[3*i + 1] << ", " << model_pts[3*i + 2] << ", "
                << rad << ",\n";
  std::cout << "\n\n";

  std::cout << "COLOR, 0.2, 0.9, 0.2,\n";
  for(size_t i = 0; i < npts; ++i)
    if(W[i]) 
      std::cout << "SPHERE, " << target_pts[3*i] << ", "
                << target_pts[3*i + 1] << ", " << target_pts[3*i + 2] << ", "
                << rad << ",\n";
  std::cout << "\n\n";
  for(size_t i = 0; i < npts; ++i)
    if(W[i]) 
      std::cout << "CYLINDER, " << model_pts[3*i] << ", "
                << model_pts[3*i + 1] << ", " << model_pts[3*i + 2] << ", "
                << target_pts[3*i] << ", "
                << target_pts[3*i + 1] << ", " << target_pts[3*i + 2] << ", "
                << rad << ", 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, \n";
  

#endif

  // Find the centroid of each set of points
  my_float_t m_centroid[3], t_centroid[3]; 
  centroid_3D(model_pts, npts, m_centroid, W);
  centroid_3D(target_pts, npts, t_centroid, W);
#if 0
  std::cout << "Target centroid: " 
	  << t_centroid[0] << " "
	  << t_centroid[1] << " "
	  << t_centroid[2] << "\n";
  std::cout << "Model centroid: " 
	  << m_centroid[0] << " "
	  << m_centroid[1] << " "
	  << m_centroid[2] << "\n";
#endif

  // Compute the matrix of sums
  my_float_t M_sums[9];
  std::fill(M_sums, M_sums + 9, 0.0);
  const my_float_t* m_pt = model_pts;
  const my_float_t* t_pt = target_pts;
  my_float_t W_sum = npts;
  if(W){
    W_sum = 0.0;
    const my_float_t* w = W;
    for(size_t p = 0; p < npts; ++p, m_pt += 3, t_pt += 3, ++w){
      W_sum += *w;
      my_float_t* M_p = M_sums;
      for(int i = 0; i < 3; ++i, M_p += 3){
        for(int j = 0; j < 3; ++j) M_p[j] += (*w) * m_pt[i] * t_pt[j]; 
      }
    }
  }else{
    for(size_t p = 0; p < npts; ++p, m_pt += 3, t_pt += 3){
      my_float_t* M_p = M_sums;
      for(int i = 0; i < 3; ++i, M_p += 3){
        for(int j = 0; j < 3; ++j) M_p[j] += m_pt[i] * t_pt[j]; 
      }
    }
  }
  // Shift by the centroids
  for(size_t i = 0; i < 3; ++i){
    for(size_t j = 0; j < 3; ++j)
      M_sums[3*i + j] = (M_sums[3*i + j] / W_sum) - m_centroid[i]*t_centroid[j];
  }
  
#if 0
  std::cout << "W_sum: " << W_sum << "\n";
  std::cout << "M_sums matrix:\n" << std::setprecision(14);
  for(size_t i = 0; i < 3; ++i){
    for(size_t j = 0; j < 3; ++j)
      std::cout << M_sums[3*i + j] << " ";
    std::cout << "\n";
  }
#endif
    
  // Build the 4x4 symmetric matrix to ??? (looking for the eigenvector 
  // associated with the largest (most positive) eigenvalue)
  my_float_t N[16];
  N[0] = M_sums[0] + M_sums[4] + M_sums[8];
  N[5] = M_sums[0] - M_sums[4] - M_sums[8];
  N[10] = M_sums[4] - M_sums[8] - M_sums[0];
  N[15] = M_sums[8] - M_sums[0] - M_sums[4];
  N[1] = N[4] = M_sums[5] - M_sums[7];
  N[2] = N[8] = M_sums[6] - M_sums[2];
  N[3] = N[12] = M_sums[1] - M_sums[3];
  N[6] = N[9] = M_sums[1] + M_sums[3];
  N[7] = N[13] = M_sums[6] + M_sums[2];
  N[11] = N[14] = M_sums[5] + M_sums[7];
#if 0
  std::cout << "N matrix:\n";
  for(size_t i = 0; i < 4; ++i){
    for(size_t j = 0; j < 4; ++j)
      std::cout << N[4*i + j] << " ";
    std::cout << "\n";
  }
#endif

  // We need to determine some conditions on N such that we would not search
  // for a rotational component.  In particular, if we are very close, then
  // the entries for the first row and first column will be very close to 
  // zero except for N[0] 
  if(N[0] + N[1] == N[0] && N[0] + N[2] == N[0] && N[0] + N[3] == N[0]){
    // Compute translation component
    *Q = Quaternion();
    std::copy(m_centroid, m_centroid + 3, T);
    my_axpy(3, -1.0, t_centroid, 1, T, 1);
    return true;
  }

  my_float_t *n = N;
  for(size_t i = 0; i < 16; ++i, ++n)
    if(-1.0*my_float_eps < *n && *n < my_float_eps) *n = 0.0;
  n = 0;

  // Use the closed form solution of the quartic to determine the rotation
  // http://www-groups.dcs.st-and.ac.uk/~history/HistTopics/Quadratic_etc_equations.html

  // Solve the characteristic equation, det(N - lambda*I) = 0, for the 
  // eigenvalues of N.  We can then write down the equation as
  // l^4 + c3*l^3 + c2*l^2 + c1*l + c0 = 0
  // c3 = 0
  // Let p = c2; c2 is -2 times the trace of M^t M
  // Let q = c1;
  // Let r = c0;
  my_float_t p = 0.0;
  for(size_t i = 0; i < 3; ++i)
    for(size_t j = 0; j < 3; ++j) p += M_sums[3*j + i]*M_sums[3*j + i]; 
  p *= -2.0;

  my_float_t q = 8.0 * (M_sums[0]*M_sums[5]*M_sums[7] +
                        M_sums[4]*M_sums[6]*M_sums[2] + 
                        M_sums[8]*M_sums[1]*M_sums[3] - 
                        M_sums[0]*M_sums[4]*M_sums[8] -
                        M_sums[5]*M_sums[6]*M_sums[1] - 
                        M_sums[7]*M_sums[3]*M_sums[2]);

  my_float_t r = (N[0]*N[5] - N[4]*N[4]) * (N[10]*N[15] - N[14]*N[14])  +
                 (N[4]*N[8] - N[0]*N[9]) * (N[9]*N[15] - N[14]*N[13])   +
                 (N[0]*N[13] - N[4]*N[12]) * (N[9]*N[14] - N[10]*N[13]) +
                 (N[4]*N[9] - N[5]*N[8]) * (N[8]*N[15] - N[14]*N[12])   +
                 (N[5]*N[12] - N[4]*N[13]) * (N[8]*N[14] - N[10]*N[12]) + 
                 (N[8]*N[13] - N[9]*N[12]) * (N[8]*N[13] - N[9]*N[12]);

  my_float_t lambda = -10000.0;

  //std::cout << "(p, q, r): " << p << " " << q << " " << r << "\n";

  // need some tolerance for M_sums to be numerically singular (i.e.
  // det(M_sums) ~ 0.0 ) 
  // I keep forgetting that 99.9 % of the time when the determinant of 
  // M_sums is ~0.0, it is because 3 or fewer points are passed in and the
  // inputs to this function are incorrect
  if(-8 * my_float_eps <= q && q <= 8 * my_float_eps){
    std::cout << "det(8*M_sums): " << q << "\n";
   
    // Solve using x^4 + px^2 + r = 0  
    my_float_t rad_arg = p*p - 4*r;
    std::cout << "radical arg: " << rad_arg << "\n";
    if(-10*my_float_eps <= rad_arg && rad_arg < 0.0) rad_arg = 0.0;
    else if(rad_arg < 0.0){
      std::cerr << "got imaginary roots\n";
      return false;
    }
    lambda = std::sqrt(0.5*(-p + std::sqrt(rad_arg)));

    std::cout << "radical arg: " << rad_arg << "\n";
    std::cout << "my_lambda: " << lambda << "\n";
    
  }else{

    my_float_t roots[4];
    quartic_roots(p, q, r, roots);
#if 0
  std::cout << "eigenvals: " << roots[0] << " "
            << roots[1] << " " << roots[2] << " " << roots[3] << "\n";
#endif

    // Need to get eigenvector corresponding to the largest (most positive) 
    // eigenvalue
    lambda = roots[0];
    for(size_t i = 1; i < 4; ++i) 
      lambda = (lambda < roots[i] ? roots[i] : lambda);
    if(lambda < 0) std::cerr << "Largest eigenvalue is negative!\n";
  }

  // To save time & space we use N to represent (N - lambdaI)
  n = N;
  for(size_t i = 0; i < 4; ++i, n +=5) *n -= lambda; 
  n = 0;
  my_float_t cofactors[16];
  cofactors_4x4(N, cofactors);

#if 0
  std::cout << "Cofactor matrix:\n";
  for(size_t i = 0; i < 4; ++i){
    for(size_t j = 0; j < 4; ++j)
      std::cout << cofactors[4*i + j] << " ";
    std::cout << "\n";
  }
#endif

  ///
  // To reduce chances of problems we average the eigenvectors
  ///

  // Find largest component of first row and its sign
  my_float_t q_vec[4];
  std::copy(cofactors, cofactors + 4, q_vec);
  my_float_t max_e_val = 0;
  size_t idx = 0;
  size_t sign = 0;
  for(size_t i = 0; i < 4; ++i){
    my_float_t tmp = (q_vec[i] > 0 ? q_vec[i] : -1.0*q_vec[i]);
    if(tmp > max_e_val){
      max_e_val = tmp;
      idx = i;
      sign = (q_vec[i] > 0 ? 1 : -1);
    }
  }

  // Add the rows and swap sign if they point in opposing directions
  my_float_t *c = cofactors + 4;
  for(size_t row = 1; row < 4; ++row, c += 4){
    size_t row_sign = (c[idx] > 0 ? 1 : -1);
    if(row_sign == sign) my_axpy(4, 1.0, c, 1, q_vec, 1);
    else my_axpy(4, -1.0, c, 1, q_vec, 1);
  }
  *Q = Quaternion(q_vec, 4);
  
  // Compute translation component -- 
  // rotate target centroid and compute translation to move rotated centroid
  // to the model centroid
  my_float_t R[9];
  Q->get_ortho_rot_mat(R);
  my_float_t rot_t_cent[3];
  my_gemm(1, 3, 3, 1.0, t_centroid, 1, R, 3, rot_t_cent, 1, 0.0);
  std::copy(m_centroid, m_centroid + 3, T);
  my_axpy(3, -1.0, rot_t_cent, 1, T, 1);
#if 0
  std::cout << "R rot cent: " << rot_t_cent[0] << " "
            << rot_t_cent[1] << " " << rot_t_cent[2] << "\n";
  std::cout << "m centroid: " << m_centroid[0] << " "
            << m_centroid[1] << " " << m_centroid[2] << "\n";
#endif
  return true;
}

void
centroid_3D(const my_float_t* pts, const size_t npts, my_float_t* C,
            const my_float_t* W)
{
  if(npts == 0){
    std::fill(C, C + 3, my_float_max);
    return;
  }
  std::fill(C, C + 3, 0);
  const my_float_t* pt = pts;
  const my_float_t* w = W;
  my_float_t W_sum = npts;
  if(W){
    W_sum = 0.0;
    for(size_t p = 0; p < npts; ++p, pt += 3, ++w){
      my_axpy(3, *w, pt, 1, C, 1);
      W_sum += *w;
    }
  }
  else for(size_t p = 0; p < npts; ++p, pt += 3) my_axpy(3, 1.0, pt, 1, C, 1);
  for(size_t i = 0; i < 3; ++i) C[i] /= W_sum;
}

bool
corresponding_point(const my_float_t* pt, const my_float_t *Vi, 
                    const my_float_t *Vj, const my_float_t *Vk, 
		    const my_float_t prev_best_d,
                    my_float_t* d, my_float_t* pt_on_plane)
{
#if 0
  std::cout.precision(6);
  std::cout << "In best point on triangle\n";
  std::cout << "pt: " << pt[0] << " " << pt[1] << " " << pt[2] << "\n";
  std::cout << "Vi: " << Vi[0] << " " << Vi[1] << " " << Vi[2] << "\n";
  std::cout << "Vj: " << Vj[0] << " " << Vj[1] << " " << Vj[2] << "\n";
  std::cout << "Vk: " << Vk[0] << " " << Vk[1] << " " << Vk[2] << "\n";
#endif

  my_float_t Eij[3], Eik[3];
  std::copy(Vj, Vj + 3, Eij);
  std::copy(Vk, Vk + 3, Eik);
  my_axpy(3, -1.0, Vi, 1, Eij, 1);
  my_axpy(3, -1.0, Vi, 1, Eik, 1);

  // unit normal to the face -- should face to the exterior of the object
  my_float_t N[3];
  cross(Eij, Eik, N);
  normalize(N);

  // Vector from Vi to pt
  my_float_t p0p1[3];
  std::copy(pt, pt + 3, p0p1);
  my_axpy(3, -1.0, Vi, 1, p0p1, 1);

  // Signed distance from pt to plane 
  *d = dot(N, p0p1);
  my_float_t abs_prev_best_d = 
    (prev_best_d > 0.0 ? prev_best_d : -1.0 * prev_best_d);
  my_float_t abs_d = (*d > 0.0 ? *d : -1.0 * (*d));
  // This is very very bad
  if(abs_d >= abs_prev_best_d) return true;

  // Projection of point onto plane (by subtracting the distance from the 
  // point in the normal direction)
  std::copy(pt, pt + 3, pt_on_plane);
  my_axpy(3, -1 *(*d), N, 1, pt_on_plane, 1);
#if 0
  std::cout << "point on plane: " << pt_on_plane[0] << " " 
            << pt_on_plane[1] << " " << pt_on_plane[2] << "\n\n";
#endif
  return pt_inside_triangle(N, Vi, Vj, Vk, pt_on_plane);
}

bool
corresponding_point(const my_float_t* pt, const my_float_t* Vi, 
                    const my_float_t* Vj, my_float_t* d, 
                    my_float_t* pt_on_line, 
                    const bool LINE_SEGMENT)
{
#if 0
  std::cout.precision(6);
  std::cout << "In best point on line\n";
  std::cout << "pt: " << pt[0] << " " << pt[1] << " " << pt[2] << "\n";
  std::cout << "Vi: " << Vi[0] << " " << Vi[1] << " " << Vi[2] << "\n";
  std::cout << "Vj: " << Vj[0] << " " << Vj[1] << " " << Vj[2] << "\n";
#endif


  // Get the length of the edge (Vi,Vj) and a unit vector in the direction of 
  // Vi to Vj
  my_float_t L[3];
  my_float_t edge_len = unit_vector(L, Vj, Vi);

  // Vector from Vi to pt
  my_float_t ViP0[3];
  std::copy(pt, pt + 3, ViP0);
  my_axpy(3, -1.0, Vi, 1, ViP0, 1);

  // compute the projection of ViP0 onto L.  Since L is a unit vector, the
  // projection of the point is the number of units in the L direction from 
  // Vi 
  my_float_t proj_pt = dot(ViP0, L);
  bool rv = true;
  
  // CASE 0:  The correspondence is the closest point on the line (not line
  // segment) containing (defined by) Vi and Vj.
  if(!LINE_SEGMENT){
    std::copy(Vi, Vi+3, pt_on_line);
    my_axpy(3, proj_pt, L, 1, pt_on_line, 1);

  // CASE 1:  The projection of the point must be contained in the line 
  // segment.  If the projection is not contained in the line segment, set
  // the corresponding point as the closer of the two end points and return
  // false to denote the correspondance is not contained in the line segment.
  }else{
    if(proj_pt < 0.0){
      std::copy(Vi, Vi+3, pt_on_line);
      rv = false;
    }else if(proj_pt > edge_len){
      std::copy(Vj, Vj+3, pt_on_line);
      rv = false;
    }else{
      std::copy(Vi, Vi+3, pt_on_line);
      my_axpy(3, proj_pt, L, 1, pt_on_line, 1);
    }
  }

#if 0
  std::cout << "point on line: " << pt_on_line[0] << " " 
            << pt_on_line[1] << " " << pt_on_line[2] << "\n\n";
#endif
  *d = dist(pt, pt_on_line);
  return rv;
}

// 3D version
// Determine if point is contained in the closed triangle 
// Use barycentric coordinates by projection onto a coordinate plane by 
// dropping the coordinate with the least variance in the plane
bool
pt_inside_triangle(const my_float_t *N, const my_float_t *Vi,
                   const my_float_t *Vj, const my_float_t *Vk, 
                   const my_float_t *pt)
{
  const my_float_t* tmp_p = N;
  my_float_t max_val = 0;
  size_t idx = 0;
  for(size_t i = 0; i < 3; ++i, ++tmp_p){
    const my_float_t n = *tmp_p;
    if(0 < n && max_val < n){
       max_val = n;
       idx = i;
    }else if(max_val < -1*n){
       max_val = -1*n;
       idx = i;
    }
  }

  my_float_t vals[8];
  my_float_t* v = vals;
  for(size_t i = 0; i < 3; ++i){
    if(idx == i) continue;
  
    *v = *(pt + i);
    ++v;
    *v = *(Vi + i);
    ++v;
    *v = *(Vj + i);
    ++v;
    *v = *(Vk + i);
    ++v;
  }
  return pt_inside_triangle(vals, vals + 4);
}

// 2D version
// From Usenet discussion on Barycentric Coordinates (Joseph O'Rourke, 1992)
// NOTE: for now, this function requires the points be ordered in a 
// counterclockwise fashion
bool
pt_inside_triangle(const my_float_t* x, const my_float_t* y)
{
  my_float_t X[3], Y[3];
  for(size_t i = 0; i < 3; ++i){
    X[i] = x[i+1] - *x;
    Y[i] = y[i+1] - *y;
  }
  my_float_t b[3];
  cross(X, Y, b);

  // Need to determine on which side the of the triangle the point lies
  for(size_t i = 0; i < 3; ++i)
    if(b[i] < 0) return false;
  return true;
}

void 
max_pair_dist(const my_float_t* pts, const size_t npts, size_t *a, size_t *b,
              my_float_t *d)
{
  const my_float_t *p_i = pts;
  const my_float_t *p_j = pts;
  my_float_t max_dist_squared = -1.0;
  size_t i_max = 0;
  size_t j_max = 0;
  for(size_t i = 0; i < npts; ++i, p_i += 3){
    p_j = p_i + 3;
    for(size_t j = i + 1; j < npts; ++j, p_j += 3){
      my_float_t dd = dist_squared(p_i, p_j);
      if(dd > max_dist_squared){
        i_max = i;
        j_max = j;
        max_dist_squared = dd;
      }
    }
  }
  *a = i_max;
  *b = j_max;
  *d = std::sqrt(max_dist_squared);
}


#if 0
// this is old code
/*
0) Get probe point and sphere center 
1) Compute cutting plane to get cap.
   a) Normal in the N-H (O-lp, ...) direction
   b) contains the point in N-H direction that is 1.75 (A) away from N
2) If the pt has positive distance to plane, then project to sphere
   else the pt gets projected to the circle

   to project to circle, first project to plane and then onto circle

   Note:  if at any point in time we are too far from what is left of the
   cap we can bail out with no close point found
*/

my_float_t 
pt2plane_dist(const my_float_t *N, const my_float_t *p0, const my_float_t *pt)
{
  my_float_t p0pt[3]; 
  std::copy(pt, pt+3, p0pt);
  my_axpy(3, -1.0, p0, 1, p0pt, 1);
  return dot(N, p0pt);
}

my_float_t 
pt2plane_proj(const my_float_t *N, const my_float_t *p0, const my_float_t *pt,
              my_float_t* pt_on_plane)
{
  std::copy(pt, pt + 3, pt_on_plane);
  my_float_t d = pt2plane_dist(N, p0, pt);
  my_axpy(3, -1.0 * d, N, 1, pt_on_plane, 1);
  return d;
}

void 
pt2plane_proj(const my_float_t *N, const my_float_t d, const my_float_t *pt,
              my_float_t* pt_on_plane)
{
  std::copy(pt, pt + 3, pt_on_plane);
  my_axpy(3, -1.0 * d, N, 1, pt_on_plane, 1);
}

void
pt2sphere_proj(const sphere &S, const my_float_t *pt, my_float_t *pt_on_sphere)
{
  my_float_t dir[3];
  unit_vector(dir, pt, S.A_c);
  std::copy(S.A_c, S.A_c + 3, pt_on_sphere);
  my_axpy(3, S.A_r, dir, 1, pt_on_sphere, 1);
}


// Plane equation in the form (N .dot. X) + (N .dot. d) = 0 -- in x,y,z form 
// ax + by + cz + d = 0
// Numerical approximation to intersection via Mathworld plane-plane page
// replace this with SVD or other linear solution to [N1,N2]^T x = -[d1;d2]
void
intersect_planes(const my_float_t *N1, const my_float_t d1, 
                 const my_float_t *N2, const my_float_t d2,
                 my_float_t *m, my_float_t *b)
{
std::cout << "Inside intersect planes\n";
std::cout << "N1: " << N1[0] << " " << N1[1] << " " << N1[2] << "\n";
std::cout << "N2: " << N2[0] << " " << N2[1] << " " << N2[2] << "\n";
std::cout << "d1: " << d1 << "\n";
std::cout << "d2: " << d2 << "\n\n";

  cross(N1, N2, m);

  my_float_t max_m_i = 0.0;
  size_t idx = 0;
  for(size_t i = 0; i < 3; ++i){
    my_float_t m_i = (m[i] < 0 ? -1.0 * m[i] : m[i]);
    if(m_i > max_m_i){
      max_m_i = m_i;
      idx = i;
    }
  }

  my_float_t A[4];
  my_float_t *a = A;
  for(size_t i = 0; i < 3; ++i){
    if(i == idx)
      b[i] = 0.0;
    else{ 
      *a = N1[i];
      *(a+2) = N2[i];
      ++a;
    }
  }
  a = 0;

  my_float_t det_A = A[0]*A[3] - A[2]*A[1];
  my_float_t A_inv[] = {A[3] / det_A,      -1.0*A[1] / det_A, 
                        -1.0*A[2] / det_A, A[0] / det_A      };
  a = A_inv;
  for(size_t i = 0; i < 3; ++i, ++b){
    if(i == idx) continue;
    *b = -1.0 * ((*a) * d1 + (*(a+1)) * d2);
    a += 2;
  }
}

bool
line_thru_sphere(const my_float_t *m, const my_float_t *b,
                 const sphere &S, my_float_t *t1, my_float_t *t2)
{
  my_float_t pt_diff[3]; 
  std::copy(b, b+3, pt_diff);
  my_axpy(3, -1.0, S.A_c, 1, pt_diff, 1);

  // quadratic equation Ax^2 + Bx + C = 0
  // A = 1
  my_float_t B = dot(m, pt_diff);
  my_float_t C = dot(pt_diff, pt_diff) - S.A_rsquared; 
  my_float_t rad_arg = B*B - C;
  // we may wish to use some tolerance around zero to set it to zero as a double
  // root
  if(rad_arg < 0) return false;
 
  my_float_t rad_val = std::sqrt(rad_arg);
  *t1 = -1.0 * B + rad_val; 
  *t2 = -1.0 * B - rad_val; 
  return true;
}

// Compute the plane containing the circle of intersection between S and Si
// If the spheres intersect at a point, we treat this as a non intersection
// If the spheres are the same, we currently consider them as a non
// intersection
bool
intersect(const sphere &A, const sphere &B, my_float_t *N, my_float_t *p0,
          my_float_t *radius)
{
  // If the spheres have the same center they are either the same spheres
  // or have no intersection -- in either case we return false
  if(A.A_c[0] == B.A_c[0] && A.A_c[1] == B.A_c[1] && A.A_c[2] == B.A_c[2])
    return false;

  std::cout << "B.A_c: " << B.A_c[0] << " " << B.A_c[1] << " " 
            << B.A_c[2] << "\n";
  std::cout << "A.A_c: " << A.A_c[0] << " " << A.A_c[1] << " " 
            << A.A_c[2] << "\n";
  std::cout << "A.A_r: " << A.A_r << "\n";
  std::cout << "B.A_r: " << B.A_r << "\n";
  my_float_t len_AB = unit_vector(N, B.A_c, A.A_c);
  std::cout << "N: " << N[0] << " " << N[1] << " " << N[2] << "\n";
  if(A.A_r + B.A_r <= len_AB) return false;

  // Note this distance will be negative if alpha > pi/2
  my_float_t d_A = (len_AB*len_AB + A.A_rsquared - B.A_rsquared) / (2*len_AB);
  std::copy(A.A_c, A.A_c + 3, p0);
  my_axpy(3, d_A, N, 1, p0, 1);
  *radius = std::sqrt(A.A_rsquared - d_A*d_A);
  if(d_A < 0) for(my_float_t *n_ptr = N; n_ptr < N + 3; ++n_ptr) *n_ptr *= -1.0;
  return true;
}

void 
pair_dist_squared(const my_float_t* pts, const size_t npts, my_float_t *d_sq)
{
  std::fill(d_sq, d_sq + npts*npts, -1.0);

  const my_float_t *p_i = pts;
  const my_float_t *p_j = pts;
  my_float_t *dd = d_sq;
  for(size_t i = 1; i < npts; ++i, p_i += 3){
    p_j = p_i + 3;
    dd += i;
    for(size_t j = i; j < npts; ++j, p_j += 3, ++dd){
      *dd = dist_squared(p_i, p_j);
    }
  }
}

bool
my_close_point(const sphere &S, const std::vector<sphere> &nbrs, 
               const my_float_t *pt, const my_float_t tolerance,
               my_float_t *close_pt)
{
  const my_float_t tol_squared = tolerance*tolerance;

  // 0) Project to sphere and return if distance is greater than tolerance
  S.proj2surf(pt, close_pt);
  std::cout << "initial projection: " << close_pt[0] << " " 
            << close_pt[1] << " " << close_pt[2] << "\n";
  if(dist_squared(pt, close_pt) > tol_squared) return false;


size_t nbr_cnt = 0;

  // 1) Iterate over each nbr sphere
  //   a) The close_pt cannot be inside any of the nbrs
  //   b) if at any time the closest point is no longer within the tolerance,
  //      return.
  std::vector<sphere>::const_iterator nbr;
  for(nbr = nbrs.begin(); nbr != nbrs.end(); ++nbr, ++nbr_cnt){
    // This nbr can be safely ignored if the closest point does not fall
    // in the current open sphere
    if(dist_squared(nbr->A_c, close_pt) < nbr->A_rsquared) break;
  }
  // If the projected point does not fall inside any neighboring sphere, the
  // projected point is the closest point
  if(nbr == nbrs.end()) return true;

std::cout << "first nbr: " << nbr_cnt << "\n";
  // Otherwise, the closest point must lie on the circle of intersection 
  // between the sphere S and the current nbr

  /////
  // project pt to plane
  /////

  // get plane normal in the S center to nbr center direction
  my_float_t N_0[3];   // "Out-facing" normal to circle of S ^ S0
  my_float_t p0[3];    // Center of circle of S ^ S0
  my_float_t radius;
  intersect(S, *nbr, N_0, p0, &radius);
  std::cout << "N_0: " << N_0[0] << " " << N_0[1] << " " << N_0[2] << "\n" ;
  std::cout << "p0 : " << p0[0] << " " << p0[1] << " " << p0[2] << "\n";
  std::cout << "radius : " << radius << "\n";
  sphere circle(p0, radius);

  // project the point to the plane:  N_0 .dot. (X - p0) = 0
  my_float_t pt_in_plane[3];
  pt2plane_proj(N_0, p0, close_pt, pt_in_plane);
  std::cout << "initial pt in plane: " << pt_in_plane[0] << " " 
            << pt_in_plane[1] << " " << pt_in_plane[2] << "\n";

  // Should compute the closest point on circle of intersection and
  // if it is too far, then return no correspondance

size_t sphere_num = 0;

  /////
  // Find all spheres that intersect (on S) with the current nbr
  /////
  std::list<circular_segment> circ_seg_list;
  std::vector<sphere>::const_iterator Si;
  for(Si = nbrs.begin(); Si != nbrs.end(); ++Si, ++sphere_num){
    // Ignore nbr itself and all spheres that do not intersect nbr
    if(Si == nbr) continue;

    if(nbr->A_r + Si->A_r <= dist(nbr->A_c, Si->A_c)) continue;
std::cout << "nbr intersects with sphere " << sphere_num << "\n";

    // Assume all spheres in the nbrs list actually do intersect nontrivially 
    // with S
    my_float_t N_i[3], p_i[3];
    intersect(S, *Si, N_i, p_i, &radius);
    std::cout << "S ^ Si radius: " << radius << "\n";

    // Compute circular segment clipped from nbr's circle on S by Si's
    // circle on S

    // Determine the intersection of the two planes containing the two circles
    // If the line of intersection does not go through the first circle, 
    //  continue
    // x = tm + b
    std::cout << "p_0: " << p0[0] << " " << p0[1] << " " 
              << p0[2] << "\n";
    std::cout << "p_i: " << p_i[0] << " " << p_i[1] << " " 
              << p_i[2] << "\n";
    my_float_t m[3], b[3];
    intersect_planes(N_0, -1.0 * dot(N_0, p0), N_i, -1.0 * dot(N_i, p_i),
                     m, b);
    std::cout << "m: " << m[0] << " " << m[1] << " " << m[2] << "\n";
    std::cout << "b: " << b[0] << " " << b[1] << " " << b[2] << "\n";
    my_float_t t0, t1;
    if(line_thru_sphere(m, b, circle, &t0, &t1) == false) continue;
std::cout << "## spheres intersect on cap\n";
    std::cout << "t0: " << t0 << "\n";
    std::cout << "t1: " << t1 << "\n";

    // Compute the 3D coordinates of the intersection points
    my_float_t x0[3], x1[3];
    std::copy(b, b + 3, x0);
    std::copy(b, b + 3, x1);
    my_axpy(3, t0, m, 1, x0, 1);
    my_axpy(3, t1, m, 1, x1, 1);
    std::cout << "x0: " << x0[0] << " " << x0[1] << " " << x0[2] << "\n";
    std::cout << "x1: " << x1[0] << " " << x1[1] << " " << x1[2] << "\n";

#if 0
    // Get the "in" direction by projecting the out facing normal of the
    // circle defined by the intersection of S_i and S to the plane containing
    // the circle defined by the intersection of nbr and S.
    // That is, compute N_i - (N_i .dot. N_0)N_0
    my_float_t in_dir[3], tmp[3];
    std::fill(tmp, tmp + 3, 0.0);
    my_axpy(3, -1.0 * dot(N_i, N_0), N_0, 1, tmp, 1);
    unit_vector(in_dir, N_i, tmp);
#endif

    my_float_t in_dir[3];
    for(size_t i = 0; i < 3; ++i) in_dir[i] = -1.0*N_i[i];
    //unit_vector(in_dir, p_i, p0); 
    //unit_vector(in_dir, p0, p_i); 
    circular_segment current_seg(x0, x1, in_dir, p0);
    //std::cout << "in_dir: " << N_i[0] << " " << N_i[1] << " " 
              //<< N_i[2] << "\n";

    // check to merge all overlapping arcs in the list (deleting if merged)
    std::list<circular_segment>::iterator csi;
    for(csi = circ_seg_list.begin(); csi != circ_seg_list.end(); ){
      if(current_seg.merge(&(*csi))) csi = circ_seg_list.erase(csi);
      else ++csi;
    }
  
    //append to the final arc to the list
    circ_seg_list.push_back(current_seg); 
  }

  std::cout << "Number of circular segments: " << circ_seg_list.size() << "\n";

  // Project pt in plane to closest pt
  circle.proj2surf(pt_in_plane, close_pt);

  std::cout << "pt_in_plane: " << pt_in_plane[0] << " " 
            << pt_in_plane[1] << " " << pt_in_plane[2] << "\n";
  std::cout << "close_pt: " << close_pt[0] << " " 
            << close_pt[1] << " " << close_pt[2] << "\n";
  
  std::list<circular_segment>::iterator csi;
  for(csi = circ_seg_list.begin(); csi != circ_seg_list.end(); ++csi)
    if(csi->contains(close_pt)){
      std::cout << "time to copy a chord end point\n";
      const my_float_t d1 = dist_squared(close_pt, csi->chord_endpt_one());
      const my_float_t d2 = dist_squared(close_pt, csi->chord_endpt_two());
      if(d1 < d2) 
        std::copy(csi->chord_endpt_one(), csi->chord_endpt_one() + 3, close_pt);
      else
        std::copy(csi->chord_endpt_two(), csi->chord_endpt_two() + 3, close_pt);
    }

 
  if(dist_squared(pt, close_pt) > tol_squared) return false;
  return true;
}

/*
    //pt2sphere_proj(circle, pt_in_plane, close_pt);
    //my_float_t v0[3];
    //unit_vector(v0, close_pt, p0);

       // !!!!!!!!!!!!!!!!
       // handle removed arcs !!!!!!!!!!!!!!!!

       // Get the arc end points and if two arcs overlap, merge them till
       // we have nonoverlapping arcs.  Idea:  if two arcs overlap the chords
       // will intersect inside the circle (determine point of intersection
       // by dropping least significant coordinate w.r.t. plane)

       // !!!!!!!!!!!!!!!!

    size_t skip_idx = 0;
    my_float_t max_val = 0.0;
    for(size_t i = 0; i < 3; ++i){
      my_float_t val = (N[i] > 0 ? N[i] : -1.0 * N[i]);
      if(val > max_val){
        max_val = val;
        skip_idx = i;
      }
    }


    std::list<circle::chord> chords;
    // cannot think of a good name for this variable
    std::vector<sphere>::const_iterator s_iter;
    for(s_iter = nbrs.begin(); s_iter != nbrs.end(); ++s_iter){
      if(s_iter == nbr) continue;
      if(nbr->A_r + s_iter->A_r <= dist(nbr->A_c, s_iter->A_c)) continue;

      my_float_t normal[3];
      //my_float_t my_dist = unit_vector(normal, s_iter->A_c, S.A_c);
      my_dist = unit_vector(normal, s_iter->A_c, S.A_c);
      //my_float_t d_S = (my_dist*my_dist + S.A_rsquared - s_iter->A_rsquared) /
      //                 2*my_dist;
      d_S = (my_dist*my_dist + S.A_rsquared - s_iter->A_rsquared) / 2*my_dist;
      my_float_t point[3];
      std::copy(S.A_c, S.A_c + 3, point);
      my_axpy(3, d_S, normal, 1, point, 1);

      // x = tm + b
      my_float_t m[3], b[3];
      intersect_planes(N, -1.0 * dot(N, p0), normal, -1.0 * dot(normal, point),
                       m, b);
      my_float_t t1, t2;
      if(line_thru_sphere(m, b, circle, &t1, &t2) == false) continue;
 
      circle::chord new_chord;
      for(size_t i = 0; i < 3; ++i){
        new_chord.x0[i] = m[i]*t1 + b[i];
        new_chord.x1[i] = m[i]*t2 + b[i];
      }
      //chords.push_front(tmp_chord);
    
      std::list<circle::chord>::iterator prev_chord = chords.begin();
      for( ; prev_chord != chords.end(); ++prev_chord){

        // drop the coordinate corresponding to skip_idx and compute
        // the intersection of the chords cur_chord and c2.  If the intersection
        // lies inside the circle  (p0, radius), then merge the two chords
        // into the current and delete c2 from the list 
        // (update the iterator c2 as well).
        if(skip_idx == 2){
              
        }

      }
      chords.push_front(new_chord);
      
*/      
/*
      // Project center of nbr & pt_in_nbr_plane to line (x = tm + b)
      my_float_t x1[3], x2[3];
      for(size_t i = 0; i < 3; ++i){
        x1[i] = m[i]*t1 + b[i];
        x2[i] = m[i]*t2 + b[i];
      }

      my_float_t c_on_line[3], pt_on_line[3];
      my_float_t c_dist, pt_dist; 
      corresponding_point(p0, x1, x2, &c_dist, c_on_line, false);
      corresponding_point(close_pt, x1, x2, &pt_dist, pt_on_line, false);
*/

     // use dot and cross to get angles for t1,t2, w.r.t. the projected point
      //my_float_t x1[3], x2[3];
      //for(size_t i = 0; i < 3; ++i){
      //  x1[i] = m[i]*t1 + b[i];
      //  x2[i] = m[i]*t2 + b[i];
     // }

      //my_float_t v1[3], v2[3];
      //unit_vector(v1, x1, p0);
      //unit_vector(v2, x2, p0);

      //my_float_t cos_theta1 = dot(v0,v1);
      //my_float_t cos_theta2 = dot(v0,v2);

      //my_float_t xprod[3];
      //cross(v0, v1, xprod);
      //my_float_t sin_theta1 = 
      //cross(v0, v2, xprod);
      //my_float_t sin_theta2 = 
     

     // pair up the angles with vectors to determine the occluded arcs
     
  
/* different method rather than plane intersect plane
      // Get the circle of s_iter in the plane of the current nbr
      my_float_t c_in_plane[3];
      my_float_t c2plane_dist = pt2plane_proj(N, p0, s_iter->A_c, c_in_plane);
      my_float_t rad_in_plane = std::sqrt(s_iter->A_rsqared - 
                                          c2plane_dist*c2plane_dist);
      if(rad_in_plane + radius <= dist(p0 c_in_plane)) continue;
*/

      // Need to determine which part of arc lies in the intersection 
      
    //} 

    // if the projected point falls in an occluded arc, push it to the closest
    // end point of the arc

      ////////////
      // Need to determine the arc of nbr in the intersection of nbr's and 
      // s_iter's circles
      ////////////

    



/*
bool
my_close_point(const my_float_t *atom_pos, const my_float_t *H_pos, 
               const my_float_t *pt, const std::vector<atom_vci> &nbrs,
               my_float_t *close_pt)
{
  // 1) Compute cutting plane to get cap.
  // a) Normal in the N-H (O-lp, ...) direction
  my_float_t N0[3];
  unit_vector(N0, H_pos, atom_pos);

  // b) contains the point in N-H direction that is 1.75 (A) away from N
  my_float_t p0[3];
  std::copy(atom_pos, atom_pos + 3, p0);
  my_axpy(3, 1.75, N0, 1, p0, 1);

  // 2a) Get signed distance from point to plane
  my_float_t d = pt2plane_dist(N0, p0, pt);

  // 2b) if the point lies "below" the plane, we need to project it to the
  // circle in the plane.  redo this as a sphere clipping so that we can 
  // feed it into the loop below
  my_float_t c0[3];
  std::copy(atom_pos, atom_pos + 3, c0);
  my_axpy(3, 0.5, N0, 1, p0, 1);
  if(d >= 0) pt2sphere_proj(c0, 2.5, pt, close_pt);
  else{
    // project point to plane
    my_float_t pt_in_plane[3];
    pt2plane_proj(N0, d, pt, pt_in_plane);  

    // project point in plane onto circle
    // The center of the circle is p0
    // The radius in our case is ??
    //pt2sphere_proj(c0, ??, pt_in_plane, close_pt);
  }

  // for each nbr, we need to check if the current point is on the 
  // "wrong" side of the plane, if so we need to reproject it.
    // Check if this sphere intersects with any of the previous spheres
    // If it does, then we need to determine where they intersect
    // Intersect planes  
}
*/
/*
bool
circle::intersects(const circle &C, my_float_t *V, my_float_t *P0,
                   my_float_t *t0, my_float_t *t1)
{
  intersect_planes(A_N, -1.0 * dot(A_N, A_C), C.A_N, -1.0 * dot(C.A_N, C.A_C),
                   V, P0);
  //sphere this_circle(A_C, A_r);
  if(line_thru_sphere(V, P0, sphere(A_C, A_r), t0, t1) == false) return false;

}
*/
#endif
