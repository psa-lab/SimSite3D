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

#ifndef QUATERION_HEADER_INCLUDED
#define QUATERION_HEADER_INCLUDED

#include <types.H>
#include <algorithm>
#include <cmath>
#include <iostream>

class Quaternion;

inline Quaternion operator*(Quaternion& R, const Quaternion& S);

class Quaternion{
public:
  //! Default constructor : identity quaternion
  Quaternion() 
  { 
    *q = 1;
    std::fill(q + 1, q + 4, 0);
    qM = NULL;
    is_unit = true;
  }

  //! Array to Quaternion constructor
  Quaternion(const my_float_t* in, const size_t len) 
  { 
    if(len == 4) std::copy(in, in + 4, q); 
    // Assume R is orthonormal
    else if(len == 9){
      const my_float_t* R = in;
      my_float_t my_q[4];
      my_q[0] = 0.5 * std::sqrt(1 + R[0] + R[4] + R[8]);
      my_q[1] = 0.5 * std::sqrt(1 + R[0] - R[4] - R[8]);
      my_q[2] = 0.5 * std::sqrt(1 - R[0] + R[4] - R[8]);
      my_q[3] = 0.5 * std::sqrt(1 - R[0] - R[4] + R[8]);
      
      size_t idx = 0;
      my_float_t max_val = my_q[0];
      for(size_t i = 1; i < 4; ++i)
        if(my_q[i] > max_val){
          max_val = my_q[i];
          idx = i;
        }

      q[idx] = my_q[idx];
      if(idx == 0){
        q[1] = 0.25 * (R[7] - R[5])/q[0];
        q[2] = 0.25 * (R[2] - R[6])/q[0];
        q[3] = 0.25 * (R[3] - R[1])/q[0];
      }else if(idx == 1){
        q[0] = 0.25 * (R[7] - R[5])/q[1];
        q[2] = 0.25 * (R[3] + R[1])/q[1];
        q[3] = 0.25 * (R[2] + R[6])/q[1];
      }else if(idx == 2){
        q[0] = 0.25 * (R[2] - R[6])/q[2];
        q[1] = 0.25 * (R[3] + R[1])/q[2];
        q[3] = 0.25 * (R[7] + R[5])/q[2];
      }else if(idx == 3){
        q[0] = 0.25 * (R[3] - R[1])/q[3];
        q[1] = 0.25 * (R[2] + R[6])/q[3];
        q[2] = 0.25 * (R[7] + R[5])/q[3];
      }

    }else{
      std::fill(q, q+4, 0.0);
    }
    qM = NULL;
    is_unit = false;
  }

  //! Rodrigues formula for quaternion -- but parameters are for the full angles
  //! and not half angles.  Use half angle identities to get the desired values.
  Quaternion(const my_float_t cos_theta, const my_float_t sin_theta, 
             const my_float_t* V)
  {
    my_float_t cos_half_theta = 1.0;
    my_float_t sin_half_theta = 0.0;
    if(cos_theta == 1.0){} // do nothing 
    else if(cos_theta == -1.0){
      cos_half_theta = 0.0;
      sin_half_theta = 1.0;
    }else{
      my_float_t cos_theta_plus_one = cos_theta + 1.0;
      cos_half_theta = std::sqrt(0.5 * cos_theta_plus_one);
      sin_half_theta = sin_theta / std::sqrt(2 * cos_theta_plus_one);
    }
    *q = cos_half_theta;
    for(uint i = 1; i < 4; ++i) q[i] = sin_half_theta * V[i - 1];
    qM = NULL;
    is_unit = false;
  }

  //! Quaternions only have 3 degrees of freedom and can be represented
  //! by 2 angles and a scalar in [0.0, 1.0]
  Quaternion(my_float_t theta0, my_float_t theta1, my_float_t s)
  {
    my_float_t sigma0 = std::sqrt(1 - s);
    my_float_t sigma1 = std::sqrt(s);
    q[0] = std::cos(theta1) * sigma1;
    q[1] = std::sin(theta0) * sigma0;
    q[2] = std::cos(theta0) * sigma0;
    q[3] = std::sin(theta1) * sigma1;
    qM = NULL;
    is_unit = false;  // really should be unit, but just check
  }

  //! Copy constructor -- only copies the quaternion
  Quaternion(const Quaternion& src)
  {
    do_copy(src);
  }

  const Quaternion& 
  operator=(const Quaternion& src)
  {
    if(this != &src) do_copy(src);
    return *this;
  }


  ~Quaternion() 
  { 
    if(qM) delete[] qM;
    qM = NULL;  
  }
 
  //! Get the equivalent 3x3 rotation matrix
  void
  get_ortho_rot_mat(my_float_t* R)
  {
    if(!is_unit) normalize();

    my_float_t q_2[4];
    const my_float_t* q_p = q;
    my_float_t* q_2_p = q_2;
    for(uint i = 0; i < 4; ++i, ++q_p, ++q_2_p) *q_2_p = *q_p*(*q_p);
    my_float_t q0qi[3];
    for(uint i = 0; i < 3; ++i) q0qi[i] = *q*q[i+1];
    my_float_t qxqy = q[1]*q[2];
    my_float_t qxqz = q[1]*q[3];
    my_float_t qyqz = q[2]*q[3];

    my_float_t *r = R;
    *r = q_2[0] + q_2[1] - q_2[2] - q_2[3];
    ++r;
    *r = 2*(qxqy - q0qi[2]);
    ++r;
    *r = 2*(qxqz + q0qi[1]);
    ++r;
    *r = 2*(qxqy + q0qi[2]);
    ++r;
    *r = q_2[0] - q_2[1] + q_2[2] - q_2[3];
    ++r;
    *r = 2*(qyqz - *q0qi);
    ++r;
    *r = 2*(qxqz - q0qi[1]);
    ++r;
    *r = 2*(qyqz + *q0qi);
    ++r;
    *r = q_2[0] - q_2[1] - q_2[2] + q_2[3];
  }

  //! Change to unit quaternion
  inline my_float_t 
  normalize()
  {
    if(is_unit) return 1.0;

    // Change the first component to be positive
    my_float_t* qp = q;
    if(*qp < 0.0)
      for(size_t i = 0; i < 4; ++i, ++qp) (*qp) *= -1.0;
    qp = q;

    // Make the quaternion a unit quaternion
    my_float_t mag = 0;
    for(uint i = 0; i < 4; ++i, ++qp) mag += *qp * (*qp);
    mag = std::sqrt(mag);
    if(mag != 1.0){
      qp = q;
      for(uint i = 0; i < 4; ++i, ++qp) *qp /= mag;
    }
    is_unit = true;
    return mag;
  }

  //! Complex conjugate
  Quaternion
  conjugate()
  {
    Quaternion rv = *this;
    my_float_t* qp = rv.q + 1;
    for(uint i = 1; i < 4; ++i, ++qp) *qp *= -1.0;
    return rv;
  }

  //! Get a string holding the quaternion written in complex coordinates
  std::string str();

  //! Get another parametrization of the quaternion
  void
  get_parameters(my_float_t *theta0, my_float_t* theta1, my_float_t *s)
  {
    *theta0 = atan2(q[1], q[2]);
    *theta1 = atan2(q[3], q[0]);
    *s = q[3]*q[3] + q[0]*q[0];
  }

  //! Generate a random quaternion as per the method of James J. Kuffner, 
  //! ICRA 2004.  
  /*! This method is computationally expensive since it calls sqrt, cos and sin.
   * If this method is used for something other than testing, it may be useful
   * to speed it up.
   */
  void
  randomize()
  {
    my_float_t s = std::rand() / static_cast<my_float_t>(RAND_MAX);
    my_float_t sigma0 = std::sqrt(1 - s);
    my_float_t sigma1 = std::sqrt(s);
    my_float_t theta0 = 2*M_PI * 
      (std::rand() / static_cast<my_float_t>(RAND_MAX));
    my_float_t theta1 = 2*M_PI * 
      (std::rand() / static_cast<my_float_t>(RAND_MAX));
    q[0] = std::cos(theta1) * sigma1;
    q[1] = std::sin(theta0) * sigma0;
    q[2] = std::cos(theta0) * sigma0;
    q[3] = std::sin(theta1) * sigma1;
    is_unit = true;
  }

  //! Generate a random quaternion with limited rotation.
  /*! Generate a random quaterion representing a rotation about an arbitrary
   * axis, but with an angle in [theta_beg, theta_end].  It is assumed that
   * 0 <= theta_beg < theta_end <= M_PI. 
   *
   * @param theta_beg smallest theta to sample (in radians)
   * @param theta_end largest theta to sample (in radians)
   */
  void 
  randomize(my_float_t theta_beg, my_float_t theta_end)
  {
    const my_float_t half_width = 0.5 * RAND_MAX;
    my_float_t v[3];
    for(size_t i = 0; i < 3; ++i) 
      v[i] = ((std::rand() - half_width) / RAND_MAX);
    
    my_float_t tmp = (std::rand() - half_width) / half_width;
    my_float_t half_theta = 0;
    if(tmp > 0)
      half_theta = theta_beg + tmp*(theta_end - theta_beg);
    else 
      half_theta = -1*theta_end - tmp*(theta_end - theta_beg);    
    
    half_theta /= 2.0;
    my_float_t cos_half_theta = std::cos(half_theta);
    my_float_t sin_half_theta = std::sin(half_theta);
    
    q[0] = cos_half_theta;
    q[1] = v[0] * sin_half_theta;
    q[2] = v[1] * sin_half_theta;
    q[3] = v[2] * sin_half_theta; 
    is_unit = false;
  }

  friend Quaternion operator*(Quaternion& r, const Quaternion& s);

  friend my_float_t distance(Quaternion& r, Quaternion& s, 
                             const my_float_t W = 1.0);
private:

  //! Fill the qM matrix
  void
  fill_qM()
  {
    if(qM) delete[] qM;
    qM = new my_float_t[16];
    *qM = *(qM + 5) = *(qM + 10) = *(qM + 15) = *q;
    for(uint i = 1; i < 4; ++i){
      qM[i] = -q[i];
      qM[4*i] = q[i];
    } 
    qM[6] = -q[3];
    qM[7] = q[2];
    qM[9] = q[3];
    qM[11] = -q[1];
    qM[13] = -q[2];
    qM[14] = q[1];
  }

  void 
  do_copy(const Quaternion& src)
  {
    std::copy(src.q, src.q + 4, q);
    is_unit = src.is_unit;
    qM = NULL;
  }

  my_float_t q[4];
  my_float_t* qM;  //!< Matrix form of the quaternion q to get q * r by matrix multiplication
  bool is_unit;
};

//! Multiply two quaternions using the "matrix" method
/*! Something is broken -- will not chain (i.e. Q1 * Q2 * Q3 or even (Q1 * (Q2*Q3))).
 */
inline Quaternion
operator*(Quaternion& R, const Quaternion& S)
{
  Quaternion RS;
  // Alwasy assume is not a unit vector just to force normalization after
  // any positve number of multiplications
  RS.is_unit = false;
 
  if(R.qM == NULL) R.fill_qM();
  // Stupid ciruclar dependance here which should go away soon
  //my_gemv(4, 4, 1.0, R.qM, 4, S.q, 1, 0, RS.q, 1);
  std::fill(RS.q, RS.q + 4, 0);
  for(uint i = 0; i < 4; ++i)
    for(uint j = 0; j < 4; ++j)
      RS.q[i] += R.qM[4*i + j] * S.q[j];
  return RS; 
}

//! Get the "distance" between two quaternions as per J. J. Kuffner's method
inline my_float_t
distance(Quaternion& r, Quaternion& s, const my_float_t W)
{
  r.normalize();
  s.normalize();
  my_float_t rv = 0.0;
  for(size_t i = 0; i < 4; ++i) rv += r.q[i] * s.q[i]; 
  if(rv < 0.0) rv *= -1.0;  
  return W * (1.0 - rv);
}

#endif
