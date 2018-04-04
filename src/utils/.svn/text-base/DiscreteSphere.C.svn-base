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
 * $Source: /psa/share/repository/pfizer_proj/src/basics/DiscreteSphere.C,v $
 * $Revision: 1.2 $
 * $Author: vanvoor4 $
 * $Date: 2008-04-17 18:36:06 $
 * 
 * $Log: not supported by cvs2svn $
 * Revision 1.1  2007/12/17 21:09:31  vanvoor4
 * Initial checkin
 *
 *
 */


#include <algorithm>
#include <mat_ops.H>
#include <DiscreteSphere.H>

using namespace SimSite3D;

const my_float_t DiscreteSphere::level_0[] = { 0.000000, 1.000000, 0.000000,
                                               0.866025, 0.500000, 0.000000,
                                               0.000000, 0.500000, 0.866025,
                                              -0.866025, 0.500000, 0.000000,
                                              -0.000000, 0.500000, -0.866025,
                                               0.000000, -1.000000, 0.000000,
                                               0.866025, -0.500000, 0.000000,
                                               0.000000, -0.500000, -0.866025,
                                              -0.866025, -0.500000, -0.000000,
                                              -0.000000, -0.500000, 0.866025};
const uint DiscreteSphere::level_0_len = 10;
const my_float_t DiscreteSphere::level_1[] = { 0.000000, 1.000000, 0.000000,
                                               0.707107, 0.707107, 0.000000,
                                               0.000000, 0.707107, 0.707107,
                                              -0.707107, 0.707107, 0.000000,
                                              -0.000000, 0.707107, -0.707107,
                                               0.707107, 0.000000, 0.707107,
                                               0.000000, 0.000000, 1.000000,
                                              -0.707107, 0.000000, 0.707107,
                                               0.000000, -1.000000, 0.000000,
                                               0.707107, -0.707107, 0.000000,
                                               0.000000, -0.707107, -0.707107,
                                              -0.707107, -0.707107, -0.000000,
                                              -0.000000, -0.707107, 0.707107,
                                               0.707107, 0.000000, -0.707107,
                                               0.000000, 0.000000, -1.000000,
                                              -0.707107, 0.000000, -0.707107,
                                               1.000000, 0.000000, 0.000000,
                                              -1.000000, 0.000000, 0.000000};
const uint DiscreteSphere::level_1_len = 18;
const my_float_t DiscreteSphere::level_2[] = { 0.000000, 1.000000, 0.000000,
                                               0.587785, 0.809017, 0.000000,
                                               0.951057, 0.309017, 0.000000,
                                               0.000000, 0.809017, 0.587785,
                                              -0.587785, 0.809017, 0.000000,
                                              -0.000000, 0.809017, -0.587785,
                                               0.672499, 0.309017, 0.672499,
                                               0.000000, 0.309017, 0.951057,
                                              -0.672499, 0.309017, 0.672499,
                                              -0.951057, 0.309017, 0.000000,
                                              -0.672499, 0.309017, -0.672499,
                                              -0.000000, 0.309017, -0.951057,
                                               0.672499, 0.309017, -0.672499,
                                               0.000000, -1.000000, 0.000000,
                                               0.587785, -0.809017, 0.000000,
                                               0.951057, -0.309017, 0.000000,
                                               0.000000, -0.809017, -0.587785,
                                              -0.587785, -0.809017, -0.000000,
                                              -0.000000, -0.809017, 0.587785,
                                               0.672499, -0.309017, -0.672499,
                                               0.000000, -0.309017, -0.951057,
                                              -0.672499, -0.309017, -0.672499,
                                              -0.951057, -0.309017, -0.000000,
                                              -0.672499, -0.309017, 0.672499,
                                              -0.000000, -0.309017, 0.951057,
                                               0.672499, -0.309017, 0.672499};
const uint DiscreteSphere::level_2_len = 26;
const my_float_t DiscreteSphere::level_3[] = { 0.000000, 1.000000, 0.000000,
                                               0.500000, 0.866025, 0.000000,
                                               0.866025, 0.500000, 0.000000,
                                               0.000000, 0.866025, 0.500000,
                                              -0.500000, 0.866025, 0.000000,
                                              -0.000000, 0.866025, -0.500000,
                                               0.612372, 0.500000, 0.612372,
                                               0.000000, 0.500000, 0.866025,
                                              -0.612372, 0.500000, 0.612372,
                                              -0.866025, 0.500000, 0.000000,
                                              -0.612372, 0.500000, -0.612372,
                                              -0.000000, 0.500000, -0.866025,
                                               0.612372, 0.500000, -0.612372,
                                               0.866025, 0.000000, 0.500000,
                                               0.500000, 0.000000, 0.866025,
                                               0.000000, 0.000000, 1.000000,
                                              -0.500000, 0.000000, 0.866025,
                                              -0.866025, 0.000000, 0.500000,
                                               0.000000, -1.000000, 0.000000,
                                               0.500000, -0.866025, 0.000000,
                                               0.866025, -0.500000, 0.000000,
                                               0.000000, -0.866025, -0.500000,
                                              -0.500000, -0.866025, -0.000000,
                                              -0.000000, -0.866025, 0.500000,
                                               0.612372, -0.500000, -0.612372,
                                               0.000000, -0.500000, -0.866025,
                                              -0.612372, -0.500000, -0.612372,
                                              -0.866025, -0.500000, -0.000000,
                                              -0.612372, -0.500000, 0.612372,
                                              -0.000000, -0.500000, 0.866025,
                                               0.612372, -0.500000, 0.612372,
                                               0.866025, 0.000000, -0.500000,
                                               0.500000, 0.000000, -0.866025,
                                               0.000000, 0.000000, -1.000000,
                                              -0.500000, 0.000000, -0.866025,
                                              -0.866025, 0.000000, -0.500000,
                                               1.000000, 0.000000, 0.000000,
                                              -1.000000, 0.000000, 0.000000};
const uint DiscreteSphere::level_3_len = 39;


DiscreteSphere::DiscreteSphere(sphere_sample_level_t lvl)
{
  switch(lvl){
  case DISCRETE_SPHERE_LEVEL_ZERO:
    A_orig_pts = level_0;
    A_num_pts = level_0_len;
    break;
  case DISCRETE_SPHERE_LEVEL_ONE:
    A_orig_pts = level_1;
    A_num_pts = level_1_len;
    break;
  case DISCRETE_SPHERE_LEVEL_TWO:
    A_orig_pts = level_2;
    A_num_pts = level_2_len;
    break;
  case DISCRETE_SPHERE_LEVEL_THREE:
    A_orig_pts = level_3;
    A_num_pts = level_3_len;
    break;
  }
  
  A_pts = 0;
  std::fill(A_center, A_center + 3, 0.0);
}

DiscreteSphere::DiscreteSphere(const DiscreteSphere& src, const my_float_t* R,  
                               const my_float_t* T)
{
  if(&src == this) return;

  A_num_pts = src.A_num_pts;
  // shouldn't be used, and src.A_orig_pts should point to one of the 
  // static arrays
  A_orig_pts = src.A_orig_pts; 

  // Transform the center
  std::copy(T, T + 3, A_center);
  my_gemm(1, 3, 3, 1.0, src.A_center, 3, R, 3, A_center, 3, 1.0);

  // Transform the sphere points
  A_pts = new my_float_t[3 * A_num_pts];
  my_float_t* pt = A_pts;
  for(uint i = 0; i < A_num_pts; ++i, pt += 3) std::copy(T, T + 3, pt);
  my_gemm(A_num_pts, 3, 3, 1.0, src.A_pts, 3, R, 3, A_pts, 3, 1.0);
}

void 
DiscreteSphere::create_copy(const my_float_t* src, my_float_t** dest)
{
  uint N = 3*A_num_pts;
  *dest = new my_float_t[N];
  std::copy(src, src + N, *dest);
}
