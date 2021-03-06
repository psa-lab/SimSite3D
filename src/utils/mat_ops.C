/******************************************************************************
 * Copyright (c) 2006-2011, Michigan State University (MSU) Board of Trustees.
 * This file is part of the SimSite3D software project.
 *
 * Authors: Jeffrey Van Voorst, jeff.vanvoorst@gmail.com
 *          Leslie Kuhn, Ph.D., KuhnL@msu.edu 
 *
 *  SimSite3D is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  SimSite3D is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc.,
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *  You may also visit http://www.gnu.org/licenses/gpl-2.0.html.
 * 
 * $Date: 2011-02-25 09:54:04 -0500 (Fri, 25 Feb 2011) $: Date of last commit
 * $Author: vanvoor4 $: Author of last commit
 * $Rev: 1617 $: svn revision of last commit
 *
 * svn file: $Id: mat_ops.C 1617 2011-02-25 14:54:04Z vanvoor4 $
 * file location: $URL: file:///psa/share/repository/SimSite3D/branches/surfaces-branch/src/utils/mat_ops.C $
 *****************************************************************************/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <mat_ops.H>
#include <string_basics.H>
#include <sstream>

void my_gemv(const int M, const int N, const my_float_t alpha,
             const my_float_t *A, const int lda, 
             const my_float_t *X, const int incX, 
             const my_float_t beta, my_float_t *Y, const int incY)
{
  const int N_diff = lda - N;
  my_float_t* y_ptr = Y;
  const my_float_t *a = A;
  if(alpha == 1.0){
    for(int m = 0; m < M; ++m, y_ptr += incY, a += N_diff){
      my_float_t& y = *y_ptr;
      y = (beta == 0.0 ? 0.0 : beta*y);
      const my_float_t* x = X;
      for(int n = 0; n < N; ++n, x += incX, ++a) y += (*a) * (*x);
    }
  }else{
    for(int m = 0; m < M; ++m, y_ptr += incY, a += N_diff){
      my_float_t& y = *y_ptr;
      y = (beta == 0.0 ? 0.0 : beta*y);
      const my_float_t* x = X;
      for(int n = 0; n < N; ++n, x += incX, ++a) y += alpha * (*a) * (*x);
    }
  }
}

void my_gemm(const int M, const int N, const int K, const my_float_t alpha,
             const my_float_t *A, const int lda, 
             const my_float_t *B, const int ldb, 
             my_float_t *C, const int ldc, const my_float_t beta)
{
  my_float_t *alpha_A = 0;
  if(alpha != 1.0){
    alpha_A = new my_float_t[M*K];
    for(int m = 0; m < M; ++m)
      for(int k = 0; k < K; ++k) alpha_A[m*K + k] = alpha * A[m*lda + k];
  }

  if(beta != 1.0){
    for(int m = 0; m < M; ++m)
      for(int n = 0; n < N; ++n)
        C[m*ldc + n] = (beta == 0.0 ? 0.0 : C[m*ldc + n] * beta);
  }

  if(alpha_A){
    for(int m = 0; m < M; ++m)
      for(int k = 0; k < K; ++k)
        for(int n = 0; n < N; ++n)
          C[m*ldc + n] += alpha_A[m*K + k] * B[k*ldb + n];
    delete [] alpha_A;
    alpha_A = 0;
  }else{
    for(int m = 0; m < M; ++m)
      for(int k = 0; k < K; ++k)
        for(int n = 0; n < N; ++n)
          C[m*ldc + n] += A[m*lda + k] * B[k*ldb + n];
  }
}

void
AAtranspose(const int M, const int N, const my_float_t *A, const int lda,
            my_float_t *B, const int ldb, const my_float_t beta)
{
  if(beta != 1.0){
    my_float_t *row_B = B;
    for(int i = 0; i < M; ++i, row_B += ldb)
      for(int j = 0; j < M; ++j) row_B[j] *= beta;
  }

  const my_float_t *row_A = A;
  my_float_t *row_B = B;
  for(int i = 0; i < M; ++i, row_A += lda, row_B += ldb){
    const my_float_t *col_Atrans = A;
    for(int j = 0; j < M; ++j, col_Atrans += lda)
      row_B[j] += row_dot(N, row_A, col_Atrans);
  }
}
              
bool 
invert_3_by_3(const my_float_t *A, my_float_t *A_inv)
{
  my_float_t D[] = {  A[4]*A[8] - A[5]*A[7],
                    -(A[1]*A[8] - A[2]*A[7]), 
                      A[1]*A[5] - A[2]*A[4],
                     
                    -(A[3]*A[8] - A[5]*A[6]),
                      A[0]*A[8] - A[2]*A[6],
                    -(A[0]*A[5] - A[2]*A[3]),

                      A[3]*A[7] - A[4]*A[6],
                    -(A[0]*A[7] - A[1]*A[6]),
                      A[0]*A[4] - A[1]*A[3]  };

 
  my_float_t det = A[0]*D[0] + A[1]*D[3]  + A[2]*D[6];
  if(det == 0){
    warn("mat_ops.C", "invert_3_by_3()", 
         "Determinant of matrix is zero -- cannot invert");
    return false;
  }else if(std::fabs(det) < 10E-15){
    std::ostringstream ostr;
    ostr << "Determinat of matrix is very small (" << det 
         << ") -- cannot invert";
    warn("mat_ops.C", "invert_3_by_3()", ostr.str());
    return false;
  }

  for(uint i = 0; i < 9; ++i) A_inv[i] = D[i] / det;
  return true;
}

void
cofactors_4x4(const my_float_t *M, my_float_t *C)
{
  my_float_t sub_mat[9];
  size_t skip_row = 0, skip_col = 0;
  for(size_t z = 0; z < 16; ++z, ++C){
    const my_float_t *m = M;
    my_float_t *s = sub_mat;
    for(size_t i = 0; i < 4; ++i){
      if(i == skip_row){
        m += 4;
        continue;
      }

      for(size_t j = 0; j < 4; ++j, ++m){
        if(j == skip_col) continue;
        *s = *m; 
        ++s;
      }
    }
    // Zero out ptrs just because
    s = 0;
    m = 0;

    *C = det_3x3(sub_mat);
    if(z % 2) *C *= -1.0;
 
    ++skip_col;
    if(skip_col > 3){
      skip_col = 0;
      ++skip_row;
    }
  }
}

my_float_t
det_3x3(const my_float_t *M)
{
  return M[0]*M[4]*M[8] + M[1]*M[5]*M[6] + M[2]*M[3]*M[7] - 
         M[2]*M[4]*M[6] - M[5]*M[7]*M[0] - M[8]*M[1]*M[3]; 
}
