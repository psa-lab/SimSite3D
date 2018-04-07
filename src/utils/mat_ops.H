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

/*
 * $Source: /psa/share/repository/pfizer_proj/src/basics/mat_ops.H,v $
 * $Revision: 1.4 $
 * $Author: vanvoor4 $
 * $Date: 2008-07-29 14:33:53 $
 * 
 * $Log: not supported by cvs2svn $
 * Revision 1.3  2007/12/17 21:05:52  vanvoor4
 * Added and axpy
 *
 * Revision 1.2  2007/09/20 20:40:35  vanvoor4
 * Changed comment to correctly indicate which is the logical dimension
 * of C
 *
 * Revision 1.1  2007/08/21 15:37:35  vanvoor4
 * Moved from gen_points dir
 *
 *
 */

#ifndef _MAT_OPTS_H
#define _MAT_OPTS_H

#include <types.H>
#include <algorithm>

// Need to have a function prototype for the LAPACK DGESDD (REAL*8 SVD) 
extern "C" {
  extern void
  dgetrf_(const int *m, const int *n, double *A, const int *lda, int *IPIV,
          int *info);

  extern void
  dgetri_(const int *n, double *A, const int *lda, int *IPIV, double *work,
          const int *lwork, int *info);

  extern void
  dgesdd_(const char *jobz, const int *m, const int *n, double *A,
          const int *lda, double *S, double *U, const int *ldu,
          double *VT, const int *lvdt, double *work, const int *lwork,
          int *iwork, int *info);
}

inline my_float_t
row_dot(const int N, const my_float_t *A, const my_float_t *B)
{
 my_float_t rv = 0.0;
  if(A == B)
    for(int i = 0; i < N; ++i, ++A) rv += (*A)*(*A);
  else
    for(int i = 0; i < N; ++i, ++A, ++B) rv += (*A)*(*B);
  return rv;
}


//! Naive AXPY -- Y := alpha * X + Y
inline void 
my_axpy(const int N, const my_float_t alpha, const my_float_t* x,
        const int incX,  my_float_t* y, const int incY)
{
  if(alpha == 1.0)
    for(int i = 0; i < N; ++i, x += incX, y += incY) *y += *x;
  else if(alpha == -1.0)
    for(int i = 0; i < N; ++i, x += incX, y += incY) *y -= *x;
  else
    for(int i = 0; i < N; ++i, x += incX, y += incY) *y += alpha * (*x);
}

inline void
vector(const int N, const my_float_t *head, const my_float_t *tail,
       my_float_t *V)
{
  std::copy(head, head + N, V);
  my_axpy(N, -1.0, tail, 1, V, 1);
}

//! Compute the matrix A(A^T)
void
AAtranspose(const int M, const int N, const my_float_t *A, const int lda,
            my_float_t *B, const int ldb, const my_float_t beta);

//! Implemetation of *GEMV with BLAS calling style. 
/*! General matrix-vector multiplication routine using the BLAS level 2 calling
 * style to allow for a move to ATLAS routines if time permits.  This routine
 * produces the result Y <-- alpha*AX + beta*Y.  Note that there is no
 * support for operations on A.
 *
 * @param M number of rows in matrix A 
 * @param N number of columns in matrix A
 * @param alpha real scalar to scale the result of AX
 * @param A pointer to the matrix A
 * @param lda should be at least max(1, N)
 * @param X pointer to the vector X
 * @param incX stride for X (number of columns in X) 
 * @param beta real scalar to scale Y before adding Y to alpha*AX
 * @param Y pointer to the vector Y
 * @param incY stride for Y (number of columns in Y)
 */
void my_gemv(const int M, const int N, const my_float_t alpha,
             const my_float_t *A, const int lda, 
             const my_float_t *X, const int incX, 
             const my_float_t beta, my_float_t *Y, const int incY);

//! Implemetation of *GEMM with BLAS calling style. 
/*! General matrix-matrix multiplication routine using the BLAS level 3 calling
 * style to allow for move to ATLAS routines if time permits.  This routine
 * produces the result C <-- alpha * AB + beta * C.  Note that there is no
 * support for operations on A and B
 * 
 * @param M the number of rows in matrices A and C
 * @param N the number of columns in matrices B and C
 * @param K the number of columns in A and rows in B
 * @param alpha real scalar to scale AB
 * @param A pointer to matrix A
 * @param lda should be K for now
 * @param B pointer to matrix B
 * @param ldb should be N for now
 * @param C pointer to matrix C
 * @param ldc should be at least max(1, N)
 * @param beta real scalar to scale C (before adding C to alpha * AB).
 */
void my_gemm(const int M, const int N, const int K, const my_float_t alpha,
             const my_float_t *A, const int lda, 
             const my_float_t *B, const int ldb, 
             my_float_t *C, const int ldc, const my_float_t beta);

//void print_mat(double matr1[DIM][MAXSIZE], double matr2[DIM][MAXSIZE], 
		               //double matr3[DIM][MAXSIZE], int index);
//
//void print_one_mat(double matr[DIM][DIM]);

//void invert(double a[DIM][DIM], double res[DIM][DIM]);

bool invert_3_by_3(const my_float_t *A, my_float_t *A_inv);

my_float_t det_3x3(const my_float_t *M);

//! Get cofactor matrix C of a 4 by 4 matrix M
void cofactors_4x4(const my_float_t *M, my_float_t *C);

inline my_float_t
sum(const int N, const my_float_t* X, const int incX = 1)
{
  my_float_t s = 0;
  const my_float_t* X_end = X + incX * N;
  for(const my_float_t* x = X; x < X_end; x += incX) s += *x; 
  return s;
}

/*
//! Wrapper to FORTRAN77 REAL*8 SVD driver routine from LAPACK
inline int
dgesdd(char jobz, int m, int n, my_float_t *A, int lda, my_float_t *S,
       my_float_t *U, int ldu, my_float_t *VT, int lvdt, my_float_t *work,
       int lwork, int *iwork)
{
  int rv;
  dgesdd_(&jobz, &m, &n, A, &lda, S, U, &ldu, VT, &lvdt, work, &lwork, iwork, 
          &rv);
  return rv;
}

*/

inline int
la_dgetrf(const int m, const int n, double *A, const int lda, int *IPIV)
{
  int rv = 0;
  dgetrf_(&m, &n, A, &lda, IPIV, &rv);
  return rv;
}

inline int
la_dgetri(const int n, double *A, const int lda, int *IPIV, double *work, 
          const int lwork)
{
  int rv = 0;
  dgetri_(&n, A, &lda, IPIV, work, &lwork, &rv); 
  return rv;
}

#endif
