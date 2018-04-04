
namespace ASCbase{
namespace optimization{

#if 1
class FirstOrderOpt{
public:

  virtual ~FirstOrderOpt();

  void move(my_float_t *X, const my_float_t *Y, const uint nrows)
  {
    // Multiply del vector to weight matrix to get gradient 
    // use stepsize to get how far to move in direction of the gradient 

    my_float_t V[] = { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    my_float_t dV[7];
    my_float_t gradient[7];
    
    const uint max_i = 100;
    h = 0.01;
    for(uint i = 0; i < max_i; ++i){

      del_vector(X, Y, nrows, V, dV);
      inner_prod(dV, V, gradient);     

    }
  }

private:
  virtual void 
  del_vector(const my_float_t *X, const my_float_t *Y, 
             const uint nrows, const my_float_t *V, my_float_t *dV) = 0;
/*
  virtual void 
  del_vector(const my_float_t *X, const my_float_t *Y, const uint nrows, 
             my_float_t **dX, uint *len_dX) = 0;
  virtual void
  inner_prod(const my_float_t* dX, const my_float_t *V, const uint len_dX, 
             my_float_t *prod) = 0;
*/

  virtual void
  inner_prod(const my_float_t* dV, const my_float_t *V, my_float_t *prod) = 0;
};
#endif

class GlobalOrientation{
public:

private:

  //! Derivative of global position at q = [1,0,0,0], t = [0, 0, 0] -- the 
  //! identity orientaion.
  /*! Assumptions:
   *  We seek the global transformation to move X to Y that minimizes the
   *  least squared error.
   *  X and Y are of dimension (nrows, 3)
   *#  V is the current transformation (length 7)
   *#    First 4 elements are for the quaterion
   *#    Last 3 elements are for the translation
   *  dV is the gradient (length 7)
   *    First 4 elements are for the quaterion
   *    Last 3 elements are for the translation
   */
  void blah(const my_float_t *X, const my_float_t *Y, const uint nrows,
            my_float_t *dV)
  {
    my_float_t *diff = new my_float_t[3*nrows];
    std::copy(X, X + 3*nrows, Y);
    my_axpy(3*nrows, -1.0, Y, 1, diff, 1);

    std::fill(dV, dV + 7, 0.0);
    const my_float_t *row_diff = diff;
    const my_float_t *x_row = X;
    for(uint i = 0; i < nrows; ++i, row_diff += 3, x_row += 3){
      my_float_t tmp[3];
      cross(x_row, row_diff, tmp);

      for(uint j = 0; j < 3; ++j){
        dV[j + 4] += row_diff[j]; 
        dV[0] += x_row[j] * row_diff[j];
        dV[j+1] += tmp[j];
      } 
    }
    for(uint j = 0; j < 3; ++j) dV[j + 4] *= 2.0;
    for(uint j = 0; j < 4; ++j) dV[j] *= 4.0;


//    const my_float_t *x = X;
//    const my_float_t *y = Y;
//    for(uint i = 0; i < nrows; ++i, x += 3, y += 3)
//      for(int j = 0; j < 3; ++j) dV[j + 4] += x[j] - y[j];
//    for(int j = 0; j < 3; ++j) dV[j + 4] *= 2.0;
     

    delete [] diff;
  }




};


#if 0
class GlobalOrientation //: public FirstOrderOpt{
public:

private:
  //! mess
  /*! Assumptions:
   *  We seek the global transformation to move X to Y that minimizes the
   *  least squared error.
   *  X and Y are of dimension (nrows, 3)
   *  V is the current transformation (length 7)
   *    First 4 elements are for the quaterion
   *    Last 3 elements are for the translation
   *  dV is the gradient (length 7)
   */
  void del_vector(const my_float_t *X, const my_float_t *Y, 
                  const uint nrows, const my_float_t *V, my_float_t *dV)
  {
    std::fill(dV, dV + 7, 0.0);
  
    // First we compute alpha = RX + tile(T) - Y
    Quaternion q(V, 4);
    my_float_t R[9];
    q.get_ortho_rot_mat(R); 
    const my_float_t *T = V + 4;

    my_float_t *alpha = new my_float_t[3*nrows];
    for(my_float_t *t = alpha; t < alpha + 3*nrow; t += 3){
      std::copy(T, T+3, t);
    }
    my_axpy(3*nrows, -1.0, Y, 1, alpha, 1);
    my_gemm(nrows, 3, 3, 1.0, X, 3, R, 3, alpha, 3, 1.0);

#if 0    
    const my_float_t M_q0[] = {  V[0], -V[3],  V[2],
                                 V[3],  V[0], -V[1],
                                -V[2],  V[1],  V[0] }; 
    const my_float_t M_q1[] = {  V[1],  V[2],  V[3],
                                 V[2], -V[1], -V[0],
                                 V[3],  V[0], -V[1] }; 
    const my_float_t M_q2[] = { -V[2],  V[1],  V[0],
                                 V[1],  V[2],  V[3],
                                -V[0],  V[3], -V[2]};
    const my_float_t M_q3[] = { -V[3], -V[0],  V[1],
                                 V[0], -V[3],  V[2],
                                 V[1],  V[2],  V[3]},
#endif
    const my_float_t M_q0[] = {  V[0],  V[3], -V[2],
                                -V[3],  V[0],  V[1],
                                 V[2], -V[1],  V[0] }; 
    const my_float_t M_q1[] = {  V[1],  V[2],  V[3],
                                 V[2], -V[1],  V[0],
                                 V[3], -V[0], -V[1] }; 
    const my_float_t M_q2[] = { -V[2],  V[1], -V[0],
                                 V[1],  V[2],  V[3],
                                 V[0],  V[3], -V[2] };
    const my_float_t M_q3[] = { -V[3],  V[0],  V[1],
                                -V[0], -V[3],  V[2],
                                 V[1],  V[2],  V[3] },

    my_float_t* U0 = new my_float_t[12*nrows];
    my_float_t* U1 = U0 + 3*nrows;
    my_float_t* U2 = U0 + 6*nrows;
    my_float_t* U3 = U0 + 9*nrows;
    my_gemm(nrows, 3, 3, 1.0, X, 3, M_q0, 3, U0, 3, 0.0);
    my_gemm(nrows, 3, 3, 1.0, X, 3, M_q1, 3, U1, 3, 0.0);
    my_gemm(nrows, 3, 3, 1.0, X, 3, M_q2, 3, U2, 3, 0.0);
    my_gemm(nrows, 3, 3, 1.0, X, 3, M_q3, 3, U3, 3, 0.0);

    for(uint i = 0; i < 3*nrows; ++i){
      dv[0] += U0[i]*alpha[i];
      dv[1] += U1[i]*alpha[i];
      dv[2] += U2[i]*alpha[i];
      dv[3] += U3[i]*alpha[i];
    }
    for(uint i = 0; i < 4; ++i) dv[i] *= 4.0;
    for(uint i = 0; i < 3*nrows; i += 3){
      dv[4] += alpha[i];
      dv[5] += alpha[i+1];
      dv[6] += alpha[i+2];
    }
    for(uint i = 0; i < 3; ++i) dv[4 + i] *= 2.0;

    delete [] U0;
    delete [] alpha;
  }

#if 0
  void del_vector(const my_float_t *X, const my_float_t *Y, 
                  const uint nrows, my_float_t **dX, uint *len_dX)
  {
    if(not nrows >= 3){
      std::cerr << "Global orientation requires at least 3 points\n"
                << std::endl;
      return;
    }

#if 0
    my_float_t M[] = { 1, -2,  2,
                       2,  1, -2,
                      -2,  2,  1};
#endif
    my_float_t M_T[] = { 1,  2, -2,
                        -2,  1,  2,
                         2, -2,  1};
    my_float_t *tmp = new my_float_t[3*nrows];
    std::copy(Y, Y + 3*nrows, tmp);
    my_gemm(nrows, 3, 3, 1.0, X, 3, M_T, 3, tmp, 3, -1.0);

    my_float_t blah[3];
    std::fill(blah, blah + 3, 0.0);
    const my_float_t *p = tmp;
    for(uint row = 0; row < nrows; ++row)
      for(int j = 0; j < 3; ++j, ++p) blah[j] += *p; 

    *dX = new my_float_t[7];
    *len_dX = 7;
    for(int j = 0; j < 3; ++j) (*dX)[j+4] = 2 * blah[j];

    my_float_t *my_dX = *dX;
    std::fill(my_dX, my_dX + 4, 0.0);
    p = tmp;
    const my_float_t *x = X;
    for(uint row = 0; row < nrows; ++row, p += 3, x += 3){
      my_dX[0] *= 4*x[0]*p[0] + 4*x[1]*p[1] + 4*x[2]*p[2];
      my_dX[1] *= -4*x[2]*p[1] + 4*x[1]*p[2];
      my_dX[2] *= 4*x[2]*p[0] -4*x[0]*p[2];
      my_dX[3] *= -4*x[1]*p[0] + 4*x[0]*p[1];
    }

    delete [] tmp;
  }
#endif

  inline void
  inner_prod(const my_float_t* dV, const my_float_t *V, my_float_t *prod)
  {
    my_gemm(1, 7, 1, 1.0, dV, 7, A_M, 7, prod, 7, 0.0);
  }

  static const my_float_t[] A_M;
};
#endif



} // end namespace optimization
} // end namespace ASCbase
