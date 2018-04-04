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
 * $Source: /psa/share/repository/pfizer_proj/src/basics/my_float_array.H,v $
 * $Revision: 1.5 $
 * $Author: vanvoor4 $
 * $Date: 2009-01-12 20:50:02 $
 * 
 * $Log: not supported by cvs2svn $
 * Revision 1.4  2007/12/17 21:03:12  vanvoor4
 * *** empty log message ***
 *
 * Revision 1.3  2007/09/20 17:37:50  vanvoor4
 * Removed underscore from define
 *
 * Revision 1.2  2007/08/29 20:14:55  vanvoor4
 * Added a method to compute the centroid of the stored points
 *
 * Revision 1.1  2007/08/21 15:44:54  vanvoor4
 * Initial checkin
 *
 *
 * 
 * 
 */

#ifndef MY_FLOAT_ARRAY_HEADER_INCLUDED
#define MY_FLOAT_ARRAY_HEADER_INCLUDED

#include <algorithm>
#include <mat_ops.H>
#include <math_basics.H>
#include <Quaternion.H>

namespace SimSite3D{

  typedef enum{
    MY_FLOAT_ARRAY_NULL,
    MY_FLOAT_ARRAY_GREW 
  }my_float_array_push_back_return_type;

  //! Note no copy constructor so copying will result in pain 
  class my_float_array{
  public:

    //! Cstor to set up floating pt array
    my_float_array(uint init_num_pts = 10, uint stride_in = 3);

    //! Copy in a reasonable manner
    my_float_array(const my_float_array& src)
    {
      init();
      do_copy(src);
    }

    //! Assignment operator -- other classes may need to take care of stale
    //! pointers
    const my_float_array&
    operator=(const my_float_array& src)
    {
      if(this != &src){
        if(orig) delete[] orig;
        if(curr) delete[] curr;
        if(scratch) delete[] scratch;
        init();

        do_copy(src);
      }
      return *this;
    }

    //! Dstor 
    ~my_float_array();

    uint 
    size() const
    { return (curr_end - curr) / _stride; }

    uint
    stride() const
    { return _stride; }
                                                                                
    uint
    reserved() const
    { return _reserved; }

    my_float_t*
    begin() const
    { return curr; }
  
    my_float_t*
    end() const
    { return curr_end; }

    //! Push back an element of length _stride onto the array
    /*! If el does not have a length of at least _stride, the operation is
     * undefined. 
     *
     * Warning:  do not use this function once the positions have been
     * transformed or the results may be different than expected
     */
    my_float_array_push_back_return_type
    push_back(const my_float_t* el)
    {
      my_float_array_push_back_return_type rv = MY_FLOAT_ARRAY_NULL;
      if(curr_end - curr >= _reserved){
        double_reserve();
        rv = MY_FLOAT_ARRAY_GREW;
      }
      std::copy(el, el + _stride, curr_end);

      curr_end += _stride;
      return rv;
    }

    //! Apply a rotation + translation current array
    /*!
     * For now a BLAS interface is more of a dream than reality and the code
     * used to implement BLAS routines would be at best the same speed as here
     * (even though this might be dastardly slow comapred with ATLAS BLAS)
     *
     * -- array version
     * for(uint z = 0; z < coords_end - coords; z += 3)
     *   for(uint i = 0; i < 3; ++i){
     *     scratch[z + i] = T[i];
     *     for(uint j = 0; j < 3; ++j)
     *       scratch[z + i] += R[3*i + j] * coords[z + j];
     *   }
     */
    void transform3D(const my_float_t* R, const my_float_t* T)
    {
      if(_stride != 3) return; // no desire to code up something other than 3x3
      int z_max = curr_end - curr;
      if(z_max > (orig_end - orig)) update_orig_pts();
    
      std::copy(curr, curr_end, scratch);
      for(my_float_t *p = curr; p < curr_end; p+=3) std::copy(T, T+3, p);
      my_gemm(z_max/_stride, 3, 3, 1.0, scratch, 3, R, 3, curr, 3, 1.0);
    }

    //! Apply a negative translation and inverse rotation
    /*! Transform the positions by -T and then rotate by R^t
     */
    void inverse_transform3D(const my_float_t* R, const my_float_t* T)
    {
      if(_stride != 3) return; // no desire to code up something other than 3x3
      int z_max = curr_end - curr;
      if(z_max > (orig_end - orig)) update_orig_pts();

      my_float_t R_transpose[9];
      for(int i = 0; i < 3; ++i)
        for(int j = 0; j < 3; ++j) R_transpose[3*j + i] = R[3*i + j];

      my_float_t T_inverse[3];
      for(int i = 0; i < 3; ++i) T_inverse[i] = -1.0 * T[i];
      my_float_t* scratch_end = scratch + z_max;
      for(my_float_t *p = scratch; p < scratch_end; p+=3) 
        std::copy(T_inverse, T_inverse+3, p);
      my_axpy(z_max, 1.0, curr, 1, scratch, 1);
    
      my_gemm(z_max/_stride, 3, 3, 1.0, scratch, 3, R_transpose, 3, curr, 3, 
              0.0);
    }
                                                                                
    //! Reset the positions of the atoms to their original positions
    void 
    revert()
    {
      // Really shouldn't ever run this here, but it is better to insure
      // that we don't core dump (or would it be better to core dump and know
      // that someone pushed additional values after a transformation was
      // applied? )
      if((curr_end - curr) > (orig_end - orig)) update_orig_pts();
      std::copy(orig, orig_end, curr);
    }

    //! Set the current values as the original values
    /*! This functionality is useful when we want to read in a set
     * of data points and apply some initial calculations and set the updated 
     * values as the original.  In this way, we will revert to the initialized 
     * values and not the values when the array was first created.
     */
    void
    set_current_values_as_original()
    {
      if(orig) delete [] orig;
      orig = 0;
      orig_end = 0;
    }

    //! Get the centroid (assumes 3D points)
    const bool 
    centroid_3D(my_float_t* C) const;

    //! Get the linear transformation to move the current points to the 
    //! original points
    /*! 
     *  Note: since a chain of transformations is in general difficult or
     *  even impossible to invert, this method is somewhat computationally
     *  expensive as it uses a least squares fit of all the points
     */
    inline void 
    get_current_inverse_3D_transform(Quaternion *Q, my_float_t *T) const
    {
      lse_3D_fit(orig, curr,  (orig_end - orig) / 3, Q, T);
    }

    //! Compute the root mean squared deviation (RMSD) between the current
    //! and orignial points
    my_float_t
    compute_RMSD() const
    {
      my_float_t my_sum = 0.0;
      const my_float_t *C = curr;
      for(const my_float_t *O = orig; O < orig_end; ++O, ++C){
        my_float_t tmp = *O - *C;
        my_sum += tmp*tmp;
      }
      size_t npts = (orig_end - orig) / 3;
      return std::sqrt(my_sum / npts) ;
    }

  private:
    void 
    init();

    //! Double the array reserve
    void double_reserve();

    //! Set the original array
    void update_orig_pts();

    //! Copy routine -- everyone else is expected to handle any stale 
    //! pointers & updating
    void do_copy(const my_float_array& src);

    uint _stride;
    int _reserved;
    my_float_t* orig;      //!< Beginning of points as added to this object
    my_float_t* orig_end;  //!< End of points as added to this object
    my_float_t* curr;      //!< Beginning of current points (after 0 or more transformations)
    my_float_t* curr_end;  //!< End of current points (after 0 or more transformations)
    my_float_t* scratch;  //!< Scratch space -- keep from frequently grabbing new mem
  };
}

#endif
