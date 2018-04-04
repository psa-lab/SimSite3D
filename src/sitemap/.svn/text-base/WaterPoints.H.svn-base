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
                                                                                
/*
 * $Source: /psa/share/repository/pfizer_proj/src/gen_points/WaterPoints.H,v $
 * $Revision: 1.1 $
 * $Author: vanvoor4 $
 * $Date: 2007-12-17 21:31:54 $
 *
 * $Log: not supported by cvs2svn $
 *
 *
 */

#ifndef WATER_POINTS_HEADER_INCLUDED
#define WATER_POINTS_HEADER_INCLUDED

#include <hbond_points.H>
#include <dir_point_storage.H>

namespace ASCbase{

  //! Setup a data class for placing water hbond points
  /*! This seems like a more straightforward way of doing things, but because
   * the hydrogen bonding method is already done for residues and it is
   * not a significant portion (in terms of computational speed) of the
   * template generation process there is no pressing need to change it
   */
  class WaterPoints{
  public: 
    //! Setup the hard coded relative positions for water hbond points
    WaterPoints();

    ~WaterPoints()
    { ; }

    //! Transform the positions and directions for the water hbond points 
    //! given the rigid transformation defined by R and T.
    void 
    transform(const my_float_t* R, const my_float_t* T)
    {
      water_pts.transform(R, T);
      water_ideal_pts.transform(R, T);
    }

    //! Transform the positions and directions for the water hbond points 
    //! by the inverse of the rigid transformation defined by R and T.
    void 
    inverse_transform(const my_float_t* R, const my_float_t* T)
    {
      water_pts.inverse_transform(R, T);
      water_ideal_pts.inverse_transform(R, T);
    }
 
    //! Revert the water points back to the positions relative to the 
    //! oxygen in H2O
    void 
    revert()
    {
      water_pts.revert();
      water_ideal_pts.revert();
    }

    hbond_ideal_pt_vec::const_iterator
    ideal_pts_begin() const
    { return water_ideal_pts.begin(); }

    hbond_ideal_pt_vec::const_iterator
    ideal_pts_end() const
    { return water_ideal_pts.end(); }

  private: 
    hbond_fit_pt_vec water_pts;
    hbond_ideal_pt_vec water_ideal_pts; 
  
    static const uint num_pts;  //!< Number of 3D points in realtive_pos
    static const my_float_t relative_pos[];  //!< Relative pos of H2O hbond pts
    static const my_float_t relative_dir[];  //!< Relative dir of H2O hbond pts
  };

}
#endif
