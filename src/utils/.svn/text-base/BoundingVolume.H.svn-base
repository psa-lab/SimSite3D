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
 * $Source: /psa/share/repository/pfizer_proj/src/basics/BoundingVolume.H,v $
 * $Revision: 1.3 $
 * $Author: vanvoor4 $
 * $Date: 2008-08-12 01:45:36 $
 * 
 * $Log: not supported by cvs2svn $
 * Revision 1.2  2008/07/28 15:25:50  vanvoor4
 * Removed old commented out code
 *
 * Revision 1.1  2008/05/14 14:26:05  vanvoor4
 * Moved from gen_points
 *
 * Revision 1.5  2007/12/17 21:30:01  vanvoor4
 * Added a distance function so that it can be changed (transparently)
 * in the future.
 *
 * Revision 1.4  2007/09/24 15:53:02  vanvoor4
 * Added the ability to change the buffer (threshold) and moved the
 * computation of point included in the volume to this file to force
 * inline
 *
 * Revision 1.3  2007/08/21 18:08:00  vanvoor4
 * Moved much of the work for bounding volumes to this class
 *
 * Revision 1.2  2007/02/07 15:23:10  vanvoor4
 * Header files changed.
 *
 * Revision 1.1  2006/08/25 17:21:09  vanvoor4
 * Initial checkin
 *
 *
 */

#ifndef _BOUNDING_VOLUME_H
#define _BOUNDING_VOLUME_H

#include <PDBStructure.H>

namespace ASCbase{

//! Simple base class to define a common interface to 3D bounding volumes.
class BoundingVolume{
public:

  typedef enum{
    NULL_VOLUME_TYPE,
    RECTANGULAR_SOLID,
    SPHERE
  }volume_t;

  //! Initialize pointers to NULL
  BoundingVolume();

  //! Basic copy cstr
  BoundingVolume(const BoundingVolume& src);

  virtual ~BoundingVolume();

  void remove_nonsurface_points(chain_const_iter _begin,
                                chain_const_iter _end);

  //! Check if any grid point is within tolerance of given point
  my_float_t
  dist_to_nearest_grid_pt(const my_float_t* given_pt, const my_float_t tol)
  {
    const my_float_t squared_tol = tol*tol;
    my_float_t min_dist_2 = my_float_max;
    for(const my_float_t* pt = grid_pts_beg; pt < grid_pts_end; pt += 3){
      my_float_t tmp = dist_squared(given_pt, pt);
      min_dist_2 = (min_dist_2 > tmp ? tmp : min_dist_2);
      if(min_dist_2 <= squared_tol) break;
    } 
    return std::sqrt(min_dist_2);
  }

  //! Is the given point inside the bounding volume?
  virtual bool contains(const my_float_t *point) const = 0;

  virtual bool BIND_vol_contains(const my_float_t *p) const = 0;

  virtual bool RAD_vol_contains(const my_float_t *p) const = 0;
   
  //! discretize the volume using the spacing for grid spacing.
  virtual size_t discretize(const my_float_t spacing, 
                          my_float_t **grid_pts) = 0;

  const my_float_t* grid_points_begin() const { return grid_pts_beg; }

  const my_float_t* grid_points_end() const { return grid_pts_end; }

  void clear();

  static const my_float_t MIN_VDW_DIST;

  virtual std::string xml_str() = 0;

protected:
  void copy_points(my_float_t *points, size_t n, size_t stride);

  void set_grid_points(my_float_t* grid_pts, const size_t npts);

  static const my_float_t MAXRADDIST;
  static const my_float_t MAXBINDDIST;
  static const std::string _fname; //!< source file name
private:
  my_float_t* grid_pts_beg;
  my_float_t* grid_pts_end;

  static const my_float_t MAX_HYDRO_DIST;
};

//! A rectangular bounding volume
class RectangularSolid : public BoundingVolume{
public:
  RectangularSolid(const my_float_t *min_corner_in,
                   const my_float_t *max_corner_in,
                   const my_float_t grid_spacing = 0.0);

  RectangularSolid(const RectangularSolid& src);

  ~RectangularSolid() {;}

  virtual bool
  contains(const my_float_t *point) const
  {
    for(size_t i = 0; i < 3; ++i)
      if(point[i] < min_corner[i] || point[i] > max_corner[i]) return false;
    return true;
  }

  virtual bool
  BIND_vol_contains(const my_float_t *p) const
  {
    for(size_t i = 0; i < 3; ++i)
      if(p[i] < BIND_min_corner[i] || p[i] > BIND_max_corner[i]) return false;
    ++p;
    return true;
  }

  virtual bool
  RAD_vol_contains(const my_float_t *p) const
  {
    for(size_t i = 0; i < 3; ++i)
      if(p[i] < RAD_min_corner[i] || p[i] > RAD_max_corner[i]) return false;
    return true;
  }

  virtual size_t discretize(const my_float_t spacing, my_float_t **grid_pts);

  virtual std::string xml_str();

private:

  my_float_t min_corner[3];
  my_float_t max_corner[3];
  my_float_t BIND_min_corner[3];
  my_float_t BIND_max_corner[3];
  my_float_t RAD_min_corner[3];
  my_float_t RAD_max_corner[3];
};

//! A spherical solid bounding volume
class Ball : public BoundingVolume{
public:
  Ball(const my_float_t *center_in, const my_float_t radius_in,
       const my_float_t grid_spacing = 0.0);

  Ball(const Ball& src);

  ~Ball() {;}

//  void
//  add_tolerance(my_float_t tol_change = 0.0)
//  { radius += tol_change; }

  virtual bool
  contains(const my_float_t *point) const
  { return dist_squared(point, center) <= rad_squared; }

  virtual bool
  BIND_vol_contains(const my_float_t *point) const
  { return dist_squared(point, center) <= BIND_rad_squared; }

  virtual bool
  RAD_vol_contains(const my_float_t *point) const
  { return dist_squared(point, center) <= RAD_rad_squared; }

  virtual size_t discretize(const my_float_t spacing, my_float_t **grid_pts);

  virtual std::string xml_str();
private:

  my_float_t center[3];
  my_float_t radius;
  my_float_t rad_squared;
  my_float_t BIND_rad_squared;
  my_float_t RAD_rad_squared;
};

//! Union of spherical solids bounding volume
class UnionOfBalls : public BoundingVolume{
public:
  UnionOfBalls(const size_t num_balls, const my_float_t *centers, 
               const my_float_t *radii);

  UnionOfBalls(const CoordFile &mol, const my_float_t tol = 3.0);

  UnionOfBalls(const UnionOfBalls& src);

  virtual ~UnionOfBalls()
  { clear(); }

  virtual bool
  contains(const my_float_t *point) const
  {
    bool rv = false;
    const my_float_t *center = A_centers_begin;
    const my_float_t *rad_squared = A_radii_squared;
    for( ; center < A_centers_end; center += 3, ++rad_squared) 
      if(dist_squared(point, center) <= *rad_squared){
        rv = true;
        break;
      }
    return rv;
  }

  virtual bool
  BIND_vol_contains(const my_float_t *point) const
  { 
    bool rv = false;
    const my_float_t *center = A_centers_begin;
    for( ; center < A_centers_end; center += 3)
      if(dist_squared(point, center) <= A_BIND_rad_squared){
        rv = true;
        break;
      }
    return rv;
  }
    
  virtual bool
  RAD_vol_contains(const my_float_t *point) const
  {
    bool rv = false;
    const my_float_t *center = A_centers_begin;
    for( ; center < A_centers_end; center += 3)
      if(dist_squared(point, center) <= A_RAD_rad_squared){
        rv = true;
        break;
      }
    return rv;
  }

  virtual size_t discretize(const my_float_t spacing, my_float_t **grid_pts)
  { std::cerr << "Unimplemented function!!!!\n"; return 0; }

  virtual std::string xml_str();
private:
  void init();
  
  void clear();

  my_float_t *A_centers_begin;
  my_float_t *A_centers_end;
  my_float_t *A_radii_squared;
  my_float_t A_BIND_rad_squared;
  my_float_t A_RAD_rad_squared;
};

}
#endif
