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
 * $Source: /psa/share/repository/pfizer_proj/src/basics/BoundingVolume.C,v $
 * $Revision: 1.2 $
 * $Author: vanvoor4 $
 * $Date: 2008/08/12 01:45:36 $
 * 
 * $Log: BoundingVolume.C,v $
 * Revision 1.2  2008/08/12 01:45:36  vanvoor4
 * Added UnionOfBalls class and uint ==> size_t
 *
 * Revision 1.1  2008/05/14 14:26:05  vanvoor4
 * Moved from gen_points
 *
 * Revision 1.5  2007/09/24 15:53:29  vanvoor4
 * Moved functions to .H to force inlining
 *
 * Revision 1.4  2007/08/21 18:07:45  vanvoor4
 * moved much of the work from GenPoints and else where to this class
 *
 * Revision 1.3  2006/11/16 20:24:00  vanvoor4
 * Added a quick, dirty method to sample a ball on a rectangular grid.
 *
 * Revision 1.2  2006/10/20 13:17:57  vanvoor4
 * Added a trace for bounding box size.
 *
 * Revision 1.1  2006/08/25 17:21:06  vanvoor4
 * Initial checkin
 *
 *
 */

#include <cstdio>
#include <cstring>
#include <sstream>
#include <BoundingVolume.H>
#include <math_basics.H>
//#define TRACE

using namespace ASCbase;

const std::string BoundingVolume::_fname = "BoundingVolume.C";
const my_float_t BoundingVolume::MIN_VDW_DIST = 2.5;
const my_float_t BoundingVolume::MAX_HYDRO_DIST = 5.2;
const my_float_t BoundingVolume::MAXRADDIST = 9.0;
const my_float_t BoundingVolume::MAXBINDDIST = 5.0;

BoundingVolume::BoundingVolume()
{
  grid_pts_beg = 0;
  grid_pts_end = 0;
}

BoundingVolume::BoundingVolume(const BoundingVolume& src)
{
  if(this == &src) return;

  grid_pts_beg = 0;
  grid_pts_end = 0;

  //std::cout << "Copy of Bounding Volume does not have grid points\n";
  grid_pts_beg = new my_float_t[src.grid_pts_end - src.grid_pts_beg];
  std::copy(src.grid_pts_beg, src.grid_pts_end, grid_pts_beg);
  grid_pts_end = src.grid_pts_end;
}

BoundingVolume::~BoundingVolume()
{
  clear();
}

void
BoundingVolume::clear()
{
  if(grid_pts_beg) delete [] grid_pts_beg;
  grid_pts_beg = 0;
  grid_pts_end = 0;
}

void 
BoundingVolume::copy_points(my_float_t *points, size_t n, size_t stride)
{
  size_t m;
  my_float_t *start = points + stride;
  for(m = 1; m <= n/2; m *= 2){
    memcpy(start, points, m * stride * my_float_size);
    start += m * stride; 
  }
  
  if(m < n) memcpy(start, points, (n-m) * stride * my_float_size); 
}

void 
BoundingVolume::set_grid_points(my_float_t* grid_pts, const size_t npts)
{
  clear();
  grid_pts_beg = grid_pts;
  grid_pts_end = grid_pts + 3*npts; 
}

void 
BoundingVolume::remove_nonsurface_points(chain_const_iter _begin,
                                         chain_const_iter _end)
{
  my_float_t* curr_pt = grid_pts_beg;
  my_float_t* surf_pts_end = grid_pts_beg;
  //size_t num_surface_pts = 0;

  for(size_t i = 0; curr_pt < grid_pts_end; i++, curr_pt += 3){
    bool too_close = false;
    bool is_surface_pt = false;

    // We cannot stop looking if a grid point is close enough untill we are
    // sure that an atom at that point will NOT have vdw overlap with any
    // protein atom.
    chain_const_iter chain = _begin;
    for( ; !too_close && chain < _end; chain++){
      residue_vci residue = chain->residues_begin;

      for( ; !too_close && residue < chain->residues_end; residue++){
        atom_vci atom = residue->atoms_begin;

        for( ; !too_close && atom != residue->atoms_end; atom++){
          my_float_t mydist = dist(curr_pt, atom->pos); 
          if((!atom->is_hbonder() && mydist < MIN_VDW_DIST)
             || (atom->is_hbonder() && mydist < MIN_HBOND_LENGTH))
            too_close = true;
          else if(mydist < MAX_HYDRO_DIST) is_surface_pt = true;
        }
      }
    }

    if(!too_close && is_surface_pt){
      if(surf_pts_end != curr_pt) std::copy(curr_pt, curr_pt +3, surf_pts_end);
      surf_pts_end += 3;
    }
  }
  
  // Precaution so that the rest of the points are less likely to be used 
  // accidentially.
  std::fill(surf_pts_end, grid_pts_end, 0);
  grid_pts_end = surf_pts_end;
}

RectangularSolid::RectangularSolid(const my_float_t *min_corner_in,
                                   const my_float_t *max_corner_in,
                                   const my_float_t grid_spacing)
 : BoundingVolume()
{
#ifdef TRACE
  fprintf(stdout, "Rectangular Bounding Box:\n");
  fprintf(stdout, "min %g %g %g\nmax %g %g %g\n", 
          min_corner_in[0], min_corner_in[1], min_corner_in[2],
          max_corner_in[0], max_corner_in[1], max_corner_in[2]);
#endif

  memcpy(min_corner, min_corner_in, 3*my_float_size);
  memcpy(max_corner, max_corner_in, 3*my_float_size);
  if(min_corner[0] >= max_corner[0] || min_corner[1] >= max_corner[1] || 
     min_corner[2] >= max_corner[2])
  err_msg(_fname, "RectangularSolid", "Bounds defining the "
          "rectangular solid are not valid");
  else{
    my_float_t *grid_points;
    size_t npts = discretize(grid_spacing, &grid_points);
    set_grid_points(grid_points, npts);
  }

  for(size_t i = 0; i < 3; ++i){
    BIND_min_corner[i] = min_corner[i] - MAXBINDDIST;
    BIND_max_corner[i] = max_corner[i] + MAXBINDDIST;
    RAD_min_corner[i] = min_corner[i] - MAXRADDIST;
    RAD_max_corner[i] = max_corner[i] + MAXRADDIST;
  }
}

RectangularSolid::RectangularSolid(const RectangularSolid& src)
 : BoundingVolume(src)
{
  if(this == &src) return;

  std::copy(src.min_corner, src.min_corner + 3, min_corner);
  std::copy(src.max_corner, src.max_corner + 3, max_corner);
  std::copy(src.BIND_min_corner, src.BIND_min_corner + 3, BIND_min_corner);
  std::copy(src.BIND_max_corner, src.BIND_max_corner + 3, BIND_max_corner);
  std::copy(src.RAD_min_corner, src.RAD_min_corner + 3, RAD_min_corner);
  std::copy(src.RAD_max_corner, src.RAD_max_corner + 3, RAD_max_corner);
}

size_t RectangularSolid::discretize(const my_float_t grid_spacing, 
                                  my_float_t **grid_points)
{
  size_t npts[3];
  for(size_t i = 0; i < 3; i++) 
    npts[i] = 
      static_cast<size_t>((max_corner[i] - min_corner[i]) / grid_spacing + 1);
 
  *grid_points = new my_float_t[3 * npts[0] * npts[1] * npts[2]];
  my_float_t *my_grid = *grid_points;
  my_float_t xval = min_corner[0];
  for(size_t i = 0; i < npts[0]; i++, xval += grid_spacing){
    my_grid[3*i] = xval;
    my_grid[3*i + 1] = min_corner[1];
    my_grid[3*i + 2] = min_corner[2]; 
  }

  const size_t stride = 3*npts[0];
  my_float_t yval = min_corner[1] + grid_spacing;
  BoundingVolume::copy_points(my_grid, npts[1], stride);
  for(size_t j = 1; j < npts[1]; j++, yval += grid_spacing)
    for(size_t i = 0; i < npts[0]; i++)
      my_grid[j*stride + 3*i + 1] = yval; 

  const size_t stride_z = 3 * npts[0] * npts[1];
  my_float_t zval = min_corner[2] + grid_spacing;
  BoundingVolume::copy_points(my_grid, npts[2], stride_z);
  for(size_t z = 1; z < npts[2]; z++, zval += grid_spacing)
    for(size_t i = 0; i < npts[0] * npts[1]; i++)
      my_grid[z*stride_z + 3*i + 2] = zval; 

  return npts[0] * npts[1] * npts[2];
}

std::string 
RectangularSolid::xml_str()
{
  my_float_t C[3], hw[3];
  for(size_t i = 0; i < 3; ++i){
    C[i] = 0.5 * (min_corner[i] + max_corner[i]);
    hw[i] = 0.5 * (max_corner[i] - min_corner[i]);
  }
  std::ostringstream rv;
  rv << "<hyperrectangle>" 
     << "<center>" << C[0] << " " << C[1] << " " << C[2] << "</center>"
     << "<half_width>" << hw[0] << " " << hw[1] << " " << hw[2]
     << "</half_width></hyperrectangle>";
  return rv.str();
}

// By no means optimal, but not a large part of computational time. 
// It may be more reasonable to discretize the sphere in shells rather than
// with a rectangular grid.

Ball::Ball(const my_float_t *center_in, const my_float_t radius_in,
           const my_float_t grid_spacing)
{
  radius = radius_in;
  if(radius < 0) radius *= -1.0;
  rad_squared = radius*radius;
  BIND_rad_squared = radius + MAXBINDDIST;
  BIND_rad_squared *= BIND_rad_squared;
  RAD_rad_squared = radius + MAXRADDIST;
  RAD_rad_squared *= RAD_rad_squared;

  memcpy(center, center_in, 3*my_float_size);
  if(grid_spacing != 0.0){
    my_float_t *grid_points;
    size_t npts = discretize(grid_spacing, &grid_points);
    set_grid_points(grid_points, npts);
  }
}

Ball::Ball(const Ball& src) : BoundingVolume(src)
{
  if(this == &src) return;

  std::copy(src.center, src.center + 3, center);
  radius = src.radius;
  rad_squared = src.rad_squared;
  BIND_rad_squared = src.BIND_rad_squared;
  RAD_rad_squared = src.RAD_rad_squared;
}

size_t 
Ball::discretize(const my_float_t grid_spacing, my_float_t **grid_points)
{
  size_t npts = static_cast<size_t>( std::ceil(2.0*radius/grid_spacing) );
  // Need to allocate mem -- cheat and use 3 * r^3? 
  *grid_points = new my_float_t[3*npts*npts*npts];
  memset(*grid_points, 0, 3*npts*npts*npts*my_float_size);

  my_float_t start[3];
  my_float_t pt[3];
  for(size_t i = 0; i < 3; i++) start[i] = center[i] - radius;
  size_t cnt = 0;
  pt[0] = start[0];
  for(size_t i = 0; i < npts; i++, pt[0] += grid_spacing){
    pt[1] = start[1];
    for(size_t j = 0; j < npts; j++, pt[1] += grid_spacing){
      pt[2] = start[2];
      for(size_t k = 0; k < npts; k++, pt[2] += grid_spacing){
        //if(this->includes_the_point(pt)){
        if(contains(pt)){
          memcpy(&((*grid_points)[3*cnt]), pt, 3*my_float_size);     
          cnt++;
        }
      }
    }
  }
      



  /* -- better distribution at some point
  // Generate slice obtained by rotating (r,0,0) to (0,r,0) 
  my_float_t *my_grid = *grid_points;
  size_t cnt = 0;
  for(my_float_t x = grid_spacing; x <= radius; x += grid_spacing){
    my_float_t ymax = std::sqrt(radius*radius - x*x);
    for(my_float_t y = grid_spacing; y <= ymax; y += grid_spacing){
      my_grid[0] = x;
      my_grid[1] = y; 
      my_grid += 3;
      cnt++;
    }
  }
  */
  
  return cnt; 
}

std::string 
Ball::xml_str()
{
  std::ostringstream rv;
  rv << "<sphere><center>" << center[0] << " " << center[1] << " " 
     << center[2] << "</center><radius>" << radius << "</radius></sphere>";
  return rv.str();
}

UnionOfBalls::UnionOfBalls(const size_t num_balls, const my_float_t *centers,
                           const my_float_t *radii)
{
  init(); 
  A_centers_begin = new my_float_t[3*num_balls];
  A_centers_end = A_centers_begin + 3*num_balls;
  A_radii_squared = new my_float_t[num_balls];
  std::copy(centers, centers + 3 * num_balls, A_centers_begin);
  for(size_t i = 0; i < num_balls; ++i) A_radii_squared[i] = radii[i]*radii[i]; 
}

UnionOfBalls::UnionOfBalls(const CoordFile &mol, const my_float_t tol)
{
  init();

  size_t sz = mol.positions_end() - mol.positions_begin();
  A_centers_begin = new my_float_t[sz];
  A_radii_squared = new my_float_t[sz/3];

  // Get the ligand atom centers -- skip the hydrogen atoms
  my_float_t* cur_pt = A_centers_begin;
  for(atom_vci a = mol.atoms_begin(); a < mol.atoms_end(); ++a){
    // Skip hydrogen atoms
    if(a->name == H || a->name == D) continue;

    std::copy(a->pos, a->pos + 3, cur_pt);
    cur_pt += 3;
  }
  A_centers_end = cur_pt;

  my_float_t *rad_squared = A_radii_squared; 
  for(size_t i = 0; i < sz/3; ++i, ++rad_squared) *rad_squared = tol*tol;
}

UnionOfBalls::UnionOfBalls(const UnionOfBalls& src)
{
  if(this == &src) return;

  init();
  size_t sz = src.A_centers_end - src.A_centers_begin;
  A_centers_begin = new my_float_t[sz];
  A_centers_end = A_centers_begin + sz;
  A_radii_squared = new my_float_t[sz/3];
  std::copy(src.A_centers_begin, src.A_centers_end, A_centers_begin);
  std::copy(src.A_radii_squared, src.A_radii_squared + sz/3, A_radii_squared); 
}

void
UnionOfBalls::init()
{
  A_centers_begin = 0;
  A_centers_end = 0;
  A_radii_squared = 0;
  A_BIND_rad_squared = MAXBINDDIST * MAXBINDDIST;
  A_RAD_rad_squared = MAXRADDIST * MAXRADDIST;
}

void
UnionOfBalls::clear()
{
  if(A_centers_begin) delete [] A_centers_begin;
  if(A_radii_squared) delete [] A_radii_squared;
  init();
}

std::string 
UnionOfBalls::xml_str()
{
  std::ostringstream rv;
  const my_float_t *center = A_centers_begin;
  const my_float_t *rad_squared = A_radii_squared;
  for( ; center < A_centers_end; center += 3, ++rad_squared)
    rv << "<sphere><center>" << center[0] << " " << center[1] << " " 
       << center[2] << "</center>"
       << "<radius>" << std::sqrt(*rad_squared) << "</radius></sphere>\n";
  return rv.str();
}
