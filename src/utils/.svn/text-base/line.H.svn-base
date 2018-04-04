
#ifndef LINE_HEADER_FILE_INCLUDED
#define LINE_HEADER_FILE_INCLUDED

#include <string_basics.H>
#include <math_basics.H>

namespace SimSite3D{
namespace geometry{

class line_t{
public:
#if 0
  line_t(const my_float_t* slope, const my_float_t* intercept,
         const my_float_t* end_pts = 0)
  {
    std::copy(slope, slope + 3, A_slope);
    std::copy(intercept, intercept + 3, A_intercept);
  }
#endif
  line_t(const my_float_t* p0, const my_float_t* p1)
  {
    A_edge_len = unit_vector(A_slope, p0, p1);
    std::copy(p0, p0 + 3, A_end_pts);
    std::copy(p1, p1 + 3, A_end_pts + 3);
  }

  line_t(const my_float_t* end_pts)
  {
    A_edge_len = unit_vector(A_slope, end_pts, end_pts + 3);
    std::copy(end_pts, end_pts + 6, A_end_pts);
  }

  ~line_t()
  { ; }

  line_t(const line_t& src)
  {
    do_copy(src);
  }

  const line_t&
  operator=(const line_t& src)
  {
    if(this != &src) do_copy(src);
    return *this;
  }


  bool
  corresponding_point(const my_float_t* pt, my_float_t* pt_on_line,
                      const bool line_segment=true)
  {
    // Get vector from first end_point to pt
    my_float_t V[3];
    vector(3, pt, A_end_pts, V);

    // Compute the projection of V onto A_slope.  Since A_slope is a unit 
    // vector, the projection of the point is the number of units in the 
    // A_slope direction from the first end point.
    my_float_t proj_pt_dist = dot(V, A_slope);
    bool rv = true;

    // CASE 0:  The correspondence is the closest point on the line (not line
    // segment) containing (defined by) Vi and Vj.
    if(!line_segment){
      std::copy(A_end_pts, A_end_pts+3, pt_on_line);
      my_axpy(3, proj_pt_dist, A_slope, 1, pt_on_line, 1);
  
    // CASE 1:  The projection of the point must be contained in the line 
    // segment.  If the projection is not contained in the line segment, set
    // the corresponding point as the closer of the two end points and return
    // false to denote the correspondance is not contained in the line segment.
    }else{
      if(proj_pt_dist < 0.0){
        std::copy(A_end_pts, A_end_pts + 3, pt_on_line);
        rv = false;
      }else if(proj_pt_dist > A_edge_len){
        std::copy(A_end_pts + 3, A_end_pts + 6, pt_on_line);
        rv = false;
      }else{
        std::copy(A_end_pts, A_end_pts+3, pt_on_line);
        my_axpy(3, proj_pt_dist, A_slope, 1, pt_on_line, 1);
      }
    }
  
    return rv;
  }

private:
  void
  do_copy(const line_t& src)
  {
    std::copy(src.A_end_pts, src.A_end_pts + 6, A_end_pts);
    std::copy(src.A_slope, src.A_slope + 6, A_slope);
    A_edge_len = src.A_edge_len;
  }
 
  my_float_t A_end_pts[6];
  my_float_t A_slope[3];
  my_float_t A_edge_len;

//  my_float_t A_intercept[3];
};

} // end namespace geometry
} // end namespace SimSite3D

#endif
