

#ifndef SPHERES_HEADER_FILE_INCLUDED
#define SPHERES_HEADER_FILE_INCLUDED

#include <string_basics.H>
#include <math_basics.H>

namespace SimSite3D{
namespace geometry{

class sphere_t{
public:
  //! Listing of the possible intersection types for two spheres
  typedef enum{
    NONE = 0,              //!< No intersection or containment
    INTERSECTION,          //!< A significant intersection -- not at 1 point
    OTHER_CONTAINS_THIS,   //!< The other sphere contains this one -- technically not an intersection
    THIS_CONTAINS_OTHER    //!< This sphere contains the other one -- technically not an intersection
  }intersectionType;

  //! Default constructor
  sphere_t()
  { 
    std::fill(A_center, A_center + 3, 0.0);
    A_sq_rad = 0.0;
    A_radius = -1.0;
  }

  //! Set the center & radius of the sphere
  sphere_t(const my_float_t* c, const my_float_t r)
  {
    std::copy(c, c+3, A_center);
    std::copy(c, c+3, A_orig_center);
    A_radius = r;
    A_sq_rad = A_radius * A_radius;
  }

  //! copy constructor
  sphere_t(const sphere_t& src)
  {
    do_copy(src);
  }

  const sphere_t&
  operator=(const sphere_t& src)
  {
    if(this != &src) do_copy(src);
    return *this;
  }

  //! Test if 3D pt is inside sphere
  const bool
  contains(const my_float_t* pt) const
  { return (A_sq_rad >= dist_squared(A_center, pt)); }

  //! Get (a pointer to) the center of the sphere
  const my_float_t*
  center() const
  { return A_center; }

  //! Get the radius of the sphere
  const my_float_t
  radius() const
  { return A_radius; }

  //! Get the squared radius of the sphere
  const my_float_t
  squared_radius() const
  { return A_sq_rad; }

  //! overload ==
  const bool
  operator==(const sphere_t& other)
  {
    if(A_radius != other.A_radius) return false;
    for(int i = 0; i < 3; ++i) 
      if(A_center[i] != other.A_center[i]) return false;
    return true;
  }

  friend std::ostream& 
  operator<<(std::ostream &out, const sphere_t& S)
  { 
    out << S.A_center[0] << " " << S.A_center[1] << " " << S.A_center[2] << " "
        << S.A_radius;
    return out;
  }

  //! Project the point to the closest point on the surface of the sphere
  void
  project_point(const my_float_t* pt, my_float_t* proj_pt) const
  {
    my_float_t V[3];
    unit_vector(V, pt, A_center);
    std::copy(A_center, A_center + 3, proj_pt);
    my_axpy(3, A_radius, V, 1, proj_pt, 1);
  }

  //! Determine the intersection of the line y = mx + b with the sphere
  /*! Assumptions:
   *    m, b are of dimension 3 
   *    pts is of at least length 6*my_float_t
   */
  int
  intersection(const my_float_t* m, const my_float_t* b,
               my_float_t* pts)
  {
    // Get vector from point on line (intercept) to the center of the sphere
    my_float_t b2center[3];
    std::copy(A_center, A_center + 3, b2center);
    my_axpy(3, -1.0, b, 1, b2center, 1); 

    // Project bc onto m
    my_float_t proj_len = dot(m, b2center);

    // Computpe the squared distance from the center of the circle to the line
    my_float_t b2center_sq_len = dot(b2center, b2center);
    my_float_t h2 = b2center_sq_len - proj_len*proj_len;

    // Line lies outside of sphere
    if(h2 >= A_sq_rad) return 0;

    // We know the radius and the "height" of the triangle, get the length of
    // the side of the triangle that is lying on the line y = mx + b
    my_float_t d = std::sqrt(A_sq_rad - h2);
    std::copy(b, b+3, pts);
    my_axpy(3, proj_len, m, 1, pts, 1);
// Need to determine how to correctly check since eps isn't always the right
// threshold
    //if(-1.0*my_float_eps < d && d < my_float_eps) return 1;
    std::copy(pts, pts+3, pts+3);
    my_axpy(3, -1.0*d, m, 1, pts, 1);
    my_axpy(3, d, m, 1, pts+3, 1);

    return 2;
  }

  //! Intersection with another sphere -- get the volume of this spherical
  //! cap that is inside the other circle
  intersectionType
  intersects(const sphere_t& other, my_float_t *vol) const
  {
    // Check if the two spheres intersect -- assume intersection at a point
    // is not a significant intersection 
    my_float_t D2 = dist_squared(A_center, other.A_center);
    if(D2 >= (A_radius + other.A_radius)*(A_radius + other.A_radius))
      return NONE;

    // Check if this sphere contains the other one
    if(contains(other)){
      *vol = 4.0/3.0 * M_PI * other.A_radius*other.A_radius*other.A_radius;
      return THIS_CONTAINS_OTHER;
    }
    if(other.contains(*this)){
      *vol = 4.0/3.0 * M_PI * A_radius*A_radius*A_radius;
      return OTHER_CONTAINS_THIS;
    }

    ///////////////////////////////////////////////////////////////////////// 
    // Determine the volume of this sphere's cap inside the other
    // V_cap = 1/3 pi h^2 (3R - h) -- mathworld.wolfram.com/SphericalCap.html
    ///////////////////////////////////////////////////////////////////////// 

    // What is the height of the cap? -- use law of cosines + algebra
    my_float_t h;
    if(A_radius >= other.A_radius){
      if(contains(other.A_center)){
        my_float_t d = -A_radius * (other.A_radius*other.A_radius - A_radius*A_radius - D2) / (-2.0 * A_radius * std::sqrt(D2));
        h = A_radius - d;
       
        // Here h is the height of the cap of the smaller sphere that lies
        // outside of the larger sphere (here the other sphere has the smaller
        // radius
        *vol = -1.0 * (1.0/3.0 * M_PI * h*h * (3*other.A_radius - h));
        *vol += 4.0/3.0 * M_PI * other.A_radius*other.A_radius*other.A_radius;
        
      }else{
        // d is the portion of the radius of this sphere that is out side of
        // the other sphere
        my_float_t d = 0.5 * ((D2 + A_radius*A_radius - 
                               other.A_radius*other.A_radius) / std::sqrt(D2));
        h = A_radius - d;
        *vol = 1.0/3.0 * M_PI * h*h * (3*A_radius - h);
      }
    }else{
      if(other.contains(A_center)){
        my_float_t d = -A_radius * (other.A_radius*other.A_radius - A_radius*A_radius - D2) / (-2.0 * A_radius * std::sqrt(D2));
        h = A_radius - d;

        // Here h is the height of the cap of the smaller sphere that lies
        // outside of the larger sphere (here this sphere has the smaller
        // radius
        *vol = -1.0 * (1.0/3.0 * M_PI * h*h * (3*A_radius - h));
        *vol += 4.0/3.0 * M_PI * A_radius*A_radius*A_radius;
        
      }else{
        // d is the portion of the radius of this sphere that is out side of
        // the other sphere
        my_float_t d = 0.5 * ((D2 + A_radius*A_radius - 
                               other.A_radius*other.A_radius) / std::sqrt(D2));
        h = A_radius - d;
        *vol = 1.0/3.0 * M_PI * h*h * (3*A_radius - h);
      } 
    }

    return INTERSECTION; 
  }

  //! Does this sphere contain the other one?
  bool
  contains(const sphere_t& other) const
  {
    if(!contains(other.A_center)) return false;

    // Get a unit vector from this sphere's center in the direction of the
    // other spheres's center
    my_float_t V[3];
    unit_vector(V, other.A_center, A_center); 

    // Project from the center of the other sphere along the unit vector's 
    // direction
    my_float_t proj_pt[3];
    std::copy(other.A_center, other.A_center + 3, proj_pt);
    my_axpy(3, other.A_radius, V, 1, proj_pt, 1);

    // If that point is inside this sphere, then it is contained in
    // this one
    return contains(proj_pt); 
  }

  //! Transform the center 
  void
  transform(const my_float_t* R, const my_float_t* T)
  { 
    my_float_t scratch[3];
    std::copy(A_center, A_center + 3, scratch);
    std::copy(T, T + 3, A_center);
    my_gemm(1, 3, 3, 1.0, scratch, 1, R, 3, A_center, 1, 1.0);
  }
  
  //! Transform the center by the inverse transformation 
  void
  inverse_transform(const my_float_t* R, const my_float_t* T)
  { 
    my_float_t scratch[3];
    for(int i = 0; i < 3; ++i) scratch[i] = A_center[i] - T[i];
    my_gemm(3, 1, 3, 1.0, R, 3, scratch, 1, A_center, 1, 0.0);
  }
  
  //! Revert the current center to the original center
  void
  revert()
  { std::copy(A_orig_center, A_orig_center + 3, A_center); }

private:
  void
  do_copy(const sphere_t& src)
  {
    std::copy(src.A_center, src.A_center + 3, A_center);
    std::copy(src.A_orig_center, src.A_orig_center + 3, A_orig_center);
    A_radius = src.A_radius;
    A_sq_rad = src.A_sq_rad;
  }

  my_float_t A_radius;     //!< Radius of the sphere
  my_float_t A_sq_rad;     //!< Squared radius 
  my_float_t A_center[3];  //!< Center of the sphere 
  my_float_t A_orig_center[3]; //!< original center of the sphere
};

//! Computational representation of an analyticaly defined 3D plane
class plane_t{
public:
  //! Default constructor
  plane_t()
  { 
    std::fill(A_N, A_N + 3, 0.0);
    std::fill(A_P0, A_P0 + 3, 0.0);
    A_P0_is_defined = false;
    A_d = 0.0;
  }

  plane_t(const my_float_t* N, const my_float_t *P0)
  {
    std::copy(N, N+3, A_N);
    std::copy(P0, P0 + 3, A_P0);
    A_P0_is_defined = true;
    A_d = -1.0 * dot(N, P0);
  }

  plane_t(const my_float_t* N, const my_float_t d)
  {
    std::copy(N, N+3, A_N);
    A_d = d;
    A_P0_is_defined = false;
  }

  plane_t(const plane_t& src)
  {
    do_copy(src);
  }

  const plane_t&
  operator=(const plane_t& src)
  {
    if(this != &src) do_copy(src);
    return *this;
  }

  ~plane_t()
  { ; }

  my_float_t 
  signed_dist(const my_float_t* pt) const 
  {
    if(!A_P0_is_defined){
      printf("ERROR, P0 is not defined (in signed_dist)!!!!!!!!!!!!\n");
      //A_P0[2] = 0.0;
       
    }
    my_float_t V[3];
    vector(3, pt, A_P0, V);
    return dot(V, A_N);
  }

  //! Plane - plane intersection
  /*! Numerical approximation to intersection via Mathworld plane-plane page
   * replace this with SVD or other linear solution to [N1,N2]^T x = -[d1;d2]
   */
  void
  intersection(const plane_t& other, my_float_t* m, my_float_t* b)
  {
    cross(A_N, other.A_N, m);

    // Determine the component of m with the largest magnitude
    my_float_t n = 0.0;
    int idx = argmaxabs<my_float_t>(m, 3, &n);

    // Solve the reduced system (ignore the coordinate with the largest normal
    // component)
    my_float_t A[4];
    my_float_t *a = A;
    for(int i = 0; i < 3; ++i){
      if(i == idx) b[i] = 0.0;
      else{
        *a = A_N[i];
        *(a+2) = other.A_N[i];
        ++a;
      }
    }
    a = 0;

    my_float_t det_A = A[0]*A[3] - A[2]*A[1];
    my_float_t A_inv[] = {A[3] / det_A,      -1.0*A[1] / det_A,
                          -1.0*A[2] / det_A, A[0] / det_A      };
    a = A_inv;
    for(int i = 0; i < 3; ++i, ++b){
      if(i == idx) continue;
      *b = -1.0 * ((*a) * A_d + (*(a+1)) * other.A_d);
      a += 2;
    }
  }

  const my_float_t*
  normal() const
  { return A_N; }

  const my_float_t*
  point() const
  { 
    if(!A_P0_is_defined){
      printf("ERROR, P0 is not defined (in point)!!!!!!!!!!!!\n");
      return 0;
    }
    return A_P0; 
  }

private:
  void
  do_copy(const plane_t& src)
  {
    std::copy(src.A_N, src.A_N + 3, A_N);
    std::copy(src.A_P0, src.A_P0 + 3, A_P0);
    A_d = src.A_d;
    A_P0_is_defined = src.A_P0_is_defined;
  }

  my_float_t A_N[3];   //!< Unit normal to the plane
  my_float_t A_P0[3];    //!< Point on the plane  -- might be full of garbage
  my_float_t A_d;        //!< The d in the plane equation NX + d = 0
  bool A_P0_is_defined;
};

//! Computational analytic representation of a 3D arc
class arc_t{
public:

  //! Initialize the arc by a sphere & arc end points & arc mid point
  arc_t(const sphere_t& S, const my_float_t* end_pt0, const my_float_t* mid_pt,
        const my_float_t* end_pt1) : A_S(S)
  {
    //A_S = S;
    std::copy(end_pt0, end_pt0 + 3, A_end_pts);
    std::copy(end_pt1, end_pt1 + 3, A_end_pts + 3);

    std::copy(mid_pt, mid_pt + 3, A_mid_pt);
    unit_vector(A_in_dir, mid_pt, S.center());

    mid_point(end_pt0, end_pt1, A_chord_mid_pt);
  }

  //! Copy constructor
  arc_t(const arc_t& src) : A_S(src.A_S)
  {
    do_copy(src);
  }

  const arc_t&
  operator=(const arc_t& src)
  {
    if(this != &src){
      do_copy(src);
      A_S = src.A_S;
    }
    return *this;
  }

  //! Is point contained in arc?
  const bool
  contains(const my_float_t* pt) const
  {
    my_float_t my_dir[3];
    std::copy(pt, pt + 3, my_dir);
    my_axpy(3, -1.0, A_chord_mid_pt, 1, my_dir, 1);
    if(dot(my_dir, A_in_dir) >= 0.0) return true;
    return false;
  }

  //! Compute the intersection of self & the other arc with out modifying self
  //! or the other arc
  /*! Assumption:
   *    Both arcs lie on the same circle 
   */
  void
  intersection(const arc_t& other, std::vector<arc_t>* I_arcs)
  {
    std::vector<bool> self_contains_other_pts, other_contains_self_pts;
    self_contains_other_pts.push_back(contains(other.A_end_pts));
    self_contains_other_pts.push_back(contains(other.A_end_pts + 3));
    other_contains_self_pts.push_back(other.contains(A_end_pts));
    other_contains_self_pts.push_back(other.contains(A_end_pts + 3));

    // Case 1) Other arc is contained in self -- intersection is other_arc
    if(self_contains_other_pts[0] && self_contains_other_pts[1] &&
       !other_contains_self_pts[0] && !other_contains_self_pts[1])
      I_arcs->push_back(other);

    // Case 2) Other arc contains self -- intersection is self
    else if(other_contains_self_pts[0] && other_contains_self_pts[1] &&
            !self_contains_other_pts[0] && !self_contains_other_pts[1])
      I_arcs->push_back(*this);

    // Case 3) Other arc has no overlap with self -- no intersection
    else if(!self_contains_other_pts[0] && !self_contains_other_pts[1] &&
            !other_contains_self_pts[0] && !other_contains_self_pts[1]);
      // Pass as there is no intersection

    // Case 4) Other arc partially overlaps self -- intersection is the overlap

    // Case 4a) all points from each arc are contained in the other arc
    else if(self_contains_other_pts[0] && self_contains_other_pts[1] &&
            other_contains_self_pts[0] && other_contains_self_pts[1]){

      // idea -- move the chord for self along the in direction until
      // we reach an endpoint of the other arc (i.e. y = mx + b)
      my_float_t sq_dists[4];
      sq_dists[0] = dist_squared(A_end_pts, other.A_end_pts);
      sq_dists[1] = dist_squared(A_end_pts, other.A_end_pts + 3);
      sq_dists[2] = dist_squared(A_end_pts + 3, other.A_end_pts);
      sq_dists[3] = dist_squared(A_end_pts + 3, other.A_end_pts + 3);

      my_float_t tmp;
      int idx = argmin(sq_dists, 4, &tmp);
      int self_idx = static_cast<int>(std::floor(idx / 2.0));
      int other_idx = idx % 2;

      // Because we are taking the two closest points (1 end point from each
      // arc), we know that the in direction is the vector from the center of
      // the circle to midpoint of the corresponding chord
      my_float_t my_dir[3];
      compute_overlap(A_end_pts + 3*self_idx, other.A_end_pts + 3*other_idx,
                      I_arcs, my_dir);

      // Next, we use the other 2 points do define another chord and compute
      // the direction from the center to the midpoint of the chord.
      int A1_self_idx = (self_idx + 1) % 2; 
      int A1_other_idx = (other_idx + 1) % 2; 
      uint my_len = I_arcs->size();
      compute_overlap(A_end_pts + 3*A1_self_idx, 
                      other.A_end_pts + 3*A1_other_idx, I_arcs, my_dir);

      // If the second arc contains one of the end points of other_arc, we need 
      // to flip the direction.
      if(I_arcs->size() > my_len){
        arc_t& A1 = I_arcs->back();
        if(A1.contains(A_end_pts + 3*self_idx) || 
           A1.contains(other.A_end_pts + 3*other_idx)){
          std::copy(A_S.center(), A_S.center() + 3, A1.A_mid_pt);
          my_axpy(3, -1.0 * A_S.radius(), my_dir, 1, A1.A_mid_pt, 1);
        }
      }

    // Case 4b) One point from each arc is contained in the other arc
    }else{
      my_float_t my_end_pts[6];
      if(self_contains_other_pts[0])
        std::copy(other.A_end_pts, other.A_end_pts + 3, my_end_pts);
      else
        std::copy(other.A_end_pts + 3, other.A_end_pts + 6, my_end_pts);
      if(other_contains_self_pts[0])
        std::copy(A_end_pts, A_end_pts + 3, my_end_pts + 3);
      else
        std::copy(A_end_pts + 3, A_end_pts + 6, my_end_pts + 3);
      
      // We will ignore very small arcs or those that are computed to be
      // a point 
      if(dist_squared(my_end_pts, my_end_pts + 3) < my_float_eps) return;

      my_float_t mid_pt[3], my_dir[3];
      arc_mid_point(my_end_pts, my_end_pts + 3, mid_pt, my_dir);
      if(!contains(mid_pt) || !other.contains(mid_pt)){
        std::copy(A_S.center(), A_S.center() + 3, mid_pt);
         my_axpy(3, -1.0 * A_S.radius(), my_dir, 1, mid_pt, 1);
      }
      
      I_arcs->push_back(arc_t(A_S, my_end_pts, mid_pt, my_end_pts + 3));
    }
  }
  
  const bool
  operator==(const arc_t& other)
  {
    if(!(A_S == other.A_S)) return false;
    for(int i = 0; i < 6; ++i)
      if(A_end_pts[i] != other.A_end_pts[i]) return false;
    for(int i = 0; i < 3; ++i)
      if(A_mid_pt[i] != other.A_mid_pt[i]) return false;
    return true;
  }

  friend std::ostream& 
  operator<<(std::ostream &out, const arc_t& A)
  { 
    out << A.A_S << " " 
        << A.A_end_pts[0] << " " << A.A_end_pts[1] << " " 
        << A.A_end_pts[2] << " " 
        << A.A_mid_pt[0] << " " << A.A_mid_pt[1] << " " << A.A_mid_pt[2] << " " 
        << A.A_end_pts[4] << " " << A.A_end_pts[5] << " " 
        << A.A_end_pts[6];
    return out;
  }

  const my_float_t*
  end_pts() const
  { return A_end_pts; }

  //! named get_mid_pt since I am too lazy to change all the local variables
  //! with the name mid_pt
  const my_float_t*
  get_mid_pt() const
  { return A_mid_pt; }

  const my_float_t*
  in_dir() const
  { return A_in_dir; }
  
private:

  //! Compute the mid point of an arc
  void
  arc_mid_point(const my_float_t* end_pt0, const my_float_t* end_pt1,
                my_float_t* mid_pt, my_float_t* mid_pt_dir)
  {
    // Compute the midpoint of the chord
    mid_point(end_pt0, end_pt1, mid_pt);

    // Find the unit vector in the direction from the center to the chord
    // mid point
    unit_vector(mid_pt_dir, mid_pt, A_S.center());

    // Compute the arc mid point by pushing the center of the circle out
    // along the unit vector by a distance of the radius
    std::copy(A_S.center(), A_S.center() + 3, mid_pt);
    my_axpy(3, A_S.radius(), mid_pt_dir, 1, mid_pt, 1);
  }

  void 
  compute_overlap(const my_float_t* self_end_pt, const my_float_t* other_end_pt,
                  std::vector<arc_t>* I_arcs, my_float_t* my_dir)
  {
      my_float_t mid_pt[3];
      arc_mid_point(self_end_pt, other_end_pt, mid_pt, my_dir);
      if(dist_squared(self_end_pt, other_end_pt) >= my_float_eps)
        I_arcs->push_back(arc_t(A_S, self_end_pt, mid_pt, other_end_pt));
  }

  void do_copy(const arc_t& src)
  {
    std::copy(src.A_end_pts, src.A_end_pts + 6, A_end_pts);
    std::copy(src.A_mid_pt, src.A_mid_pt + 3, A_mid_pt);
    std::copy(src.A_in_dir, src.A_in_dir + 3, A_in_dir);
    std::copy(src.A_chord_mid_pt, src.A_chord_mid_pt + 3, A_chord_mid_pt);
  }
         
  sphere_t A_S;             //!< Sphere on which the arc lies
  my_float_t A_end_pts[6];  //!< End points of the arc
  my_float_t A_mid_pt[3];   //!< Mid point of the arc
  my_float_t A_in_dir[3];   //!< vector; S.center -> mid_pt
  my_float_t A_chord_mid_pt[3];  //!< Mid point of chord defined by end_pts
};


//! Computational representation of an analytical intersection circle
/*! NOTE: the contains function is from the sphere_t class
 *  
 *  NOTE: some functions from the sphere_t class might not be applicable to 
 *  the iCircle idea -- be careful.  As an example, project_point is not
 *  overloaded and will project to sphere not necessarily on the circle.
 */
class iCircle : public sphere_t{
public:
  iCircle(const my_float_t* C, my_float_t R, const my_float_t* N, 
          const std::vector<arc_t>& initial_arcs) 
    : sphere_t(C, R), A_plane(N, C)
  {
    if(initial_arcs.size())
    {
      A_full_circle = false;
      A_final_arcs = 
        std::vector<arc_t>(initial_arcs.begin(), initial_arcs.end());
    }else A_full_circle = true;
    
  }

  //!< Simple copy constructor
  iCircle(const iCircle& src) : sphere_t(src), A_plane(src.A_plane)
  {
    do_copy(src);
  }

  //!< Simple assignment function (operator)
  const iCircle& 
  operator=(const iCircle& src)
  {
    if(this != &src){
      do_copy(src);
      A_plane = src.A_plane;
      sphere_t::operator= (src);
    }
    return *this;
  }

  //!< Remove the part of self that overlaps with the other iCircle
  void
  remove_overlap(const iCircle& other)
  {
    // If the two iCircles do not intersect, we have nothing to do
    my_float_t D2 = dist_squared(center(), other.center());
    my_float_t center_dist = radius() + other.radius();
    if(D2 >= center_dist * center_dist) return;

    // Determine if any of self's arcs are full contained in other's volume
    std::vector<arc_t> tmp_arcs;
    for(uint i = 0; i < A_final_arcs.size(); ++i){
      arc_t& AA = A_final_arcs[i];
      // We test the mid point to rule out the cases where the arc angle is
      // greater than pi
      if(!other.contains(AA.end_pts()) || !other.contains(AA.end_pts() + 3) ||
         !other.contains(AA.get_mid_pt()))
        tmp_arcs.push_back(AA);
    }
    if(!A_full_circle and tmp_arcs.size() < 1) return;

    // Check if self is fully contained in the other sphere
    if(other.contains(center())){
      my_float_t proj_pt[3];
      other.sphere_t::project_point(center(), proj_pt);
      if(!contains(proj_pt)){
        A_full_circle = false;
        A_final_arcs.clear();
        return;
      }
    }

    // Determine the line (chord) where the two circles intersect
    // Plane intersection is the line Y = mX + b
    my_float_t m[3], b[3];
    A_plane.intersection(other.A_plane, m, b);
    //intersect_planes(A_N, A_d, other.A_N, other.A_d, m, b);

    // Determine intersection between line defined by m & b and self
    my_float_t pts[6];
    int npts = intersection(m, b, pts); 
    if(!npts) return;

    // Get the 2 arcs for self -- both the remaining or kept arc and the
    // removed arc or arc that falls inside the other sphere
    my_float_t mid_pt[3], my_dir[3], mid_pts[6], sq_dists[2];
    mid_point(pts, pts + 3, mid_pt);
    unit_vector(my_dir, mid_pt, center());
    std::copy(center(), center() + 3, mid_pts);
    std::copy(center(), center() + 3, mid_pts + 3);
    my_axpy(3, radius(), my_dir, 1, mid_pts, 1);
    my_axpy(3, -1.0 * radius(), my_dir, 1, mid_pts + 3, 1);
    sq_dists[0] = dist_squared(other.center(), mid_pts);
    sq_dists[1] = dist_squared(other.center(), mid_pts + 3);

    my_float_t* keep_mid_pt = 0;
    my_float_t* rm_mid_pt = 0;
    if(sq_dists[0] < other.squared_radius() &&
       sq_dists[1] < other.squared_radius()){
      err_msg("sphere.H", "iCircle.remove_overlap()",
              "Entire circle is contained in other circle");
      return;
    }else if(sq_dists[0] < other.squared_radius()){
      keep_mid_pt = mid_pts + 3;
      rm_mid_pt = mid_pts;
    }else if(sq_dists[1] < other.squared_radius()){
      keep_mid_pt = mid_pts;
      rm_mid_pt = mid_pts + 3;
    }else{
      err_msg("sphere.H", "iCircle.remove_overlap()",
              "Circles do not intersect");
      return;
    }

    if(dist_squared(pts, pts + 3) < 1E-14){
      // Keep arc is very small, the contribution of self to total difference is
      // negligible
      if(dist_squared(keep_mid_pt, pts) < 1E-14){
        A_full_circle = false;
        A_final_arcs.clear();
      // Remove arc is very small, the contribution of other to the total
      // difference is negligible
      }else if(dist_squared(rm_mid_pt, pts) < 1E-14) return;
    }

    // In direction of overlap is the direction from the center to midpoint of
    // chord
    arc_t keep_arc(*this, pts, keep_mid_pt, pts + 3);
   
    // Update the arc vector for self
    if(A_full_circle){
      A_full_circle = false;
      A_final_arcs.push_back(keep_arc);
    }else{
      std::vector<arc_t> new_final_arcs;
      for(uint i = 0; i < A_final_arcs.size(); ++i)
        keep_arc.intersection(A_final_arcs[i], &new_final_arcs);
      A_final_arcs = new_final_arcs;
    }
  }

  //! Project point to the intersection circle
  bool
  project_point(const my_float_t* pt, my_float_t* proj_pt, my_float_t *sq_dist,
                const my_float_t tol=1.0)
  {
    my_float_t sq_tol = tol*tol;
    my_float_t dist_to_plane = A_plane.signed_dist(pt);

    // 1st, project pt to plane
    //my_float_t proj_pt[3];
    std::copy(pt, pt + 3, proj_pt);
    my_axpy(3, -1.0 * dist_to_plane, A_plane.normal(), 1, proj_pt, 1);
 
    // 2nd, project the pt in plane to the circle (sphere restricted to the
    // plane
    my_float_t in_plane_dir[3];
    unit_vector(in_plane_dir, proj_pt, center());
    std::copy(center(), center() + 3, proj_pt);
    my_axpy(3, radius(), in_plane_dir, 1, proj_pt, 1);
   
    *sq_dist = dist_squared(proj_pt, pt);
    if(*sq_dist > sq_tol) return false;
     
    if(!A_full_circle){
      point_type pt_on_A;
      std::vector<point_type> arc_pts;
      std::vector<my_float_t> pts_sq_dists;
      std::vector<geometry::arc_t>::const_iterator A;
      for(A = A_final_arcs.begin(); A < A_final_arcs.end(); ++A){
/*
        std::cout << "  arc end points " << A->end_pts()[0] << " "
                  << A->end_pts()[1] << " " << A->end_pts()[2] << "   "
                  << A->end_pts()[3] << " " << A->end_pts()[4] << " "
                  << A->end_pts()[5] << "\n";
        std::cout << "  arc in dir: " << A->in_dir()[0] << " "
                  << A->in_dir()[1] << " " << A->in_dir()[2] << "\n";
*/
        // If the projected point is in the arc, keep it
        if(A->contains(proj_pt)){
//          std::cout << "  Arc contains pt\n";
          std::copy(proj_pt, proj_pt + 3, pt_on_A.pt);
          arc_pts.push_back(pt_on_A);
          pts_sq_dists.push_back(*sq_dist);
        // Otherwise choose the closest end point of the arc
        }else{
//          std::cout << "  Arc does not contain pt\n";
          my_float_t d2z[] = { dist_squared(A->end_pts(), pt),
                               dist_squared(A->end_pts() + 3, pt) };
          my_float_t min_d2;
          int idx = argmin<my_float_t>(d2z, 2, &min_d2);
          if(min_d2 <= sq_tol){
            std::copy(A->end_pts() + 3*idx, A->end_pts() + 3*idx + 3,
                      pt_on_A.pt);
            arc_pts.push_back(pt_on_A);
            pts_sq_dists.push_back(min_d2);
          }
        }
      }

      // Get the closest point of all computed points
      if(pts_sq_dists.size()){
        my_float_t min_dist = pts_sq_dists[0];
        uint idx = 0;
        for(uint i = 1; i < pts_sq_dists.size(); ++i){
          if(pts_sq_dists[i] < min_dist){
            min_dist = pts_sq_dists[i];
            idx = i;
          }
        }
        std::copy(arc_pts[idx].pt, arc_pts[idx].pt + 3, proj_pt);
        *sq_dist = min_dist; 
      }else return false;
    }
    return true;
  }

  //! Returns true if this intersection circle does not intersect with the 
  //! cutting plane or any other intersection circle
  bool
  full_circle() const
  { return A_full_circle; }

  //! Get the number of final arcs
  uint
  num_final_arcs() const
  { return A_final_arcs.size(); }

  //! Get a const iterator to the first final arc
  std::vector<arc_t>::const_iterator
  final_arcs_begin() const
  { return A_final_arcs.begin(); }

  //! Get a const iterator to one past the last final arc
  std::vector<arc_t>::const_iterator
  final_arcs_end() const
  { return A_final_arcs.end(); }

  //! Structure to allow storage of 3D pts in a vector
  struct point_type{
    my_float_t pt[3];
  };

private:
  //!< Simple copy function
  void
  do_copy(const iCircle& src)
  {
    A_full_circle = src.A_full_circle;
    A_final_arcs = src.A_final_arcs;
  }

  bool A_full_circle; 
  plane_t A_plane;
  std::vector<arc_t> A_final_arcs;
};



} // end namespace utils
} // end namespace SimSite3D


#endif
