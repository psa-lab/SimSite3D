#ifndef TRIANGLE_MESH_SPHERE_HEADER_INCLUDED
#define TRIANGLE_MESH_SPHERE_HEADER_INCLUDED

#include <types.H>
#include <line.H>
#include <HbondSurfaces.H>

namespace SimSite3D{
namespace geometry{

//! forward declaration of vertex_t for half_edge_t
class vertex_t;
typedef dir_point_storage<vertex_t>::iterator vertex_vi;
typedef dir_point_storage<vertex_t>::const_iterator vertex_vci;

//! Simple half edge structure
struct half_edge_t{
  std::vector<half_edge_t>::iterator next;
  vertex_vci tail;
};
typedef std::vector<half_edge_t> half_edge_vec;
typedef half_edge_vec::iterator half_edge_vi;
typedef half_edge_vec::const_iterator half_edge_vci;

//! Might be better named "node"?
class vertex_t : public dir_point_t{
public:
  vertex_t(alloc_t a = ALLOC_POSITION) : dir_point_t(a)
  {

  }

  vertex_t(const vertex_t& v) : dir_point_t(v)
  {

  }

  const vertex_t&
  operator=(const vertex_t& v)
  {
    if(&v != this){
      dir_point_t::operator= (v);
    }
    return *this;
  }

  ~vertex_t()
  { ; }


  void
  setup_one_ring(half_edge_vci beg_in, half_edge_vci end_in)
  {

    // Get this list sorted to see if we are missing a boundary?
    for(half_edge_vci hedge = beg_in; hedge != end_in; ++hedge)
      if(this == &(*(hedge->tail))){

        std::list<half_edge_vci>::iterator ring_iter;
        for(ring_iter = A_one_ring_edges.begin(); 
            ring_iter != A_one_ring_edges.end(); ++ ring_iter){
          if(hedge->next->next->tail == (*ring_iter)->tail){
            break; 
          }
        }
        A_one_ring_edges.insert(ring_iter, hedge->next);
      }
    std::cout << "num edges in 1 ring of (" << pos[0] << " " << pos[1] << " "
              << pos[2] << "): " << A_one_ring_edges.size() << std::endl;

    std::list<half_edge_vci>::iterator ring_iter;
    for(ring_iter = A_one_ring_edges.begin(); 
        ring_iter != A_one_ring_edges.end(); ++ ring_iter){
      std::cout << (*ring_iter)->tail->pos[0] << " "
                << (*ring_iter)->tail->pos[1] << " "
                << (*ring_iter)->tail->pos[2] << "\n";
    }
  }

  bool
  inside_one_ring(const my_float_t* pt)
  {
    // intersection of line segments by closest points
    std::list<half_edge_vci>::const_iterator edge_list_iter;
    for(edge_list_iter = A_one_ring_edges.begin();
        edge_list_iter != A_one_ring_edges.end(); ++edge_list_iter){
      const half_edge_vci& edge = *edge_list_iter;
      geometry::line_t(edge->tail->pos, edge->next->tail->pos);
    }
    // function is not finished -- good way to test the triangle & 
    // one ring code -- they should produce same results (up to numerical
    // accuracy)
    return false;
  }

private:
  std::list<half_edge_vci> A_one_ring_edges;
};

//! Requires vertex_t class and should be independent of vertex placement
/*!
 *  Since this class is supposed to be independent of the vertices positions
 *  and normals, it is currently the users reponsibility to insure that the
 *  vertices and triangle are consistent. 
 */
class triangle_t{
public:
  triangle_t(vertex_vci a, vertex_vci b, vertex_vci c)
  {
    A_vertices.push_back(a);
    A_vertices.push_back(b);
    A_vertices.push_back(c);
    my_float_t N[3];
    A_area = face_normal(N);

    // We want the order of the vertices to be right handed with respect to the
    // out facing normal -- this could be troublesome if we have a messy mesh 
//    std::cout << "dot product is: " << dot(a->dir, N) << std::endl;
    if(dot(a->dir, N) < 0.0){
      A_vertices[1] = c;
      A_vertices[2] = b;
    }
  }

  ~triangle_t()
  { ; }

  triangle_t(const triangle_t& src)
  {
    do_copy(src);
  }

  const triangle_t&
  operator=(const triangle_t& src)
  {
    if(this != &src) do_copy(src);
    return *this;
  }

  //! 3D version that drops the coordinate with the least variance and calls
  //! the 2D method
  bool
  contains(const my_float_t* pt) const
  {
#if 0
    std::cout << "triangle: \n";
    for(uint i = 0; i < 3; ++i)
      std::cout << A_vertices[i]->pos[0] << " " << A_vertices[i]->pos[1] << " "
                << A_vertices[i]->pos[2] << "\n";
    std::cout << "point: " << pt[0] << " " << pt[1] << " " << pt[2] << "\n";
#endif

    my_float_t U[3], V[3], N[3];
    vector(3, A_vertices[1]->pos, A_vertices[0]->pos, U);
    vector(3, A_vertices[2]->pos, A_vertices[0]->pos, V);
    cross(U, V, N);

    my_float_t maxval;
    const int drop_idx = argmaxabs<my_float_t>(N, 3, &maxval);

//    std::cout << "3D normal: " << N[0] << " " << N[1] << " " << N[2] << "\n";
//    U[drop_idx] = 0.0;
//    V[drop_idx] = 0.0;
//    cross(U, V, N);
//    std::cout << "2D normal: " << N[0] << " " << N[1] << " " << N[2] << "\n";

    my_float_t vals[8];
    my_float_t* v = vals;
    for(int i = 0; i < 3; ++i){
      if(drop_idx == i) continue;

      *v = pt[i]; 
      ++v;
      // For right handed coordinate systems we need to handle the Y direction
      // opposite from the others (i.e. [1,0,0] X [0,0,1] = [0, -1, 0]
      if((drop_idx != 1 && N[drop_idx] > 0.0) || 
         (drop_idx == 1 && N[drop_idx] < 0.0))
        for(int j = 0; j < 3; ++j, ++v) *v = A_vertices[j]->pos[i];
      else{
//        std::cout << "swapping point order\n";
        *v = A_vertices[0]->pos[i];
        ++v;
        *v = A_vertices[2]->pos[i];
        ++v;
        *v = A_vertices[1]->pos[i];
        ++v;
      }
    }
    return contains(vals, vals + 4);
  }

  //! Compute the normal to the out facing triangle face
  my_float_t
  face_normal(my_float_t* N)
  {
    my_float_t U[3], V[3];
    vector(3, A_vertices[1]->pos, A_vertices[0]->pos, U);
    vector(3, A_vertices[2]->pos, A_vertices[0]->pos, V);
    cross(U, V, N);
//      std::cout << U[0] << " " << U[1] << " " << U[2] << "\n"
//                << V[0] << " " << V[1] << " " << V[2] << "\n"
//                << N[0] << " " << N[1] << " " << N[2] << std::endl;
    my_float_t rv = normalize(N);
//    std::cout << "normal magnitude: " << rv << std::endl;
    return rv;
  }

  std::vector<vertex_vci>::const_iterator
  vertices_begin() const
  { return A_vertices.begin(); }

  std::vector<vertex_vci>::const_iterator
  vertices_end() const
  { return A_vertices.end(); }

  my_float_t
  area() const
  { return A_area; }

private:
  //! 2D version
  /*! From Usenet discussion on Barycentric Coordinates (Joseph O'Rourke,
   * 1992)
   * NOTE: for now, this function requires the points be ordered in a 
   * counterclockwise fashion
   */ 
  bool
  contains(const my_float_t *x, const my_float_t *y) const
  {
    my_float_t X[3], Y[3];
//    std::cout << *x << " " << *y << std::endl;
    for(size_t i = 0; i < 3; ++i){
      X[i] = x[i+1] - *x;
      Y[i] = y[i+1] - *y;
//      std::cout << x[i+1] << " " << y[i+1] << std::endl;
    }
    my_float_t b[3];
    cross(X, Y, b);

//    std::cout << "b: " << b[0] << " " << b[1] << " " << b[2] << std::endl;

    // Need to determine on which side the of the triangle the point lies
    for(size_t i = 0; i < 3; ++i)
      if(b[i] < 0) return false;
    return true;
  }

  void
  do_copy(const triangle_t& src)
  {
    for(uint i = 0; i < src.A_vertices.size(); ++i)
      A_vertices.push_back(src.A_vertices[i]);
    A_area = src.A_area;
  }

  std::vector<vertex_vci> A_vertices;     //!< 3 vertices of the triangle
  my_float_t A_area;           //!< Area of the triangle

#if 0
  // what if the vertices move?
  my_float_t A_face_normal[3]; //!< Out facing normal to the face of the triangle
  my_float_t A_xvals[3];       //!< 2D X coordinates of the triangle vertices
  my_float_t A_yvals[3];       //!< 2D Y coordinates of the triangle vertices
#endif
};

typedef std::vector<triangle_t> triangle_vec;
typedef triangle_vec::iterator triangle_vi;
typedef triangle_vec::const_iterator triangle_vci;

//! Straightforward class to work with a mesh sphere that has portions of it
//! occluded
class TriMeshSphere{
public:

  //! More descriptive method to denote initialization type for spherical mesh
  typedef enum sphere_init_type{
    FULL_SPHERE,
    SPHERICAL_CAP
  }sphere_init_t;

  //! Default constructor -- allows class to be an object in another class
  //! before one knows the values for the sphere
  TriMeshSphere()
  {
    init();
  }

  //! Constructor to generate a mesh representation of a sphere or part of 
  //! sphere depending on the init_obj variable.
  TriMeshSphere(const my_float_t* R, const my_float_t* T, const my_float_t* dir,
                const sphere_init_t init_obj, const my_float_t radius = 3.0,
                const uint level = 2);

  TriMeshSphere(const TriMeshSphere& src)
  {
    init();
    do_copy(src);
  }

  TriMeshSphere& operator=(const TriMeshSphere& src)
  {
    if(this != &src)
    {
      do_copy(src);
    }
    return *this;
  }
  
  ~TriMeshSphere()
  { init(); }

  void
  transform(const my_float_t* R, const my_float_t* T)
  { A_nodes.transform(R, T); }

  void
  inverse_transform(const my_float_t* R, const my_float_t* T)
  { A_nodes.inverse_transform(R, T); }

  void
  revert()
  { A_nodes.revert(); }

  void
  adjust_points(hbond_surface_t& surf)
  {
    A_adjust_points(surf, A_triangles, A_num_tris);
  }

  my_float_t
  complementary_surface_area(const std::map<vertex_vci, bool>& kept_vertices)
  const;

  //! Get a constant iterator to the first vertex in the vector
  const vertex_vci
  verts_begin() const
  { return A_nodes.begin(); }

  //! Get a constant iterator to the end of the vertex vector
  const vertex_vci
  verts_end() const
  { return A_nodes.end(); }

  //! Get the number of vertices in the surface
  const uint
  num_verts() const
  { return A_nodes.size(); }

  //! Write out the triangles in the mesh
  void
  write(std::ostream &out, const interactionType act_type, const char delim='|')
  const;

  void
  write_msms_headers(std::ostream &vert_out, std::ostream &face_out) const;

  void
  write_msms_cap(std::ostream &vert_out, std::ostream &face_out, 
                 const interactionType act_type) const;

private:
  void
  init()
  {
    A_local_points = 0;
    A_triangles = 0;
    A_npts = 0;
    A_num_tris = 0;
  }

  //! The function that does all the work for copying the class
  void
  do_copy(const TriMeshSphere& src);

  void
  set_global_pos(const my_float_t* R, const my_float_t* T, 
                 const my_float_t* bond_dir)
  {
//    std::cout << "Time to set the global pos of the cap" << std::endl;
    my_float_t local_bond_dir[3];
    my_gemv(3, 3, 1.0, R, 3, bond_dir, 1, 0.0, local_bond_dir, 1);
    my_float_t A[] = {1.0, 0.0, 0.0,
                      0.0, 1.0, 0.0,
                      0.0, 0.0, 1.0};
    my_float_t B[9], AB[9], ABR[9];
//    std::cout << "local bond dir: " << local_bond_dir[0] << " " 
//              << local_bond_dir[1] << " " << local_bond_dir[2] << std::endl;

    // Check if we need to rotate about Y first
    my_float_t tmp_dir[3];
    if(local_bond_dir[2] < -1E-03 || local_bond_dir[2] > 1E-03){
      my_float_t sin_t = local_bond_dir[2];
      my_float_t cos_t = std::sqrt(1.0 - sin_t*sin_t);
      my_float_t X_dir[] = {1.0, 0.0, 0.0};

      A[0] = cos_t;
      A[1] = 0.0;
      A[2] = sin_t;
      A[3] = 0.0;
      A[4] = 1.0;
      A[5] = 0.0;
      A[6] = -1.0*sin_t;
      A[7] = 0.0;
      A[8] = cos_t;
      my_gemm(1, 3, 3, 1.0, X_dir, 3, A, 3, tmp_dir, 3, 0.0);
    }else{
      tmp_dir[0] = 1.0;
      tmp_dir[1] = 0.0;
      tmp_dir[2] = 0.0;
    }

    // Compute the angle of the projection of both vectors in the Z=0 plane
//    std::cout << "tmp_dir[X] = " << tmp_dir[0] << std::endl;
    my_float_t cos_t = local_bond_dir[0] / tmp_dir[0];
//    std::cout << "cos(t): " << cos_t << std::endl;
    my_float_t sin2_t = 1.0 - cos_t*cos_t;
    my_float_t sin_t = 0.0;
    if(sin2_t > 1E-8) sin_t = std::sqrt(sin2_t);
    if(local_bond_dir[1] < 0.0) sin_t *= -1.0;

    B[0] = cos_t;
    B[1] = sin_t;
    B[2] = 0.0;
    B[3] = -1.0*sin_t;
    B[4] = cos_t;
    B[5] = 0.0;
    B[6] = 0.0;
    B[7] = 0.0;
    B[8] = 1.0;
    my_gemm(3, 3, 3, 1.0, A, 3, B, 3, AB, 3, 0.0);
//    std::cout << "AB: " << AB[0] << " " << AB[1] << " " << AB[2] << "\n"
//              << "     " << AB[3] << " " << AB[4] << " " << AB[5] << "\n"
//              << "     " << AB[6] << " " << AB[7] << " " << AB[8] << "\n"
//              << std::endl;
    my_gemm(3, 3, 3, 1.0, AB, 3, R, 3, ABR, 3, 0.0);

//    std::cout << "ABR: " << ABR[0] << " " << ABR[1] << " " << ABR[2] << "\n"
//              << "     " << ABR[3] << " " << ABR[4] << " " << ABR[5] << "\n"
//              << "     " << ABR[6] << " " << ABR[7] << " " << ABR[8] << "\n"
//              << std::endl;


    A_nodes.transform(ABR, T);
    A_nodes.set_current_positions_and_directions_as_original();
  }

  void 
  init_storages(const my_float_t* pts, const uint npts, const uint* triangles,
                const uint ntris, const my_float_t radius);

  void 
  A_adjust_points(hbond_surface_t& surf, const uint* tri_vert_idz, 
                  const uint num_tris);


  void
  add_triangle(const uint* idz);

  void
  setup_vert_to_deltas_map();

  const my_float_t* A_local_points;
  uint A_npts;
  const uint* A_triangles;
  uint A_num_tris;

  dir_point_storage<vertex_t> A_nodes;
  triangle_vec A_deltas;
  std::map<vertex_vci, std::vector<triangle_vci> > A_verts_to_deltas_map;
  half_edge_vec A_half_edges;

  // Points for a 1.0 (A) sphere.  Unlike the hbond caps, the radius of the
  // sphere will vary depending on which metal we are modeling
  static const my_float_t A_level1_sphere_points[];
  static const uint A_level1_sphere_npts;
  static const uint A_level1_sphere_triangles[];
  static const uint A_level1_num_sphere_tris;

  static const my_float_t A_level2_sphere_points[];
  static const uint A_level2_sphere_npts;
  static const uint A_level2_sphere_triangles[];
  static const uint A_level2_num_sphere_tris;

  // Points for a cap on a 1.0 (A) sphere
  static const my_float_t A_level2_cap_points[];
  static const uint A_level2_cap_npts;
  static const uint A_level2_cap_triangles[];
  static const uint A_level2_num_cap_tris;
};



} // end namespace geometry
} // end namespace SimSite3D

#endif
