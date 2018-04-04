#ifndef TRANSFORMABLE_TRIMESH_HEADER_INCLUDED
#define TRANSFORMABLE_TRIMESH_HEADER_INCLUDED

#include <SimpleTrimeshTwo.H>

namespace SimSite3D{
namespace geometry{

class TransformableTrimesh : public SimpleTrimeshTwo{
public:

  //! This mesh may have the vertex positions moved (and other edits as well)
  /*! Unlike ImmovableTrimesh, this mesh data structure does not make any
   *  assumptions that are not transformation invariant.
   */
  TransformableTrimesh(const std::string fname,
                       const bool check_vertex_order=true)
    : SimpleTrimeshTwo(fname)
  { 
    if(!read_face_file(fname, check_vertex_order, false)){
      set_fail_flag();
      return;
    }
  }
  
  ~TransformableTrimesh()
  {
  
  }

  //! Get a non-constant pointer to the beginning of the vertices array
  /*! I am unsure of the "correct" way to tell the g++ compiler which of
   * the vertices begin functions we want -- it always chooses the unconstant
   * one with the same name
   */
  my_float_t*
  mutable_vertices_begin()
  { return SimpleTrimeshTwo::mutable_vertices_begin(); }

#if 0
  bool
  compare(const my_float_t *other_vert_begin, const uint num_other_vert,
          my_float_t *closest_pts, my_float_t *dists,
          const my_float_t max_dist = my_float_max) const
  {
    std::cout << "Using TransformableTrimesh::compare\n";
    std::fill(dists, dists + num_other_vert, my_float_max);
    std::fill(closest_pts, closest_pts + 3*num_other_vert, my_float_max);

    my_float_t *cp = closest_pts;
    my_float_t *d = dists;
    const my_float_t *other_vert_end = other_vert_begin + 3 * num_other_vert;
    for(const my_float_t *v = other_vert_begin; v < other_vert_end; v += 3){
      face_vci face;
      for(face = faces_begin(); face < faces_end(); ++face){
        compute_close_point(v, face->vertices[0], face->vertices[1], 
                            face->vertices[2], d, cp);
      }

      d += 1;
      cp += 3; 
    } 
    return true;
  }
#endif

};

#if 0
class TransformableTrimesh : public SimpleTrimeshTwo{
public:

  typedef struct{
    const my_float_t *vertices[3];
    const my_float_t *normals[3];
    my_float_t area;
  }triangle_t;

  typedef std::vector<triangle_t> face_vec;
  typedef face_vec::iterator face_vi;
  typedef face_vec::const_iterator face_vci;
  typedef std::map<const my_float_t*, std::vector<face_vci> > vert_face_map;
  typedef vert_face_map::iterator vert_face_mi;
  typedef vert_face_map::const_iterator vert_face_mci;

  //! This mesh may have the vertex positions moved (and other edits as well)
  /*! Unlike ImmovableTrimesh, this mesh data structure does not make any
   *  assumptions that are not transformation invariant.
   */
  TransformableTrimesh(const std::string fname)
    : SimpleTrimeshTwo(fname)
  { 
    if(!read_face_file(fname, false, false)){
      set_fail_flag();
      return;
    }
    compute_areas(); 
  }
  
  ~TransformableTrimesh()
  {
  
  }

  //! Find closest point on surface mesh for the given point (pt)
  /*! This method is inefficient since it checks each triangle in the mesh
   *  to find the closest point for each input vertex.  This function exists
   *  to help with prototyping and as a fallback, simpler method to check
   *  that the more advanced methods still compute the closest point.
   *  In particular, we do the most straightforward and naive check of each
   *  triangle without first checking vertex distances, etc.
   */
  bool
  compare(const my_float_t *other_vert_begin, const uint num_other_vert,
          my_float_t *closest_pts, my_float_t *dists,
          const my_float_t max_dist = my_float_max) const
  {
    std::fill(dists, dists + num_other_vert, my_float_max);
    std::fill(closest_pts, closest_pts + 3*num_other_vert, my_float_max);
   
    my_float_t *cp = closest_pts;
    my_float_t *d = dists;
    const my_float_t *other_vert_end = other_vert_begin + 3 * num_other_vert;
    for(const my_float_t *v = other_vert_begin; v < other_vert_end; v += 3){
      face_vci face;
      for(face = A_triangles.begin(); face < A_triangles.end(); ++face){
        compute_close_point(v, face->vertices[0], face->vertices[1], 
                            face->vertices[2], d, cp);
      }
      d += 1;
      cp += 3; 
    } 
    return true;
  }

private:
  //! Initialize the triangle storage
  void
  init_triangle_data_structure(const uint num_vert, const uint num_tri)
  {
    A_triangles.reserve(num_tri);
  }

  //! Add a triangle to the triangle storage
  void
  update_triangle_data_structure(const uint tri_idx, const uint *vert_idz,
                                 const my_float_t **vertices,
                                 const my_float_t **normals)
  {
    triangle_t tmp;
    std::copy(vertices, vertices + 3, tmp.vertices);
    std::copy(normals, normals + 3, tmp.normals);
    A_triangles.push_back(tmp);  
  }

  //! Compute the area of the triangles
  void
  compute_areas()
  {
    my_float_t U[3], V[3], U_x_V[3];
    for(face_vi face = A_triangles.begin(); face < A_triangles.end(); ++face){
      std::copy(face->vertices[0], face->vertices[0] + 3, U);
      std::copy(face->vertices[2], face->vertices[2] + 3, V);
      my_axpy(3, -1.0, face->vertices[1], 1, U, 1);
      my_axpy(3, -1.0, face->vertices[1], 1, V, 1);
      cross(U, V, U_x_V);
      face->area = 0.5 * normalize(U_x_V);
    }
  }

  face_vec A_triangles;              //!< Triangle storage
  vert_face_map A_map_vert_to_faces; //!< Map of each vertex to all its triangles
};
#endif

} // end namespace geometry
} // end namespace SimSite3D
#endif
