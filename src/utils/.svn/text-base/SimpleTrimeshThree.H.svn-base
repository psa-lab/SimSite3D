#ifndef SIMPLE_TRIMESH_THREE_HEADER_INCLUDED
#define SIMPLE_TRIMESH_THREE_HEADER_INCLUDED

#include <FaceAttrib.H>

namespace SimSite3D{
namespace geometry{

class Edge{
public:
  Edge()
  { init(); }

  Edge(VertAttrib::vi V0, VertAttrib::vi V1, FaceAttrib::vi F0, 
       FaceAttrib::vi F1)
  {
    init();

    // We want the order to be similar to that in the faces --
    // If we start rotation about V0 with face F0, we want F1 to be
    // the last face in the (circular linked) list before F1
    A_V[0] = V0;
    A_V[1] = V1;
    if(F0 != FaceAttrib::NULL_VI){
      if(F0->next(V0) == F1){
        A_nbrs[0] = F1;
        A_nbrs[1] = F0;
      }else{
        A_nbrs[0] = F0;
        A_nbrs[1] = F1;
      }
    }else if(F1 != FaceAttrib::NULL_VI){
      if(F1->next(V1) == F0){
        A_nbrs[0] = F1;
        A_nbrs[1] = F0;
      }else{
        A_nbrs[0] = F0;
        A_nbrs[1] = F1;
      }
    }
  }

  Edge(const Edge& E)
  {
    init();
    do_copy(E);
  }

  const Edge&
  operator=(const Edge& E)
  {
    if(this != &E) do_copy(E);
    return *this;
  }

  ~Edge()
  { ; }

  void
  get_verts(VertAttrib::vi *A, VertAttrib::vi *B) const
  {
    *A = A_V[0];
    *B = A_V[1];
  }

  void
  get_faces(FaceAttrib::vi *A, FaceAttrib::vi *B) const
  {
    *A = A_nbrs[0];
    *B = A_nbrs[1];
  }

private:
  void 
  init();
#if 0
  {
    std::fill(A_V, A_V + 2, VertAttrib::NULL_VI);
    std::fill(A_nbrs, A_nbrs + 2, FaceAttrib::NULL_VI);
  }
#endif

  void
  do_copy(const Edge& E)
  {
    std::copy(E.A_V, E.A_V + 2, A_V);
    std::copy(E.A_nbrs, E.A_nbrs + 2, A_nbrs);
  }

  VertAttrib::vi A_V[2];      //!< End points of the line segment (edge)
  FaceAttrib::vi A_nbrs[2];  //!< Faces that share the edge
};


//! Simple surface class built to contain basic features of a triangle mesh
class SimpleTrimeshThree{
public:

  //! Default constructor -- initializes some class variables
  SimpleTrimeshThree();

  ~SimpleTrimeshThree()
  { ; }

  bool
  read_vert_file(const std::string fname, const bool sample_surface = false);

  bool
  read_face_file(const std::string fname, std::vector<Edge> *edges,
                 const bool check_vertex_order, 
                 const bool compute_max_dist_to_centroid);

  //! Check if class has failed -- if it has, this class's use at best will be 
  //! undefined
  const bool
  fail() const
  { return A_fail; }

  const size_t
  num_verts() const
  { return A_vertices.size(); }

  // This is public for now -- testing
  //! Compute the energy at each vertex and add to a priority queue based
  //! on the estimated change in energy
  void
  initialize_opt(std::vector<Edge> &edges);

  const my_float_t
  max_dist_to_centroid() const
  { return A_max_dist_to_cent; }

  VertAttrib::vci
  vertices_begin() const
  { return A_vertices.begin(); }

  VertAttrib::vci
  vertices_end() const
  { return A_vertices.end(); }

  FaceAttrib::vci
  faces_begin() const
  { return A_faces.begin(); }

  FaceAttrib::vci
  faces_end() const
  { return A_faces.end(); }

protected:
  void
  set_fail_flag()
  { A_fail = true; }

private:
  //! Get all the faces for the given vertex
  /*! This method does not add the starting face 
   * First this method will check all the faces that can be found
   * via the face->next links.  If a "null" face is found,
   * it will check all the faces via the face->prev links.
   *
   * This method will fail if the mesh is of poor quality and all of the
   * faces cannot be found by both next and prev traversal
   */
  bool
  get_faces(const VertAttrib::vi &central_vert,
            const FaceAttrib::vi &start_face,
            std::vector<FaceAttrib::vi> *faces_out,
	    const bool skip_first_prev = false)
  {
    if(start_face == FaceAttrib::NULL_VI) return false;
    FaceAttrib::vi F = start_face->next(central_vert);
    for( ; F != start_face && F != FaceAttrib::NULL_VI;  
        F = F->next(central_vert))
      faces_out->push_back(F);

    // Try the prev links if there is a "gap"
    if(F == FaceAttrib::NULL_VI){
      F = start_face->prev(central_vert); 
      if(F != FaceAttrib::NULL_VI && skip_first_prev) F = F->prev(central_vert);
      for( ; F != start_face && F != FaceAttrib::NULL_VI;  
          F = F->prev(central_vert))
        faces_out->push_back(F);
    }
  }

  //! This method will fail in the case of a "bad" mesh
  void
  get_faces(const Edge &E, std::vector<FaceAttrib::vi> *faces_out)
  {
    VertAttrib::vi end_pts[2];
    FaceAttrib::vi edge_faces[2];
    E.get_verts(&end_pts[0], &end_pts[1]);
    E.get_faces(&edge_faces[0], &edge_faces[1]);

    // get_faces does not add the starting face to the vector
    // this way faces do not get counted twice
    get_faces(end_pts[0], edge_faces[0], faces_out);
    get_faces(end_pts[1], edge_faces[1], faces_out);
  }

  bool
  get_faces(FaceAttrib::vi central_face, 
            std::vector<FaceAttrib::vi> *faces_out)
  {
    faces_out->push_back(central_face);

    VertAttrib::vi central_verts[3];
    central_face->get_vertices(central_verts);

    // The nice idea gets messier here when we try to deal with 
    // special cases.  Ideally, we would traverse each circular array
    // starting with the next face after the central face and not add
    // the starting face for that traversal as it is covered by the 
    // (next, previous) traversel.  When we hit boundaries, we use
    // prev as well to see how far we can get, if the next face does not
    // exist then we need to start with the central face and we do not
    // want get_faces to add central face
    for(int i = 0; i < 3; ++i){
      FaceAttrib::vi next_face = central_face->next(central_verts[i]);
      if(next_face != FaceAttrib::NULL_VI)
        get_faces(central_verts[i], next_face, faces_out, true);
      else
        get_faces(central_verts[i], central_face, faces_out, false);
    }
  }


  inline void
  face_centroid(const VertAttrib::vi *V, my_float_t *centroid)
  {
    std::copy(V[0]->pos, V[0]->pos + 3, centroid);
    my_axpy(3, 1.0, V[1]->pos, 1, centroid, 1);
    my_axpy(3, 1.0, V[2]->pos, 1, centroid, 1);
    for(size_t j = 0; j < 3; ++j) centroid[j] /= 3.0;
  }


  dir_point_storage<VertAttrib> A_vertices;
  std::vector<FaceAttrib> A_faces;
  point_storage<point_t> A_sample_pts; //!< Points 
  my_float_t A_max_dist_to_cent;
  bool A_fail;

};




} // end namespace geometry
} // end namespace SimSite3D
#endif
