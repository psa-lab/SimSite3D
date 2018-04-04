#ifndef SIMPLE_TRIMESH_TWO_HEADER_INCLUDED
#define SIMPLE_TRIMESH_TWO_HEADER_INCLUDED

#include <dir_point.H>
#include <dir_point_storage.H>


namespace ASCbase{
namespace geometry{

//! Simple surface class built to contain basic features of a triangle mesh
class SimpleTrimeshTwo{
public:

  //! Simple struct to contain the basic features of a triangle 
  typedef struct{
    const my_float_t *vertices[3];
    const my_float_t *normals[3];
    my_float_t area;
  }triangle_t;

  typedef std::vector<triangle_t> face_vec;  //!< Vector of faces (triangles)
  typedef face_vec::iterator face_vi;
  typedef face_vec::const_iterator face_vci;
  typedef std::map<const my_float_t*, std::vector<face_vi> > vert_face_map;
  typedef vert_face_map::iterator vert_face_mi;
  typedef vert_face_map::const_iterator vert_face_mci;

  //! Load MSMS vert file using fname + ".vert"
  SimpleTrimeshTwo(const std::string fname, 
                   const bool load_mutable_faces=false);

  //! Basic destruction
  virtual
  ~SimpleTrimeshTwo();
  
  //! Get a pointer to the beginning of the vertices array
  const my_float_t*
  vertices_begin() const
  { return A_vertices; }

  //! Number of vertices in the array (array size is 3*number_of_vertices)
  size_t
  number_of_vertices() const
  { return A_num_verts; }

  //! Get a pointer to the beginning of the normals array; number of normals
  //! is the same as number of vertices
  const my_float_t*
  normals_begin() const
  { return A_normals; }

  face_vci
  faces_begin() const
  { return A_triangles.begin(); }

  face_vci
  faces_end() const
  { return A_triangles.end(); }

  const size_t
  number_of_faces() const
  { return A_triangles.size(); }

  //! Apply the transform VR + T to the vertices V
  virtual void
  transform(const my_float_t* R, const my_float_t* T);

  //! Apply the transform (V-T)R.t to the vertices V, where R.t is R transposed
  virtual void
  inverse_transform(const my_float_t* R, const my_float_t* T);

  //! Revert the vertices to the saved positions
  void
  revert()
  {
    if(A_original_vertices){
      size_t arr_len = 3*A_num_verts;
      std::copy(A_original_vertices, A_original_vertices + arr_len, A_vertices);
      std::copy(A_original_normals, A_original_normals + arr_len, A_normals);
    }
  }

  //! Get transform to move current vertices to original position of
  //! the vertices
  /*! The reason we find the transform from the current back to original
   *  is SimSite3D moves the query to the db reference frame.  The user 
   *  typically wants to align all dbas hits to the query reference frame.
   *  Thus we compute current -> original positions.
   */
  inline virtual void
  get_current_inverse_3D_transform(Quaternion *Q, my_float_t *T) const
  {
    if(A_original_vertices)
      lse_3D_fit(A_original_vertices, A_vertices, A_num_verts, Q, T);
    else{
      my_float_t q[] = {1.0, 0.0, 0.0, 0.0};
      *Q = Quaternion(q, 4);
      std::fill(T, T+3, 0.0);
    }
  }

  virtual bool
  closest_point(const my_float_t *pt, my_float_t *dist, my_float_t *cp,
                const my_float_t max_dist = 2.0) const
  {
    std::cerr << "Unimplemented function: SimpleTrimeshTwo::closest_point\n"
              << "Please use ImmovableMesh::closest_point instead\n";
    return false;
  }

  //! Find closest point on surface mesh for the given point (pt)
  /*! This method is inefficient since it checks each triangle in the mesh
   *  to find the closest point for each input vertex.  This function exists
   *  to help with prototyping and as a fallback, simpler method to check
   *  that the more advanced methods still compute the closest point.
   *  In particular, we do the most straightforward and naive check of each
   *  triangle without first checking vertex distances, etc.
   */
  virtual bool
  compare(const my_float_t *other_vert_begin, const uint num_other_vert,
          my_float_t *closest_pts, my_float_t *dists,
          const my_float_t *other_normals_begin = 0, 
          my_float_t *best_Ns = 0, my_float_t *N_dists = 0,
          const my_float_t max_dist = my_float_max) const;

  //! compute features for an alignment
  /*!
   *  Assumes the array dists is of the same length as the number of vertices
   */
  void
  compute_features(const my_float_t *dists, const my_float_t max_dist,
                   const my_float_t *norm_dots, my_float_t *vertex_count, 
                   my_float_t *face_count, my_float_t *RMSE, 
                   my_float_t *RMS_norm_err, my_float_t *area, 
                   my_float_t *max_q_pt_dist_p) const;

  //! Get the maximum allowed distance from centroid
  my_float_t
  max_dist_to_centroid() const
  { return A_max_dist_to_cent; }

  //! Returns true if something went wrong, otherwise false
  bool
  fail() const
  { return A_fail; }

  //! Wrapper to call corresponding_point and update closest_pt and best dist
  /*! 
   *  This function is not required to be public, but it doesn't matter since
   *  it doesn't modify any class variables.
   */
  void
  compute_close_point(const my_float_t *pt, const my_float_t *v0,
                      const my_float_t *v1, const my_float_t *v2,
                      const my_float_t prev_best_d,
                      my_float_t *best_d, my_float_t *closest_pt) const;

  //! Write the vert & face files
  bool
  write(const std::string ofname) const;

  bool
  remove_small_triangles();

  const my_float_t
  get_total_SA() const
  {
    my_float_t rv = 0;
    for(face_vci face = A_triangles.begin(); face < A_triangles.end(); ++face)
      rv += face->area; 
    return rv;
  }

protected:

  void 
  set_fail_flag()
  { A_fail = true; }

  //! Read an MSMS face file and build the triangles.  If check_vertex_order,
  //! make sure the vertices have the correct order or swap them.
  bool
  read_face_file(const std::string fname, bool check_vertex_order,
                 bool compute_max_dist_to_centroid);

  //! Initialize the triangle storage
  virtual void
  init_triangle_data_structure(const uint num_vert, const uint num_tri)
  {
    A_triangles.reserve(num_tri);
  }

  //! Add a triangle to the triangle storage
  virtual void
  update_triangle_data_structure(const uint tri_idx, const uint *vert_idz,
                                 const my_float_t **vertices,
                                 const my_float_t **normals)
  {
    triangle_t tmp;
    std::copy(vertices, vertices + 3, tmp.vertices);
    std::copy(normals, normals + 3, tmp.normals);
    A_triangles.push_back(tmp);  
  }

  //! Simple function to update closest point and distance
  /*! 
   *  This function is not required to be protected rather than public, 
   *  but it is rather dangerous in that it uses pointers.
   */
  inline bool
  update_closest_pt(my_float_t *best_d, my_float_t *closest_pt,
                    const my_float_t d, const my_float_t *cp) const
  {
#if 0
    my_float_t abs_d = (d < 0.0 ? -1.0 * d : d);
    my_float_t abs_saved_d = (*best_d < 0.0 ? -1.0 * (*best_d) : (*best_d));
#endif
    my_float_t abs_d = d;
    if(abs_d < 0.0) abs_d *= -1.0;
    my_float_t abs_saved_d = *best_d ;
    if(abs_saved_d < 0.0) abs_saved_d *= -1.0;
    if(abs_d < abs_saved_d){
      *best_d = d;
      std::copy(cp, cp + 3, closest_pt);
      return true;
    }
    return false;
  }

  //! Get a non-constant pointer to the beginning of the vertices array
  /*! I am unsure of the "correct" way to tell the g++ compiler which of
   * the vertices begin functions we want -- it always chooses the unconstant
   * one with the same name
   */
  my_float_t*
  mutable_vertices_begin()
  { return A_vertices; }

private:

  //! Initialize class variables
  void
  init();

  //! Read an MSMS vert input file and store the vertices and normals
  bool
  read_vert_file(const std::string fname);

  //! Compute the area of the triangles
  void
  compute_areas();

  //! Load the faces in a mutable data structure so that we can modify
  //! the vertices, faces, and associated pointers during mesh adjustments
  bool
  load_mutable_faces(const std::string faces_fname);

  size_t A_num_verts;      //!< Number of vertices (and normals)
  my_float_t *A_vertices;  //!< Vertices array -- 3*A_num_verts
  my_float_t *A_normals;   //!< Normals array -- 3*A_num_verts
  my_float_t *A_original_vertices;
  my_float_t *A_original_normals;
  my_float_t A_max_dist_to_cent;  //!< Maximum distance from a triangle centroid to one of its vertices
  face_vec A_triangles;              //!< Triangle storage
#if 0
  vert_face_map A_map_vert_to_faces; //!< Map of each vertex to all its triangles
#endif

  bool A_fail;

  //! Simple struct to contain the basic features of a triangle 
  typedef struct{
    my_float_t *vertices[3];
    my_float_t *normals[3];
    my_float_t area;
  }mutable_face_t;

  typedef std::vector<mutable_face_t> mutable_face_vec;  //!< Vector of mutable faces (triangles)
  typedef mutable_face_vec::iterator mutable_face_vi;
  typedef mutable_face_vec::const_iterator mutable_face_vci;
  typedef std::map<const my_float_t*, std::vector<mutable_face_vi> > mutable_vert_face_map;
  typedef mutable_vert_face_map::iterator mutable_vert_face_mi;
  typedef mutable_vert_face_map::const_iterator mutable_vert_face_mci;

  // I am lazy -- it is easier for me to store pointers to the normals
  // rather than use pointer arithmetic to look them up
  typedef struct{
    my_float_t **vert0_p;
    my_float_t **vert1_p;
    my_float_t **normal0_p;
    my_float_t **normal1_p;
  }local_edge_ptrs_t;

  mutable_face_vec A_mutable_faces;

  
};
} // end namespace geometry
} // end namespace ASCbase
#endif
