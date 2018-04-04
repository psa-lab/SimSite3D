#ifndef VERT_ATTRIB_HEADER_FILE_INCLUDED
#define VERT_ATTRIB_HEADER_FILE_INCLUDED

#include <VertAttrib.H>

namespace SimSite3D{
namespace geometry{

//! Class to hold the attributes associated with a face
/*!
 *  Note: inside this class I tended to use FaceAttrib::vi so that it is 
 *  clear which vector iterator is being used
 */
class FaceAttrib{
public:

  typedef std::vector<FaceAttrib>::iterator vi;  //!< Iterator in a vector of faces
  typedef std::vector<FaceAttrib>::const_iterator vci; //!< Const iterater in a vector of faces

  class point_data_t{
  public:

    point_data_t()
    { init(); }

    point_data_t(const my_float_t *sample_pt)
    {
      init();
      x = sample_pt;
    }

    point_data_t(const point_data_t &src)
    {
      init();
      do_copy(src);
    }

    const point_data_t&
    operator=(const point_data_t &src)
    {
      if(this != &src) do_copy(src);
      return *this;
    }

    ~point_data_t()
    { ; }

    const my_float_t *x;  //!< Sample point 
    FaceAttrib::vi face;  //!< Face holding the projection of x
    my_float_t dist;      //!< Distance between x and its projection
    my_float_t B[3];      //!< Barycentric cooridinates in face for x

  private:

    //! Initialize class variables
    void
    init();

    void
    do_copy(const point_data_t &src)
    {
      x = src.x;
      face = src.face;
      dist = src.dist;
      std::copy(src.B, src.B + 3, B);
    }
  };

  //! Default constructor
  FaceAttrib()
  {
    initialize();
  }

  //! Construct using 3 vertices
  /*!
   * Assumption: The vertices are in counter clockwise order so that 
   * barycentric method to determine if the projection of a point lies inside 
   * the face places the centroid of the face inside.  Besides, the right
   * hand rule has been typically followed up to this point
   *
   * NEED TO THINK ABOUT THIS MORE -- after rotation we could mess up 
   * the cross product?  If it flips to be clockwise with respect to
   * the right hand rule we should be ok since the resulting vector of the
   * cross products will still face the same direction for all 3 vectors
   */
  FaceAttrib(VertAttrib::vi V0, VertAttrib::vi V1, VertAttrib::vi V2)
  {
    initialize();
    A_V[0] = V0;
    A_V[1] = V1;
    A_V[2] = V2;
    my_float_t U[3], V[3];
    vector(3, A_V[0]->pos, A_V[1]->pos, U);
    vector(3, A_V[2]->pos, A_V[1]->pos, V);
    cross(U, V, A_normal);
    A_area = normalize(A_normal);
  }

  //! Copy constructor
  FaceAttrib(const FaceAttrib& F)
  {
    initialize();
    do_copy(F);
  }

  //! Assignment operator
  const FaceAttrib&
  operator=(const FaceAttrib& F)
  {
    if(&F != this) do_copy(F);
    return *this;
  }
    
  virtual ~FaceAttrib()
  {
  }
  
  //! assign the 3 vertices
  /*!
   * Assumption: The vertices are in counter clockwise order so that 
   * barycentric method to determine if the projection of a point lies inside 
   * the face places the centroid of the face inside.  Besides, the right
   * hand rule has been typically followed up to this point
   */
  void
  set_vertices(VertAttrib::vi V0, VertAttrib::vi V1, VertAttrib::vi V2)
  {
    A_V[0] = V0;
    A_V[1] = V1;
    A_V[2] = V2;
    my_float_t U[3], V[3];
    vector(3, A_V[0]->pos, A_V[1]->pos, U);
    vector(3, A_V[2]->pos, A_V[1]->pos, V);
    cross(U, V, A_normal);
    A_area = normalize(A_normal);
  }

  //! Copy constructor
  void
  add_sample_pt(const my_float_t *pt)
  { A_X.push_back(pt); }
    
  void
  get_other_verts(const VertAttrib::vci V, VertAttrib::vi A, VertAttrib::vi B)
  { 
    if(A_V[0] == V){
      A = A_V[1];
      B = A_V[2];
    }else if(A_V[1] == V){
      A = A_V[2];
      B = A_V[0];
    }else if(A_V[2] == V){
      A = A_V[0];
      B = A_V[1];
    }else{
      A = VertAttrib::NULL_VI;
      B = VertAttrib::NULL_VI;
    }
  }

  //! Get the point on face that is the closest to P
  /*! ASSUMPTION: A_normal is valid!
   */
  void
  closest_point(const my_float_t *P, const my_float_t prev_best_d,
                my_float_t *d, my_float_t *closest_pt,
                my_float_t *corr_N) const;

  //! Use the area method to compute barycentric coordinates for p
  /*!
   * Assumption: p lies in the plane defined by the 3 vertices of this face
   * This function constrains the barycentric coordinates to be non-negative,
   * but it does return false if the point is outside of the face.
   * Note that the "outside" is subject to numerical errors and it
   * could potentially be "inside", but that is irrelavent with respect to
   * the way we currently intend to use this function.
   *
   * @return True if p is inside this face, false otherwise
   */
  bool
  barycentric_coords(const my_float_t *p, my_float_t *b) const;

  //! Return the next face after the current face given a vertex V of the
  //! current face -- use right hand rule for next (i.e. counter-clockwise
  //! rotation)
  FaceAttrib::vi 
  next(VertAttrib::vi V) const
  {
    if(V == A_V[0]) return A_fnei[2];
    else if(V == A_V[1]) return A_fnei[0];
    else if(V == A_V[2]) return A_fnei[1];
    else{
      std::cerr << "The vertex for FaceAttrib::next() is not a vertex of "
                << "the face\n";
      return FaceAttrib::NULL_VI;
    }
  }

  //! Return the previous face before the current face given a vertex V of the
  //! current face -- use right hand rule for previou (i.e. clockwise
  //! rotation)
  FaceAttrib::vi
  prev(VertAttrib::vi V) const
  {
    if(V == A_V[0]) return A_fnei[0];
    else if(V == A_V[1]) return A_fnei[1];
    else if(V == A_V[2]) return A_fnei[2];
    else{ 
      std::cerr << "The vertex for FaceAttrib::prev() is not a vertex of "
                << "the face\n";
      return FaceAttrib::NULL_VI;
    }
  } 

  bool
  set_face_neighbor(VertAttrib::vi V0, VertAttrib::vi V1, FaceAttrib::vi nbr)
  {
    if((V0 == A_V[0] && V1 == A_V[1]) || (V0 == A_V[1] && V1 == A_V[0]))
      A_fnei[0] = nbr;
    else if((V0 == A_V[1] && V1 == A_V[2]) || (V0 == A_V[2] && V1 == A_V[1]))
      A_fnei[1] = nbr;
    else if((V0 == A_V[2] && V1 == A_V[0]) || (V0 == A_V[0] && V1 == A_V[2]))
      A_fnei[2] = nbr;
    else{
    
//      std::cerr << "One or more of the vertices in set_face_neigbor() is not "
//                << "in the current face\n";
      return false;
    }
    return true;
  }
  
  void
  get_vertices(VertAttrib::vi *V) const
  {
    std::copy(A_V, A_V + 3, V);
  }

  static void
  test_removing_edge(VertAttrib::vi *end_pts, FaceAttrib::vi *edge_faces,
                     VertAttrib::vi verts_begin)
  {
/*
    std::cout << "End points (vertices): " 
              << end_pts[0] - verts_begin << " " 
              << end_pts[1] - verts_begin << "\n";

    std::cout << "faces:\n"
              << "\tvertices: "
              << edge_faces[0]->A_V[0] - verts_begin << " "
              << edge_faces[0]->A_V[1] - verts_begin << " "
              << edge_faces[0]->A_V[2] - verts_begin << "\n"
              << edge_faces[1]->A_V[0] - verts_begin << " "
              << edge_faces[1]->A_V[1] - verts_begin << " "
              << edge_faces[1]->A_V[2] - verts_begin << "\n\n";
*/

    // Get a copy of the top right (TR), bottom right (SR), top left (TL)
    // and bottom left (SL) when the edge is oriented such that v_t is the
    // top and v_s is the bottom vertex
    std::vector<FaceAttrib> temp_faces(4);
    temp_faces[0] = *(edge_faces[0]->next(end_pts[0]));
    temp_faces[1] = *(edge_faces[1]->prev(end_pts[0]));
    temp_faces[2] = *(edge_faces[0]->prev(end_pts[1]));
    temp_faces[3] = *(edge_faces[1]->next(end_pts[1]));
    FaceAttrib::vi TR = temp_faces.begin();
    FaceAttrib::vi TL = temp_faces.begin() + 1;
    FaceAttrib::vi SR = temp_faces.begin() + 2;
    FaceAttrib::vi SL = temp_faces.begin() + 3;
/*
    std::cout << "TR:  "
              << temp_faces[0].A_V[0] - verts_begin << " "
              << temp_faces[0].A_V[1] - verts_begin << " "
              << temp_faces[0].A_V[2] - verts_begin << "\n\n"
              << "TL:  "
              << temp_faces[1].A_V[0] - verts_begin << " "
              << temp_faces[1].A_V[1] - verts_begin << " "
              << temp_faces[1].A_V[2] - verts_begin << "\n\n";
*/

    // Change the face neighbors for TR, SR, TL and SL as they would be
    // if the given edge was removed
    TR->A_fnei[edge_faces[0]->prev_fnei_idx(end_pts[0])] = SR;
    SR->A_fnei[edge_faces[0]->next_fnei_idx(end_pts[0])] = TR;
    TL->A_fnei[edge_faces[1]->prev_fnei_idx(end_pts[1])] = SL;
    SL->A_fnei[edge_faces[1]->next_fnei_idx(end_pts[1])] = TL;
  
    // Test with 3 different starting locations
    std::vector<VertAttrib> new_vert(3);  
    new_vert[0] = *(end_pts[0]);
    new_vert[1] = *(end_pts[1]);
    for(int i = 0; i < 3; ++i) 
      new_vert[2].pos[i] = 0.5*(end_pts[0]->pos[i] + end_pts[1]->pos[1]);
    // don't forget to update normal at the end

    // We now need to replace the iterator to the VertAttribs corresponding
    // to the end points of the edge to remove in TR, SR, TL, and SL
    TR->replace_vert(end_pts[0], new_vert.begin());
    TL->replace_vert(end_pts[0], new_vert.begin());
    SR->replace_vert(end_pts[1], new_vert.begin());
    SL->replace_vert(end_pts[1], new_vert.begin());

    // We need to be careful as we cannot rotate about vertices to get
    // all of the faces :( using the copies since we are not actually
    // modifying the mesh
    // Thus for each triangle that has one of the end_points as a vertex,
    // we need to get all of it's neighbor triangles
    
    std::vector<FaceAttrib::vi> faces_of_edge;
    // get_faces does not add the starting face to the vector
    // this way faces do not get counted twice
    get_faces(end_pts[0], edge_faces[0], &faces_of_edge);
    get_faces(end_pts[1], edge_faces[1], &faces_of_edge);
    
    // For each face in the faces the share exactly one vertex with
    // the edge, get all their neighbor faces.  
    std::vector< std::vector<FaceAttrib::vi> > nbrs_of_edge_faces;
    for(size_t i = 0; i < faces_of_edge.size(); ++i){
      // Ignore the two faces that share the edge as they are getting optimized
      // out
      if(faces_of_edge[i] == edge_faces[0] || faces_of_edge[i] == edge_faces[1])
        continue;

      std::vector<FaceAttrib::vi> tmp_vec;
      nbrs_of_edge_faces.push_back(tmp_vec);
      get_faces(faces_of_edge[i],  &nbrs_of_edge_faces.back());
  
      // If a face is TR, TL, SR, or SL get all the faces for it and
      // then replace it with the modified local copy of that face (it is the
      // first face in the most recently added vector)
      // TR
      if(faces_of_edge[i] == edge_faces[0]->next(end_pts[0])){
        nbrs_of_edge_faces.back()[0] = TR;
      // SR
      }else if(faces_of_edge[i] == edge_faces[0]->prev(end_pts[0])){
        nbrs_of_edge_faces.back()[0] = SR;
      // TL
      }else if(faces_of_edge[i] == edge_faces[1]->next(end_pts[1])){
        nbrs_of_edge_faces.back()[0] = TL;
      // SL
      }else if(faces_of_edge[i] == edge_faces[1]->prev(end_pts[1])){
        nbrs_of_edge_faces.back()[0] = SL;
      }
    }

    // For each vector in nbrs_of_edge_faces, the first face is the
    // one for which we need to check the point distances (A_X) in that face
    // to all the faces in the current vector
    // Cheat and use a map for now -- too tired to think straight 
//    std::map< FaceAttrib::vi, std::vector<const my_float_t*> >new_correspondences; 

    std::vector<point_data_t> correspondences;

    for(size_t i = 0; i < nbrs_of_edge_faces.size(); ++i){
      std::vector<FaceAttrib::vi> &nbrs = nbrs_of_edge_faces[i];
      FaceAttrib::vi curr_face = nbrs[0];
      std::vector<const my_float_t*>::const_iterator x;
      for(x = curr_face->A_X.begin(); x < curr_face->A_X.end(); ++x){
        // Check distance of *x to each face in nbrs_of_edg_faces[i] 
        // Keep the lowest distance + barycenter coords for the projection

        // If the current face is TR or SR also check the points for edge_faces[0]
        // If the current face is TL or SL also check the points for edge_faces[1]

        point_data_t best_proj(*x);
        std::vector<FaceAttrib::vi>::iterator face_iter;
        for(face_iter = nbrs.begin(); face_iter < nbrs.end(); ++face_iter){
          FaceAttrib::vi face = *face_iter;

          // project point to current face
          my_float_t x_in_face_plane[3];
          my_float_t d = face->proj_pt_to_plane(*x, x_in_face_plane);

          // get barycentric coords -- if not if face, see if we can use
          my_float_t B[3];
          //bool pt_in_face = face->barycentric_coords(x_in_face_plane, B);
          face->barycentric_coords(x_in_face_plane, B);
              std::cout << "\nbarycentric coords for face: ";
              std::cout << "B[0]: " << B[0] ;
              std::cout << "   B[1]: " << B[1] ;
              std::cout << "   B[2]: " << B[2]  << "\n";

          if(d < best_proj.dist){
            best_proj.dist = d; 
            std::copy(B, B + 3, best_proj.B);
            best_proj.face = face;
          } 
        }

        if(best_proj.face == FaceAttrib::NULL_VI){
          std::cerr << "Error in gettin face for a sample point\n\n";
        }else{
          correspondences.push_back(best_proj);  
        }
      }
    }
  }

  my_float_t
  proj_pt_to_plane(const my_float_t *pt, my_float_t *proj_pt) const
  {
    // Get a vector from one of the vertices to the point
    my_float_t V[3]; 
    vector(3, pt, A_V[0]->pos, V);
    std::cout << "The vector V[0] -> pt: " << V[0] << " " << V[1] << " " 
              << V[2] << "\n";
    std::cout << "normal to face: " << A_normal[0] << " " << A_normal[1] << " "
              << A_normal[2] << "\n";

    // signed distance from pt to plane
    my_float_t rv = dot(A_normal, V);
    std::cout << "Signed dist to plane: " << rv << "\n";

    // Project the point by adding the negative of rv times the normal to the 
    // given point
    std::copy(pt, pt + 3, proj_pt);
    my_axpy(3, -1.0 * rv, A_normal, 1, proj_pt, 1);
    return rv;
  }

  friend std::ostream& 
  operator<<(std::ostream &out, const FaceAttrib& F)
  { 
    for(int i = 0; i < 3; ++i)
      out << F.A_V[i]->pos[0] << " " << F.A_V[i]->pos[1] << " " 
          << F.A_V[i]->pos[2] << " ";

    return out;
  }

  static std::vector<FaceAttrib> NULL_FACES_STORAGE;
  static const FaceAttrib::vi NULL_VI;

private:
  //! Initialize class variables
  void
  initialize();

  void
  do_copy(const FaceAttrib& F)
  {
    std::copy(F.A_V, F.A_V + 3, A_V);
    std::copy(F.A_fnei, F.A_fnei + 3, A_fnei);
    A_X.resize(F.A_X.size());
    std::copy(F.A_X.begin(), F.A_X.end(), A_X.begin());
    std::copy(F.A_normal, F.A_normal + 3, A_normal);
    A_area = F.A_area;
  }

  //! Return the next face after the current face given a vertex V of the
  //! current face -- use right hand rule for next (i.e. counter-clockwise
  //! rotation)
  size_t
  next_fnei_idx(VertAttrib::vi V) const
  {
    if(V == A_V[0]) return 2;
    else if(V == A_V[1]) return 0;
    else if(V == A_V[2]) return 1;
    else{
      std::cerr << "The vertex for FaceAttrib::next() is not a vertex of "
                << "the face\n";
      return 100000;
    }
  }

  //! Return the previous face before the current face given a vertex V of the
  //! current face -- use right hand rule for previous (i.e. clockwise
  //! rotation)
  size_t
  prev_fnei_idx(VertAttrib::vi V) const
  {
    if(V == A_V[0]) return 1;
    else if(V == A_V[1]) return 2;
    else if(V == A_V[2]) return 0;
    else{
      std::cerr << "The vertex for FaceAttrib::prev() is not a vertex of "
                << "the face\n";
      return 100000;
    }
  }

  //! Replace the old vertex with the new vertex
  bool
  replace_vert(VertAttrib::vi V_old, VertAttrib::vi V_new)
  {
    if(V_old == A_V[0]) A_V[0] = V_new;
    else if(V_old == A_V[1]) A_V[1] = V_new;
    else if(V_old == A_V[2]) A_V[2] = V_new;
    else{
      std::cerr << "The old vertex for FaceAttrib::replace_vert() is not a "
                << "vertex of the face\n";
      return false;
    }
    return true;
  } 
  
  //! Get all the faces for the given vertex
  /*! This method does not add the starting face 
   * First this method will check all the faces that can be found
   * via the face->next links.  If a "null" face is found,
   * it will check all the faces via the face->prev links.
   *
   * This method will fail if the mesh is of poor quality and all of the
   * faces cannot be found by both next and prev traversal
   */
  static bool 
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
    
  static bool
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
    return true;
  }

  //! Example "sharp" function -- always returns false
  bool
  not_sharp(VertAttrib::vi A, VertAttrib::vi B)
  { return false; }

  my_float_t
  E_disc(VertAttrib::vi V_s, VertAttrib::vi V_t, VertAttrib::vi V_l,
         VertAttrib::vi V_r,
         bool(*sharpFunc)(VertAttrib::vi A, VertAttrib::vi B))
  {
    return 0.0;
  }

  VertAttrib::vi A_V[3];       //!< Vertices of this face
  FaceAttrib::vi A_fnei[3];    //!< Neighbor faces -- first corresponds to A_V[0] -> A_V[1], ...
  // should store these elsewhere ?
  std::vector<const my_float_t*> A_X;
  //point_storage<point_t> A_X; //!< Points of X that could project to this face

  my_float_t A_normal[3]; //!< Normal to the face (V[0] - V[1]) x (V[2] - V[1])
  my_float_t A_area;      //!< Area of the face

};


} // end namespace geometry
} // end namespace SimSite3D



#endif
