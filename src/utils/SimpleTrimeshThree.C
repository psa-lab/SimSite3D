#include <SimpleTrimeshThree.H>
#include <stream_basics.H>

using namespace SimSite3D::geometry;

void 
Edge::init()
{
  std::fill(A_V, A_V + 2, VertAttrib::NULL_VI);
  std::fill(A_nbrs, A_nbrs + 2, FaceAttrib::NULL_VI);
}

SimpleTrimeshThree::SimpleTrimeshThree()
{
  A_fail = 0;
  A_max_dist_to_cent = -1;
}

bool
SimpleTrimeshThree::read_vert_file(const std::string fname, 
                                   const bool sample_surface)
{
  std::ifstream vert_file;
  if(!open_ifstream(vert_file, fname + ".vert")) return false;
  
  // Eat the first two lines
  std::string line;
  std::getline(vert_file, line);
  std::getline(vert_file, line);
  
  size_t num_spheres, n_verts;
  my_float_t density, probe_radius;
  vert_file >> n_verts >> num_spheres >> density >> probe_radius;
  
  // Notice that the getline is used only to eat up the remaining portion 
  // of the line and that it is ok to eat it here since we have a part of 
  // the third line left at the beginning of the loop
  for(size_t i = 0; i < n_verts && std::getline(vert_file, line); ++i){
    VertAttrib tmp_VA;
    vert_file >> tmp_VA.pos[0] >> tmp_VA.pos[1] >> tmp_VA.pos[2]
              >> tmp_VA.dir[0] >> tmp_VA.dir[1] >> tmp_VA.dir[2];
    A_vertices.push_back(tmp_VA);

    if(sample_surface){
      point_t tmp;
      std::copy(tmp_VA.pos, tmp_VA.pos + 3, tmp.pos);
      A_sample_pts.push_back(tmp);
    }
  }

  // Simple check
  if(n_verts > A_vertices.size()) return false;
  return true;
}

bool
SimpleTrimeshThree::read_face_file(const std::string fname, 
                                   std::vector<Edge> *edges,
                                   const bool check_vertex_order, 
                                   const bool compute_max_dist_to_centroid)
{
  std::ifstream face_file;
  if(!open_ifstream(face_file, fname + ".face")) return false;

  // Eat the first two lines
  std::string line;
  std::getline(face_file, line);
  std::getline(face_file, line);

  size_t num_faces, num_spheres;
  my_float_t density, probe_radius;
  face_file >> num_faces >> num_spheres >> density >> probe_radius;

  // Use a 2D vector (vector of vectors) to index edges to initialize the
  // fnei (face neighbors) array in each face.  The idea is if the
  // value of an index is false, set it to the current face index and we have 
  // 1 triangle on the edge.  If the index is not negative, then we know the
  // index of the "other" face and can set fnei in both faces.
  std::vector<int> edges_to_face_idz(
    (A_vertices.size() - 1) * A_vertices.size() / 2);

  if(edges){
    std::fill(edges_to_face_idz.begin(), edges_to_face_idz.end(), -1);
    edges->clear();
    edges->reserve(3 * A_vertices.size());
  }

  if(A_faces.size()) A_faces.clear();
  A_faces.resize(num_faces);
  uint v_idz[3];
  uint v_idz_in[3];
  for(size_t i = 0; i < num_faces && std::getline(face_file, line); ++i){
    face_file >> v_idz_in[0] >> v_idz_in[1] >> v_idz_in[2];

    // MSMS vertex indices are 1 based -- we need zero based
    dir_point_storage<VertAttrib>::iterator V[3];
    for(size_t j = 0; j < 3; ++j){
      v_idz[j] = v_idz_in[j] - 1; 
      V[j] = A_vertices.begin() + v_idz[j];
    }

    // Since this check is expected to be used only at template creation time,
    // in the interest of having less code we use the corresponding point 
    // function although the centroid is already in the plane.
    my_float_t centroid[3];
    if(check_vertex_order){
      face_centroid(V, centroid);

      my_float_t d;
      my_float_t pt[3];
      // Need the static casts or the compiler gets confused on which
      // corresponding point function to call -- technically we only need
      // the third one cast
      if(!corresponding_point(centroid, 
                              //static_cast<const my_float_t*>(V[0]->pos), 
                              //static_cast<const my_float_t*>(V[1]->pos), 
                              //static_cast<const my_float_t*>(V[2]->pos), 
			      V[0]->pos, V[1]->pos, V[2]->pos,
			      my_float_max,
			      &d, pt)){
        std::swap(V[1], V[2]);
	std::swap(v_idz[1], v_idz[2]);
      }
    }

    // This should be computed at template generation time and stored in the
    // .csv file rather than computing at screening time -- we want
    // to make screening as fast as possible
    if(compute_max_dist_to_centroid){
      if(!check_vertex_order) face_centroid(V, centroid);
      
      for(size_t j = 0; j < 3; ++j){
         my_float_t d2 = dist_squared(centroid, V[j]->pos);
         A_max_dist_to_cent = 
           (A_max_dist_to_cent < d2 ? d2 : A_max_dist_to_cent);
      }
    }

    A_faces[i].set_vertices(V[0], V[1], V[2]);
    // Check for neighbors
    for(int j = 0; j < 3; ++j){
      int next_j = (j+1) % 3;
      std::vector<int>::iterator idx = edges_to_face_idz.begin();
      if(V[j] > V[next_j]){
        idx += (v_idz[j]*(v_idz[j] - 1)/2 + v_idz[next_j]);
      }else{
        idx += (v_idz[next_j]*(v_idz[next_j] - 1)/2 + v_idz[j]);
      }

      // Set the index for the edge to that of the current face
      if(*idx < 0) *idx = i;
      // We have both faces that share this edge -- update their data
      // structures (fnei array)
      else{
        FaceAttrib::vi curr_face = A_faces.begin() + i;
        FaceAttrib::vi other_face = A_faces.begin() + *idx;
        curr_face->set_face_neighbor(V[j], V[next_j], other_face);
	other_face->set_face_neighbor(V[j], V[next_j], curr_face);

        if(edges)
	  edges->push_back(Edge(V[j], V[next_j], curr_face, other_face));
      }

    // add the triangle -- add by assignment so that we can store an iterator
    // to the face at this point rather than after the loop
      if(edges)
        A_faces[i].add_sample_pt((A_sample_pts.begin() + v_idz[j])->pos);
    }
  }

  A_max_dist_to_cent = sqrt(A_max_dist_to_cent);
  // compute the area of each face

  return true;
}

void
SimpleTrimeshThree::initialize_opt(std::vector<Edge> &edges)
{ // Check for all the edges that we could collapse
  //   std::vector< std::vector<FaceAttrib>::iterator > faces;
  // spin through each face and take current vert to the other two
  // ignoring those for which the idx of the other vert is less than current
  // vert
  //

  // Do this step wise -- adding each term after we get the previous term 
  // implemented.
  std::vector<Edge>::const_iterator E;
  for(E = edges.begin(); E < edges.end(); ++E){
    // Do the cyclic optimization here -- at most 10 steps -- 5 may be 
    // sufficient
    
    // 1) For each sample point find its best projection -- because
    // this is initial optimization and we use the vertices as given 
    // we don't actually do the projections (they are not needed) and
    // the barycentric coordinates can be determined from index for each vertex

    // 2) Determine the input to LSQR and call LSQR -- one of the papers 
    // has a short section on the form of the optimization matrix

    
    VertAttrib::vi V[2];
    FaceAttrib::vi F[2];
    E->get_verts(&V[0], &V[1]);
    E->get_faces(&F[0], &F[1]);
    if((A_vertices.begin() + 4 == V[0] && A_vertices.begin() + 5 == V[1]) ||
       (A_vertices.begin() + 4 == V[1] && A_vertices.begin() + 5 == V[0])){
      std::cout << "found it\n";

      // Get the faces for the current edge
//      std::vector<FaceAttrib::vi> faces;
//      get_faces(*E, &faces);
//      std::cout << "num faces: " << faces.size() << std::endl;

      // Copy the faces in *(v_s, v_t)
      // This is reasonable since we will not be modifying vertex
      // attributes or neighbors 
      // Area of the face is borked --

      FaceAttrib::test_removing_edge(V, F, A_vertices.begin());
  





    }




  }
  

}


