#include <TriMeshSphere.H>
#include <iomanip>

using namespace ASCbase::geometry;

#if 0
const my_float_t TriMeshSphere::A_level2_points[] =
  { 0.866025, 0.154508, 0.475528,
    0.866025, 0.500000, 0.000000,
    0.580892, 0.658524, 0.478446,
    1.000000, 0.000000, 0.000000,
    0.500000, 0.866025, 0.000000,
    0.500000, 0.267617, 0.823639,
    0.866025, -0.404508, 0.293893,
    0.580892, -0.251534, 0.774141,
    0.500000, -0.700629, 0.509037,
    0.866025, -0.404508, -0.293893,
    0.580892, -0.813980, 0.000000,
    0.500000, -0.700629, -0.509037,
    0.866025, 0.154508, -0.475528,
    0.580892, -0.251534, -0.774141,
    0.500000, 0.267617, -0.823639,
    0.866025, 0.500000, -0.000000,
    0.580892, 0.658524, -0.478446,
    0.500000, 0.866025, -0.000000,
  };
#endif

#if 0
const my_float_t TriMeshSphere::A_level2_points[] =
  { 2.5980762,0.4635255,1.4265848,
    2.5980762,1.5000000,0.0000000,
    1.7426771,1.9755719,1.4353370,
    3.0000000,0.0000000,0.0000000,
    1.5000000,2.5980762,0.0000000,
    1.5000000,0.8028497,2.4709173,
    2.5980762,-1.2135255,0.8816779,
    1.7426771,-0.7546013,2.3224240,
    1.5000000,-2.1018878,1.5271109,
    2.5980762,-1.2135255,-0.8816779,
    1.7426771,-2.4419411,0.0000000,
    1.5000000,-2.1018878,-1.5271109,
    2.5980762,0.4635255,-1.4265848,
    1.7426771,-0.7546013,-2.3224240,
    1.5000000,0.8028497,-2.4709173,
    1.7426771,1.9755719,-1.4353370
  };
#endif

const my_float_t TriMeshSphere::A_level1_sphere_points[] =
  { 0.707107,0.218508,0.672498,
    0.707107,0.707107,0.000000,
    0.000000,0.809017,0.587785,
    1.000000,0.000000,0.000000,
    0.000000,1.000000,0.000000,
    0.000000,0.309017,0.951056,
    -0.707107,0.218508,0.672498,
    -0.707107,0.707107,0.000000,
    -1.000000,0.000000,0.000000,
    0.707107,-0.572061,0.415627,
    0.000000,-0.309017,0.951056,
    0.000000,-0.809017,0.587785,
    -0.707107,-0.572061,0.415627,
    0.707107,-0.572061,-0.415627,
    0.000000,-1.000000,0.000000,
    0.000000,-0.809017,-0.587785,
    -0.707107,-0.572061,-0.415627,
    0.707107,0.218508,-0.672498,
    0.000000,-0.309017,-0.951056,
    0.000000,0.309017,-0.951056,
    -0.707107,0.218508,-0.672498,
    0.000000,0.809017,-0.587785
  };
const uint TriMeshSphere::A_level1_sphere_npts = 22;

const uint TriMeshSphere::A_level1_sphere_triangles[] = 
  { 0,1,2,
    3,1,0,
    4,2,1,
    5,0,2,
    6,7,2,
    8,7,6,
    4,2,7,
    5,6,2,
    9,0,10,
    3,0,9,
    5,10,0,
    11,9,10,
    12,6,10,
    8,6,12,
    5,10,6,
    11,12,10,
    13,9,14,
    3,9,13,
    11,14,9,
    15,13,14,
    16,12,14,
    8,12,16,
    11,14,12,
    15,16,14,
    17,13,18,
    3,13,17,
    15,18,13,
    19,17,18,
    20,16,18,
    8,16,20,
    15,18,16,
    19,20,18,
    1,17,21,
    3,17,1,
    19,21,17,
    4,1,21,
    7,20,21,
    8,20,7,
    19,21,20,
    4,7,21,
  }; 
const uint TriMeshSphere::A_level1_num_sphere_tris = 40;

const my_float_t TriMeshSphere::A_level2_sphere_points[] =
  { 0.398782,0.579486,0.710753,
    0.777438,0.508840,0.369694,
    0.398782,0.855037,0.331489,
    0.707107,0.218508,0.672498,
    0.707107,0.707107,0.000000,
    0.000000,0.809017,0.587785,
    0.923879,0.118256,0.363954,
    0.923879,0.382683,0.000000,
    1.000000,0.000000,0.000000,
    0.382683,0.923879,0.000000,
    0.000000,0.951056,0.309017,
    0.000000,1.000000,0.000000,
    0.000000,0.587785,0.809017,
    0.382683,0.285494,0.878662,
    0.000000,0.309017,0.951056,
    -0.398782,0.579486,0.710753,
    -0.777438,0.508840,0.369694,
    -0.398782,0.855037,0.331489,
    -0.707107,0.218508,0.672498,
    -0.707107,0.707107,0.000000,
    -0.923879,0.118256,0.363954,
    -0.923879,0.382683,0.000000,
    -1.000000,0.000000,0.000000,
    -0.382683,0.923879,0.000000,
    -0.382683,0.285494,0.878662,
    0.398782,-0.496895,0.770758,
    0.777438,-0.194359,0.598177,
    0.398782,-0.051044,0.915624,
    0.707107,-0.572061,0.415627,
    0.000000,-0.309017,0.951056,
    0.923879,-0.309597,0.224936,
    0.000000,0.000000,1.000000,
    0.000000,-0.587785,0.809017,
    0.382683,-0.747434,0.543043,
    0.000000,-0.809017,0.587785,
    -0.398782,-0.496895,0.770758,
    -0.777438,-0.194359,0.598177,
    -0.398782,-0.051044,0.915624,
    -0.707107,-0.572061,0.415627,
    -0.923879,-0.309597,0.224936,
    -0.382683,-0.747434,0.543043,
    0.398782,-0.886584,-0.234398,
    0.777438,-0.628960,0.000000,
    0.398782,-0.886584,0.234398,
    0.707107,-0.572061,-0.415627,
    0.000000,-1.000000,0.000000,
    0.923879,-0.309597,-0.224936,
    0.000000,-0.951056,0.309017,
    0.000000,-0.951056,-0.309017,
    0.382683,-0.747434,-0.543043,
    0.000000,-0.809017,-0.587785,
    -0.398782,-0.886584,-0.234398,
    -0.777438,-0.628960,0.000000,
    -0.398782,-0.886584,0.234398,
    -0.707107,-0.572061,-0.415627,
    -0.923879,-0.309597,-0.224936,
    -0.382683,-0.747434,-0.543043,
    0.398782,-0.051044,-0.915624,
    0.777438,-0.194359,-0.598177,
    0.398782,-0.496895,-0.770758,
    0.707107,0.218508,-0.672498,
    0.000000,-0.309017,-0.951056,
    0.923879,0.118256,-0.363954,
    0.000000,-0.587785,-0.809017,
    0.000000,0.000000,-1.000000,
    0.382683,0.285494,-0.878662,
    0.000000,0.309017,-0.951056,
    -0.398782,-0.051044,-0.915624,
    -0.777438,-0.194359,-0.598177,
    -0.398782,-0.496895,-0.770758,
    -0.707107,0.218508,-0.672498,
    -0.923879,0.118256,-0.363954,
    -0.382683,0.285494,-0.878662,
    0.398782,0.855037,-0.331489,
    0.777438,0.508840,-0.369694,
    0.398782,0.579486,-0.710753,
    0.000000,0.809017,-0.587785,
    0.000000,0.587785,-0.809017,
    0.000000,0.951056,-0.309017,
    -0.398782,0.855037,-0.331489,
    -0.777438,0.508840,-0.369694,
    -0.398782,0.579486,-0.710753,
  }; 
const uint TriMeshSphere::A_level2_sphere_npts = 82;

const uint TriMeshSphere::A_level2_sphere_triangles[] = 
  { 0,1,2,
    3,1,0,
    4,2,1,
    5,0,2,
    6,7,1,
    8,7,6,
    4,1,7,
    3,6,1,
    9,10,2,
    11,10,9,
    5,2,10,
    4,9,2,
    12,13,0,
    14,13,12,
    3,0,13,
    5,12,0,
    15,16,17,
    18,16,15,
    19,17,16,
    5,15,17,
    20,21,16,
    22,21,20,
    19,16,21,
    18,20,16,
    23,10,17,
    11,10,23,
    5,17,10,
    19,23,17,
    12,24,15,
    14,24,12,
    18,15,24,
    5,12,15,
    25,26,27,
    28,26,25,
    3,27,26,
    29,25,27,
    30,6,26,
    8,6,30,
    3,26,6,
    28,30,26,
    13,31,27,
    14,31,13,
    29,27,31,
    3,13,27,
    32,33,25,
    34,33,32,
    28,25,33,
    29,32,25,
    35,36,37,
    38,36,35,
    18,37,36,
    29,35,37,
    39,20,36,
    22,20,39,
    18,36,20,
    38,39,36,
    24,31,37,
    14,31,24,
    29,37,31,
    18,24,37,
    32,40,35,
    34,40,32,
    38,35,40,
    29,32,35,
    41,42,43,
    44,42,41,
    28,43,42,
    45,41,43,
    46,30,42,
    8,30,46,
    28,42,30,
    44,46,42,
    33,47,43,
    34,47,33,
    45,43,47,
    28,33,43,
    48,49,41,
    50,49,48,
    44,41,49,
    45,48,41,
    51,52,53,
    54,52,51,
    38,53,52,
    45,51,53,
    55,39,52,
    22,39,55,
    38,52,39,
    54,55,52,
    40,47,53,
    34,47,40,
    45,53,47,
    38,40,53,
    48,56,51,
    50,56,48,
    54,51,56,
    45,48,51,
    57,58,59,
    60,58,57,
    44,59,58,
    61,57,59,
    62,46,58,
    8,46,62,
    44,58,46,
    60,62,58,
    49,63,59,
    50,63,49,
    61,59,63,
    44,49,59,
    64,65,57,
    66,65,64,
    60,57,65,
    61,64,57,
    67,68,69,
    70,68,67,
    54,69,68,
    61,67,69,
    71,55,68,
    22,55,71,
    54,68,55,
    70,71,68,
    56,63,69,
    50,63,56,
    61,69,63,
    54,56,69,
    64,72,67,
    66,72,64,
    70,67,72,
    61,64,67,
    73,74,75,
    4,74,73,
    60,75,74,
    76,73,75,
    7,62,74,
    8,62,7,
    60,74,62,
    4,7,74,
    65,77,75,
    66,77,65,
    76,75,77,
    60,65,75,
    78,9,73,
    11,9,78,
    4,73,9,
    76,78,73,
    79,80,81,
    19,80,79,
    70,81,80,
    76,79,81,
    21,71,80,
    22,71,21,
    70,80,71,
    19,21,80,
    72,77,81,
    66,77,72,
    76,81,77,
    70,72,81,
    78,23,79,
    11,23,78,
    19,79,23,
    76,78,79,
  };
const uint TriMeshSphere::A_level2_num_sphere_tris = 160;
    
const my_float_t TriMeshSphere::A_level2_cap_points[] =
  { 0.866025,0.154508,0.475528,
    0.866025,0.500000,0.000000,
    0.500000,0.700629,0.509037,
    1.000000,0.000000,0.000000,
    0.500000,0.866025,0.000000,
    0.500000,0.267617,0.823639,
    0.866025,-0.404508,0.293893,
    0.500000,-0.267617,0.823639,
    0.500000,-0.700629,0.509037,
    0.866025,-0.404508,-0.293893,
    0.500000,-0.866025,0.000000,
    0.500000,-0.700629,-0.509037,
    0.866025,0.154508,-0.475528,
    0.500000,-0.267617,-0.823639,
    0.500000,0.267617,-0.823639,
    0.500000,0.700629,-0.509037
  };
    
#if 0
  { 2.5980762,0.4635255,1.4265848,
    2.5980762,1.5000000,0.0000000,
    1.500000,2.101888,1.527111,
    3.0000000,0.0000000,0.0000000,
    1.500000,2.598076,0.000000,
    1.500000,0.802850,2.470917,
    2.5980762,-1.2135255,0.8816779,
    1.500000,-0.802850,2.470917,
    1.500000,-2.101888,1.527111,
    2.5980762,-1.2135255,-0.8816779,
    1.500000,-2.598076,0.000000,
    1.500000,-2.101888,-1.527111,
    2.5980762,0.4635255,-1.4265848,
    1.500000,-0.802850,-2.470917,
    1.500000,0.802850,-2.470917,
    1.500000,2.101888,-1.527111,
  };
#endif
const uint TriMeshSphere::A_level2_cap_npts = 16;

const uint TriMeshSphere::A_level2_cap_triangles[] = 
  { 0,1,2,
    3,1,0,
    4,2,1,
    5,0,2,
    6,0,7,
    3,0,6,
    5,7,0,
    8,6,7,
    9,6,10,
    3,6,9,
    8,10,6,
    11,9,10,
    12,9,13,
    3,9,12,
    11,13,9,
    14,12,13,
    1,12,15,
    3,12,1,
    14,15,12,
    4,1,15,
  };
const uint TriMeshSphere::A_level2_num_cap_tris = 20;

TriMeshSphere::TriMeshSphere(const my_float_t* R, const my_float_t* T, 
                             const my_float_t* dir, 
                             const sphere_init_t init_obj,
                             const my_float_t radius, const uint level)
{
  init();

  if(init_obj == FULL_SPHERE){
    std::cout << "Creating sphere with radius " << radius << " (A)\n";

    if(level == 2){
      A_local_points = A_level2_sphere_points;
      A_npts = A_level2_sphere_npts;
      A_num_tris = A_level2_num_sphere_tris;
      A_triangles = A_level2_sphere_triangles;
    }else{
      std::cerr << "Unknown sphere sampling level in TriMeshSphere::Cstr" 
                << std::endl;
      return;
    }
  }else if(init_obj == SPHERICAL_CAP){
    std::cout << "Creating spherical cap with radius " << radius << " (A)\n";

    if(level == 2){
      A_local_points = A_level2_cap_points;
      A_npts = A_level2_cap_npts;
      A_num_tris = A_level2_num_cap_tris;
      A_triangles = A_level2_cap_triangles;
    }else{
      std::cerr << "Unknown cap sampling level in TriMeshSphere::Cstr" 
                << std::endl;
      return;
    }
  }

  init_storages(A_local_points, A_npts, A_triangles, A_num_tris, radius);
  set_global_pos(R, T, dir);
}

void
TriMeshSphere::init_storages(const my_float_t* pts, const uint npts, 
                             const uint* triangles, const uint ntris,
                             const my_float_t radius)
{
  const my_float_t* pts_end = pts + 3*npts;
  for(const my_float_t* p = pts; p < pts_end; p += 3){
    vertex_t tmp;
    std::copy(p, p + 3, tmp.pos);
    for(int i = 0; i < 3; ++i) tmp.pos[i] *= radius; 
    std::copy(p, p + 3, tmp.dir);
//    normalize(tmp.dir);
    A_nodes.push_back(tmp);
  }

  ///////
  // NOTE: from now on do not modify the vector A_nodes
  ///////

  // Set the vertex iterators -- we cannot set the half_edge iterators till
  // the vector has its final shape
  const uint* tris_end = triangles + 3*ntris;
  for(const uint* idx = triangles; idx < tris_end; idx += 3) add_triangle(idx);
  setup_vert_to_deltas_map();

  ///////
  // NOTE: from now on do not modify the vectors A_nodes or A_half_edges
  ///////
  // Set the "next" iterator
#if 0
  half_edge_vi hedge;
  for(hedge = A_half_edges.begin(); hedge < A_half_edges.end(); ){
    hedge->next = hedge + 1;
    hedge = hedge->next;
    hedge->next = hedge + 1;
    hedge = hedge->next;
    hedge->next = hedge - 2;
    ++hedge;
  }

  dir_point_storage<vertex_t>::iterator node;
  for(node = A_nodes.begin(); node < A_nodes.end(); ++node){
    node->setup_one_ring(A_half_edges.begin(), A_half_edges.end());
  }
#endif
}

void
TriMeshSphere::A_adjust_points(hbond_surface_t& surf, 
                               const uint* tri_vert_idz, const uint num_tris)
{
//  std::cout << "Adjusting points for atom: " << surf.atom->res_num << " " << surf.atom->atom_num << std::endl;


// Skip the 1 ring stuff for now since it is more complicated using my current
// coding of 1 ring, than it is to check the triangles containing the vertex
#if 0
  // For each point, check the 1 ring & if needed, the boundary edges and 
  // only keep those points that did not move over the 1 ring (or boundary 
  // edges)
#endif

  // Spin through each vertex to see if the closest point is outside the 
  // triangles containing the current vertex (similar to 1 ring)
  std::vector<bool> keep_vec(A_nodes.size());
  std::fill(keep_vec.begin(), keep_vec.end(), false);
  std::vector<bool>::iterator keep = keep_vec.begin();
  my_float_t* close_pts = new my_float_t[3 * A_nodes.size()];
  my_float_t* cp = close_pts;
  std::map<vertex_vci, std::vector<triangle_vci> >::iterator my_map_iter;
  for(my_map_iter = A_verts_to_deltas_map.begin();
      my_map_iter != A_verts_to_deltas_map.end(); ++my_map_iter, cp += 3){
    vertex_vci V = my_map_iter->first;
    surf.closest_point(V->pos, cp, 3.0);

    std::vector<triangle_vci>& my_deltas = my_map_iter->second;
    std::vector<triangle_vci>::iterator d_iter;
    my_float_t d2 = dist_squared(V->pos, cp);
    if(d2 < 1E-06) *keep = true;
    else if(d2 < my_float_max){
      for(d_iter = my_deltas.begin(); d_iter != my_deltas.end(); ++d_iter){

#if 0
        my_float_t my_cent[] = {0.0, 0.0, 0.0};
        std::vector<vertex_vci>::const_iterator tri_vert;
        for(tri_vert = (*d_iter)->vertices_begin();
            tri_vert < (*d_iter)->vertices_end(); ++tri_vert)
          for(int zz = 0; zz < 3; ++zz) my_cent[zz] += (*tri_vert)->pos[zz];
        for(int zz = 0; zz < 3; ++zz) my_cent[zz] /= 3.0;
        if(!(*d_iter)->contains(my_cent))        
          std::cout << "triangle does not contain its centroid" << std::endl;
#endif
 
        if((*d_iter)->contains(cp)){
//          std::cout << "point is kept\n";
          *keep = true;
          break;
        }//else std::cout << "point is NOT kept\n";
      }
    }
    ++keep;
  }
  cp = 0;

  // Remove the trimmed vertices
  uint new_idx = 0;
  uint old_idx = 0;
  std::map<uint, uint> new_idz;
  dir_point_storage<vertex_t> kept_nodes;
  cp = close_pts;
  keep = keep_vec.begin();
  for(vertex_vci V = A_nodes.begin(); V != A_nodes.end(); ++V, ++keep, cp += 3){
    if(*keep){
      kept_nodes.push_back(*V);
//      std::cout << "Keeping V(" << old_idx << "): " 
//                << V->pos[0] << " " << V->pos[1] << " "
//                << V->pos[2] << "\n"
//                << "  as cp " << cp[0] << " "
//                << cp[1] << " " << cp[2] << std::endl;
      std::copy(cp, cp + 3, kept_nodes.back().pos);
      new_idz[old_idx] = new_idx;
      ++new_idx;
    }
    ++old_idx;
  }

  A_nodes = kept_nodes;
  A_nodes.set_current_positions_and_directions_as_original();

  // Fill the triangle vectors & other useful vectors using the "new" indices
  A_deltas.clear(); 
  A_half_edges.clear();
  A_verts_to_deltas_map.clear();

  const uint* tri_vert_idz_end = tri_vert_idz + 3*num_tris;
  for(const uint* t = tri_vert_idz; t < tri_vert_idz_end; t += 3){
    if(keep_vec[t[0]] && keep_vec[t[1]] && keep_vec[t[2]]){
      uint idz[] = {new_idz[t[0]], new_idz[t[1]], new_idz[t[2]]};
      add_triangle(idz);
    }
  }
  setup_vert_to_deltas_map();

  // Update the normals
  for(vertex_vi V = A_nodes.begin(); V != A_nodes.end(); ++V)
    unit_vector(V->dir, V->pos, surf.atom->pos);

  delete[] close_pts;
}

my_float_t
TriMeshSphere::complementary_surface_area(const std::map<vertex_vci, bool>&
                                          kept_vertices) const
{
  my_float_t total_SA = 0.0;
  for(triangle_vci tri = A_deltas.begin(); tri < A_deltas.end(); ++tri){
    std::vector<vertex_vci>::const_iterator v_iter = tri->vertices_begin();
    std::map<vertex_vci, bool>::const_iterator kept;
    for( ; v_iter < tri->vertices_end(); ++v_iter){
      kept = kept_vertices.find(*v_iter);
      if(kept == kept_vertices.end()) break;
    }
    if(kept != kept_vertices.end()) total_SA += tri->area();
  }
  return total_SA;
}


void
TriMeshSphere::write(std::ostream &out, const interactionType act_type, 
                     const char delim) const
{
  out.setf(std::ios_base::fixed, std::ios_base::floatfield);
  out.precision(3);

  triangle_vci tri;
  for(tri = A_deltas.begin(); tri < A_deltas.end(); ++tri){
    if(act_type == ACCEPTOR) out << "ACCEPTOR" << delim;
    else if(act_type == DONOR) out << "DONOR" << delim;
    else if(act_type == DONEPTOR) out << "DONEPTOR" << delim;
    else if(act_type == METAL_1 || act_type == METAL_2) out << "METAL" << delim;

    std::vector<vertex_vci>::const_iterator vvci;
    for(vvci = tri->vertices_begin(); vvci < tri->vertices_end(); ++vvci){
      vertex_vci vert = *vvci;
      out << vert->pos[0] << " " << vert->pos[1] << " " << vert->pos[2] << delim
          << vert->dir[0] << " " << vert->dir[1] << " " << vert->dir[2] 
          << delim;
    }
    out << "\n";
  }
}

void
TriMeshSphere::write_msms_headers(std::ostream &vert_out, 
                                  std::ostream &face_out) const
{
  vert_out << "# SimSite3D hbond surfaces vertices ... \n"
           << "#faces  #sphere density probe_r -- #sphere is unknown \n"
           << std::setw(7) << 0 << std::setw(8) << 0 << "  0.00  0.00";

  face_out << "# SimSite3D hbond surfaces faces ... \n"
           << "#faces  #sphere density probe_r -- #sphere is unknown \n"
           << std::setw(7) << 0 << std::setw(8) << 0 << "  0.00  0.00";
} 


void
TriMeshSphere::write_msms_cap(std::ostream &vert_out, std::ostream &face_out,
                              const interactionType act_type) const
{
  vert_out.setf(std::ios_base::fixed, std::ios_base::floatfield);
  vert_out.precision(3);

  for(vertex_vci V = verts_begin(); V != verts_end(); ++V){
    vert_out << "\n";
    for(int i = 0; i < 3; ++i) vert_out << std::setw(9) << V->pos[i] << " ";
    for(int i = 0; i < 3; ++i) vert_out << std::setw(9) << V->dir[i] << " ";
    vert_out << "      0       0 ";

    // Abuse the last position by writing out the point "color"
    if(act_type == ACCEPTOR) vert_out << " 0";
    else if(act_type == DONOR) vert_out << " 1";
    else if(act_type == DONEPTOR) vert_out << " 2";
    else if(act_type == METAL_1 || act_type == METAL_2) vert_out << " 3";
  }

  triangle_vci tri;
  for(tri = A_deltas.begin(); tri < A_deltas.end(); ++tri){
    std::vector<vertex_vci>::const_iterator vvci;
    face_out << "\n";
    for(vvci = tri->vertices_begin(); vvci < tri->vertices_end(); ++vvci)
      face_out << std::setw(7) << 1 + (*vvci - verts_begin()) << " ";
    face_out << "0    0"; 
  }
}


void
TriMeshSphere::do_copy(const TriMeshSphere& src)
{
  A_nodes = src.A_nodes;
  A_local_points = src.A_local_points;
  A_triangles = src.A_triangles;
  A_npts = src.A_npts;
  A_num_tris = src.A_num_tris;

  // Rather than copying the vectors we must replicate them since
  // the iterators must point to nodes in A_nodes and not in src.A_nodes
  // That is, we need to compute the offset (index) for each vertex in
  // each triangle in src
  triangle_vci src_t;
  for(src_t = src.A_deltas.begin(); src_t != src.A_deltas.end(); ++src_t){
    std::vector<vertex_vci>::const_iterator src_t_v_iter;
    uint idz[3];
    uint i = 0;
    for(src_t_v_iter = src_t->vertices_begin();
        src_t_v_iter != src_t->vertices_end(); ++src_t_v_iter, ++i)
      idz[i] = *src_t_v_iter - src.A_nodes.begin();
    add_triangle(idz);
  }
  setup_vert_to_deltas_map();
} 

void
TriMeshSphere::add_triangle(const uint* idz)
{
  triangle_t delta(A_nodes.begin() + idz[0], A_nodes.begin() + idz[1],
                   A_nodes.begin() + idz[2]);
  A_deltas.push_back(delta);
  std::vector<vertex_vci>::const_iterator vp;
  for(vp = delta.vertices_begin(); vp < delta.vertices_end(); ++vp){
    half_edge_t tmp_edge;
    tmp_edge.tail = *vp;
    A_half_edges.push_back(tmp_edge);
  }
}

void
TriMeshSphere::setup_vert_to_deltas_map()
{ 
  triangle_vci tri;
  for(tri = A_deltas.begin(); tri < A_deltas.end(); ++tri){
    std::vector<vertex_vci>::const_iterator vv;
    for(vv = tri->vertices_begin(); vv < tri->vertices_end(); ++vv)
      A_verts_to_deltas_map[*vv].push_back(tri);
  }
}
