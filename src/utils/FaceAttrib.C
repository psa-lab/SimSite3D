#include <FaceAttrib.H>

using namespace ASCbase;
using namespace ASCbase::geometry;

std::vector<FaceAttrib> FaceAttrib::NULL_FACES_STORAGE;
const FaceAttrib::vi FaceAttrib::NULL_VI = FaceAttrib::NULL_FACES_STORAGE.end();

void
FaceAttrib::initialize()
{
  std::fill(A_V, A_V + 3, VertAttrib::NULL_VI);
  std::fill(A_fnei, A_fnei + 3, FaceAttrib::NULL_VI);
  A_area = 0;
}

void
FaceAttrib::closest_point(const my_float_t *P, const my_float_t prev_best_d,
                          my_float_t *d, 
                          my_float_t *closest_pt, my_float_t *corr_N) const
{
  // 1) Project P onto the face plane
  // ASSUMPTION: A_normal is valid!
  
  // Vector from V[0] to P
  my_float_t V0_P[3];
  vector(3, P, A_V[0]->pos, V0_P);
  
  // Signed distance from pt to plane
  *d = dot(A_normal, V0_P);

  // Check for early termination
  my_float_t abs_prev_best_d = 
    (prev_best_d > 0.0 ? prev_best_d : -1.0 * prev_best_d);
  my_float_t abs_d = (*d > 0.0 ? *d : -1.0 * (*d));
  if(abs_d >= abs_prev_best_d) return;

  // Projection of point onto plane (by subtracting the distance from the 
  // point in the normal direction)
  my_float_t *proj_pt = closest_pt;
  std::copy(P, P + 3, proj_pt);
  my_axpy(3, -1 *(*d), A_normal, 1, proj_pt, 1);

  // 2) Get the barycentric coordinates of the projection of P
  my_float_t b[3];
  bool pt_in_face = barycentric_coords(proj_pt, b);

  // 3) If the projected point is not in the face, get the closest point
  // on the perimeter of the face to the projected point
  if(not pt_in_face){
    std::fill(proj_pt, proj_pt + 3, 0.0);
    for(int i = 0; i < 3; ++i)
      if(b[i] != 0.0) my_axpy(3, b[i], A_V[i]->pos, 1, proj_pt, 1);
    *d = dist(proj_pt, P);
  }

  // 4) Do a simple interpolation for the normal
  if(corr_N){
    std::fill(corr_N, corr_N + 3, 0.0);
    for(int i = 0; i < 3; ++i)
      if(b[i] != 0.0) my_axpy(3, b[i], A_V[i]->dir, 1, corr_N, 1);
  }
}

#if 0
// hmm thought maybe that we could use closest point for faster computation,
// but then I forgot how to compute the normal ...
//
//! Assumes that p is in the face plane
void
FaceAttrib::closest_point(const my_float_t *p, my_float_t *cp) const
{
  // Compute vectors corresponding to the edges -- note that
  // vectors have an associated direction but edges do not
  my_float_t V1_V2[3], V0_V2[3], V0_V1[3];
  vector(3, A_V[2]->pos, A_V[1]->pos, V1_V2);
  vector(3, A_V[2]->pos, A_V[0]->pos, V0_V2);
  vector(3, A_V[1]->pos, A_V[0]->pos, V0_V1);

  // Use cross product to compute twice the area of the triangle
  // U denotes cross products (vectors)
  my_float_t U[3];
  cross(V0_V1, V0_V2, U);
  
  // Compute 3 cross products that will give us the area of the 
  // "smaller" triangles that correspond to the barycentric coordinates
  my_float_t V0_p[3], V1_p[3];
  vector(3, p, A_V[0]->pos, V0_p);
  vector(3, p, A_V[1]->pos, V1_p);
  my_float_t U0[3], U1[3], U2[3];
  cross(V1_V2, V1_p, U0);
  cross(V0_p, V0_V2, U1); 
  cross(V0_V1, V0_p, U2); 
  
  // Using the dot product we can tell if the "areas" are positive or
  // negative.  By using pairs of dot products we can determine the
  // location of the point with respect to the face.  Because the
  // actual dot products correspond to parallel or antiparallel vectors,
  // numerical issues shouldn't be an issue
  my_float_t sign_U[3];
  sign_U[0] = dot(U, U0);
  sign_U[1] = dot(U, U1); 
  sign_U[2] = dot(U, U2);
  
  // NOTE: we cannot have all 3 signs be negative
  // If two signs are negative "project" to the "remaining" vertex
  if(sign_U[0] < 0.0 && sign_U[1] < 0.0)
    std::copy(A_V[2]->pos, A_V[2]->pos + 3, cp);
  else if(sign_U[0] < 0.0 && sign_U[2] < 0.0)
    std::copy(A_V[1]->pos, A_V[1]->pos + 3, cp);
  else if(sign_U[1] < 0.0 && sign_U[2] < 0.0)
    std::copy(A_V[0]->pos, A_V[0]->pos + 3, cp);
  else{

    // If we get here we have at most one sign that is negative
    // and the areas get shared between the two remaining triangles
    // Ugghhh -- all these methods seem to require sqrt which is 
    // quite expensive
    my_float_t A[3];
    if(sign_U[0] < 0.0){
      A[1] = std::sqrt(dot(U1,U1));
      A[2] = std::sqrt(dot(U2,U2));
      b[1] = A[1] / (A[1] + A[2]);
      b[2] = A[2] / (A[1] + A[2]);
    }else if(sign_U[1] < 0.0){
      A[0] = std::sqrt(dot(U0,U0));
      A[2] = std::sqrt(dot(U2,U2));
      b[0] = A[0] / (A[0] + A[2]);
      b[2] = A[2] / (A[0] + A[2]);
    }else if(sign_U[2] < 0.0){
      A[0] = std::sqrt(dot(U0,U0));
      A[1] = std::sqrt(dot(U1,U1));
      b[0] = A[0] / (A[0] + A[1]);
      b[1] = A[1] / (A[0] + A[1]);
    // If we get here then no signs are negative and the point is inside
    // the triangle
    }else{
      std::copy(p, p + 3, cp);
    }
  }
}
#endif

bool
FaceAttrib::barycentric_coords(const my_float_t *p, my_float_t *b) const
{
  // Compute vectors corresponding to the edges -- note that
  // vectors have an associated direction but edges do not
  my_float_t V1_V2[3], V0_V2[3], V0_V1[3];
  vector(3, A_V[2]->pos, A_V[1]->pos, V1_V2);
  vector(3, A_V[2]->pos, A_V[0]->pos, V0_V2);
  vector(3, A_V[1]->pos, A_V[0]->pos, V0_V1);
/*
  std::cout << "P: " << p[0] << " " << p[1] << " " << p[2] << "\n";
  for(int i = 0; i < 3; ++i)
    std::cout << "V[" << i << "]: " << A_V[i]->pos[0] << " "
              << A_V[i]->pos[1] << " " << A_V[i]->pos[2] << "\n";
*/

  // Use cross product to compute twice the area of the triangle
  // U denotes cross products (vectors)
  my_float_t U[3];
  cross(V0_V1, V0_V2, U);
  
  // Compute 3 cross products that will give us the area of the 
  // "smaller" triangles that correspond to the barycentric coordinates
  my_float_t V0_p[3], V1_p[3];
  vector(3, p, A_V[0]->pos, V0_p);
  vector(3, p, A_V[1]->pos, V1_p);
  my_float_t U0[3], U1[3], U2[3];
  cross(V1_V2, V1_p, U0);
  cross(V0_p, V0_V2, U1); 
  cross(V0_V1, V0_p, U2); 
  
  // Using the dot product we can tell if the "areas" are positive or
  // negative.  By using pairs of dot products we can determine the
  // location of the point with respect to the face.  Because the
  // actual dot products correspond to parallel or antiparallel vectors,
  // numerical issues shouldn't be an issue
  my_float_t sign_U[3];
  sign_U[0] = dot(U, U0);
  sign_U[1] = dot(U, U1); 
  sign_U[2] = dot(U, U2);
  //std::cout << "Sign of Uz: " << sign_U[0] << " " << sign_U[1] << " "  << sign_U[2] << "\n";
  
  // NOTE: we cannot have all 3 signs be negative
  // If two signs are negative "project" to the "remaining" vertex
  std::fill(b, b + 3, 0.0);
  if(sign_U[0] < 0.0 && sign_U[1] < 0.0){
    b[2] = 1.0;
    return false;
  }else if(sign_U[0] < 0.0 && sign_U[2] < 0.0){
    b[1] = 1.0;
    return false;
  }else if(sign_U[1] < 0.0 && sign_U[2] < 0.0){
    b[0] = 1.0;
    return false;
  }

  // If we get here we have at most one sign that is negative
  // and the areas get shared between the two remaining triangles
  // My idea about ratio of areas was either wrong or must take into
  // account the negative area of the remaining triangle 
  my_float_t A[3];
  if(sign_U[0] < 0.0){
    my_float_t edge_len = normalize(V1_V2);
    // Projection of V1_p onto the unit vector V1_V2
    my_float_t proj = dot(V1_p, V1_V2);
    if(proj <= 0.0) b[1] = 1.0;
    else if(proj >= edge_len) b[2] = 1.0;
    else{
      my_float_t t = proj / edge_len;
      b[2] = t;
      b[1] = 1.0 - t;
    }
  }else if(sign_U[1] < 0.0){
    my_float_t edge_len = normalize(V0_V2);
    // Projection of V0_p onto the unit vector V0_V2
    my_float_t proj = dot(V0_p, V0_V2);
    if(proj <= 0.0) b[0] = 1.0;
    else if(proj >= edge_len) b[2] = 1.0;
    else{
      my_float_t t = proj / edge_len;
      b[2] = t;
      b[0] = 1.0 - t;
    }
  }else if(sign_U[2] < 0.0){
    my_float_t edge_len = normalize(V0_V1);
    // Projection of V0_p onto the unit vector V0_V1
    my_float_t proj = dot(V0_p, V0_V1);
    if(proj <= 0.0) b[0] = 1.0;
    else if(proj >= edge_len) b[1] = 1.0;
    else{
      my_float_t t = proj / edge_len;
      b[1] = t;
      b[0] = 1.0 - t;
    }

  // If we get here then no signs are negative and the point is inside
  // the triangle
  }else{
    my_float_t face_area = dot(U,U);
    A[0] = dot(U0,U0);
    A[1] = dot(U1,U1);
    A[2] = dot(U2,U2);
    //std::cout << "(2 * total area)^2: " << face_area << "\n";
    //std::cout << "(2 * areas)^2: " << A[0] << " "<< A[1] << " " << A[2] << "\n";
    for(int i = 0; i < 3; ++i) b[i] = std::sqrt(A[i] / face_area);
    return true;
  }
  return false;
}

void
FaceAttrib::point_data_t::init()
{
  x = 0;
  face = FaceAttrib::NULL_VI;
  dist = my_float_max;
  std::fill(B, B+3, 0.0);
}

