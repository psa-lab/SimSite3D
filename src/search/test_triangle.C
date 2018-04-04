#include "TriMeshCap.H"

using namespace ASCbase;
using namespace ASCbase::geometry;

int main(){

  dir_point_storage<vertex_t> nodes;
  vertex_t tmp;
  my_float_t C[] = {0.0, 0.0, 0.0};
  tmp.pos[0] = 0.0;
  tmp.pos[1] = 0.0;
  tmp.pos[2] = 0.0;
  for(int i = 0; i < 3; ++i) C[i] += tmp.pos[i];
  nodes.push_back(tmp);
   
  tmp.pos[0] = 0.0;
  tmp.pos[1] = 1.0;
  tmp.pos[2] = 0.5;
  for(int i = 0; i < 3; ++i) C[i] += tmp.pos[i];
  nodes.push_back(tmp);

  tmp.pos[0] = 1.0;
  tmp.pos[1] = 0.0;
  tmp.pos[2] = 0.5;
  for(int i = 0; i < 3; ++i) C[i] += tmp.pos[i];
  nodes.push_back(tmp);

  for(int i = 0; i < 3; ++i) C[i] /= 3.0;
  
  triangle_t my_tri(nodes.begin(), nodes.begin() + 1, nodes.begin() + 2);
  if(my_tri.contains(C))
    std::cout << "Triangle contains its centroid\n";
  else
    std::cout << "Triangle does NOT contain its centroid\n";

  my_float_t U[3], V[3];
  vector(3, (nodes.begin()+1)->pos, nodes.begin()->pos, U);
  vector(3, (nodes.begin()+2)->pos, nodes.begin()->pos, V);
  my_float_t N[3];
  cross(U, V, N);
  std::cout << "N: " << N[0] << " " << N[1] << " " << N[2] << std::endl;
  cross(V, U, N);
  std::cout << "N: " << N[0] << " " << N[1] << " " << N[2] << std::endl;
  
}
