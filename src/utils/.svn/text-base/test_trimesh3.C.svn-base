#include <ImmovableTrimesh.H>
#include <ImmovableTrimeshThree.H>

using namespace SimSite3D::geometry;

int main(int argc, char **argv)
{
#if 1
  ImmovableTrimeshThree my_mesh3("/psa/results/SimSite3D_datasets/testing/adenines/ade_pockets/1ua2_A_atp_surf");
  ImmovableTrimesh my_mesh("/psa/results/SimSite3D_datasets/testing/adenines/ade_pockets/1ua2_A_atp_surf");
//  std::vector<Edge> edges;
//  my_mesh.read_face_file("simple_test", &edges, true, true);
  std::cout << "number of vertices: " << my_mesh3.num_verts() << "\n";
//  my_mesh.initialize_opt(edges);
#endif

  for(VertAttrib::vci V = my_mesh3.vertices_begin(); V != my_mesh3.vertices_end(); ++V){
    FaceBins::bin_vci bin = my_mesh3.get_bin(V->pos);
    std::cout << "\nnext V\n";

    typedef std::vector< FaceAttrib::vci >::const_iterator face_vec_vci;

    for(face_vec_vci F_i = bin->begin(); F_i < bin->end(); ++F_i){
      std::cout << "Face index: " << *F_i - my_mesh3.faces_begin() << "\n";
      std::cout << "V: " << V->pos[0] << " " 
               << V->pos[1] << " " << V->pos[2] << "\n" ;
      my_float_t B[3], pt_in_plane[3];
      my_float_t d = (*F_i)->proj_pt_to_plane(V->pos, pt_in_plane);
      std::cout << "Dist to plane: " << d << "\n";
      if((*F_i)->barycentric_coords(pt_in_plane, B)){
        std::cout << "B: " << B[0] << " " << B[1] << " " << B[2] << "\n";
      }
    }

  }














}
