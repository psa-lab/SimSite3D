#include <SurfDepsJoints.H>
#include <point_and_surf_score.H>
using namespace SimSite3D;

int main(int argc, char **argv)
{
  PDBStructure prot("/psa/results/SimSite3D_datasets/testing/adenines/proteins/1ua2_A_atp_p.pdb");

  ModelSitemap model("/psa/results/SimSite3D_datasets/testing/adenines/ade_pockets/1ua2_A_atp_s.csv");
  DbaseSitemap dbase("/psa/results/SimSite3D_datasets/testing/adenines/dbase/1b38_atp_s.csv");

  geometry::TransformableTrimesh prot_surf("/psa/results/SimSite3D_datasets/testing/adenines/ade_pockets/1ua2_A_atp_surf");
  geometry::ImmovableTrimesh db_surf("/psa/results/SimSite3D_datasets/testing/adenines/dbase/1b38_atp_surf");

  prot_joint_dep the_joints;
  
  SurfDepsJoints my_surf_to_joints(the_joints, &prot, &prot_surf);

  // For each block, compute ... just a check to see if it works
  my_surf_to_joints.compute_dihedrals_and_verts_J();

  // don't forget to update axis of rotations after each iteration
  const my_float_t *J = my_surf_to_joints.J();
  const int nrows = my_surf_to_joints.nrows_J();
  const int ncols = my_surf_to_joints.ncols_J();

  std::cout << "J is an " << nrows << " X " << ncols << " matrix\n";
  const my_float_t *jacob = J;
  for(int i = 0; i < nrows; ++i){
    std::cout << "\n"; 
    for(int j = 0; j < ncols; ++j, ++jacob){
      if(j > 0) std::cout << ", ";
      //std::cout << jacob[i*ncols + j];
      std::cout << *jacob;
    }
  }

  // Compute the gradient of the objective function
  my_float_t *grad_F = new my_float_t[ncols];
  std::fill(grad_F, grad_F + ncols, 0.0);

  // 1) Minimize the distance between corresponding surface points
  dbase_surf.compare(model_surf.vertices_begin(),
                     model_surf.number_of_vertices(), A_surf_closest_pts,
                     A_dists, A_max_surf_pt_dist);
  // 2) Minimize the overlap penalty for atoms
  // 3) Minimize angle deviations -- need to think about this one more



}
