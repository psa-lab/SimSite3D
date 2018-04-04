#include <prot_joint_dep.H>
using namespace ASCbase;

int main(int argc, char **argv)
{
  PDBStructure prot(argv[1]);
//  geometry::TransformableTrimesh prot_surf("/psa/results/SimSite3D_datasets/testing/adenines/ade_pockets/1ua2_A_atp_surf");

  prot_joint_dep joints_info;

  // Print a nice header  
 
  prot_joint_dep stuff;
  for(residue_vci R = prot.residues_begin(); R != prot.residues_end(); ++R){
    residue_joints joints(&prot, R, stuff.get_joints_info(R->name));
    const std::vector<my_float_t>& angles = joints.orig_joint_angles();
    std::vector<my_float_t>::const_iterator a;
    for(a = angles.begin(); a != angles.end(); ++a) std::cout << *a << " ";
    std::cout<< "\n";
    
  }
  // For each residue in prot: 
  //   const joint_deps_t j_info = joints_info.get_joints_info(res->resName);
  //   residue_joints res_joints(&prot, res, j_info);
  //   std::vector<my_float_t> angles;
  //   res_joints.compute_dihedral_angles(&angles);

  //   print the residue info and angles out

}
