//#include <stream_basics.H>
#include <PDBStructure.H>
#include <math_basics.H>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>

using namespace SimSite3D;

void
get_MC_positions(my_float_t *pos, const std::string fname, const int *res_nums, 
                 const int num_res);

//! align 2qx0 to each snapshot in md_idz and write out PH2 in that orientation
void
get_2qx0_ligs(const std::string snapshot_dir, 
              const int *md_idz, const int num_md_idz, 
              const int *res_nums, const int num_res);
//! Align the snapshots to each other -- writeout alignment in files to use
//! as input to the search method
void
get_transforms(const std::string snapshot_dir, 
               const std::string aligns_file_dir,
               const int *md_idz, const int num_md_idz, 
               const int *res_nums, const int num_res);

int main(int argc, char **argv)
{
  std::string snapshot_dir = "/psa/cukier/suli/YP/pdb_5to7ns_YP_wATP_normal_fittedTotheCore";

  int res_nums[] = {43, 44, 45, 46, 54, 55, 56, 96, 98, 122, 123, 124, 125};
  int num_res = 13; 

  int md_idz[] = {1452, 1490, 1542, 1602, 2850, 1539, 1600, 1471, 1468, 1536};

  std::string aligns_dir = "/psa/results/SimSite3D_datasets/testing";
  aligns_dir += "/new_YP_HPPK_MD/cluster0/align_files";

  get_transforms(snapshot_dir, aligns_dir, md_idz, 10, res_nums, 13);


#if 0
  for(int i = 1; i < 3000; ++i){
    std::stringstream fname1;
    fname1 << snapshot_dir << "/pdb_" << i << "ps.pdb";
    PDBStructure pdb_1(fname1.str());

    residue_vci res = pdb_1.residues_begin();
    int res_num_idx = 0;
    int pos_idx = 0;
    for( ; res < pdb_1.residues_end() && res_num_idx < num_res; ++res){
      if(res->number < res_nums[res_num_idx]) continue;

      atom_vci a = res->atoms_begin;
      for( ; a < res->atoms_end; ++a){
        if(a->name == N || a->name == CA || a->name == C || a->name == O){
          std::copy(a->pos, a->pos +3, &pos_one[pos_idx]);
          pos_idx += 3;       
        }
      }

      ++res_num_idx;
    }

    for(int j = i+1; j < 3000; ++j){
      std::stringstream fname2;
      fname2 << snapshot_dir << "/pdb_" << j << "ps.pdb";
      PDBStructure pdb_2(fname2.str());

      res = pdb_2.residues_begin();
      res_num_idx = 0;
      pos_idx = 0;
      for( ; res < pdb_2.residues_end() && res_num_idx < num_res; ++res){
        if(res->number < res_nums[res_num_idx]) continue;
  
        atom_vci a = res->atoms_begin;
        for( ; a < res->atoms_end; ++a){
          if(a->name == N || a->name == CA || a->name == C || a->name == O){
            std::copy(a->pos, a->pos +3, &pos_two[pos_idx]);
            pos_idx += 3;       
          }
        }
  
        ++res_num_idx;
      }

      Quaternion Q;
      my_float_t T[3];
      lse_3D_fit(pos_one, pos_two, 4*num_res, &Q, T, 0);
      my_float_t R[9];
      Q.get_ortho_rot_mat(R);

      std::copy(pos_two, pos_two + 12*num_res, scratch);
      for(my_float_t *p = pos_two; p < pos_two + 12*num_res; p+=3) 
        std::copy(T, T+3, p);
      my_gemm(4*num_res, 3, 3, 1.0, scratch, 3, R, 3, pos_two, 3, 1.0);
      //my_gemm

      my_float_t rmsd = 0.0;
      for(int zz = 0; zz < 12*num_res; ++zz){
        my_float_t tmp = pos_one[zz] - pos_two[zz];
        rmsd += tmp*tmp;
        //std::cout << "stuff: " << rmsd << "\n";
      }
      rmsd /= (4.0*num_res);
      rmsd = std::sqrt(rmsd);

      std::cout << "MYRMSD: " << i << "," << j << "," << rmsd << ",\n";
    }
  }
#endif
}

void
get_MC_positions(my_float_t *pos, const std::string fname, const int *res_nums, 
                 const int num_res)
{
  PDBStructure pdb(fname);

  residue_vci res = pdb.residues_begin();
  int res_num_idx = 0;
  int pos_idx = 0;
  for( ; res < pdb.residues_end() && res_num_idx < num_res; ++res){
    if(res->number < res_nums[res_num_idx]) continue;

    atom_vci a = res->atoms_begin;
    for( ; a < res->atoms_end; ++a){
      if(a->name == N || a->name == CA || a->name == C || a->name == O){
        std::copy(a->pos, a->pos +3, &pos[pos_idx]);
        pos_idx += 3;       
      }
    }
    ++res_num_idx;
  }
}

void
get_2qx0_ligs(const std::string snapshot_dir, 
              const int *md_idz, const int num_md_idz, 
              const int *res_nums, const int num_res)
{
  my_float_t *pos_one = new my_float_t[num_res *4 *3];
  my_float_t *pos_two = new my_float_t[num_res *4 *3];
  my_float_t *scratch = new my_float_t[num_res *4 *3];

  PDBStructure my_2qx0("/psa/results/SimSite3D_datasets/testing/new_YP_HPPK_MD/2qx0_A.pdb");
  get_MC_positions(pos_one, "/psa/results/SimSite3D_datasets/testing/new_YP_HPPK_MD/2qx0_A.pdb", res_nums, num_res);

  for(int i = 0; i < num_md_idz; ++i){
    std::stringstream fname1;
    fname1 << snapshot_dir << "/pdb_" << md_idz[i] << "ps.pdb";
    get_MC_positions(pos_two, fname1.str(), res_nums, num_res);

    Quaternion Q;
    my_float_t T[3];
    lse_3D_fit(pos_two, pos_one, 4*num_res, &Q, T, 0);
    my_float_t R[9];
    Q.get_ortho_rot_mat(R);
    my_2qx0.transform(R,T);

    std::stringstream lig_fname;
    lig_fname << "/psa/results/SimSite3D_datasets/testing/new_YP_HPPK_MD/"
              << "cluster0/ligands/pdb_" << md_idz[i] << "_PH2_l.pdb";
    std::ofstream out;
    open_ofstream(out, lig_fname.str());

    atom_vci a;
    for(a = my_2qx0.atoms_begin(); a < my_2qx0.atoms_end(); ++a){
      if(a->res_str != "PH2") continue;

      out.setf(std::ios_base::fixed, std::ios_base::floatfield);
      out.precision(3);
      out << "HETATM" << std::setw(5) << a->atom_num << " " << a->name_str
          << a->altLoc << a->res_str << " "
          << a->chainID << std::setw(4) << a->res_num << a->iCode << "   ";
      for(uint j = 0; j < 3; ++j) out << std::setw(8) << a->pos[j];
      out.precision(2);
      out << std::setw(6) << a->occupancy
          << std::setw(6) << a->tempFactor << "\n";
    }
    my_2qx0.revert();
  }
}

void
get_transforms(const std::string snapshot_dir, 
               const std::string aligns_file_dir,
               const int *md_idz, const int num_md_idz, 
               const int *res_nums, const int num_res)
{
  my_float_t *pos_one = new my_float_t[num_res *4 *3];
  my_float_t *pos_two = new my_float_t[num_res *4 *3];
  my_float_t *scratch = new my_float_t[num_res *4 *3];

  //PDBStructure my_2qx0("/psa/results/SimSite3D_datasets/testing/new_YP_HPPK_MD/2qx0_A.pdb");
  //get_MC_positions(pos_one, "/psa/results/SimSite3D_datasets/testing/new_YP_HPPK_MD/2qx0_A.pdb", res_nums, num_res);

  for(int i = 0; i < num_md_idz; ++i){
    std::stringstream fname1;
    fname1 << snapshot_dir << "/pdb_" << md_idz[i] << "ps.pdb";
    get_MC_positions(pos_one, fname1.str(), res_nums, num_res);

    std::ofstream out;
    std::stringstream ostr;
    ostr << aligns_file_dir << "/pdb_" << md_idz[i] << "_PH2_xforms.out";
    open_ofstream(out, ostr.str());
    out << "# Aligns file -- for HPPK MD\n";

    for(int j = 0; j < num_md_idz; ++j){
      std::stringstream fname2;
      fname2 << snapshot_dir << "/pdb_" << md_idz[j] << "ps.pdb";
      get_MC_positions(pos_two, fname2.str(), res_nums, num_res);

      Quaternion Q;
      my_float_t T[3];
      lse_3D_fit(pos_one, pos_two, 4*num_res, &Q, T, 0);
      my_float_t R[9];
      Q.get_ortho_rot_mat(R);

      out << "pdb_" << md_idz[j] << "_PH2|100.0|";
      // Our rotation matrices are to be muliplied from the left by the 
      // positions, but the output matrices are for post multiplication by 
      // vectors
      out << std::fixed << std::setprecision(14);
      std::string follower = "        |";
      for(size_t iii = 0; iii < 3; ++iii)
        for(size_t jjj = 0; jjj < 3; ++jjj){
          out << R[3*jjj + iii] << follower[3*iii + jjj];
        }
      for(uint iii = 0; iii < 3; ++iii) out << T[iii] << follower[6 + iii];
      out << "\n";
    }
  }
}
