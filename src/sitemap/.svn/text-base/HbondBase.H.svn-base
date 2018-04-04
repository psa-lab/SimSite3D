#ifndef HBOND_BASE_HEADER_FILE_INCLUDED
#define HBOND_BASE_HEADER_FILE_INCLUDED

#include <hbond_triads.H>

namespace ASCbase{

class HbondBase{
public:

  HbondBase()
  {
    build_info_map();
  }

  virtual ~HbondBase()
  {

  }

  void build_info_map()
  {
    if(A_local_hbond_info.size()) return;

    for(uint i = 0; i < num_new_hbond_triads; ++i){
      if(new_hbond_triads[i].residue == MAIN_CHAIN){
        for(uint j = 0; j < A_handled_residues.size(); ++j){
          if(!new_hbond_triads[i].B == N && !A_handled_residues[j] == PRO)
            A_local_hbond_info[A_handled_residues[j]][new_hbond_triads[i].B].push_back(new_hbond_triads[i]);
        }
      }else
        A_local_hbond_info[new_hbond_triads[i].residue][new_hbond_triads[i].B].push_back(new_hbond_triads[i]);
    }
  }

  static std::vector<residue_type>&
  handled_residues()
  {
    if(A_handled_residues.size() == 0){
      A_handled_residues.push_back(ALA);
      A_handled_residues.push_back(ARG);
      A_handled_residues.push_back(ASN);
      A_handled_residues.push_back(ASP);
      A_handled_residues.push_back(CYS);
      A_handled_residues.push_back(GLN);
      A_handled_residues.push_back(GLU);
      A_handled_residues.push_back(GLY);
      A_handled_residues.push_back(HIS);
      A_handled_residues.push_back(ILE);
      A_handled_residues.push_back(LEU);
      A_handled_residues.push_back(LYS);
      A_handled_residues.push_back(MET);
      A_handled_residues.push_back(PHE);
      A_handled_residues.push_back(PRO);
      A_handled_residues.push_back(SER);
      A_handled_residues.push_back(THR);
      A_handled_residues.push_back(TRP);
      A_handled_residues.push_back(TYR);
      A_handled_residues.push_back(VAL);
    }
    return A_handled_residues;
  }

  const std::vector<hbond_triad_t>&
  get_hbond_info(const residue_type res, const atom_type atom)
  { return A_local_hbond_info[res][atom]; }


protected:
  std::map<residue_type, std::map<atom_type, std::vector<hbond_triad_t> > > 
    A_local_hbond_info;

private:
  static std::vector<residue_type> A_handled_residues;
};

} // end namespace ASCbase

#endif
