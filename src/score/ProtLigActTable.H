/******************************************************************************
 * Copyright (c) 2006,2007, Michigan State University (MSU) Board of Trustees.
 *   All rights reserved.
 *
 * This file is part of the SimSite3D Software project.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 * Authors: Jeffrey Van Voorst, jeff.vanvoorst@gmail.com
 *          Leslie Kuhn, Ph.D., KuhnL@msu.edu 
 *****************************************************************************/

#include <ProtLigScore.H>
#include <ProtLigScoreParameters.H>

namespace SimSite3D{

class ProtLigActTable{

public:
  typedef struct{
    std::string lig_id;
    my_float_t affi;
    my_float_t efficiency;
    my_float_t orient;
    std::vector<int> atom_idz;
  }score_t;

  // -- NOTE::::::::
  // the ligand iterators will be invalidated anyhow -- we should only keep
  // protein atom info -- is protein atom indexes faster than protein atom
  // iterators?
  
  ProtLigActTable(PDBStructure &prot)
  {
    A_prot = &prot;
    build_hphob_SC_to_atom_map(prot);
    A_col_sums.resize(prot.atoms_end() - prot.atoms_begin());
    std::fill(A_col_sums.begin(), A_col_sums.end(), 0);
  }   

  ~ProtLigActTable() { ; }

  void
  add_row(ProtLigScore &s_info)
  {
    // Assume we only have prot-lig salt bridges and hbonds and atom A is 
    // always from the protein and atom B is always a ligand atom
     
    std::string tmp_lig_name = s_info.lig_name();
    size_t slash_pos = tmp_lig_name.rfind("/");
    if(slash_pos != std::string::npos) 
      tmp_lig_name = tmp_lig_name.substr(slash_pos + 1);
    size_t dot_pos = tmp_lig_name.rfind(".");
    if(dot_pos != std::string::npos)
      tmp_lig_name = tmp_lig_name.substr(0, dot_pos);
     
    score_t tmp;
    tmp.lig_id = tmp_lig_name;
    tmp.affi = s_info.affiscore();
    tmp.efficiency = s_info.ligand_efficiency();
    tmp.orient = s_info.orientscore();

    // Salt bridges
    std::vector<ProtLigScore::salt_bridge_t>::const_iterator sb_i;
    sb_i = s_info.salt_bridges_begin();
    for( ; sb_i < s_info.salt_bridges_end(); ++sb_i){
      size_t idx = sb_i->A - A_prot->atoms_begin();
      ++(*(A_col_sums.begin() + idx));
      tmp.atom_idz.push_back(idx); 
    }

    // Hbonds
    std::vector<ProtLigScore::hbond_t>::const_iterator hb_i;
    hb_i = s_info.hbonds_begin();
    for( ; hb_i < s_info.hbonds_end(); ++hb_i){
      size_t idx = hb_i->A - A_prot->atoms_begin();
      ++(*(A_col_sums.begin() + idx));
      tmp.atom_idz.push_back(idx); 
    }

    // Hphob interactions
    std::map<residue_vci, bool>::const_iterator map_i;
    map_i = s_info.hphob_SC_begin();
    for( ; map_i != s_info.hphob_SC_end(); ++map_i){
      atom_vci SC_CB = hphob_SC_to_atom(map_i->first);
      if(SC_CB == atom_t::NULL_ATOM_VCI) continue;

      size_t idx = SC_CB - A_prot->atoms_begin();
      ++(*(A_col_sums.begin() + idx));
      tmp.atom_idz.push_back(idx); 
    }
    A_scores.push_back(tmp);
  }

  void
  print_table(std::ostream &out, char delim = '|')
  {
    // Print table headings
    out << "ligand ID" << delim << "AffiScore" << delim << "ligand efficiency"
        << delim << "OrientScore" << delim;

    std::vector<int> column_atom_idz;
    std::vector<int>::const_iterator col_sum;
    for(col_sum = A_col_sums.begin(); col_sum < A_col_sums.end(); ++col_sum){
      if(*col_sum == 0) continue;

      int idx = col_sum - A_col_sums.begin();
      column_atom_idz.push_back(idx);
      atom_vci atom = A_prot->atoms_begin() + idx;
      
      if(atom->chainID != ' ') out << "(" << atom->chainID << ")";
      out << PDB_residues::residue_to_string(atom->res) << atom->res_num;
      if(atom->iCode != ' ') out << atom->iCode;

      // CB is placeholder atom for hydrophobic sidechain
      if(atom->name == CB) out << " SC" << delim;
      // Polar interaction
      else out << " " << PDB_residues::atom_to_string(atom->name) << delim;
    }
    out << "\n";

    // Print out the rows
    out.precision(3);  // 3 significant digits seems reasonable
    std::vector<score_t>::iterator score_iter = A_scores.begin();
    for( ; score_iter < A_scores.end(); ++score_iter){
      out << score_iter->lig_id << delim << score_iter->affi << delim
          << score_iter->efficiency << delim << score_iter->orient << delim;

      // Sort the hit atoms & side chains CBs and get a unique set of the
      // sorted indices
      std::sort(score_iter->atom_idz.begin(), score_iter->atom_idz.end());
      std::vector<int>::const_iterator hit_atom_idx, hit_atoms_end;
      hit_atom_idx = score_iter->atom_idz.begin();
      hit_atoms_end = 
        std::unique(score_iter->atom_idz.begin(), score_iter->atom_idz.end());

      std::vector<int>::const_iterator col_idx = column_atom_idz.begin();
      for( ; col_idx < column_atom_idz.end(); ++col_idx){
        if(*col_idx == *hit_atom_idx){
          if(hit_atom_idx >= score_iter->atom_idz.end()){
            std::cerr << "Too many hit atom indices \n";
            return;
          }
          out << "1" << delim;
          ++hit_atom_idx; 
        }else if(*col_idx > *hit_atom_idx && 
                 hit_atom_idx < hit_atoms_end){
          std::cerr << "Column indices do not match hit atom indices\n";
          return;
        }else out << "0" << delim;
      }
      out << "\n";
    }
  }

  atom_vci
  hphob_SC_to_atom(const residue_vci res)
  {
    std::map<residue_vci, atom_vci>::const_iterator my_iter;
    my_iter = A_hphob_SC_to_atom.find(res);
    if(my_iter == A_hphob_SC_to_atom.end()){
      if(res->name != GLY && res->name != SER)
        std::cerr << "Unknown residue in hphob_SC_to_atom\n";
      return atom_t::NULL_ATOM_VCI;
    }
    return my_iter->second;
  }
private:
  void
  build_hphob_SC_to_atom_map(PDBStructure &prot);

  const PDBStructure *A_prot;
  std::vector<int> A_col_sums;
  std::map<residue_vci, atom_vci> A_hphob_SC_to_atom;
  std::vector<score_t> A_scores;
};
#if 0
int main(const int argc, const char **argv)
{
  std::cout << "\n" << argv[0] << " (" << PACKAGE_NAME << ") " 
            << PACKAGE_VERSION << "\n\n";

  // Do not return -1 here since system, fork, etc return -1 on failure, and we
  // wish to distinguish between system and program failure
  ProtLigScoreParameters my_params(argc, argv);
  BaseParameters::status_t status = my_params.status();
  if(status == BaseParameters::DISPLAY_HELP_ONLY) return 0;
  else if(status != BaseParameters::READY){
    std::cerr << "\n" << argv[0]
              << " *FAILED* \n\tCould not initialize parameters\n";
    return 1;
  }
  
  PDBStructure prot(my_params.prot_fname);
  if(my_params.lig_fname.length()){
    mol2File lig(my_params.lig_fname);
    lig.calc_charge_sums();
    ProtLigScore my_score(prot, lig);
    my_score.report_scores(std::cout);
  }else{
    std::ifstream lig_list_file;
    open_ifstream(lig_list_file, my_params.lig_list_fname);
    if(lig_list_file.fail()) return 1;

    std::vector<score_t> scores;
    typedef std::map<residue_vci, std::vector<bool> > hphob_map_t;
    typedef std::map<atom_vci, std::vector<bool> > hbond_map_t;
    std::map<residue_vci, std::vector<bool> > prot_hphob_SCs;
    std::map<atom_vci, std::vector<bool> > prot_hbond_atoms;
    for(std::string line; std::getline(lig_list_file, line); ){
      mol2File lig(line);
      lig.calc_charge_sums();
      ProtLigScore my_score(prot, lig);
      my_score.report_scores(std::cout);
      if(my_params.build_interact_tbl){
        score_t s;
        s.lig_id = line;
        s.affi = my_score.affiscore();
        s.affieff = my_score.affiefficientscore();
        s.orient = my_score.orientscore(); 
        scores.push_back(s);
        // Assumption -- atom A in hbonds and salt bridges are always from the
        // protein and atom B is always from the ligand
        std::vector<ProtLigScore::salt_bridge_t>::const_iterator sb_iter;
        sb_iter = my_score.salt_bridges_begin();
        for( ; sb_iter < my_score.salt_bridges_end(); ++sb_iter){
          hbond_map_t::iterator hb = prot_hbond_atoms.lower_bound(sb_iter->A);
        }
       
      }
    }
    if(my_params.build_interact_tbl){
    
    }
  }
  return 0; 
}
#endif
} // end namespace SimSite3D
