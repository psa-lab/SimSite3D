#ifndef ALIGN_DATA_WITH_PRE_IK_DATA_HEADER_FILE_INCLUDED
#define ALIGN_DATA_WITH_PRE_IK_DATA_HEADER_FILE_INCLUDED

#include <AlignDataWithRigidData.H>

namespace ASCbase{

class align_w_preIK_data_t : public align_w_rigid_data_t{
public:
  align_w_preIK_data_t() : align_w_rigid_data_t()
  { init(); }

  //! Basic copy constructor
  align_w_preIK_data_t(const align_w_preIK_data_t &other) 
   : align_w_rigid_data_t(other)
  {
    init();
    do_copy(other);
  }

  //! Basic assignment operator
  const align_w_preIK_data_t &
  operator=(const align_w_preIK_data_t &other)
  {
    if(&other != this){
      do_copy(other);
      align_w_rigid_data_t::operator= (other);
    }
    return *this;
  }

  //! Do nothing destructor
  virtual ~align_w_preIK_data_t()
  { ; }

  virtual void
  save_preIK_vals()
  {
    preIK_score = score;
    std::copy(R, R + 9, preIK_R);
    std::copy(T, T + 3, preIK_T);
    preIK_terms.resize(terms.size());
    std::copy(terms.begin(), terms.end(), preIK_terms.begin());
  }

  virtual void
  get_score_field_labels(std::vector<std::string> *fields,
                         const bool normalize_score) const
  {
    rigid_align_t::get_score_field_labels(fields, normalize_score);
    fields->push_back("Initial site atomic RMSD");
    fields->push_back("ICP site atomic RMSD");
    fields->push_back("Final site atomic RMSD");
  }


  virtual void
  write_score_fields(std::ostream& out, const uint orient_num,
                     const bool wrote_ligs, const std::string& ext_SF_id_in,
                     const std::string& struct_id,
                     const std::string& lig_id) const
  {
    // First write the base class's score lines
    align_w_rigid_data_t::write_score_fields(out, orient_num, wrote_ligs,
                                             ext_SF_id_in, struct_id, lig_id);
    out << "\n";
  
    out << "preIK_" << struct_id << "|";

    // Score is likely accurate to at best 2 or 3 decimal places 
    out << std::fixed << std::setprecision(3) << std::left << std::setfill(' ')
        << preIK_score << "|";
  
    // Transformation (rotation followed by translation) 
    out << std::fixed << std::setprecision(14);
    // Our rotation matrices are to be muliplied from the left by the 
    // positions, but the output matrices are for post multiplication by 
    // vectors
    std::string follower = "        |";
    for(size_t i = 0; i < 3; ++i)
      for(size_t j = 0; j < 3; ++j){
        out << preIK_R[3*j + i] << follower[3*i + j];
      }
    for(uint i = 0; i < 3; ++i) out << preIK_T[i] << follower[6 + i];
  
    out << "||";
  
    if(preIK_terms.size()){
      out << std::fixed << std::setprecision(7) << preIK_terms[0];
      for(size_t i = 1; i < preIK_terms.size(); ++i) 
        out << " " << preIK_terms[i];
    }
    out << "|" << std::fixed << std::setprecision(3);

    if(A_init_rmsd < my_float_max) out << A_init_rmsd;
    out << "|";
    if(A_ICP_rmsd < my_float_max) out << A_ICP_rmsd;
    out << "|";
    if(A_final_rmsd < my_float_max) out << A_final_rmsd;
    out << "|" << std::fixed << std::setprecision(7);
  }

  //! Compute the atomic rmsd for the query and dataset residues 
  /*! WARNING! at the present this is only valid for proteins with the
   * same residues ...
   */
  virtual void
  set_site_atomic_rmsd(const my_float_t rmsd, const size_t rmsd_type)
  {
    if(rmsd_type == 0) A_init_rmsd = rmsd;
    else if(rmsd_type == 1) A_ICP_rmsd = rmsd;
    else if(rmsd_type == 2) A_final_rmsd = rmsd; 
    else{
      std::cerr << "Error: unknown rmsd type -- should be 0-2\n";
    }
  }


#if 0
  virtual void
  get_score_rigid_field_labels(std::vector<std::string> *fields,
                               const bool normalize_score) const
  {
    std::string tmp;
    if(normalize_score)
      tmp = "Normalized ASCbase initial alignment score of target to query";
    else tmp = "Raw ASCbase initial alignment score of target to query";
    fields->push_back(tmp);
  
    fields->push_back("Initial alignment rotation matrix (before ICP/IK");
    fields->push_back("Initial alignment translation vector (before ICP/IK");
    tmp = "Match print of the query's sitemap points satisfied by ";
    tmp += "sitemap points\n#     in the database hit";
    fields->push_back(tmp);
    tmp = "Ligand fragment binary string:  1 or 0 in nth position ";
    tmp += "implies that the\n#     nth mol2 ligand atom is or is not in ";
    tmp += "the mol2 ligand fragment (resp.)";
    fields->push_back(tmp);
    fields->push_back("Initial alignment features--terms (before ICP/IK");
  }
#endif

  my_float_t preIK_score;
  my_float_t preIK_R[9];
  my_float_t preIK_T[3];
  std::vector<my_float_t> preIK_terms;

private:
  my_float_t A_init_rmsd;
  my_float_t A_ICP_rmsd;
  my_float_t A_final_rmsd; 

  void
  init()
  {
    preIK_score = my_float_max;
    std::fill(preIK_R, preIK_R + 9, 0.0);
    std::fill(preIK_T, preIK_T + 3, 0.0);
    A_init_rmsd = my_float_max;
    A_ICP_rmsd = my_float_max;
    A_final_rmsd = my_float_max;
  }

  void
  do_copy(const align_w_preIK_data_t &other)
  {
    preIK_score = other.preIK_score;
    std::copy(other.preIK_R, other.preIK_R + 9, preIK_R);
    std::copy(other.preIK_T, other.preIK_T + 3, preIK_T);
    preIK_terms.resize(other.preIK_terms.size());
    std::copy(other.preIK_terms.begin(), other.preIK_terms.end(),
              preIK_terms.begin());
    A_init_rmsd = other.A_init_rmsd;
    A_ICP_rmsd = other.A_ICP_rmsd;
    A_final_rmsd = other.A_final_rmsd;
  }
};

} // end namespace ASCbase

#endif
