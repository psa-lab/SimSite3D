#ifndef ALIGN_DATA_WITH_RIGID_DATA_HEADER_FILE_INCLUDED
#define ALIGN_DATA_WITH_RIGID_DATA_HEADER_FILE_INCLUDED

#include <ScoreMapBase.H>

namespace ASCbase{

class align_w_rigid_data_t : public rigid_align_t{
public:
  align_w_rigid_data_t() : rigid_align_t()
  { init(); }

  //! Basic copy constructor
  align_w_rigid_data_t(const align_w_rigid_data_t &other) : rigid_align_t(other)
  {
    init();
    do_copy(other);
  }

  //! Basic assignment operator
  const align_w_rigid_data_t &
  operator=(const align_w_rigid_data_t &other)
  {
    if(&other != this){
      do_copy(other);
      rigid_align_t::operator= (other);
    }
    return *this;
  }

  //! Do nothing destructor
  virtual ~align_w_rigid_data_t()
  { ; }

  virtual void
  save_rigid_alignment_vals()
  {
    rigid_score = score;
    std::copy(R, R + 9, rigid_R);
    std::copy(T, T + 3, rigid_T);
    rigid_terms.resize(terms.size());
    std::copy(terms.begin(), terms.end(), rigid_terms.begin());
  }

  virtual void
  write_score_fields(std::ostream& out, const uint orient_num,
                     const bool wrote_ligs, const std::string& ext_SF_id_in,
                     const std::string& struct_id, 
                     const std::string& lig_id) const
  {
    // First write the base class's score line
    rigid_align_t::write_score_fields(out, orient_num, wrote_ligs, ext_SF_id_in,
                                      struct_id, lig_id);
    out << "\n";

    out << "init_" << struct_id << "|";

    // Score is likely accurate to at best 2 or 3 decimal places 
    out << std::fixed << std::setprecision(3) << std::left << std::setfill(' ')
        << rigid_score << "|";
  
    // Transformation (rotation followed by translation) 
    out << std::fixed << std::setprecision(14);
    // Our rotation matrices are to be muliplied from the left by the 
    // positions, but the output matrices are for post multiplication by 
    // vectors
    std::string follower = "        |";
    for(size_t i = 0; i < 3; ++i)
      for(size_t j = 0; j < 3; ++j){
        out << rigid_R[3*j + i] << follower[3*i + j];
      }
    for(uint i = 0; i < 3; ++i) out << rigid_T[i] << follower[6 + i];
  
    out << "||";
  
    if(rigid_terms.size()){
      out << std::fixed << std::setprecision(7) << rigid_terms[0];
      for(size_t i = 1; i < rigid_terms.size(); ++i) out << " " << rigid_terms[i];
    }
    out << "|" << std::fixed << std::setprecision(7);
  }

/*
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
*/

private:
  my_float_t rigid_score;
  my_float_t rigid_R[9];
  my_float_t rigid_T[3];
  std::vector<my_float_t> rigid_terms;

  void
  init()
  {
    rigid_score = my_float_max;
    std::fill(rigid_R, rigid_R + 9, 0.0);
    std::fill(rigid_T, rigid_T + 3, 0.0);
  }

  void
  do_copy(const align_w_rigid_data_t &other)
  {
    rigid_score = other.rigid_score;
    std::copy(other.rigid_R, other.rigid_R + 9, rigid_R);
    std::copy(other.rigid_T, other.rigid_T + 3, rigid_T);
    rigid_terms.resize(other.rigid_terms.size());
    std::copy(other.rigid_terms.begin(), other.rigid_terms.end(),
              rigid_terms.begin());
  }
};

} // end namespace ASCbase

#endif
