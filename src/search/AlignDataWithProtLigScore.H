#ifndef ALIGN_DATA_WITH_PROT_LIG_SCORE_HEADER_FILE_INCLUDED
#define ALIGN_DATA_WITH_PROT_LIG_SCORE_HEADER_FILE_INCLUDED

#include <ScoreMapBase.H>

namespace SimSite3D{

class align_w_pl_score_t : public rigid_align_t{
public:
  align_w_pl_score_t() : rigid_align_t()
  { init(); }

  align_w_pl_score_t(const align_w_pl_score_t &other) : rigid_align_t(other)
  {
    init();
    do_copy(other);
  }

  const align_w_pl_score_t&
  operator=(const align_w_pl_score_t &other)
  {
    if(&other != this){
      do_copy(other);
      rigid_align_t::operator=(other);
    }
    return *this;
  }

  virtual ~align_w_pl_score_t()
  { ; }

  virtual void
  write_score_fields(std::ostream& out, const uint orient_num,
                     const bool wrote_ligs, const std::string& ext_SF_id_in,
                     const std::string& struct_id, 
                     const std::string& lig_id) const
  {
    // First write the base class's score line
    rigid_align_t::write_score_fields(out, orient_num, wrote_ligs, ext_SF_id_in,
                                      struct_id, lig_id);
    out << A_orientscore << " " << A_affiscore << " " 
        << A_ligand_efficiency << "|";
  }


  virtual void
  write_score_fields(std::ostream &out) const
  { 
    rigid_align_t::write_score_fields(out);
    out << A_orientscore << " " << A_affiscore << " " 
        << A_ligand_efficiency << "|";
  }

  virtual void
  get_score_field_labels(std::vector<std::string> *fields,
                         const bool normalize_score) const
  {
    rigid_align_t::get_score_field_labels(fields, normalize_score);
    fields->push_back("PL Score: Orientscore, Affiscore, and Ligand "
                      "Efficiency");
  }

  virtual void
  set_prot_lig_scores(const my_float_t orientscore, const my_float_t affiscore,
                      const my_float_t ligand_efficiency)
  {
    A_orientscore = orientscore;
    A_affiscore = affiscore;
    A_ligand_efficiency = ligand_efficiency;
  }

  virtual bool
  compute_prot_lig_score()
  { return true; }

private:
  my_float_t A_affiscore;
  my_float_t A_orientscore;
  my_float_t A_ligand_efficiency;

  void
  init()
  {
    A_affiscore = my_float_max;
    A_orientscore = my_float_max;
    A_ligand_efficiency = my_float_max;
  }

  void
  do_copy(const align_w_pl_score_t &other)
  {
    A_affiscore = other.A_affiscore;
    A_orientscore = other.A_orientscore;
    A_ligand_efficiency = other.A_ligand_efficiency;
  }

};
} // end namespace SimSite3D
#endif
