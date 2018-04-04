#ifndef ALIGN_W_TRI_PARAMS_HEADER_FILE_INCLUDED
#define ALIGN_W_TRI_PARAMS_HEADER_FILE_INCLUDED

#include <ScoreMapBase.H>

namespace ASCbase{

//! addition of triangle parameters to saved alignments
class align_w_tri_params_t : public rigid_align_t{
public:
  align_w_tri_params_t() : rigid_align_t()
  { init(); }

  //! Basic copy constructor
  align_w_tri_params_t(const align_w_tri_params_t &other) : rigid_align_t(other)
  {
    init();
    do_copy(other);
  }

  //! Basic assignment operator
  const align_w_tri_params_t &
  operator=(const align_w_tri_params_t &other)
  {
    if(&other != this){
      do_copy(other);
      rigid_align_t::operator= (other);
    }
    return *this;
  }

  //! Do nothing destructor
  virtual ~align_w_tri_params_t()
  { ; }

  virtual void
  write_score_fields(std::ostream& out, const uint orient_num,
                     const bool wrote_ligs, const std::string& ext_SF_id_in,
                     const std::string& struct_id, 
                     const std::string& lig_id) const 
  {
    rigid_align_t::write_score_fields(out, orient_num, wrote_ligs, ext_SF_id_in,
                                      struct_id, lig_id);
  }

  virtual void
  write_score_fields(std::ostream &out) const
  {
    rigid_align_t::write_score_fields(out);
    out << A_perimeter << " " << A_long_len << " " << A_short_len << "|"
        << A_num_orients << "|";
  }

  virtual void
  get_score_field_labels(std::vector<std::string> *fields, 
                         const bool normalize_score) const
  {
    rigid_align_t::get_score_field_labels(fields, normalize_score);
    fields->push_back("Query triangle perimeter, longest side, and shortest sitein Angstroms");
    fields->push_back("# of orientations");
  }

  virtual void
  set_triangle_params(const my_float_t perimeter, const my_float_t long_len,
                      const my_float_t short_len)
  {
    A_perimeter = perimeter;
    A_long_len = long_len;
    A_short_len = short_len;
  }

  virtual void
  set_number_of_orientations(const size_t num)
  {
    A_num_orients = num;
  }

private:
  my_float_t A_perimeter;
  my_float_t A_long_len;
  my_float_t A_short_len;
  size_t A_num_orients;

  void
  init()
  { 
    A_perimeter = 0.0;
    A_long_len = 0.0;
    A_short_len = 0.0;
    A_num_orients = 0;
  }

  void
  do_copy(const align_w_tri_params_t &other)
  {
    A_perimeter = other.A_perimeter;
    A_long_len = other.A_long_len;
    A_short_len = other.A_short_len;
    A_num_orients = other.A_num_orients;
  }
};

} // End namespace ASCbase
#endif

