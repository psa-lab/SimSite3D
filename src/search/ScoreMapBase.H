/******************************************************************************
 * Copyright (c) 2006,2007, Michigan State University (MSU) Board of Trustees.
 *   All rights reserved.
 *
 * This file is part of the ASCbase Software project.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 * Authors: Jeffrey Van Voorst, vanvoor4@msu.edu
 *          Leslie Kuhn, Ph.D., KuhnL@msu.edu 
 *****************************************************************************/
#ifndef SCORE_MAP_BASE_HEADER_FILE
#define SCORE_MAP_BASE_HEADER_FILE

#include <mol2File.H>
#include <DrugScoreInterface.H>
#include <SFCscoreInterface.H>
#include <SearchParameters.H>
#include <sstream>
#include <iomanip>

namespace ASCbase{

//! Data class: Holds the numerous features of an alignment
/*! This class was supposed to make it easier to swap things in and out.
 *  However, it was put together relatively quickly and its major design
 *  flaw is that one cannot pick and choose score fields.  It might be
 *  useful for future development if one can pick and choose pieces for
 *  scoring and reporting.
 */
class rigid_align_t{
public:

  //! Initialize the features to zero or very large values
  rigid_align_t()
  { init(); }

  //! Basic copy constructor
  rigid_align_t(const rigid_align_t &other)
  { 
    init();
    do_copy(other); 
  }

  //! Basic assignment operator
  const rigid_align_t &
  operator=(const rigid_align_t &other)
  {
    if(&other != this){
      do_copy(other);
    }
    return *this;
  }

  //! Do nothing destructor
  virtual ~rigid_align_t()
  { ; }

  //! Two alignments are equal if their score is the same and they have
  //! the same transformation
  inline bool 
  operator== (const rigid_align_t &b) const
  {
    if(score != b.score) return false;
    for(size_t i = 0; i < 9; ++i)
      if(R[i] != b.R[i]) return false;
    for(size_t i = 0; i < 3; ++i)
      if(T[i] != b.T[i]) return false;
    return true;
  }

  //! Two alignments are not equal if their score is not the same or they do
  //! not have the same transformation
  inline bool 
  operator!= (const rigid_align_t &b) const
  {
    if(score != b.score) return true;
    for(size_t i = 0; i < 9; ++i)
      if(R[i] != b.R[i]) return true;
    for(size_t i = 0; i < 3; ++i)
      if(T[i] != b.T[i]) return true;
    return false;
  }

  virtual void
  write_score_fields(std::ostream& out, const uint orient_num,
                     const bool wrote_ligs, const std::string& ext_SF_id_in,
                     const std::string& struct_id, 
                     const std::string& lig_id) const ;

  virtual void
  write_score_fields(std::ostream& out) const;

/*
  virtual void
  write_score_rigid_fields(std::ostream &out) const
  {
    std::cerr << "rigid_align_t does not support writing rigid scores\n";
  }

  virtual void
  write_preIK_fields(std::ostream &out) const
  {
    std::cerr << "rigid_align_t does not support writing pre-IK fields\n";
  }
*/

  virtual void
  get_score_field_labels(std::vector<std::string> *fields, 
                         const bool normalize_score) const;

  virtual void
  set_prot_lig_scores(const my_float_t affiscore, const my_float_t orientscore,
                      const my_float_t ligand_efficiency)
  {
    std::cerr << "rigid_align_t does not support saving prot-lig scores\n";
  }

  virtual void
  set_triangle_params(const my_float_t perimeter, const my_float_t long_len,
                      const my_float_t short_len)
  {
  //std::cerr << "rigid_align_t does not support setting triangle params\n";
  }

  virtual void
  set_number_of_orientations(const size_t num)
  {

  }

  virtual void
  save_rigid_alignment_vals()
  {
  //std::cerr << "rigid_align_t does not support saving rigid alignment vals\n";
  }

  virtual void
  save_preIK_vals()
  {
  //std::cerr << "rigid_align_t does not support saving pre-IK vals\n";
  }

  virtual bool
  compute_prot_lig_score()
  { return false; }

  virtual void
  set_site_atomic_rmsd(const my_float_t rmsd, const size_t rmsd_type)
  {
  }
  
  
  my_float_t tier1_score;      //!< ASCbase fast orientation score (sieve)
  my_float_t score;            //!< Final ASCbase score -- tier2
  std::vector<my_float_t> terms;
  std::vector<my_float_t> ext_scores;  //!< External score(s) 
  std::vector<bool> frag_atoms_flags;  
  my_float_t R[9];       //!< Rotation matrix for search to query
  my_float_t T[3];       //!< Translation vector for search to query
  Quaternion Q;
  mol2File *frag_file;   //!< Pointer to the ligand fragment in query pocket
  std::vector<bool> match_print;

  // Added for prototype of hbond caps - used in only ASCbasePy as this point
  my_float_t hb_caps_score;
  std::string hb_caps_match_print;

private:

  void
  init()
  {
    score = my_float_max;
    tier1_score = my_float_max;
    frag_file = 0;
    hb_caps_score = my_float_max;
  }

  //! A straightforward copy method
  void
  do_copy(const rigid_align_t &other)
  {
    tier1_score = other.tier1_score;
    score = other.score;
    terms = other.terms;
    ext_scores = other.ext_scores;
    frag_atoms_flags = other.frag_atoms_flags;
    std::copy(other.R, other.R + 9, R);
    std::copy(other.T, other.T + 3, T);
    Q = other.Q;
    // It can be very dangerous to copy the pointer for frag file, at the
    // present this is only done in the scoring section and the local copy
    // of this class is tossed after insertion into a container
    frag_file = other.frag_file;
    match_print = other.match_print;

    hb_caps_score = other.hb_caps_score;
    hb_caps_match_print = other.hb_caps_match_print; 
  }

};


//! Base class for classes which use multimaps to facility the scoring and 
//! and storing of alignments of target sitemaps to a given model (query).
/*! Because of the gcc model, we must have all template definitions 
 * included in all the source files which use the functions, classes,
 * etc that are prototyped and defined using the template keyword.
 */
class ScoreMapBase{
public:

  //! Default constructor for scoring methods which use multimaps to
  //! store the alignments.
  /*!
   * @param params Reference to the search's runtime parameters
   */
  ScoreMapBase(const SearchParameters& params)
  {
    A_external_SF = 0;
    A_write_ligands = params.write_ligands;
    A_scratch_dir = params.scratch_dir;
    A_prot_lig_score = params.do_internal_prot_lig_score;
  }

  //! basic destruction
  virtual ~ScoreMapBase()
  {
    if(A_external_SF) delete(A_external_SF);
    A_external_SF = 0;
    ext_SF_id = "";
    A_scratch_dir = "";
  }

  /*! Setup the external scoring function.  This function is in this class
   * since it acts as a storage class for the alignments
   *
   * @param sf_name Path to the external scoring function file
   * @param method_name Name of the external scoring function (as stored in
   *        the external_scoring_functions.txt file)
   * @return True if method was found, else false
   */
  bool 
  set_external_scoring_method(const std::string sf_name,
                              std::string method_name)
  {
    std::string src_fname = "ScoreMapBase.H";

    // Read system scoring functions file
    std::ifstream in;
    if(!open_ifstream(in, sf_name)) return false;
    std::string str;
    bool first_method = true;
    while(std::getline(in, str)){
      if(str[0] == '%' || str[0] == '#') continue;
  
      strip_trailing_comments(&str, "%#");
      std::vector<std::string> toks; 
      string_tok(str, &toks, ',');
      if(toks.size() != 3){
        bad_file_line(src_fname, "set_external_scoring_method", sf_name, str);
        continue; 
      }

      // The default method is the one listed first in the external scoring 
      // function file
      if(first_method && method_name == "default"){
        method_name = toks[1];
        first_method = false;
      }
  
      // toks[0]: The user specified prot-lig sf identifier (on commandline)
      // toks[1]: The internal scoring function type -- class to construct
      // toks[2]: The prot-lig sf cmdline with variables for prot and lig file
      //          names
      if(toks[1] == method_name){
        ext_SF_id = toks[0];
        if(method_name == "NONE"){
          A_external_SF = 0;
          return true;
        }else if(toks[1] == "DrugScore"){
          A_external_SF = new DrugScoreInterface(toks[2]);
          if(!A_external_SF->fail()) return true;
        }else if(toks[1] == "SFCscore"){
          A_external_SF = new SFCscoreInterface(toks[2]);
          if(!A_external_SF->fail()) return true;
        }else{
          std::cout << "Warning:  Unknown external scoring function \""
                    << toks[0] << "\" skipping . . . \n";
        }
      }
    }
  
    if(!A_write_ligands && A_external_SF){
      err_msg(src_fname, "set_external_scoring_method",
              "Using a prot-lig scoring function requires writing the ligand files");
      return false;
    }
  
    if(A_external_SF && A_external_SF->fail()){
      delete(A_external_SF);
      A_external_SF = 0;
      std::ostringstream ostr;
      ostr << "Unable to initialize the external scoring function named "
           << method_name;
      err_msg(src_fname, "set_external_scoring_method", ostr.str(), std::cerr);
      return false;
    }
  
    std::ostringstream ostr;
    ostr << "Unable to find an external scoring function in " << sf_name
         << " named " << method_name;
    err_msg(src_fname, "set_external_scoring_method", ostr.str(), std::cerr);
    return false;
  }
  
  //! Run the external scoring function over all the saved alignments
  /*! For a given ASCbase run, the model sitemap is held constant.  Thus,
   * we need only one protein_file.  For each saved alignment, score the 
   * protein and the saved ligand fragment.
   *
   * @param protein_file Protein PDB corresponding to the model sitemap
   * @param output_dir Directory holding the ligand fragments
   * @return True if the scoring function was initialized, else false
   */
  template<class align_mmap_iter> bool
  run_external_scoring_method(align_mmap_iter a_beg, align_mmap_iter a_end,
                              const std::string &prot_fname,  
                              const std::string &struct_id, 
                              const std::string &output_dir)
  {
    if(!A_external_SF) return false;

    uint cnt = 1;
    std::string prefix = output_dir + "/ligand_fragments/" + struct_id + "_";

    for(align_mmap_iter align = a_beg; align != a_end; ++align, ++cnt){
      std::ostringstream tmp;
      tmp << std::setfill('0') << std::setw(5) << cnt << "_f.mol2";
      A_external_SF->score(prot_fname, prefix + tmp.str(),
                           &(align->second->ext_scores));
    }
    return true;
  }

  //! Get the name(s) of the external scoring function(s)
  void 
  get_ext_SF_names(std::vector<std::string>* names)
  { A_external_SF->sf_names(names); }

protected:
  std::string A_scratch_dir;
  bool A_write_ligands;  //!< True => write ligs and frags, else do not write
  bool A_prot_lig_score; //!< True ==> compute query prot -- dbase lig frag affi and orient scores
  std::string ext_SF_id;  //!< Label given to the ext. prot-lig scoring fcn

private:
  ExternalScoringFunction* A_external_SF; //!< external prot-lig score fcn ptr

};
} // end namespace ASCbase
#endif
