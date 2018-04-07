/******************************************************************************
 * Copyright (c) 2010, Michigan State University (MSU) Board of Trustees.
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

#ifndef ALIGNMENTS_FROM_RESULTS_FILE_HEADER_FILE
#define ALIGNMENTS_FROM_RESULTS_FILE_HEADER_FILE

namespace SimSite3D{

//! Class to score alignments in a given SimSite3D results file
/*! NOTE: this method is inefficient in that it will scan the entire results
 * file each time it is called -- this is to avoid the problem of very long
 * results files that cannot fit in memory.
 */
class AlignmentsFromResultsFile{
public:
  //! Cstor 
  AlignmentsFromResultsFile(const std::string res_fname)
  { 
    A_res_fname = res_fname; 
  }

  //! basic destruction
  virtual ~AlignmentsFromResultsFile()
  { ; }

  //! Required alignment method
  template<typename align_T> bool
  align(ModelSitemap *query_site, Sitemap& dset_site, 
        std::vector<align_T> *alignments)
  {
    std::ifstream res_file;
    if(!open_ifstream(res_file, A_res_fname)) return false;

    std::string dset_struct_id = dset_site.get_struct_id();
    for(std::string line; std::getline(res_file, line); ){
      if(line[0] == '#' || line.length() == 0) continue;

      std::vector<std::string> toks;
      string_tok(line, &toks, '|');
     
      std::string line_struct_id = toks[0];
      size_t pos = toks[0].rfind("_f.mol2");
      if(pos != std::string::npos) line_struct_id = toks[0].substr(0, pos);
      if(line_struct_id != dset_struct_id) continue;

      // NOTE: the matrices written in the results file are for multiplying
      // each position from the left.  HOWEVER, in SimSite3D the matrices are 
      // the transpose of those written in the results file as the rotation 
      // matrices are premultiplied by the coordinates (in other words the
      // coordinates are multiplied by rotation matrices from the right inside 
      // of SimSite3D).

      if(toks.size() < 4){
        std::cerr << "Too few fields in the alignment from " << A_res_fname 
                  << "\n" << "offending line: \"" << line << "\"" << std::endl;
        continue;
      }

      std::vector<std::string> R_toks, T_toks;
      string_tok(toks[2], &R_toks, ' ');
      string_tok(toks[3], &T_toks, ' ');
      if(R_toks.size() != 9 || T_toks.size() != 3){
        std::cerr << "Could not read the transformation from " << A_res_fname 
                  << "\n" << "offending line: \"" << line << "\"" << std::endl;
        continue;
      }

      align_T tmp_align;
      bool failed = false;
      for(int i = 0; i < 3 && !failed; ++i)
        for(int j = 0; j < 3 && !failed; ++j)
          if(!my_strtod(R_toks[3*i + j], &(tmp_align.R[3*j+i]))){
            std::cerr << "Could not read the rotation from " 
                      << A_res_fname << "\n" << "offending line: \"" << line 
                      << "\"" << std::endl;
            failed = true;
          }

      for(int i = 0; i < 3 && !failed; ++i)
        if(!my_strtod(T_toks[i], &(tmp_align.T[i]))){
          std::cerr << "Could not read the translation from " 
                    << A_res_fname << "\n" << "offending line: \"" << line 
                    << "\"" << std::endl;
          failed = true;
        }

      if(!failed) alignments->push_back(tmp_align);
    }
    return true;
  }

private:
  std::string A_res_fname;
};
} // end namespace SimSite3D
#endif
