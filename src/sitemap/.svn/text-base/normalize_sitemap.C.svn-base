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

/*
 * $Source: /psa/share/repository/pfizer_proj/src/gen_points/normalize_sitemap.C,v $
 * $Revision: 1.4 $
 * $Author: vanvoor4 $
 * $Date: 2007-11-01 19:09:05 $
 * 
 * $Log: not supported by cvs2svn $
 * Revision 1.3  2007/11/01 16:11:50  vanvoor4
 * Updated to reflect the location of the diverse sitemaps database
 * and change to handling parameters
 *
 * Revision 1.2  2007/10/18 04:27:55  vanvoor4
 * Added slash
 *
 * Revision 1.1  2007/08/21 18:35:43  vanvoor4
 * Moved from ../basics
 *
 *
 *
 */

#include <sstream>
#include <errno.h>
#include <cmath>
#include <iomanip>
#include <iterator>
#include <unistd.h>
#include <cstring>
#include <BaseParameters.H>
#include <basics.H>

bool calc_mean_and_stdev(std::string fname, my_float_t *mean, 
                         my_float_t* stdev);

bool write_stats(std::string fname, my_float_t mean, my_float_t stdev,
                 const bool hydrophobic_query);


bool 
normalize_sitemap(const std::string site_fname, 
                  const ASCbase::BaseParameters& params,
                  const bool hydrophobic_query)
{
  std::cout << "Calculating normalization stats for " << site_fname << "\n";
  
  // Make a temporary directory to store the "hits" against the diverse set
  // Lazily convert string to cstring so that mkdtemp can modify the cstring
  std::string tmp_dir = params.scratch_dir +  "/.ASCbase_tmp_XXXXXX";
  char* output_dir = new char[tmp_dir.length() + 1];
  std::strcpy(output_dir, tmp_dir.c_str());
  mkdtemp(output_dir);
  sleep(1);
  if(!dir_exists(output_dir)) return false;

  // Under the current scoring method we probably do not want a score
  // threshold for the normalization step.
  std::string search_prog = params.install_dir + "/bin/search_sitemaps";
  std::string res_fname = std::string(output_dir) + "/diverse_blast.out";
  std::string cmd = search_prog + " --dbase_sites " + params.diverse_sites;
  cmd += " --dbase_ligs " + params.diverse_ligs + " --score_threshold 100.0 ";
  cmd += std::string(" --proj_output ") + output_dir + " --no_normalization";
  if(params.require_min_npts == false) cmd += " --allow_small_site_maps ";
  if(hydrophobic_query) cmd += " --hydrophobic_query ";
  cmd += " --prot_lig_score NONE -o " + res_fname + " " + site_fname;

  // system returns -1 on error and search-sitemaps should return 0 upon
  // success
  int rv = system(cmd.c_str());
  if(rv != 0){
    int errsv = errno;
    std::ostringstream msg;
    msg << "Unable to run the normalization of the sitemap " << site_fname;

    // Assume fork failed
    if(rv == -1) msg << ": " << strerror(errsv);
    // Assume search_sitemaps failed
    else msg << ": search_sitemaps failed";

    err_msg("normalize_sitemap.C", "normalize_sitemap", msg.str(), std::cerr);
    return false;  
  }
  
  my_float_t mean, stdev;
  if(!calc_mean_and_stdev(res_fname, &mean, &stdev)) return false;
  if(!write_stats(site_fname, mean, stdev, hydrophobic_query)) return false;

  // cleanup  -- the "easy" way
  cmd = std::string("/bin/rm -rf ") + output_dir;
  system(cmd.c_str());
 
  return true;
}

bool 
calc_mean_and_stdev(std::string fname, my_float_t *mean, my_float_t* stdev)
{
  *mean = 0;
  *stdev = 0;

  std::ifstream in; 
  if(!open_ifstream(in, fname)) return false;

  my_float_t sum = 0, sum_of_squares = 0;
  uint n = 0;

  std::string line;
  while(std::getline(in, line)){
    if(strip_trailing_comments(&line, "%#") == 0) continue;

    std::vector<std::string> toks;
    string_tok(line, &toks, '|');

    // Want the ASCbase score field
    my_float_t score;
    if(my_strtof(toks[1], &score)){
      ++n;
      sum += score;
      sum_of_squares += score*score;  
    }
  }
  if(n == 0){
    warn("normalize_sitemap.C", "calc_mean_and_stdev",
         "Unable to calculate normalization statistics; search_sitemaps failed to produce any hits to the diverse sitemaps database", std::cerr);
    return false; 
  }
  
  *mean = sum / n;
  *stdev = std::sqrt(sum_of_squares / n - (*mean) * (*mean));
  
  return true;
}

bool 
write_stats(std::string fname, my_float_t mean, my_float_t stdev,
            const bool hydrophobic_query)
{
  std::ifstream in;
  if(!open_ifstream(in, fname)) return false;

  std::ostringstream norm_line;
  norm_line << std::scientific << std::setprecision(8);
  norm_line << "  " << mean << " " << stdev;

  std::vector<std::string> lines;
  std::string line;
  while(std::getline(in, line)){
    if(!hydrophobic_query && 
       line[0] == '<' && line.substr(0,21) == "<normalization_stats>"){
      lines.push_back(line);
      std::getline(in, line);
      line = norm_line.str(); 
    }else if(hydrophobic_query &&
             line[0] == '<' && line.substr(0,18) == "<hphob_norm_stats>"){
      lines.push_back(line);
      std::getline(in, line);
      line = norm_line.str(); 
    }
    lines.push_back(line);
  }
  in.close();

  // Write lines
  std::ofstream out;
  if(!open_ofstream(out, fname)) return false;
  std::copy(lines.begin(), lines.end(), 
            std::ostream_iterator<std::string>(out, "\n"));
  out.close();

  return true;
}
