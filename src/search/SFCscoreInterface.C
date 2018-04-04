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
 * $Source: /psa/share/repository/pfizer_proj/src/search/SFCscoreInterface.C,v $
 * $Revision: 1.5 $
 * $Author: vanvoor4 $
 * $Date: 2007-11-01 16:39:32 $
 * 
 * $Log: not supported by cvs2svn $
 * Revision 1.4  2007/08/21 19:31:01  vanvoor4
 * changed license header
 *
 * Revision 1.3  2007/06/29 13:46:52  vanvoor4
 * unlink() requires <unistd.h>
 *
 * Revision 1.2  2007/06/13 17:21:07  vanvoor4
 * Modified to reflect the requirement of .spf file to follow the
 * protein and ligand specifications
 *
 * Revision 1.1  2007/06/07 20:54:13  vanvoor4
 * initial checkin
 *
 * 
 * 
 */

#include <SFCscoreInterface.H>
#include <basics.H>
#include <sstream>
#include <unistd.h>

using namespace ASCbase;
const std::string SFCscoreInterface::_fname = "SFCscoreInterface.C";

SFCscoreInterface::SFCscoreInterface(std::string cmdline)
 : ExternalScoringFunction(cmdline)
{
  // Look through the post_opts string for the .spf file
  std::string::size_type pos = cmdline.find(".spf");
  if(pos == std::string::npos){
    std::ostringstream out;
    out << "Unable to find a .spf file specified in the command line string\n"
        << "(column 2 of $ASCBASE_SOFTWARE_DIR/ASCbaseSoftParams/"
        << "external_scoring_functions.txt)";
    err_msg(_fname, "SFCscoreInterface", out.str());
  }else{
    std::string::size_type beg = cmdline.find_last_of(" ", pos) + 1;
    get_scoring_functions(cmdline.substr(beg, pos + 1 - beg)); 
  }
}

bool 
SFCscoreInterface::score(std::string prot_path, std::string lig_path, 
                         std::vector<my_float_t> *scores)
{
  // Not sure what we should use here -- depends on SFCscore
  scores->clear();
  if(!run(prot_path, lig_path)) return false;

  // Get the score from the SFCscore .scores file and clean up 
  size_t pos = prot_path.rfind("/") + 1;
  std::string prot_name = prot_path.substr(pos);

  for(uint i = 0; i < sfc_names.size(); ++i){
    std::string pref = prot_name + "." + sfc_names[i];
    std::ifstream in;
    if(!open_ifstream(in, pref + ".scores")) return false;

    std::string line;
    while(std::getline(in, line)){
      // Ignore comments
      pos = line.find("#");
      if(pos != std::string::npos) line.erase(pos);
      if(line.length() == 0) continue; 
   
      std::string prot, lig;
      my_float_t tmp; 
      std::istringstream lstr(line);
      lstr >> prot >> lig >> tmp;
      // Assume only one score since we are operating under the assumption 
      // that we are scoring one protein ligand pair.
      if(!lstr.fail()){
        scores->push_back(tmp);
        break;
      }
    }

    in.close();
    // Clean up the score files
    unlink((pref + ".scores").c_str());
    unlink((pref + ".contrib").c_str());
  }

  // Clean up the description files
  unlink((prot_name + ".desc").c_str());
  unlink((prot_name + ".raw_desc").c_str());

  return true;
}

void 
SFCscoreInterface::sf_names(std::vector<std::string> *names)
{
  names->resize(sfc_names.size());
  std::copy(sfc_names.begin(), sfc_names.end(), names->begin());
}

bool
SFCscoreInterface::get_scoring_functions(std::string spf_fname)
{
  std::ifstream spf_fstr;
  if(!open_ifstream(spf_fstr, spf_fname)) return false;

  std::string line;
  while(std::getline(spf_fstr, line)){
    // Ignore comments
    size_t pos = line.find("#");
    if(pos != std::string::npos) line.erase(pos);

    // look for the score line(s)
    std::istringstream lstr(line);
    std::string tok;
    lstr >> tok;
    if(tok != "SCORE") continue;
    lstr >> tok;
    while(!lstr.fail()){
      sfc_names.push_back(tok);       
      lstr >> tok;
    }
  }

  spf_fstr.close();
  return true;
}
