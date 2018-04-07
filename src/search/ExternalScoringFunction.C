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

/*
 * $Source: /psa/share/repository/pfizer_proj/src/search/ExternalScoringFunction.C,v $
 * $Revision: 1.4 $
 * $Author: vanvoor4 $
 * $Date: 2007-11-01 16:41:32 $
 * 
 * $Log: not supported by cvs2svn $
 * Revision 1.3  2007/08/21 19:21:37  vanvoor4
 * Changed license header
 *
 * Revision 1.2  2007/06/13 17:22:22  vanvoor4
 * Added support for pre and post options with respect to ligand and
 * protein specifications on scoring function command line
 *
 * Revision 1.1  2007/06/07 20:53:38  vanvoor4
 * initial checkin
 *
 * 
 * 
 */

#include <sstream>
#include <cstdlib>
#include <stream_basics.H>
#include <ExternalScoringFunction.H>
#include <string_basics.H>

using namespace SimSite3D;

ExternalScoringFunction::ExternalScoringFunction(std::string cmdline)
{
  a_fail = true;
  int lig_pos = -1, prot_pos = -1;

  for(int i = 0; i < 2; ++i){
    size_t pos = cmdline.find_first_of("$");
    if(pos == std::string::npos) return;
    cmdline_toks.push_back(cmdline.substr(0, pos));
    cmdline.erase(0, pos);

    if(lig_pos == -1 && cmdline.substr(0, 7) == "$LIGAND"){
      cmdline_toks.push_back("");
      lig_pos = cmdline_toks.size() - 1;
      cmdline.erase(0, 7);
    }else if(prot_pos == -1 && cmdline.substr(0, 8) == "$PROTEIN"){
      cmdline_toks.push_back("");
      prot_pos = cmdline_toks.size() - 1;
      cmdline.erase(0, 8);
    }
  }
  if(prot_pos == -1 || lig_pos == -1) return;
  if(cmdline.length()) cmdline_toks.push_back(cmdline);
  prot_str_iter = cmdline_toks.begin() + prot_pos;
  lig_str_iter = cmdline_toks.begin() + lig_pos;

  a_fail = false;
}

bool
ExternalScoringFunction::run(std::string prot_path, std::string lig_path)
{
  if(a_fail) return false;
  if(!normal_file_exists(prot_path) || !normal_file_exists(lig_path)) 
    return false; 

  std::ostringstream sys_str;
  std::vector<std::string>::const_iterator s_iter;
  for(s_iter = cmdline_toks.begin(); s_iter != cmdline_toks.end(); ++s_iter){
    if(s_iter == prot_str_iter) sys_str << prot_path;
    else if(s_iter == lig_str_iter) sys_str << lig_path;
    else sys_str << *s_iter;
  }
  system(sys_str.str().c_str()); 
  return true;
}
