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
 * $Source: /psa/share/repository/pfizer_proj/src/search/DrugScoreInterface.C,v $
 * $Revision: 1.5 $
 * $Author: vanvoor4 $
 * $Date: 2007-11-01 16:39:44 $
 * 
 * $Log: not supported by cvs2svn $
 * Revision 1.4  2007/08/21 19:20:38  vanvoor4
 * Changed license header
 *
 * Revision 1.3  2007/06/29 13:46:36  vanvoor4
 * unlink() requires <unistd.h>
 *
 * Revision 1.2  2007/06/13 17:22:54  vanvoor4
 * Added support for pre and post options with respect to ligand and
 * protein specifications on scoring function command line
 *
 * Revision 1.1  2007/06/07 20:53:10  vanvoor4
 * initial checkin
 *
 * 
 * 
 */

#include <DrugScoreInterface.H>
#include <basics.H>
#include <unistd.h>

using namespace SimSite3D;
const std::string DrugScoreInterface::_fname = "DrugScoreInterface.C";

DrugScoreInterface::DrugScoreInterface(std::string cmdline)
 : ExternalScoringFunction(cmdline)
{
}

bool 
DrugScoreInterface::score(std::string prot_path, std::string lig_path, 
                          std::vector<my_float_t> *scores)
{
  // operating under the assumption that most scoring functions like to 
  // mimic change in energy (more negative is better score).
  scores->clear();
  if(!run(prot_path, lig_path)) return false;

  // Get the score from the DrugScore .scr file and clean up 
  size_t pos = lig_path.rfind("/") + 1;
  std::string lig_name = lig_path.substr(pos);
  lig_name = lig_name.substr(0, lig_name.length() - 5);
  std::string ds_prefix = "PAIR_10_" + lig_name;

  std::ifstream in;
  if(!open_ifstream(in, ds_prefix + ".scr")) return false;
  my_float_t tmp;
  in >> tmp;
  scores->push_back(tmp);
  in.close();

  unlink((ds_prefix + ".scr").c_str());
  unlink((ds_prefix + ".cor").c_str());
  unlink((ds_prefix + ".log").c_str());

  return true;
}

void 
DrugScoreInterface::sf_names(std::vector<std::string> *names)
{
  names->clear();
  names->push_back("DrugScore");
}
