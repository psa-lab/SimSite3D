/******************************************************************************
 * Copyright (c) 2007, Michigan State University (MSU) Board of Trustees.
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
 * Authors: Jeffrey Van Voorst, vanvoor4@msu.edu
 *          Leslie Kuhn, Ph.D., KuhnL@msu.edu 
 *****************************************************************************/

/*
 * $Source: /psa/share/repository/pfizer_proj/src/basics/PDB_metals.C,v $
 * $Revision: 1.3 $
 * $Author: vanvoor4 $
 * $Date: 2008-02-26 18:20:16 $ 
 * 
 * $Log: not supported by cvs2svn $
 * Revision 1.2  2008/01/04 20:54:13  vanvoor4
 * Added a check so that we can use this to query if an atom is a metal.
 *
 * Revision 1.1  2007/12/17 21:09:16  vanvoor4
 * initial checkin
 * 
 * 
 */

#include <PDB_metals.H>
#include <string_basics.H>
#include <sstream>

using namespace SimSite3D;

const std::string PDB_metals::A_fname = "PDB_metals.C";

const my_float_t PDB_metals::MAX_METAL_1_SQUARED_LENGTH = 2.9 * 2.9;
const my_float_t PDB_metals::MAX_METAL_2_SQUARED_LENGTH = 2.6 * 2.6;
const my_float_t PDB_metals::MIN_METAL_1_SQUARED_LENGTH = 2.0 * 2.0;
const my_float_t PDB_metals::MIN_METAL_2_SQUARED_LENGTH = 1.7 * 1.7;

// These distances are the same as those used in SLIDE 3.0
const uint PDB_metals::A_metal_array_size = 11;
const pdb_metal_info_t PDB_metals::A_metal_array[] = {
  {"CA  ", CALCIUM, METAL_1, 2.4, 2.0, 2.9 },
  {"CO  ", CO, METAL_2, 1.9, 1.7, 2.6 },
  {"CU  ", CU, METAL_2, 2.1, 1.7, 2.6 },
  {"CD  ", CADMIUM, METAL_2, 2.2, 1.7, 2.6 },
  {"FE  ", FE, METAL_2, 2.2, 1.7, 2.6 },
  {" K  ", K,  METAL_1, 2.4, 2.0, 2.9 },
  {"MG  ", MG, METAL_2, 2.1, 1.7, 2.6 },
  {"MN  ", MN, METAL_2, 2.2, 1.7, 2.6 },
  {"NA  ", NA, METAL_1, 2.4, 2.0, 2.9 },
  {"NI  ", NI, METAL_2, 2.2, 1.7, 2.6 },
  {"ZN  ", ZN, METAL_2, 2.1, 1.7, 2.6 }
};

std::map<atom_type, std::vector<std::string> > PDB_metals::A_metal_res_names;


const pdb_metal_info_t* 
PDB_metals::lookup(const std::string metal_name, const std::string metal_res,
                   metal_warn_type msg_type)
{
  if(metal_name.length() != 4){
    err_msg(A_fname, "lookup(str)", 
            "PDB metals names must contain exactly 4 characters");
    return 0;
  }else if(metal_res.length() != 3){
    err_msg(A_fname, "lookup(str)", 
            "PDB metals residue's id must contain exactly 3 characters");
    return 0;
  }

  // Check if the residue identifier for this metal atom is one of the
  // residue ids for metal ions -- we do not want to include metal atoms
  // which are part of a larger molecule as part of the protein.
  if(A_metal_res_names.size() == 0) set_metal_res_names();
  const pdb_metal_info_t* m = A_metal_array; 
  for(uint i = 0; i < A_metal_array_size; ++i, ++m) 
    if(m->pdb_metal_str == metal_name){
      const std::vector<std::string> &valid_res_names = 
        A_metal_res_names[m->metal_name];
      std::vector<std::string>::const_iterator res_name;
      res_name = valid_res_names.begin();
      bool found = false;
      for( ; res_name < valid_res_names.end() && !found; ++res_name)
        if(*res_name == metal_res) found = true;
      if(!found){
        std::ostringstream msg;
        msg << "The residue id (" << metal_res << ") for the metal "
            << metal_name << " does not correspond to a metal ion.\n"
            << "Assuming that this metal atom is part of a larger molecule\n";
        warn(A_fname, "lookup(str, res)", msg.str());
        return 0;
      }
      return m;    
    }

  std::ostringstream msg;
  if(msg_type == WARN_UKNOWN_METALS){
    msg << "Unknown PDB metal name: \"" << metal_name << "\""; 
    warn(A_fname, "lookup(str)", msg.str());
  }
  return 0;
}

const pdb_metal_info_t* 
PDB_metals::lookup(const atom_type metal, metal_warn_type msg_type)
{
  if(A_metal_res_names.size() == 0) set_metal_res_names();

  const pdb_metal_info_t* m = A_metal_array; 
  for(uint i = 0; i < A_metal_array_size; ++i, ++m) 
    if(m->metal_name == metal) return m;    

  std::ostringstream msg;
  if(msg_type == WARN_UKNOWN_METALS){
    msg << "Unknown PDB metal type: \" AtomType enum " << metal << "\""; 
    warn(A_fname, "lookup(metal_type)", msg.str());
  }
  return 0;
}

void
PDB_metals::set_metal_res_names()
{
  std::vector<std::string> res_names;
  res_names.push_back(" CA");
  A_metal_res_names[CALCIUM] = res_names;

  res_names.clear();
  res_names.push_back(" CO");
  res_names.push_back("3CO");
  A_metal_res_names[CO] = res_names;

  res_names.clear();
  res_names.push_back("CU1");
  res_names.push_back(" CU");
  res_names.push_back("CU3");
  A_metal_res_names[CU] = res_names;
  
  res_names.clear();
  res_names.push_back(" CD");
  A_metal_res_names[CADMIUM] = res_names;

  res_names.clear();
  res_names.push_back("FE2");
  res_names.push_back(" FE");
  A_metal_res_names[FE] = res_names;

  res_names.clear();
  res_names.push_back(" K ");
  A_metal_res_names[K] = res_names;

  res_names.clear();
  res_names.push_back(" MG");
  A_metal_res_names[MG] = res_names;

  res_names.clear();
  res_names.push_back(" MN");
  res_names.push_back("MN3");
  A_metal_res_names[MN] = res_names;

  res_names.clear();
  res_names.push_back(" NA");
  A_metal_res_names[NA] = res_names;

  res_names.clear();
  res_names.push_back(" NI");
  res_names.push_back("3NI");
  A_metal_res_names[NI] = res_names;

  res_names.clear();
  res_names.push_back(" ZN");
  A_metal_res_names[ZN] = res_names;
}
