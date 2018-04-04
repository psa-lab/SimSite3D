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
 * $Source: /psa/share/repository/pfizer_proj/src/gen_points/SitemapPointsFile.C,v $
 * $Revision: 1.6 $
 * $Author: vanvoor4 $
 * $Date: 2009-01-12 21:12:10 $
 * 
 * $Log: not supported by cvs2svn $
 * Revision 1.5  2008/07/30 16:48:25  vanvoor4
 * Occupancy is now at 1.00 for sitemap points rather than taking
 * the default occupancy (currently 0.00 is the default).
 *
 * Revision 1.4  2008/03/31 17:55:31  vanvoor4
 * Removed the explicit call to base constructor.
 *
 * Revision 1.3  2007/12/17 21:23:44  vanvoor4
 * PDB sitemap points need to be labeled as hetero to be written
 * correctly
 *
 * Revision 1.2  2007/08/29 20:23:23  vanvoor4
 * Changed handling of hbond point types
 *
 * Revision 1.1  2007/08/21 18:43:07  vanvoor4
 * initial checkin
 *
 *  
 *  
 *  
 */ 
    
#include <SitemapPointsFile.H>
#include <iomanip>

using namespace ASCbase;
const std::string SitemapPointsFile::_fname = "SitemapPointsFile.C";

SitemapPointsFile::SitemapPointsFile()
{
  counts.resize(4);
  std::fill(counts.begin(), counts.end(), 0);
}

SitemapPointsFile::SitemapPointsFile(const std::string filename, 
                                     const verbose_level_t verbosity)
 : PDBBase(filename, verbosity)
{
  counts.resize(4);
  std::fill(counts.begin(), counts.end(), 0);

  if(fail()) return;
  std::ifstream pdb_file;
  if(open_ifstream(pdb_file, filename)){
//    read_header(pdb_file);
    if(!read_data(pdb_file)) red_light();
    pdb_file.close();

    for(atom_vci a = atoms_begin(); a != atoms_end(); ++a){
      if(a->tempFactor == 0.0) ++(counts[0]);
      else if(a->tempFactor == 50.0) ++(counts[1]);
      else if(a->tempFactor == 25.0) ++(counts[2]);
      else if(a->tempFactor == 100.0) ++(counts[3]);
    }

  }else red_light();
}


SitemapPointsFile::~SitemapPointsFile()
{
}

#if 0
bool 
SitemapPointsFile::read_header(std::ifstream &fin)
{
  fin.seekg(0);
  bool rv = false;
  std::string line;
  std::string method_str;
  method_str = "REMARK  10 Sitemap Generation method: ball (point and radius)";
  while(std::getline(fin, line) && line.substr(0,6) == "REMARK")
    if(line.substr(0,62) == method_str && std::getline(fin, line)){
      std::string s = line.substr(line.find("Params:" + 7))
      eat_leading_whitespace(&s);
      std::vector<std::string> toks;
      string_tok(s, &toks, ' ');
      for(size_t i = 0; i < params.size(); ++i)

      rv = true;
      break;
    }

  return rv; 
}
#endif

bool 
SitemapPointsFile::add_points(const hbond_fit_pt_vci _begin, 
                              const hbond_fit_pt_vci _end)
{
  uint serial_adds[] = {1000, 2000, 3000};
  for(hbond_fit_pt_vci pt = _begin; pt < _end; ){
    uint* serial = 0;
    uint serial_add = 0;
    my_float_t tempFactor = 0.0;
    interactionType prev_act_type = pt->act_type;
    if(pt->act_type == ACCEPTOR){ 
      tempFactor = 0.0;
      serial = &(counts[0]);
      serial_add = serial_adds[0];
    }else if(pt->act_type == DONOR){ 
      tempFactor = 50.0;
      serial = &(counts[1]);
      serial_add = serial_adds[1];
    }else if(pt->act_type == DONEPTOR){ 
      tempFactor = 25.0;
      serial = &(counts[2]);
      serial_add = serial_adds[2];
    }else{
      warn(_fname, "add_points()", "Unknown hydrogen bond point type");
      return false;
    }

    atom_t atom;
    for( ; pt < _end && pt->act_type == prev_act_type; ++pt){
      (*serial)++;
      atom.atom_num = *serial + serial_add;
      atom.res_num = *serial + serial_add;
      atom.name = O;
      atom.res = HOH;
      atom.name_str = " O  ";
      atom.res_str = "HOH";
      atom.occupancy = 1.0;
      atom.tempFactor = tempFactor;
      atom.is_hetero = true;
      std::copy(pt->pos, pt->pos + 3, atom.pos);
      append_atom(atom);
    }
  }
  return true;
}

bool
SitemapPointsFile::add_points(const hphob_point_lci _begin,
                              const hphob_point_lci _end)
{
  atom_t atom;
  for(hphob_point_lci pt = _begin; pt != _end; ++pt){
    std::copy(pt->pos, pt->pos + 3, atom.pos);
    ++(counts[3]);
    atom.atom_num = counts[3];
    atom.res_num = counts[3];
    atom.occupancy = 1.0;
    atom.tempFactor = 100.0;
    atom.name = O;
    atom.res = HOH;
    atom.name_str = " O  ";
    atom.res_str = "HOH";
    atom.is_hetero = true;
    append_atom(atom);
  } 
  return true;
}
