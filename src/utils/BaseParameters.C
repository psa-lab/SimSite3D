/******************************************************************************
 * Copyright (c) 2006-2011, Michigan State University (MSU) Board of Trustees.
 * This file is part of the SimSite3D software project.
 *
 * Authors: Jeffrey Van Voorst, jeff.vanvoorst@gmail.com
 *          Leslie Kuhn, Ph.D., KuhnL@msu.edu 
 *
 *  SimSite3D is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  SimSite3D is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc.,
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *  You may also visit http://www.gnu.org/licenses/gpl-2.0.html.
 * 
 * $Date: 2011-02-25 09:54:04 -0500 (Fri, 25 Feb 2011) $: Date of last commit
 * $Author: vanvoor4 $: Author of last commit
 * $Rev: 1617 $: svn revision of last commit
 *
 * svn file: $Id: BaseParameters.C 1617 2011-02-25 14:54:04Z vanvoor4 $
 * file location: $URL: file:///psa/share/repository/SimSite3D/branches/surfaces-branch/src/utils/BaseParameters.C $
 *****************************************************************************/
#include <BaseParameters.H>
#include <string_basics.H>
#include <stream_basics.H>
#include <cstdlib>

using namespace SimSite3D;
typedef std::pair<std::string, std::string*> str_var_pair;
const std::string BaseParameters::A_fname = "BaseParameters.C";

BaseParameters::BaseParameters()
{
  A_status = INITIALIZING;
  init_str_to_var_map();
 
  // 0) This variable is needed till we decide how to handle the input data
  //    I recommend folding the static data into header files.
  if(!get_env_var("SIMSITE3D_INSTALL_DIR", &install_dir)){
    err_msg(A_fname, "BaseParameters()", "$SIMSITE3D_INSTALL_DIR is not set.  Please set $SIMSITE3D_INSTALL_DIR before\ncontinuing");
    A_status = FATAL_ERROR;
    return;
  }

  // 1) Attempt to load /etc/simsite3d/simsite3d.conf
  load_conf_file("/etc/simsite3d/simsite3d.conf");

  // 2) Attempt to load ${HOME}/.simsite3d/simsite3d.conf
  std::string home;
  if(get_env_var("HOME", &home)) 
    load_conf_file(home + "/.simsite3d/simsite3d.conf");
  else warn(A_fname, "cstr()", "The variable $HOME is not set.");

  // 3) Environment variables
  load_environment();

  // 4) values specified on cmd line -- must be set by derived classes

  // This may be introducing unintended path assignments
  if(proj_output.length() == 0 && scratch_dir.length() == 0)
    proj_output = scratch_dir = ".";
  else if(proj_output.length() == 0) proj_output = ".";
  else if(scratch_dir.length() == 0) scratch_dir = proj_output;

  // for testing
  load_surf_files = false;

  require_min_npts = true;
}

bool
BaseParameters::load_conf_file(const std::string conf_fname)
{
  std::ifstream conf_file;
  if(!normal_file_exists(conf_fname, false)) return false;
  if(!open_ifstream(conf_file, conf_fname)) return false;
  
  for(std::string line; std::getline(conf_file, line); ){
    if(line[0] == '%' || line[0] == '#') continue;

    strip_trailing_comments(&line, "#%");
    eat_leading_whitespace(&line);
    std::string::size_type pos = line.find_first_of(" \t");
    if(pos == std::string::npos) continue;

    std::string var = line.substr(0, pos);
    line.erase(0, pos);
    eat_leading_whitespace(&line);
    std::string val = line;
    pos = line.find_first_of(" \t");
    if(pos != std::string::npos) val.erase(pos);
       
    std::map<std::string, std::string*>::iterator m_iter;
    m_iter = A_str_to_var.find(var);
    if(m_iter == A_str_to_var.end()){
      std::string msg = "Unknown option in the parameters file: ";
      msg += conf_fname + "\n\t" + var + "\n"; 
      warn(A_fname, "load_conf_file", msg);
    }else *(m_iter->second) = val;
  }
  return true;
}

void
BaseParameters::print_version(std::ostream &out, std::string prog_name)
{
  out << "Copyright (C) 2006-2011, Michigan State University "
      << "(MSU) Board of Trustees.\n"
      << "\nWritten by Jeffrey R. Van Voorst and Leslie A. Kuhn\n"
      << "SimSite3D comes with ABSOLUTELY NO WARRANTY;\n"
      << prog_name 
      << " is part of SimSite3D, and is free software;\n"
      << "and you are welcome to redistribute it and/or "
      << "modifiy it under the terms of\nthe GNU General Public License as "
      << "published by the Free Software Foundation;\neither version 2 of the "
      << "License or (at your option) any later version.\n";
}

void
BaseParameters::load_environment()
{
  std::map<std::string, std::string*>::iterator m_iter;
  for(m_iter = A_str_to_var.begin(); m_iter != A_str_to_var.end(); ++m_iter)
    get_env_var(m_iter->first, m_iter->second); 
}

void
BaseParameters::init_str_to_var_map()
{
  A_str_to_var.insert(str_var_pair("SIMSITE3D_DBASE_LIGS", &dbase_ligs)); 
  A_str_to_var.insert(str_var_pair("SIMSITE3D_DBASE_PROTS", &dbase_prots)); 
  A_str_to_var.insert(str_var_pair("SIMSITE3D_DBASE_SITES", &dbase_sites)); 
  A_str_to_var.insert(str_var_pair("SIMSITE3D_DIVERSE_LIGS", &diverse_ligs)); 
  A_str_to_var.insert(str_var_pair("SIMSITE3D_DIVERSE_SITES", &diverse_sites)); 
  A_str_to_var.insert(str_var_pair("SIMSITE3D_SCRATCH_DIR", &scratch_dir));
  A_str_to_var.insert(str_var_pair("SIMSITE3D_PROJ_OUTPUT", &proj_output));
}

bool
BaseParameters::get_env_var(const std::string var, std::string *val)
{
  const char *tmp = std::getenv(var.c_str());
  if(!tmp) return false;
  *val = tmp;
  return true;
}
