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
 * $Source: /psa/share/repository/pfizer_proj/src/basics/param_tools.C,v $
 * $Revision: 1.3 $
 * $Author: vanvoor4 $
 * $Date: 2007-11-01 15:57:08 $
 * 
 * $Log: not supported by cvs2svn $
 * Revision 1.2  2007/08/21 15:43:25  vanvoor4
 * Changed warn call so that it used the default value for out
 *
 * Revision 1.1  2007/02/07 16:07:15  vanvoor4
 * Wrappers for popt and getenv functions.
 *
 * 
 */

#include <cstring>
#include <cstdlib>
#include <param_tools.H>
#include <string_basics.H>

static const std::string _fname = "param_tools.C";

bool
get_popt_arg(poptContext& optCon, char** arg)
{   
  const char* rv = poptGetArg(optCon);
  if(rv == NULL) return false; 
  size_t sz = strlen(rv) + 1;
  *arg = (char*) my_malloc(sz * sizeof(char));
  if(!*arg) return false;
  std::copy(rv, rv + sz, *arg);
  return true; 
} 

bool
check_file_name(char** fname, const std::string dir)
{
  if(!normal_file_exists(*fname, false)){
    if(dir.length() == 0) return false;
    std::string tmp = dir + "/" + *fname;
    if(!normal_file_exists(tmp)) return false;

    size_t sz = tmp.length() + 1;
    *fname = (char*) my_realloc(*fname, sz);
    if(!*fname) return false;
    std::strcpy(*fname, tmp.c_str());
  }
  return true;
}

bool
get_env_var(const std::string var, char** str)
{
  const char *tmp = getenv(var.c_str());
  if(tmp){
    uint sz = strlen(tmp) + 1;
    *str = (char*) my_realloc(*str, sz);
    if(*str == NULL) return false;
    std::copy(tmp, tmp + sz, *str);
    return true;
  }
  std::string msg = "Unable to read the environment variable ";
  warn(_fname, "get_env_var", msg + var);
  return false;
}
