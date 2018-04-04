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
 * $Source: /psa/share/repository/pfizer_proj/src/basics/string_basics.C,v $
 * $Revision: 1.11 $
 * $Author: vanvoor4 $
 * $Date: 2009-01-21 01:40:59 $
 * 
 * $Log: not supported by cvs2svn $
 * Revision 1.10  2007/11/07 14:59:04  vanvoor4
 * Moved some string fcns from stream_basics to here.
 *
 * Revision 1.9  2007/11/01 15:54:53  vanvoor4
 * Fixed a bug in the get path and struct id function.
 * Added an ostream parameter to bad_file_line so that we
 * can make global changes in the header file or local changes
 * to the chosen ostream for messages and errors
 *
 * Revision 1.8  2007/10/18 04:27:36  vanvoor4
 * The file extentions for the coordinate files need not have all
 * have the same number of characters
 *
 * Revision 1.7  2007/10/11 16:01:58  vanvoor4
 * Added function to write messages to a given ostream -- defaults to
 * std::err
 *
 * Revision 1.6  2007/08/21 16:39:39  vanvoor4
 * Actually removed them this time
 *
 * Revision 1.5  2007/08/21 16:38:49  vanvoor4
 * Moved the atom types checking to atom_t
 *
 * Revision 1.4  2007/07/05 16:39:30  vanvoor4
 * Added is_nitrogen() and is_oxygen()
 *
 * Revision 1.3  2007/06/06 19:52:40  vanvoor4
 * Changed string_tok to a more reentrant version
 *
 * Revision 1.2  2007/05/09 14:48:29  vanvoor4
 * Added string_tok()
 *
 * Revision 1.1  2007/02/07 16:01:43  vanvoor4
 * The string functions from utils
 *
 *
 */

#include <cstring>
#include <cstdlib>
#include <string_basics.H>

const std::string _fname = "string_basics.C";

void* 
my_malloc(size_t size)
{
  void *loc = std::malloc(size);
  if(loc) return loc;
  err_msg(_fname, "malloc", "memory allocation failed", std::cerr);
  return NULL;
}

void* 
my_realloc(void* loc, size_t size)
{
  loc = std::realloc(loc, size);
  if(loc) return loc;
  err_msg(_fname, "realloc", "memory allocation failed", std::cerr);
  return NULL;
}

bool
my_strtof(const std::string str, my_float_t* rv)
{
  char* endptr;
  const char* nptr = str.c_str();
  *rv = std::strtof(nptr, &endptr);
  if(endptr == nptr){
     warn(_fname, "my_strtof", "Non numeric value in coversion to float",
          std::cerr);
     return false;
  }
  return true;
}

bool
my_strtod(const std::string str, my_float_t* rv)
{
  char* endptr;
  const char* nptr = str.c_str();
  *rv = std::strtod(nptr, &endptr);
  if(endptr == nptr){
     warn(_fname, "my_strtod", "Non numeric value in coversion to float",
          std::cerr);
     return false;
  }
  return true;
}

bool 
my_strtoui(const std::string str, uint* rv)
{
  char* endptr;
  const char* nptr = str.c_str();
  *rv = static_cast<uint>(strtoul(nptr, &endptr, 10));
  if(endptr == nptr){
     warn(_fname, "my_strtoui", "Non numeric value in coversion to uint",
          std::cerr);
     return false;
  }
  return true;
}

bool
my_strtoi(const std::string str, int *rv)
{
  char* endptr;
  const char* nptr = str.c_str();
  *rv = static_cast<int>(strtol(nptr, &endptr, 10));
  if(endptr == nptr){
     warn(_fname, "my_strtoi", "Non numeric value in coversion to int",
          std::cerr);
     return false;
  }
  return true;
}

void 
err_msg(std::string file, std::string fcn, std::string msg, std::ostream& out)
{
  // use std::endl because we always want to write out errors
  out << "\n\t*****   ERROR   *****\n"
      << "File:  " << file << "\n"
      << "Function:  " << fcn << "\n"
      << msg << std::endl;
}

void 
warn(std::string file, std::string fcn, std::string msg, std::ostream& out)
{
  out << "\n\tWarning!     "
      << "File:  " << file << "   Function:  " << fcn << "\n"
      << msg << "\n";
}

void 
message(std::string msg, std::ostream& out)
{ out << msg << "\n"; }

void 
bad_file_line(std::string source, std::string function, std::string offender, 
              std::string line, std::ostream& out)
{
  std::string msg = "Bad line in ";
  msg += offender + std::string("\n\"") + line + "\"";
  warn(source, function, msg, out);
}

void
string_tok(std::string str, std::vector<std::string>* toks, char delim)
{
  std::string::const_iterator prev;
  std::string::const_iterator curr;
  for(curr = str.begin(); curr < str.end(); ++curr){
    for(prev = curr; curr < str.end() && *curr != delim; ++curr);
    std::string tmp(prev, curr);
    toks->push_back(tmp);
  }
}

bool
get_path_and_struct_id(const std::string fname, std::string* path,
                       std::string* struct_id)
{
  std::string tmp = fname;
  size_t pos = tmp.rfind('/');
  *struct_id = tmp.substr(pos + 1);

  // If a file extension was given, then we need to remove the _Y.XXX
  // portion of the file name to get the structure id -- where _Y is the 
  // the file characterization (_p: protein, _l: ligand, etc) and XXX is the
  // file extension.
  size_t dot_pos = struct_id->rfind('.');
  if(dot_pos != std::string::npos)
    struct_id->erase(struct_id->begin() + (dot_pos - 2), struct_id->end());
  if(pos != std::string::npos) *path = tmp.substr(0, pos+1);
  else *path = ".";
  return true;
}

int
strip_trailing_comments(std::string *s, std::string chars)
{
  size_t pos = s->find_first_of(chars); 
  if(pos == std::string::npos) return s->length();
  s->erase(pos);
  return pos;
}

void
eat_leading_whitespace(std::string *s)
{
  std::string::iterator end_pt;
  for(end_pt = s->begin(); end_pt < s->end() && isspace(*end_pt); end_pt++);
  s->erase(s->begin(), end_pt);
}
