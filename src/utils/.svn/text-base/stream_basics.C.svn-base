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
 * Authors: Jeffrey Van Voorst, vanvoor4@msu.edu
 *          Leslie Kuhn, Ph.D., KuhnL@msu.edu 
 *****************************************************************************/

/*
 * $Source: /psa/share/repository/pfizer_proj/src/basics/stream_basics.C,v $
 * $Revision: 1.8 $
 * $Author: vanvoor4 $
 * $Date: 2009-01-12 21:04:43 $
 * 
 * $Log: not supported by cvs2svn $
 * Revision 1.7  2007/11/14 16:15:00  vanvoor4
 * Changed the output for checking existence of files and for opening
 * streams for reading and writing.
 *
 * Revision 1.6  2007/11/07 14:59:28  vanvoor4
 * Moved some string fcns from here to stream_basics.
 *
 * Revision 1.5  2007/11/01 15:58:10  vanvoor4
 * Removed a commented out line -- wasn't needed
 *
 * Revision 1.4  2007/08/21 15:58:17  vanvoor4
 * Added eat leading whitespace and changed the assignment of a substr
 * to using erase from pos to end.
 *
 * Revision 1.3  2007/06/06 19:54:01  vanvoor4
 * Added strip_trailing_comments -- modifies the given string.
 *
 * Revision 1.2  2007/05/09 14:46:38  vanvoor4
 * Added dir_exists() and my_mkdir().  Changed the dir_exists()
 * and normal_file_exists() to be more restrictive.  Unfortunately
 * such checks might fail on links, and I am not sure if we would
 * like to use pipes or sockets in the future.
 *
 * Revision 1.1  2007/02/07 16:00:58  vanvoor4
 * The stream functions from utils.
 *
 *
 */

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include <iostream>
#include <sstream>
#include <cstring>
#include <algorithm>
#include <stream_basics.H>
#include <string_basics.H>

bool 
normal_file_exists(std::string name, bool warn)
{
  struct stat buf;
  // stat returns 0 on success
  if(stat(name.c_str(), &buf)){ 
    int errsv = errno;
    if(warn){
      std::ostringstream ostr;
      if(errsv == ENOENT) ostr << "Unable to find the file: " << name << "\n";
      else
        ostr << "Unable to access the file: " << name << ":\n" 
             << "\t" << std::strerror(errsv) << "\n";
      err_msg("stream_basics.C", "normal_file_exists()", ostr.str());  
    }
    return false;
  }

  // It is quite possible that this check needs to be modified based on
  // the types of files to be used.  For example, we may wish to use a pipe
  if(!S_ISREG(buf.st_mode) && warn){
    std::ostringstream ostr;
    ostr << "The file " << name << " is not a regular file, expected"
         << " a file name.\n";
    err_msg("stream_basics.C", "normal_file_exists()", ostr.str());  
    return false;
  }

  return true;
}

bool 
dir_exists(std::string name, bool warn)
{
  struct stat buf;
  // stat returns 0 on success
  if(stat(name.c_str(), &buf)){ 
    int errsv = errno;
    if(warn && errsv == ENOENT)
      std::cerr << "Unable to find the directory " << name << "\n" ;
    else if(warn)
      std::cerr << "Unable to open the directory " << name << ":\n" 
                << "\t" << std::strerror(errsv) << "\n";
    return false;
  }

  // It is quite possible that this check needs to be modified based on
  // the types of files to be used.
  if(!S_ISDIR(buf.st_mode) && warn){
    std::cerr << "The file " << name << " is a not directory, expected"
              << " a directory.\n";
    return false;
  }

  return true;
}

bool
my_mkdir(std::string dir, mode_t mode)
{
   // Checking the return of mkdir fails under NFS since mkdir expects the
   // creation of the directory to be immediate 
   mkdir(dir.c_str(), mode);
   sleep(1);
   return dir_exists(dir);
}

bool 
open_ifstream(std::ifstream& infile, std::string name)
{
  infile.open(name.c_str());
  if(infile.fail()){
    normal_file_exists(name);
    return false;
  }
  return true;
}

bool 
open_ofstream(std::ofstream& outfile, std::string name)
{
  outfile.open(name.c_str(), std::ios_base::out|std::ios_base::trunc);
  if(outfile.fail()){
    if(normal_file_exists(name)){
      std::ostringstream ostr;
      ostr << "Unable to open the file " << name << " for writing.\n";
      err_msg("stream_basics.C", "open_ofstream()", ostr.str());  
    }  
    return false;
  }
  return true;
}

// Submitted by Chris Theis 2005-12-14, 7:02 pm on www.codecomments.com
int 
count_lines(std::ifstream& In)
{
  int rv = static_cast<int>(std::count(std::istreambuf_iterator<char>(In),
				       std::istreambuf_iterator<char>(),'\n'));
  In.seekg(0);
  return rv;
}

bool
find_section(std::ifstream &in, std::string tag)
{
  for(std::string line; std::getline(in, line); )
    if(line.find(tag) != std::string::npos) return true;
  return false;
}
