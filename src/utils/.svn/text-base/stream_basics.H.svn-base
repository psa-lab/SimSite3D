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
 * $Source: /psa/share/repository/pfizer_proj/src/basics/stream_basics.H,v $
 * $Revision: 1.6 $
 * $Author: vanvoor4 $
 * $Date: 2009-01-12 21:04:43 $
 * 
 * $Log: not supported by cvs2svn $
 * Revision 1.5  2007/11/07 14:59:23  vanvoor4
 * Moved some string fcns from here to stream_basics.
 *
 * Revision 1.4  2007/08/21 15:58:39  vanvoor4
 * Added eat_leading_whitespace
 *
 * Revision 1.3  2007/06/06 19:53:30  vanvoor4
 * Added strip_trailing_comments -- modifies the given string.
 *
 * Revision 1.2  2007/05/09 14:44:56  vanvoor4
 * Added dir_exists() and my_mkdir()
 *
 * Revision 1.1  2007/02/07 16:01:06  vanvoor4
 * The stream functions from utils
 *
 *
 * 
 * 
 */

#ifndef _STREAM_BASICS_HEADER_INCLUDED
#define _STREAM_BASICS_HEADER_INCLUDED

#include <fstream>
#include <sys/types.h>

  //! Check if the file name exists and is a normal file
  /*! Check if the file whose path is stored in name exists and is a normal 
   * file.  The second parameter allows the suppresion of the warning messages.
   * This is typically useful when checking for existence in multiple 
   * directories/paths
   *
   * @param name Name of the file to check for existence
   * @param warn False implies the suppression of warning messages
   * @return True if file exists and is normal, otherwise false.
   */
  bool normal_file_exists(std::string name, bool warn = true);

  //! Check if the file name exists and is a directory
  /*! Check if the file whose path is stored in name exists and is a directory
   * The second parameter allows the suppresion of the warning messages.
   * This is typically useful when checking for existence in multiple 
   * directories/paths
   *
   * @param name Name of the directory to check for existence
   * @param warn False implies the suppression of warning messages
   * @return True if file exists and is a directory, otherwise false.
   */
  bool dir_exists(std::string name, bool warn = true);

  //! Wrapper for making directories -- paritially to handle NFS 
  /*! Wrapper for the C mkdir function.  Under NFS (and possibly other 
   * distributed filesystems) the mkdir call does not return a 0 upon 
   * sucessfully creating a directory as specified in the manpages for mkdir.
   * The problem is mkdir checks too quickly.  my_mkdir sleeps for 1 second
   * after calling mkdir and then checks using the dir_exists() function.
   *
   * @param dir Name of the directory to create; all directories in "dir" 
   *            except for the last must exist.
   * @param mode Permissions on the newly created directory (defaults to 0700).
   * @return value of dir_exits(dir)
   */
  bool my_mkdir(std::string dir, mode_t mode = 0700);
  
  //! Opens the input file stream with the given name.
  bool open_ifstream(std::ifstream& infile, std::string name);
  
  //! Opens the output file stream with the given name.
  bool open_ofstream(std::ofstream& outfile, std::string name);
  
  //! Counts the number of '\n' characters in a file
  int count_lines(std::ifstream& In);

  bool find_section(std::ifstream &in, std::string tag);

#endif
