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
 * $Source: /psa/share/repository/pfizer_proj/src/basics/string_basics.H,v $
 * $Revision: 1.8 $
 * $Author: vanvoor4 $
 * $Date: 2007-11-07 14:59:00 $
 * 
 * $Log: not supported by cvs2svn $
 * Revision 1.7  2007/11/01 15:54:58  vanvoor4
 * Fixed a bug in the get path and struct id function.
 * Added an ostream parameter to bad_file_line so that we
 * can make global changes in the header file or local changes
 * to the chosen ostream for messages and errors
 *
 * Revision 1.6  2007/10/11 16:02:21  vanvoor4
 * Added a function to send messages to an ostream -- defaults to stderr
 *
 * Revision 1.5  2007/08/21 16:41:22  vanvoor4
 * Removed the atom checking routines and added the struct id function
 *
 * Revision 1.4  2007/07/05 16:39:37  vanvoor4
 * Added is_nitrogen() and is_oxygen()
 *
 * Revision 1.3  2007/06/06 19:52:18  vanvoor4
 * Changed string_tok to a more reentrant version and added a comment.
 *
 * Revision 1.2  2007/05/09 14:48:22  vanvoor4
 * Added string_tok()
 *
 * Revision 1.1  2007/02/07 16:01:34  vanvoor4
 * The string functions from utils
 *
 *
 * 
 * 
 */

#ifndef _STRING_BASICS_HEADER_INCLUDED
#define _STRING_BASICS_HEADER_INCLUDED

#include <string>
#include <iostream>
#include <cstdlib>
#include <types.H>
#include <defs.H> // go away in a bit

  //! Wrapper for malloc -- used mainly for popt
  void* my_malloc(size_t size);

  //! Wrapper for realloc -- used mainly for popt
  void* my_realloc(void* loc, size_t size);

  bool my_strtof(const std::string str, my_float_t* rv);

  bool
  my_strtod(const std::string str, my_float_t* rv);

  bool my_strtoui(const std::string str, uint* rv);

  bool my_strtoi(const std::string str, int *rv);

  //! Wrapper for error messages -- typically fatal errors.  out should be
  //! left to the default value so that the default value need be changed only
  //! here.
  void err_msg(std::string file, std::string fcn, std::string msg,
               std::ostream& out = std::cerr);
  
  //! Wrapper for warnings and nonfatal error messages.  out should be
  //! left to the default value so that the default value need be changed only
  //! here.
  void warn(std::string file, std::string fcn, std::string msg,
            std::ostream& out = std::cerr);

  //! Send an error/warning message with out printing file info header.
  /*! The intended use is for calling methods to send additional information
   */
  void message(std::string msg, std::ostream& out = std::cerr);

  void bad_file_line(std::string source, std::string function,
                     std::string filename, std::string line,
                     std::ostream& out = std::cerr);

  //! A quick and dirty std::string tokenizer
  /*! NOTE:  Unlike strtok, this method does not modify the input string nor 
   * does it treat consecutive delimiters as one.  The reason is, when parsing
   * .csv files, pipe delimited files, etc. we may wish to denote an empty 
   * field or missing value by leaving no space between the delimiters.
   * 
   * @param str string to tokenize
   * @param toks Pointer to vector of string tokens
   * @param delim The delimiter character
   */
  void string_tok(std::string str, std::vector<std::string>* toks, char delim);

  //! Given a path to a structure, ligand, or sitemap file get the leading 
  //! path (directory structure) and the structure id
  bool get_path_and_struct_id(const std::string fname, std::string* path,
                              std::string* struct_id);

  //! Strip comments from string
  /*! NOTE: modifies s
   *
   * @param s Pointer to the string to strip
   * @param chars Characters which indicate the rest of the line is a comment
   * @return Remaining length of the string (0 if first char of string is a 
   *         comment char)
   */
  int strip_trailing_comments(std::string *s, std::string chars);

  //! Remove the whitespace from the string
  /*! NOTE: modifies s
   *
   */
  void eat_leading_whitespace(std::string *s);

#endif
