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
 * $Source: /psa/share/repository/pfizer_proj/src/basics/param_tools.H,v $
 * $Revision: 1.3 $
 * $Author: vanvoor4 $
 * $Date: 2007-11-01 15:57:21 $
 * 
 * $Log: not supported by cvs2svn $
 * Revision 1.2  2007/08/21 15:43:46  vanvoor4
 * Changed license header
 *
 * Revision 1.1  2007/02/07 16:07:21  vanvoor4
 * Wrappers for popt and getenv functions.
 *
 *
 * 
 * 
 */

#ifndef _PARAM_TOOLS_HEADER_INCLUDED
#define _PARAM_TOOLS_HEADER_INCLUDED

#include <stream_basics.H>
#include <popt.h>

//! Wrapper for the poptGetArg function
/*!
 * @param optCon Reference to the current poptContext
 * @param arg Pointer to the cstring to get the argument
 * @return false if poptGetArg fails or malloc fails, else true
 */
bool get_popt_arg(poptContext& optCon, char** arg);

//! Check if the filename is for a normal file or if it is in the dbase dir
/*! First check if the filename (as given by the user) is the name of a 
 * normal file.  If it is, return true.  Otherwise, check if the filename
 * (as given) exists in the given directory (as a normal file).
 *
 * @param fname Pointer to a cstring holding the name of a file
 * @param dir String holding the directory to check for existence of the file
 * @return True if the sitemap file is a normal file or is in the dbase dir
 */
bool check_file_name(char** fname, const std::string dir);

//! getenv wrapper for strings
/*!
 * @param var Name of the environment variable
 * @param str Pointer to the cstring to hold the environment var
 * @return True if variable exists and has a value, otherwise false
 */
//bool get_env_var(const std::string var, char** str);

#endif
