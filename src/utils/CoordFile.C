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
 * $Source: /psa/share/repository/pfizer_proj/src/basics/CoordFile.C,v $
 * $Revision: 1.5 $
 * $Author: vanvoor4 $
 * $Date: 2009-01-12 20:42:35 $
 * 
 * $Log: not supported by cvs2svn $
 * Revision 1.4  2008/07/28 15:25:00  vanvoor4
 * Some variable name changes
 *
 * Revision 1.3  2008/02/26 17:52:39  vanvoor4
 * Added the constant A_fname
 *
 * Revision 1.2  2007/09/20 17:43:23  vanvoor4
 * Changed bump parameters to mirror those used by other PSA software 
 *
 * Revision 1.1  2007/08/29 20:20:45  vanvoor4
 * Initial checkin
 *
 * 
 */

#include <CoordFile.H>
#include <iostream>
#include <iomanip>

using namespace SimSite3D;

const my_float_t CoordFile::ALLOWED_OVERLAP = 0.5;
const my_float_t CoordFile::min_hbond_dist = 2.5;
const std::string CoordFile::A_fname = "CoordFile.H(C)"; 

CoordFile::CoordFile() 
{ 
  A_fail = false; 
}

CoordFile::CoordFile(const std::string filename_in, 
                     const verbose_level_t verbosity)
{ 
  name_of_file = filename_in; 
  if(verbosity != VERBOSE_SILENT)
    std::cout << "Reading the file: " << filename_in << "\n";

  A_fail = true; 
  if(filename_in.length() == 0){
    err_msg(A_fname, "Cstr()", "filename has zero length");
  }
  else A_fail = false;
}

CoordFile::CoordFile(const CoordFile& src)
{
  if(&src == this) return;
  name_of_file = src.name_of_file;
}

CoordFile::~CoordFile() 
{
}

void
CoordFile::write_xyzr_line(std::ostream& out, const my_float_t* pos,
                           const my_float_t vdw_radius)
{
  out.precision(3);
  out << std::setw(8) << pos[0] << " " << std::setw(8) << pos[1] << " "
      << std::setw(8) << pos[2] << " ";
  out.precision(2);
  out << vdw_radius << "\n";
}
