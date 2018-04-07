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
 * $Source: /psa/share/repository/pfizer_proj/src/basics/Quaternion.C,v $
 * $Revision: 1.1 $
 * $Author: vanvoor4 $
 * $Date: 2007-10-04 18:07:09 $
 * 
 * $Log: not supported by cvs2svn $
 *
 */

#include <Quaternion.H>
#include <sstream>

std::string
Quaternion::str()
{
  std::ostringstream ostr;
  ostr.setf(std::ios_base::fixed, std::ios_base::floatfield);
  ostr.precision(4);
  ostr << q[0] << " + i" << q[1] << " + j" << q[2] << " + k" << q[3];
  return ostr.str();
}

