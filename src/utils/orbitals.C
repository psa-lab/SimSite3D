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
 * $Source: /psa/share/repository/pfizer_proj/src/basics/orbitals.C,v $
 * $Revision: 1.1 $
 * $Author: vanvoor4 $
 * $Date: 2007-10-11 16:01:15 $
 * 
 * $Log: not supported by cvs2svn $
 *
 *
 */


#include <orbitals.H>

using namespace SimSite3D;

orbit_defs::orbit_defs()
{
  //std::map<std::string, orbit_type> str_to_orbital;
  str_to_orbital[".err"] = UNKNOWN_ORBIT;
  str_to_orbital[".1"] = SP1;
  str_to_orbital[".2"] = SP2;
  str_to_orbital[".3"] = SP3;
  str_to_orbital[".4"] = SP4;
  str_to_orbital[".am"] = AMIDE;
  str_to_orbital[".ambig"] = AMBIG;
  str_to_orbital[".co2"] = CO2;
  str_to_orbital[".cat"] = CAT;
  str_to_orbital[".pl3"] = PL3;
  str_to_orbital[".ar"] = AR;
  str_to_orbital[".O"] = O_ORBITAL;
  str_to_orbital[".O2"] = O2;
  str_to_orbital[""] = DEFAULT_ORBIT;
  
  //static std::map<orbit_type, std::string> orbital_to_str;
  orbital_to_str[UNKNOWN_ORBIT] = ".err";
  orbital_to_str[SP1] = ".1";
  orbital_to_str[SP2] = ".2";
  orbital_to_str[SP3] = ".3";
  orbital_to_str[SP4] = ".4";
  orbital_to_str[AMIDE] = ".am";
  orbital_to_str[AMBIG] = ".ambig";
  orbital_to_str[CO2] = ".co2";
  orbital_to_str[CAT] = ".cat";
  orbital_to_str[PL3] = ".pl3";
  orbital_to_str[AR] = ".ar";
  orbital_to_str[O_ORBITAL] = ".O";
  orbital_to_str[O2] = ".O2";
  orbital_to_str[DEFAULT_ORBIT] = "";
  orbital_to_str[ORBIT_ENUM_END] = ".err";
}
