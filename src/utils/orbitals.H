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
 * $Source: /psa/share/repository/pfizer_proj/src/basics/orbitals.H,v $
 * $Revision: 1.1 $
 * $Author: vanvoor4 $
 * $Date: 2007-10-11 16:01:26 $
 * 
 * $Log: not supported by cvs2svn $
 *
 *
 */

#ifndef ORBITALS_HEADER_INCLUDED
#define ORBITALS_HEADER_INCLUDED

#include <string>
#include <map>

namespace SimSite3D{

  //! atom orbital types 
  typedef enum orbit_enum{
    UNKNOWN_ORBIT  =  0,
    SP1,       //!< sp1 orbital
    SP2,       //!< sp2 orbital
    SP3,       //!< sp3 orbital
    SP4,       //!< sp4 orbital
    AMIDE,     //!< Amide nitrogen
    AMBIG,     //!< ambiguous -- typically unresolved sidechain O and N
    CO2,
    CAT,
    PL3,
    AR,         //!< aromatic
    O_ORBITAL,
    O2,
    DEFAULT_ORBIT,   //!< No orbit was specified
    ORBIT_ENUM_END   //!< End of enum marker
  }orbit_type;

  //! This class was added in a partial braindead moment -- is currently kept
  //! in the event of it being useful at a later date.  In particular, if we
  //! wish to convert orbit_type ==> string.
  class orbit_defs{
  public:
    //! Setup the maps
    orbit_defs();

    ~orbit_defs() { ; }

    orbit_type 
    operator[](std::string s_in)
    { 
      std::map<std::string, orbit_type>::const_iterator pos;
      pos = str_to_orbital.find(s_in);
      return (pos == str_to_orbital.end() ? UNKNOWN_ORBIT : pos->second);
    }

    std::string 
    operator[](orbit_type orb)
    { 
      std::map<orbit_type, std::string>::const_iterator pos;
      pos = orbital_to_str.find(orb);
      return (pos == orbital_to_str.end() ? ".err" : pos->second);
    }

  private:
    std::map<std::string, orbit_type> str_to_orbital;  
    std::map<orbit_type, std::string> orbital_to_str;
  };

}
#endif
