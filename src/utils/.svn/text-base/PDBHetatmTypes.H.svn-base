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
 * $Source: /psa/share/repository/pfizer_proj/src/basics/PDBHetatmTypes.H,v $
 * $Revision: 1.2 $
 * $Author: vanvoor4 $
 * $Date: 2008-07-29 16:56:00 $
 * 
 * $Log: not supported by cvs2svn $
 * Revision 1.1  2008/05/13 17:13:54  vanvoor4
 * Initial checkin --
 * vdw radii are from table XI in Li & Nussinov
 *
 * 
 */

#ifndef PDB_HETATM_TYPES_HEADER_INCLUDED
#define PDB_HETATM_TYPES_HEADER_INCLUDED

#include <cctype>
#include <AtomTypes.H>
#include <iostream>

namespace SimSite3D{

  class PDBHetatmTypes : public AtomTypes{
  public:
    //! A somewhat arbtrary method of figuring out element type from PDB 
    //! Hetatm 4 character atom name string
    static const atom_type
    string_to_atom(const std::string atom_str)
    {
      if(A_string_to_atom_type.empty()) build_the_map();

      std::string new_str = atom_str;
      // old skool hydrogens, in particular, have an leading number
      if(isdigit(new_str[0])) new_str[0] = ' ';
      // remediated hydrogens are typically H[A-Z][0-9]*
      if(new_str[0] == 'H' && isalpha(new_str[1]) && isdigit(new_str[2]))
        new_str = " H  ";
      // assume that we don't know how to handle atom names such as " O1P"
      new_str[2] = new_str[3] = ' ';
      std::map<std::string, atom_type>::const_iterator my_iter;
      my_iter = A_string_to_atom_type.find(new_str);
      if(my_iter != A_string_to_atom_type.end()) return my_iter->second;
      return UNKNOWN_ATOM;
    }

  private:
    static void build_the_map();   

    static std::map<std::string, atom_type> A_string_to_atom_type;

    static const uint A_num_atom_strings;
    static const atom_conv_type A_atom_strings[];
  };
}
#endif
