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
 * $Source: /psa/share/repository/pfizer_proj/src/gen_points/interact_point_t.H,v $
 * $Revision: 1.2 $
 * $Author: vanvoor4 $
 * $Date: 2008/04/17 18:50:58 $
 * 
 * $Log: interact_point_t.H,v $
 * Revision 1.2  2008/04/17 18:50:58  vanvoor4
 * Added some typedefs to make the code more human readable.
 *
 * Revision 1.1  2007/08/21 18:33:38  vanvoor4
 * Initial checkin
 *
 *
 * 
 * 
 */ 

#ifndef _INTERACT_POINT_HEADER_INCLUDED
#define _INTERACT_POINT_HEADER_INCLUDED

#include <dir_point.H>
#include <defs.H>
#include <map>

namespace SimSite3D{

  //! Interactions points
  class interact_point_t : public dir_point_t{
  public:
    interact_point_t(alloc_t a = ALLOC_POSITION) : dir_point_t(a)
    {
      act_type = UNKNOWN_INTERACTION;
      // problem with iterators is cannot initialize atom to anything 
      // safe/usefull
    }

    interact_point_t(const interact_point_t& p) : dir_point_t(p)
    {
      do_copy(p); 
    }

    const interact_point_t& 
    operator=(const interact_point_t& p)
    {   
      if(&p != this){
        do_copy(p);
        dir_point_t::operator= (p);
      }
      return *this;
    }

    ~interact_point_t() { ; }

    static bool
    cmp(const interact_point_t &A, const interact_point_t &B)
    { return A.act_type < B.act_type; }

    interactionType act_type; //!< Interaction type of the points
  private:
    void
    do_copy(const interact_point_t& p)
    {
      act_type = p.act_type;
    }
  };

  typedef struct{
    chain_const_iter chain;
    residue_vci residue;
    atom_vci atom;
  }interact_atom_t;
  typedef std::vector<interact_atom_t> interact_atoms_vec;
  typedef interact_atoms_vec::const_iterator interact_atoms_vci;

}

#endif
