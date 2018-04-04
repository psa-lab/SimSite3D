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
 * $Source: /psa/share/repository/pfizer_proj/src/gen_points/hphob_point_t.H,v $
 * $Revision: 1.2 $
 * $Author: vanvoor4 $
 * $Date: 2007-09-24 15:35:20 $
 * 
 * $Log: not supported by cvs2svn $
 * Revision 1.1  2007/08/21 18:33:17  vanvoor4
 * Initial checkin
 *
 *
 * 
 * 
 */ 

#ifndef HPHOB_POINT_HEADER_INCLUDED
#define HPHOB_POINT_HEADER_INCLUDED

#include <point.H>
#include <defs.H>
#include <map>

namespace SimSite3D{

  //! Use the base classes position as the computed center "ideal" template
  //! points.
  class hphob_point_t : public point_t{
  public:
    hphob_point_t(alloc_t a = ALLOC_POSITION) : point_t(a)
    {
    }

    hphob_point_t(const hphob_point_t& p) : point_t(p)
    {
       if(&p == this) return;

       for(std::list<atom_vci>::const_iterator ai = p.atoms.begin();
           ai != p.atoms.end(); ++ai)
         atoms.push_back(*ai);
    }

    ~hphob_point_t() { ; }

    std::list<atom_vci> atoms;  //!< Protein heavy atoms that could interact
  private:
  };

}

#endif
