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
 * $Source: /psa/share/repository/pfizer_proj/src/gen_points/fit_point_t.H,v $
 * $Revision: 1.1 $
 * $Author: vanvoor4 $
 * $Date: 2007-08-21 18:09:41 $
 *
 * $Log: not supported by cvs2svn $
 *
 *
 *
 */
                                                                                
#ifndef _HBOND_POINT_HEADER_INCLUDED
#define _HBOND_POINT_HEADER_INCLUDED
                                                                                
#include <point.H>
                                                                                
namespace SimSite3D{

  //! Used to "harbor" the points used to align sitemaps.
  //! For now the direction is handled in a separate "array"
  class fit_point_t : public point_t{
  public:
    //! Basic cstr -- calls point_t::cstr to handle the position
    fit_point_t(alloc_t a = ALLOC_POSITION) : point_t(a)
    {
      act_type = UNKNOWN_INTERACTION; 
    }  

    //! Basic copy cstr -- required so that we can call the copy cstr for the 
    //! point_t class.
    fit_point_t(const fit_point_t& p) : point_t(p)
    {
      do_copy(p);
    }

    const fit_point_t&
    operator=(const fit_point_t& p) 
    {
      if(&p != this){
        do_copy(p);
        point_t::operator= (p);
      }
      return *this;
    }
  private:
    void
    do_copy(const fit_point_t& p)
    {
      act_type = p.act_type;
      atoms = p.atoms;
    }

    interactionType act_type;    //!< Acceptor, Donor or Hydrophobic
    std::vector<atom_vci> atoms; //!< Atom or atoms making the interaction at this point
  };

}
