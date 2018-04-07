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
 * $Source: /psa/share/repository/pfizer_proj/src/gen_points/hbond_points.H,v $
 * $Revision: 1.4 $
 * $Author: vanvoor4 $
 * $Date: 2008-05-15 17:30:47 $
 * 
 * $Log: not supported by cvs2svn $
 * Revision 1.3  2008/03/31 17:50:40  vanvoor4
 * added an iterator to point to the carbon neighbor of the
 * protein atom making the hbond -- useful for acceptors.
 *
 * Revision 1.2  2007/12/17 21:23:00  vanvoor4
 * Added atom for hbond_points
 *
 * Revision 1.1  2007/08/29 20:22:38  vanvoor4
 * Initial checkin
 *
 * 
 */ 

#ifndef HBOND_POINT_HEADER_INCLUDED
#define HBOND_POINT_HEADER_INCLUDED

#include <map>
#include <dir_point.H>
#include <dir_point_storage.H>
#include <defs.H>
#include <atom.H>
#include <BoundingVolume.H>


namespace SimSite3D{

  //! Forward declaration of "ideal" point class
  class hbond_ideal_point_t;
  typedef dir_point_storage<hbond_ideal_point_t>	hbond_ideal_pt_vec;
  typedef hbond_ideal_pt_vec::iterator			hbond_ideal_pt_vi;
  typedef hbond_ideal_pt_vec::const_iterator		hbond_ideal_pt_vci;

  //! A "fit" hbonding point
  class hbond_point_t : public dir_point_t{
  public:
    hbond_point_t(alloc_t a = ALLOC_POSITION) : dir_point_t(a)
    {
      act_type = UNKNOWN_INTERACTION;
      atom = atom_t::NULL_ATOM_VCI;
    }

    hbond_point_t(const hbond_point_t& p) : dir_point_t(p)
    {
      do_copy(p); 
    }

    const hbond_point_t& 
    operator=(const hbond_point_t& p)
    {
      if(&p != this){
        do_copy(p);
        dir_point_t::operator= (p);
      }
      return *this;
    }
    ~hbond_point_t() { ; }

    static bool
    cmp(const hbond_point_t &A, const hbond_point_t &B)
    { return A.act_type < B.act_type; }

    const SimSite3D::atom_t&
    __get_atom() const
    { return *atom; }

    const SimSite3D::hbond_ideal_point_t&
    __get_ideal_pt() const
    { return *ideal_pt; }

    interactionType act_type; //!< Interaction type of the points
    hbond_ideal_pt_vci ideal_pt;  //!< Iterator to the associated ideal point
    atom_vci atom;            //!< Heavy atom giving rise to the points 
  private:
    void
    do_copy(const hbond_point_t& p)
    {
      act_type = p.act_type;
      ideal_pt = p.ideal_pt;
      atom = p.atom;
    }
  };

  typedef dir_point_storage<hbond_point_t>	hbond_fit_pt_vec;
  typedef hbond_fit_pt_vec::iterator		hbond_fit_pt_vi;
  typedef hbond_fit_pt_vec::const_iterator	hbond_fit_pt_vci;

  //! An "ideal" hbond point base class
  class hbond_ideal_pt_base : public dir_point_t{
  public:
    hbond_ideal_pt_base(alloc_t a = ALLOC_POSITION) : dir_point_t(a)
    {
      pt_num = 0;
      act_type = UNKNOWN_INTERACTION;
      atom = atom_t::NULL_ATOM_VCI;
      carbon_nbr = atom_t::NULL_ATOM_VCI;
      second_nbr = atom_t::NULL_ATOM_VCI;
    }

    hbond_ideal_pt_base(const atom_vci hbond_atom, const atom_vci C_nbr_atom, 
                        const atom_vci second_nbr_atom, 
                        const BoundingVolume& site_vol, const int cap_number, 
                        const bool include_metals = false, 
                        const alloc_t a = ALLOC_POSITION) : dir_point_t(a)
    {
      pt_num = cap_number;
      //act_type = UNKNOWN_INTERACTION;
      atom = hbond_atom;
      carbon_nbr = C_nbr_atom;
      second_nbr = second_nbr_atom;
    }

    hbond_ideal_pt_base(const hbond_ideal_pt_base& p) : dir_point_t(p)
    {
      do_copy(p);
    }

    const hbond_ideal_pt_base&
    operator=(const hbond_ideal_pt_base& p)
    {
      if(&p != this){
        do_copy(p);
        dir_point_t::operator= (p);
      }
      return *this;
    }

    virtual ~hbond_ideal_pt_base() { ; }

    static bool
    cmp(const hbond_ideal_pt_base &A, const hbond_ideal_pt_base &B)
    { return A.act_type < B.act_type; }

    const SimSite3D::atom_t&
    __get_atom() const
    { return *atom; }

    const SimSite3D::atom_t&
     __get_carbon_nbr() const
     { return *carbon_nbr; }

    uint pt_num;              //!< Number given to the ideal point
    interactionType act_type; //!< Interaction type of the points
    atom_vci atom;            //!< Protein heavy atom giving rise to the points
    atom_vci carbon_nbr;      //!< Carbon atom used to define the acceptor/donor plane
    atom_vci second_nbr;      //!< Third atom used to define the plane

  private:

    void
    do_copy(const hbond_ideal_pt_base& p)
    {
      pt_num = p.pt_num;
      act_type = p.act_type;
      atom = p.atom;
      carbon_nbr = p.carbon_nbr;
      second_nbr = p.second_nbr;
    }

  };
  

  //! An "ideal" hbond point
  class hbond_ideal_point_t : public hbond_ideal_pt_base{
  public:
    hbond_ideal_point_t(alloc_t a = ALLOC_POSITION) : hbond_ideal_pt_base(a)
    {
      pt_num = 0;
      act_type = UNKNOWN_INTERACTION;
      // problem with iterators is cannot initialize atom to anything 
      // safe/usefull
    }

    hbond_ideal_point_t(const hbond_ideal_point_t& p) : hbond_ideal_pt_base(p)
    {
      do_copy(p);
    }

    const hbond_ideal_point_t&
    operator=(const hbond_ideal_point_t& p)
    {
      if(this != &p){
        do_copy(p);
        hbond_ideal_pt_base::operator= (p);
      }
      return *this;
    }

    ~hbond_ideal_point_t() { ; }

    hbond_fit_pt_vci fit_pts_beg;  //!< Const iter to the first fit point corresponding to this ideal point
    hbond_fit_pt_vci fit_pts_end;  //!< Const iter to 1 past the end of the fit points corresponding to this ideal point
  private:
    void
    do_copy(const hbond_ideal_point_t& p)
    {
      fit_pts_beg = p.fit_pts_beg;
      fit_pts_end = p.fit_pts_end;
    }
  };


}

#endif
