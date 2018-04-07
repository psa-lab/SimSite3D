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
 * $Source: /psa/share/repository/pfizer_proj/src/basics/point.H,v $
 * $Revision: 1.3 $
 * $Author: vanvoor4 $
 * $Date: 2009-01-12 20:50:42 $
 * 
 * $Log: not supported by cvs2svn $
 * Revision 1.2  2007/12/17 20:55:06  vanvoor4
 * Added missing header files.
 *
 * Revision 1.1  2007/08/21 15:54:43  vanvoor4
 * initial checkin
 *
 *
 * 
 * 
 */ 

#ifndef _POINT_HEADER_INCLUDED
#define _POINT_HEADER_INCLUDED

#include <types.H>
#include <algorithm>
#include <vector>

namespace SimSite3D{

  //! Listing of the possible interaction types for a given vertex
  typedef enum{
    UNKNOWN_INTERACTION =    0,
    ACCEPTOR            =    1,
    DONOR               =    4,
    HYDROPHOB           =   16,
    DONEPTOR            =   64,
    DONOR_HYDROGEN      =  256,
    METAL_1             = 1024,
    METAL_2             = 4096,
    HBVEC,
    NOTHING
  }interactionType;

  typedef enum{
    NO_ALLOCATION,
    ALLOC_POSITION
  }alloc_t;

  //! Base class for a point that supports local allocation of position or
  //! allocation within a class such as the point_storage class.
  class point_t{
  public:
    point_t(alloc_t a = ALLOC_POSITION){
      A_mem_allocated = a;
      if(a == ALLOC_POSITION) pos = new my_float_t[3];
      else pos = 0;
    } 

    point_t(const point_t& p)
    {
      do_copy(p);
    }

    const point_t&
    operator=(const point_t& p)
    {
      if(&p != this) do_copy(p);
      return *this;
    }

    ~point_t() 
    {
      delete_pos();
    }

    void delete_pos()
    {
      if(A_mem_allocated == ALLOC_POSITION && pos) delete [] pos;
      pos = 0;
      A_mem_allocated = NO_ALLOCATION;
    }

    my_float_t* pos; 
  private:
    void
    do_copy(const point_t& p)
    {
      A_mem_allocated = p.A_mem_allocated;
      if(p.A_mem_allocated == ALLOC_POSITION){
        pos = new my_float_t[3];
        std::copy(p.pos, p.pos+3, pos);
      }else pos = p.pos;
    }

    alloc_t A_mem_allocated;
  };

  typedef std::vector<point_t> point_vec_t;
  typedef point_vec_t::iterator point_vi;
  typedef point_vec_t::const_iterator point_vci;

}

#endif
