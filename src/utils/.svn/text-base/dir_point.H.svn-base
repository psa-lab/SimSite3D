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

#ifndef _DIRECTED_POINT_HEADER_INCLUDED
#define _DIRECTED_POINT_HEADER_INCLUDED

#include <point.H>

namespace SimSite3D{
  
  //! A point and associated direction
  class dir_point_t : public point_t{
  public:
    dir_point_t(alloc_t a = ALLOC_POSITION) : point_t(a)
    {
      _dir_allocated = a;
      if(a == ALLOC_POSITION) dir = new my_float_t[3];
      else dir = 0;
    } 

    dir_point_t(const dir_point_t& p) : point_t(p)
    {
      do_copy(p);
    }

    const dir_point_t&
    operator=(const dir_point_t& p)
    {
      if(&p != this){
        do_copy(p);
        point_t::operator= (p);
      }
      return *this;
    }

    ~dir_point_t() 
    {
      delete_dir();
    }

    void delete_dir()
    {
      if(_dir_allocated == ALLOC_POSITION && dir) delete [] dir;
      dir = 0;
      _dir_allocated = NO_ALLOCATION;
    }

    my_float_t* dir; 
  private:
    void
    do_copy(const dir_point_t& p)
    {
      _dir_allocated = p._dir_allocated;
      if(p._dir_allocated == ALLOC_POSITION){
        dir = new my_float_t[3];
        std::copy(p.dir, p.dir+3, dir);
      }else dir = p.dir;
    }

    alloc_t _dir_allocated;
  };

  typedef std::vector<dir_point_t> dir_point_vec_t;
  typedef dir_point_vec_t::iterator dir_point_vi;
  typedef dir_point_vec_t::const_iterator dir_point_vci;

}

#endif
