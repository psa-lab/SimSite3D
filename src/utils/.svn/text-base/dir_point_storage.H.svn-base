/******************************************************************************
 * Copyright (c) 2006,2007, Michigan State University (MSU) Board of Trustees.
 *   All rights reserved.
 *
 * This file is part of the ASCbase Software project.
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
 * $Source: /psa/share/repository/pfizer_proj/src/basics/dir_point_storage.H,v $
 * $Revision: 1.3 $
 * $Author: vanvoor4 $
 * $Date: 2007-08-29 20:17:00 $
 * 
 * $Log: not supported by cvs2svn $
 * Revision 1.2  2007/08/21 20:58:59  vanvoor4
 * Added __GLIBCXX__ macro to support older versions of gcc
 *
 * Revision 1.1  2007/08/21 15:32:06  vanvoor4
 * Initial checkin
 *
 *
 * 
 * 
 */ 

#ifndef DIR_POINT_STORAGE_HEADER_INCLUDED
#define DIR_POINT_STORAGE_HEADER_INCLUDED

#include <point_storage.H>

namespace ASCbase{

  //! Storage class for a directed point (point and associated normal direction)
  template<typename _Tp, typename _Sequence = std::vector<_Tp> >
  class dir_point_storage : public point_storage<_Tp, _Sequence>{
    // Concept requirements
    typedef typename _Sequence::value_type _Sequence_value_type;
#ifdef __GLIBCXX__
    __glibcxx_class_requires(_Tp, _SGIAssignableConcept);
    __glibcxx_class_requires(_Sequence, _BackInsertionSequenceConcept);
    __glibcxx_class_requires2(_Tp, _Sequence_value_type, _SameTypeConcept);
#else
    __glibcpp_class_requires(_Tp, _SGIAssignableConcept);
    __glibcpp_class_requires(_Sequence, _BackInsertionSequenceConcept);
    __glibcpp_class_requires2(_Tp, _Sequence_value_type, _SameTypeConcept);
#endif

    typedef point_storage<_Tp, _Sequence>		_Base;

  public:
    typedef typename _Sequence::value_type		value_type;
    typedef typename _Sequence::reference		reference;
    typedef typename _Sequence::const_reference		const_reference;
    typedef typename _Sequence::iterator		iterator;
    typedef typename _Sequence::const_iterator		const_iterator;
    typedef typename _Sequence::size_type		size_type;
    typedef          _Sequence				container_type;

  protected:
    my_float_array directions;

  public:
    dir_point_storage(uint init_num_pts = 150, uint stride_in = 3)
      : _Base(init_num_pts, stride_in), directions(init_num_pts, stride_in)
    {
    }

    dir_point_storage(const dir_point_storage& src)
      : _Base(src), directions(src.directions)
    {
      do_copy(src);
    }

    const dir_point_storage&
    operator=(const dir_point_storage& src)
    {
      if(this != &src){
        _Base::operator= (src);
        directions = src.directions;
        do_copy(src);
      }
      return *this;
    }

    ~dir_point_storage()
    {
    }

    const my_float_t*
    dirs_begin() const
    { return directions.begin(); }

    const my_float_t*
    dirs_end() const
    { return directions.end(); }

    void 
    push_back(const value_type& __x)
    {
      _Base::push_back(__x);
      value_type& b = _Base::back();
      b.delete_dir();
      if(directions.push_back(__x.dir) == MY_FLOAT_ARRAY_GREW){
        my_float_t* d = directions.begin();
        for(iterator s = _Base::Container.begin(); 
            d < directions.end() && s != _Base::Container.end(); ++s){
          s->dir = d;
          d += directions.stride(); 
        }
      }else b.dir = directions.end() - directions.stride();
    }

    void 
    transform(const my_float_t* R, const my_float_t* T)
    { 
      const my_float_t zero[] = {0, 0, 0};
      directions.transform3D(R, zero); 
      _Base::transform(R, T); 
    }

    void 
    inverse_transform(const my_float_t* R, const my_float_t* T)
    { 
      const my_float_t zero[] = {0, 0, 0};
      directions.inverse_transform3D(R, zero); 
      _Base::inverse_transform(R, T); 
    }

    void 
    revert()
    { 
      directions.revert(); 
      _Base::revert();
    }

    void
    set_current_positions_and_directions_as_original()
    {
      directions.set_current_values_as_original();
      _Base::set_current_positions_as_original();
    }

/*
    _Base::const_iterator
    corresponding_point(const my_float_t* P, const my_float_t* D)
    {
      const my_float_t* d = directions.begin();
      const my_float_t* p = positions.begin();

      my_float_t min_dist = my_float_max; 
      uint best_idx = 0;
      uint cnt = 0;
      my_float_t* b_pos;
      for(b_pos = positions.begin(); b_pos < positions.end(); b_pos += 3){
        my_float_t d = dist(a_pos, b_pos);
	if(d < min_dist){
          best_idx = cnt;
	  min_dist = d;
	}
	++cnt;
      }
      return Container.begin() + best_idx;
    }
*/

  private:
  
    void
    do_copy(const dir_point_storage& src)
    {
      my_float_t* d = directions.begin();
      for(iterator s = _Base::Container.begin(); 
          d < directions.end() && s != _Base::Container.end(); ++s){
        s->dir = d;
        d += directions.stride(); 
      }
    }

  };
}

#endif
