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
 * $Source: /psa/share/repository/pfizer_proj/src/basics/point_storage.H,v $
 * $Revision: 1.7 $
 * $Author: vanvoor4 $
 * $Date: 2009-01-12 20:50:42 $
 * 
 * $Log: not supported by cvs2svn $
 * Revision 1.6  2008/07/28 15:14:24  vanvoor4
 * Added function to get all points within some radius of another point
 * Added a number of Doxygen comments
 *
 * Revision 1.5  2008/02/26 19:15:22  vanvoor4
 * Removed the old octree stuff, and started with the new, but not
 * implemented as yet.
 *
 * Revision 1.4  2007/09/20 17:37:27  vanvoor4
 * Added support for Octree
 *
 * Revision 1.3  2007/08/29 20:13:49  vanvoor4
 * Added wrapper to centroid
 *
 * Revision 1.2  2007/08/21 20:59:11  vanvoor4
 * Added __GLIBCXX__ macro to support older versions of gcc
 *
 * Revision 1.1  2007/08/21 15:55:03  vanvoor4
 * Initial checkin
 *
 *
 * 
 * 
 */ 

#ifndef POINT_STORAGE_HEADER_INCLUDED
#define POINT_STORAGE_HEADER_INCLUDED

#include <iostream>
#include <iterator>
#include <map>
#include <point.H>
#include <my_float_array.H>
#include <math_basics.H>
#include <point_bins.H>

namespace ASCbase{

  //! A wrapper class to an STL container class for the point_t and derived
  //! datatypes (classes)
  template<typename _Tp, typename _Sequence = std::vector<_Tp> >
  class point_storage{
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

  public:
    typedef typename _Sequence::value_type		value_type;
    typedef typename _Sequence::reference		reference;
    typedef typename _Sequence::const_reference		const_reference;
    typedef typename _Sequence::iterator		iterator;
    typedef typename _Sequence::const_iterator		const_iterator;
    typedef typename _Sequence::size_type		size_type;
    typedef typename _Sequence::iterator::iterator_category iterator_category;
    typedef          _Sequence				container_type;

    typedef typename point_bins<const_iterator>::bin_vci	bins_vci_type;
    typedef          point_bins<const_iterator>			bins_type;

    typedef std::pair<my_float_t, const_iterator> float_const_iter_pair;
    typedef std::multimap<my_float_t, const_iterator> float_const_iter_map;
  protected:
    _Sequence Container;
    my_float_array positions;

  public:
    point_storage(uint init_num_pts = 150, uint stride_in = 3)
      : positions(init_num_pts, stride_in)
    {
      A_has_been_transformed = false;
    }

    point_storage(const point_storage& src) : positions(src.positions)
    {
      do_copy(src);
    }

    const point_storage&
    operator=(const point_storage& src)
    {
      if(this != &src){
        positions = src.positions;
        do_copy(src);
      }
      return *this;
    }

    virtual ~point_storage()
    {
    }

    //! Constant pointer to the beginning of the points array
    const my_float_t*
    points_begin() const
    { return positions.begin(); }

    //! Constant pointer to one past the end of the points array
    const my_float_t*
    points_end() const
    { return positions.end(); }

    //! Constant iterator to the beginning of the container
    const_iterator
    begin() const 
    { return Container.begin(); }

    //! Constant iterator to 1 past the end of the container
    const_iterator
    end() const 
    { return Container.end(); }

    //! Iterator to the beginning of the container
    iterator
    begin()
    { return Container.begin(); }

    //! Iterator to 1 past the end of the container
    iterator
    end()
    { return Container.end(); }

    //! Get the size of the container
    size_type
    size() const
    { return Container.size(); }

    //! Push an element on the back of the container -- updating the position
    //! pointers if the points array was resized
    void 
    push_back(const value_type& __x)
    {
      Container.push_back(__x);
      value_type& b = Container.back();
      if(positions.push_back(__x.pos) == MY_FLOAT_ARRAY_GREW){
        b.delete_pos();
        my_float_t* p = positions.begin();
        for(iterator s = Container.begin(); 
            p < positions.end() && s != Container.end(); ++s){
          s->pos = p;
          p += positions.stride(); 
        }
      }else{
        b.delete_pos();
        b.pos = positions.end() - positions.stride();
      }
    }

    //! Erase the element pointed to by the iterator
    iterator 
    erase(iterator __i)
    { return Container.erase(__i); }

    // Don't forget to transform A_otree if it is not null
    //! Apply a rigid body transformation to the points (current positions)
    void 
    transform(const my_float_t* R, const my_float_t* T)
    { 
      A_has_been_transformed = true;
      positions.transform3D(R, T); 
      A_bins.transform_bins(R, T);
    }

    //! Apply a rigid body transformation to the points (current positions)
    void 
    inverse_transform(const my_float_t* R, const my_float_t* T)
    { 
      A_has_been_transformed = true;
      positions.inverse_transform3D(R, T); 
      A_bins.inverse_transform_bins(R, T);
    }

    //! Revert the positions to the original positions
    void 
    revert() 
    { 
      positions.revert(); 
      A_bins.revert();
      A_has_been_transformed = false;
    }

    //! Set the current positions as those to revert to when this.revert() is
    //! called
    void
    set_current_positions_as_original()
    {
      positions.set_current_values_as_original();
      if(A_bins.bins_are_setup()){
        std::cerr << "Need to rebin the points!!\n";
      }
    }

    //! Get transform to move current positions to their original positions
    inline void 
    get_current_inverse_3D_transform(Quaternion *Q, my_float_t *T) const
    {
      positions.get_current_inverse_3D_transform(Q, T);
    }

    //! Compute the root mean squared deviation (RMSD) between the current
    //! and orignial positions
    inline my_float_t
    compute_RMSD() const
    { return positions.compute_RMSD(); }


    //! Compute the centroid of the positions
    const bool
    centroid_3D(my_float_t* C) const
    { return positions.centroid_3D(C); }

    //! Get the closest point to the given point a and the distance between a
    //! and the closest point.
    /*!
     * @param a_pos Constant pointer to a position 
     * @param d_out Pointer to a float to store the distance between the two
     *              points
     * @return Constant iterator to the element in the container whose
     *         position was the closest to a
     */
    const_iterator
    closest_point(const my_float_t *a_pos, my_float_t *d_out) const
    { 
      if(A_bins.bins_are_setup()){
        *d_out = my_float_max;
        my_float_t min_dist_sq = my_float_max;
        bins_vci_type bin, moved_bin; 
        A_bins.get_bin(a_pos, &bin, &moved_bin);
	if(bin >= A_bins.bins_end()) return Container.end();

        size_t best_idx = 0;
        for(size_t i = 0; i < bin->positions.size(); ++i){
          my_float_t dist_sq = dist_squared(bin->positions[i], a_pos);
          if(dist_sq < min_dist_sq){
            min_dist_sq = dist_sq;
            best_idx = i;
          }
        }
        const_iterator rv =  bin->assoc_data[best_idx]; 

        bool closer_to_moved_pt = false;
        for(size_t i = 0; i < moved_bin->positions.size(); ++i){
          my_float_t dist_sq = dist_squared(moved_bin->positions[i], a_pos);
          if(dist_sq < min_dist_sq){
            min_dist_sq = dist_sq;
            best_idx = i;
            closer_to_moved_pt = true;
          }
        }
	if(closer_to_moved_pt) rv = moved_bin->assoc_data[best_idx]; 

        *d_out = std::sqrt(min_dist_sq);
        return rv;
      }else
        return _closest_point(a_pos, d_out, iterator_category()); 
    }

    //! Random access iterator version of closest_point
    const_iterator
    _closest_point(const my_float_t *a_pos, my_float_t *d_out, 
                   std::random_access_iterator_tag tag) const
    {
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
      *d_out = min_dist;
      return Container.begin() + best_idx;
    }

    //! Forward iterator version of closest_point
    const_iterator
    _closest_point(const my_float_t *a_pos, my_float_t* d_out, 
                   std::forward_iterator_tag tag) const
    {
      my_float_t min_dist = my_float_max; 
      const_iterator best;
      const_iterator curr = Container.begin();
      my_float_t* b_pos;
      for(b_pos = positions.begin(); b_pos < positions.end(); b_pos += 3){
        my_float_t d = dist(a_pos, b_pos);
	      if(d < min_dist){
          best = curr;
	        min_dist = d;
	      }
	      ++curr;
      }
      *d_out = min_dist;
      return best;
    }

    //! Find all points within some radius of a
    void
    close_points(const my_float_t *a_pos, const my_float_t radius,
                 float_const_iter_map* pts_map) const
    {
      const my_float_t radius_squared = radius*radius;

      if(A_bins.bins_are_setup()){
        bins_vci_type bin, moved_bin; 
        A_bins.get_bin(a_pos, &bin, &moved_bin);
        for(size_t i = 0; i < bin->positions.size(); ++i){
          my_float_t d_squared = dist_squared(bin->positions[i], a_pos);
          if(d_squared <= radius_squared)
            pts_map->insert(float_const_iter_pair(std::sqrt(d_squared), 
                                                  bin->assoc_data[i]));
        }
        for(size_t i = 0; i < moved_bin->positions.size(); ++i){
          my_float_t d_squared = dist_squared(moved_bin->positions[i], a_pos);
          if(d_squared <= radius_squared)
            pts_map->insert(float_const_iter_pair(std::sqrt(d_squared), 
                                                  moved_bin->assoc_data[i]));
        }
      }else{
        const_iterator iter = Container.begin();
        my_float_t* b_pos = positions.begin();
        for( ; iter != Container.end(); ++iter, b_pos += 3){
          my_float_t d_squared = dist_squared(a_pos, b_pos);
          if(radius_squared > d_squared)
            pts_map->insert(float_const_iter_pair(std::sqrt(d_squared), iter));
        }
      }
    }

    //! Reference to the front of the container (first point)
    reference
    front()
    { return Container.front(); }

    //! Constant reference to the front of the container (first point)
    const reference
    front() const
    { return Container.front(); }

    //! Reference to the back of the container (last point)
    reference
    back()
    { return Container.back(); }

    //! Constant reference to the back of the container (last point)
    const reference
    back() const 
    { return Container.back(); }


#if 0
    // Removed the const function "descriptor" since the first time through
    // we will need to build the octree -- octree is not built first since
    // we do not want the over head if the otree will not be used
    const_iterator
    closest_point_test(const my_float_t *a_pos, my_float_t *_d)
    { return _closest_point_test(a_pos, _d, iterator_category()); }

    const_iterator
    _closest_point_test(const my_float_t *a_pos, my_float_t *_d, 
                        std::random_access_iterator_tag tag)
    {
      if(A_otree == 0) build_otree();
    }
#endif

    //! Set up the grid/bins that can be used for neighbor checks of no more
    //! than max_dist
    /*! Only call this function after the container is complete (no more objects
     *  will be added) and before closest_point, etc.
     */
    void
    bin_points(const my_float_t max_dist)
    {
      if(A_has_been_transformed)  
        std::cerr << "Cannot bin points after a global transform has been "
                  << "applied\n"; 
      else if(!A_bins.bins_are_setup())
        A_bins.setup_grid(Container.begin(), Container.end(), max_dist);
    }

    inline void
    get_bin(const my_float_t *pt, bins_vci_type *nbrs, 
            bins_vci_type *moved_nbrs) const
    { return A_bins.get_bin(pt, nbrs, moved_nbrs); }

    inline void
    write_bin_bounds(std::ostream &out) const
    { A_bins.write_bounds(out); }

  private:

    //! This function has not been verified
    void
    do_copy(const point_storage& src)
    {
      A_has_been_transformed = src.A_has_been_transformed;
      Container.clear();
      Container.reserve(src.Container.size());
      const_iterator c;
      for(c = src.Container.begin(); c != src.Container.end(); ++c)
        Container.push_back(*c);

      //std::copy(src.Container.begin(), src.Container.end(), Container);

      my_float_t* p = positions.begin();
      for(iterator s = Container.begin(); 
          p < positions.end() && s != Container.end(); ++s){
        s->pos = p;
        p += positions.stride(); 
      }
    }

    // We are prototyping methods that can make use of an axis aligned
    // grid.  In the future, if we want a datastructure that can handle
    // global transformations, we could consider spheres centered at
    // the center of each grid cell and with a radius of the distance
    // from a grid cell center to one of that grid cell's end corners
    
    bool A_has_been_transformed;
    point_bins<const_iterator> A_bins;
  
  };

}

#endif
