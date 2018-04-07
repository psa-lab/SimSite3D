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
 * $Source: /psa/share/repository/pfizer_proj/src/basics/Octree.H,v $
 * $Revision: 1.1 $
 * $Author: vanvoor4 $
 * $Date: 2007-09-20 18:17:39 $
 * 
 * $Log: not supported by cvs2svn $
 * 
 * 
 */ 

#ifndef OCTREE_HEADER_INCLUDED
#define OCTREE_HEADER_INCLUDED

#include <types.H>
#include <cstring>
#include <vector>
#include <string_basics.H>

template<typename _Tp, typename _Sequence = std::vector<_Tp> >
class octree_data_t{
public:
  //typedef typename _Sequence::value_type              value_type;
  //typedef typename _Sequence::reference               reference;
  //typedef typename _Sequence::const_reference         const_reference;
  //typedef typename _Sequence::iterator                iterator;
  //typedef typename _Sequence::const_iterator          const_iterator;
  //typedef typename _Sequence::size_type               size_type;
  //typedef typename _Sequence::iterator::iterator_category iterator_category;
  typedef          _Sequence                          container_type;

  octree_data_t()
  {
    a_head = 0;
    std::memset(a_children, 0, 8 * sizeof(octree_data_t*));
    std::memset(a_box_min, 0, 3*my_float_size);
    std::memset(a_box_max, 0, 3*my_float_size);
    std::memset(a_centroid, 0, 3*my_float_size);
  }

  //! Straight copy -- of course both objects will have same data
  //! Current assumption is the value_type and container_type both have
  //! valid copy constructors
  octree_data_t(octree_data_t& src)
  {
    if(this == &src) return;

    a_head = src.a_head;
    std::copy(src.a_children, src.a_children + 8, a_children);
    std::copy(src.a_box_min, src.a_box_min + 3, a_box_min);
    std::copy(src.a_box_max, src.a_box_max + 3, a_box_max);
    std::copy(src.a_centroid, src.a_centroid + 3, a_centroid);
    a_data = src.a_data; 
  }

  octree_data_t* a_head;
  octree_data_t* a_children[8];
  my_float_t a_box_min[3];    //!< Min corner of box
  my_float_t a_box_max[3];    //!< Max corner of box
  my_float_t a_centroid[3];

  container_type a_data; 
private:
};

//! Simple Octree for any datatype with an associated 3D position
/*! It is reasonably efficient when compared with brute force, but a few more
 * checks and the number of returned points, etc would be reduced further.
 */
template<typename _Tp, typename _Sequence = std::vector<_Tp> >
class Octree{
public:
  typedef typename _Sequence::value_type              value_type;
  //typedef typename _Sequence::reference               reference;
  typedef typename _Sequence::const_reference         const_reference;
  //typedef typename _Sequence::iterator                iterator;
  typedef typename _Sequence::const_iterator          const_iterator;
  //typedef typename _Sequence::size_type               size_type;
  //typedef typename _Sequence::iterator::iterator_category iterator_category;
  typedef          _Sequence                          container_type;

  Octree(/*const uint max_nsubdiv,*/ const my_float_t min_side_len,
         const uint max_numel = 10)
  {
    //a_max_nsubdiv = max_nsubdiv;
    a_min_side_len = min_side_len;
    a_max_numel = max_numel;
  }


  //! The positions are used only in building of the tree and may or may not
  //! be stored in the data -- this is entirely up to the user of the class
  //! The assumption is that the number of positions is the same as the number
  //! of the items in data.  Furthermore, the data must be stored in the
  //! container type used to construct the Octree class
  void
  build(const my_float_t* positions_begin, const my_float_t* postions_end,
        const my_float_t* centroid, const my_float_t root_side_len,
        const container_type& data)
  {
    a_root_side_len = root_side_len; 
    std::copy(centroid, centroid + 3, a_root.a_centroid);
    std::vector<uint> idz(data.size());
    uint i = 0;
    for(const_iterator r = data.begin(); r != data.end(); ++r, ++i){   
      a_root.a_data.push_back(r);
      idz[i] = i;
    }
    subdivide(&a_root, idz, positions_begin, a_root_side_len); 
  }

  void 
  subdivide(octree_data_t<const_iterator>* cnode, 
            const std::vector<uint>& idz, const my_float_t* positions_begin,
            const my_float_t side_len)
  {
    std::cout << "side_len is " << side_len << std::endl;
    std::cout << "number of elements in cube is " << cnode->a_data.size() << std::endl;

    my_float_t h_width = 0.5 * side_len;
    if(h_width <= a_min_side_len) return;
 
    std::vector< std::vector<uint> > idz_vec(8);
    my_float_t* centroid = cnode->a_centroid;
  
    // partition based on the position of each point 
    uint cube_idx;
    for(uint i = 0; i < cnode->a_data.size(); ++i){
      const my_float_t* pos = positions_begin + 3*idz[i];
      cube_idx = (pos[0] > centroid[0] ? 0 : 1); 
      cube_idx += (pos[1] > centroid[1] ? 0 : 2); 
      cube_idx += (pos[2] > centroid[2] ? 0 : 4); 
  
      // add to smaller cube
      octree_data_t<const_iterator>* child = *(cnode->a_children + cube_idx);
      if(!child){
        *(cnode->a_children + cube_idx) = new octree_data_t<const_iterator>;
        child = *(cnode->a_children + cube_idx);
        child->a_head = cnode;
        uint mask = 1; 
        for(uint j = 0; j < 3; ++j, mask <<= 1)
          child->a_centroid[j] = (cube_idx & mask ? centroid[j] - h_width : 
                                                    centroid[j] + h_width);
      }
      idz_vec[cube_idx].push_back(idz[i]);
      child->a_data.push_back(cnode->a_data[i]); 
    }

    // Parent is not leaf node so clear its data
    cnode->a_data.clear();
    for(uint i = 0; i < 8; ++i)
      if(idz_vec[i].size() > a_max_numel) 
        subdivide(cnode->a_children[i], idz_vec[i], positions_begin, h_width);
  }

  void 
  near_by_points(const my_float_t* pos, const my_float_t max_dist, 
                 std::vector<const_iterator>* pts)
  {
    uint cube_idx;
    const my_float_t* centroid = a_root.a_centroid;
    my_float_t side_len = a_root_side_len;
    octree_data_t<const_iterator>* cnode = &a_root;
    while(side_len > max_dist){
      cube_idx = (pos[0] > centroid[0] ? 0 : 1);
      cube_idx += (pos[1] > centroid[1] ? 0 : 2);
      cube_idx += (pos[2] > centroid[2] ? 0 : 4);
 
      if(cnode->a_children[cube_idx] == 0) break;
      cnode = cnode->a_children[cube_idx];
            
        
      side_len *= 0.5;
    }

    pts->clear();
    get_pts(cnode, pts);
  }

private:
  void
  get_pts(const octree_data_t<const_iterator>* cnode, 
          std::vector<const_iterator>* pts)
  {
    bool leaf = true;
    for(uint i = 0; i < 8; ++i)
      if(cnode->a_children[i]){
        leaf = false;
        get_pts(cnode->a_children[i], pts);
      }     

    if(!leaf) return;
    
    for(uint i = 0; i < cnode->a_data.size(); ++i)
      pts->push_back(cnode->a_data[i]);
  } 

  octree_data_t<const_iterator> a_root;
  my_float_t a_root_side_len;
  uint a_max_numel; //!< Max number of elements in a cube at a level less than the max level
//  uint a_max_nsubdiv;  //!< Maximum level of splitting 
  
  //my_float_t a_box_min[3];   //!< Min corner of root cube
  //my_float_t a_box_max[3];   //!< Max corner of root cube
  my_float_t a_min_side_len; //!< Smallest cube should be about this size

  static const std::string a_fname;
};

// not sure how do do this
//const std::string Octree<_Tp>::a_fname = "Octree.C";

#endif
