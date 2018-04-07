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
 * $Source: /psa/share/repository/pfizer_proj/src/basics/SimpleOctree.H,v $
 * $Revision: 1.3 $
 * $Author: vanvoor4 $
 * $Date: 2009-01-12 21:05:24 $
 * 
 * $Log: not supported by cvs2svn $
 * Revision 1.2  2008/05/13 15:50:02  vanvoor4
 * Not sure why I had previously decided to copy the data in a
 * node container rather than returning a pointer to the data
 * container.
 *
 * Revision 1.1  2008/02/26 19:14:30  vanvoor4
 * Initial checkin of a "simple" Octree.  Just need something that works
 * for now and is faster than brute force.  The actual tree method is
 * no different from the previous implementation, but the data supported
 * is not quite as general.
 *
 * 
 */ 

#ifndef SIMPLE_OCTREE_HEADER_INCLUDED
#define SIMPLE_OCTREE_HEADER_INCLUDED

#include <vector>
#include <queue>
#include <algorithm>
#include <iostream>
#include <cstring>
#include <types.H>

template<typename _Tp>
class octree_node_t{
  typedef _Tp					value_type;

public:

  octree_node_t()
  {
    parent = 0;
    half_width = 0;
    std::memset(children, 0, 8 * sizeof(octree_node_t*));
    std::fill(centroid, centroid + 3, 0);
  }

  //! Straight copy -- of course both objects will have same data
  //! Current assumption is the value_type and container_type both have
  //! valid copy constructors
  octree_node_t(octree_node_t& src)
  {
    if(this == &src) return;

    parent = src.parent;
    half_width = src.half_width;
    std::copy(src.children, src.children + 8, children);
    std::copy(src.centroid, src.centroid + 3, centroid);
    data = src.data; 
  }

  octree_node_t* parent;
  octree_node_t* children[8];
  my_float_t centroid[3];
  my_float_t half_width;  

  std::vector<_Tp> data;   //!< Container holding the data
};

//! Simple Octree for any datatype with an associated 3D position
/*! It is reasonably efficient when compared with brute force, but a few more
 * checks and the number of returned points, etc would be reduced further.
 */
template<typename _Tp>
class SimpleOctree{
  typedef _Tp					value_type;

public:

  typedef struct{
    const std::vector<_Tp>* bin_data;
    std::vector<const my_float_t*> positions;
  }bin_type;


  SimpleOctree(const my_float_t min_half_width, const uint max_numel = 1)
  {
    A_min_half_width = min_half_width;
    A_max_numel = max_numel;
  }

  ~SimpleOctree()
  {
    // Need to deallocate the tree here
  }

  void
  build(_Tp pts_begin, _Tp pts_end)
  {
    if(pts_end <= pts_begin) return;

    // Spin through the positions to find the min and max values
    my_float_t max_pt[3], min_pt[3];
    std::copy(pts_begin->pos, pts_begin->pos + 3, max_pt);
    std::copy(pts_begin->pos, pts_begin->pos + 3, min_pt);
    for(_Tp p = pts_begin; p < pts_end; ++p){
      for(uint i = 0; i < 3; ++i){
        if(max_pt[i] < p->pos[i]) max_pt[i] = p->pos[i];
        if(min_pt[i] > p->pos[i]) min_pt[i] = p->pos[i];
      }
    }

    // Compute the centroid and add the tolerance to each coordinate of min 
    // and max to obtain the root side len.
    my_float_t max_side_len = 0;
    for(uint i = 0; i < 3; ++i){
      max_pt[i] += A_min_half_width;
      min_pt[i] -= A_min_half_width; 
      my_float_t len = max_pt[i] - min_pt[i];
      if(len > max_side_len) max_side_len = len;
    }

    my_float_t root_side_len = A_min_half_width;
    while(max_side_len > A_min_half_width){
      max_side_len *= 0.5;
      root_side_len *= 2;
    }

    for(uint i = 0; i < 3; ++i)
      A_root.centroid[i] = 0.5 * (max_pt[i] + min_pt[i]);
    A_root.half_width = 0.5 * root_side_len;

    // Copy data to root node
    A_root.data.reserve(pts_end - pts_begin);
    for(_Tp p = pts_begin; p < pts_end; ++p) A_root.data.push_back(p);
    
    subdivide_with_overlap(&A_root);
  }


  // Should make non overlapping bins  -- this makes the 
  // "search" easier at the cost of storing more pointers.  However, to me,
  // this seems like a much easier method to implement and a bit faster 
  // than searching for the nodes that are near a query point.
  void
  subdivide_with_overlap(octree_node_t<_Tp>* cnode)
  {
    if(cnode->data.size() <= A_max_numel) return;
    my_float_t h_width = 0.5 * cnode->half_width;
    if(h_width < A_min_half_width) return;

    // partition based on the position of each point 
    const my_float_t* centroid = cnode->centroid;

    std::vector<_Tp>& data = cnode->data;
    for(size_t pt_idx = 0; pt_idx < data.size(); ++pt_idx){
      const _Tp& pt = data[pt_idx];
   
      std::vector<bool> above(3, false), below(3, false);
      for(uint i = 0; i < 3; ++i){
        if((pt->pos[i] - A_min_half_width) <= centroid[i]) below[i] = true;
        if((pt->pos[i] + A_min_half_width) >= centroid[i]) above[i] = true;
      }

      // Explicitly write out the 8 cases -- likely to be easier to read
      std::vector<bool> pt_in_bins(8, false);
      pt_in_bins[0] = (below[0] && below[1] && below[2]);
      pt_in_bins[1] = (above[0] && below[1] && below[2]);
      pt_in_bins[2] = (below[0] && above[1] && below[2]);
      pt_in_bins[3] = (above[0] && above[1] && below[2]);
      pt_in_bins[4] = (below[0] && below[1] && above[2]);
      pt_in_bins[5] = (above[0] && below[1] && above[2]);
      pt_in_bins[6] = (below[0] && above[1] && above[2]);
      pt_in_bins[7] = (above[0] && above[1] && above[2]);
       
      for(uint i = 0; i < 8; ++i){
        if(!pt_in_bins[i]) continue;

        octree_node_t<_Tp>* child = cnode->children[i];
        if(!child){
          cnode->children[i] = new octree_node_t<_Tp>;
	  child = cnode->children[i];
	  child->parent = cnode;
	  child->half_width = h_width;
          uint dim_mask = 1; 
          for(uint j = 0; j < 3; ++j, dim_mask <<= 1)
	    if(i & dim_mask) child->centroid[j] = centroid[j] + h_width;
	    else child->centroid[j] = centroid[j] - h_width;
        }
        child->data.push_back(pt);
      }
    }

    // Recursively subdivide children
    for(uint i = 0; i < 8; ++i)
      if(cnode->children[i]){
        subdivide_with_overlap(cnode->children[i]);
        cnode->data.clear();
      }
  }

  const std::vector<_Tp>*
  near_by_points(const my_float_t* pos, const my_float_t max_dist)
  {
    const octree_node_t<_Tp>* cnode = &A_root;
    my_float_t h_width;
    const my_float_t h_max_dist = 0.5 * max_dist;
    for(h_width = A_root.half_width; h_width > h_max_dist; h_width *= 0.5){
      uint cube_idx = 0;
      if(pos[0] > cnode->centroid[0]) cube_idx += 1;
      if(pos[1] > cnode->centroid[1]) cube_idx += 2;
      if(pos[2] > cnode->centroid[2]) cube_idx += 4;

      if(cnode->children[cube_idx] == 0) break;
      cnode = cnode->children[cube_idx];
    }
   
    return &(cnode->data);
  }

  //! Use a depth first search to bin the given positions -- of course,
  //! positions falling outside of the root node or in a an empty bin are
  //! discarded.
  /*! This function was not tested enough and it should also be timed to
   * check if it is really any faster than near_by_points for each point
   */
  void
  bin_positions(const std::vector<const my_float_t*> &pos_vec, 
                const my_float_t max_dist,
                std::vector<bin_type> *bins)
  {
    bin_positions(A_root, max_dist, pos_vec, bins);
  }

private:

  // Needs more testing
  void
  bin_positions(octree_node_t<_Tp> &cnode, const my_float_t max_dist,
                const std::vector<const my_float_t*> &pos_vec,
                std::vector<bin_type> *bins)
  {
    // If data is stored in the current it node, the assumption is the current
    // node is a leaf node.
    if(cnode.data.size()){
      bin_type bin;
      bin.bin_data = &(cnode.data);
      bin.positions.resize(pos_vec.size()); 
      std::vector<const my_float_t*>::const_iterator p_iter;
      std::vector<const my_float_t*>::iterator dest = bin.positions.begin();
      for(p_iter = pos_vec.begin(); p_iter < pos_vec.end(); ++p_iter, ++dest)
        *dest = *p_iter;
      bins->push_back(bin);
      return;
    }
    if(cnode.half_width <= max_dist){
      std::cerr << "Unhandled condition in bin_positions!" << std::endl;
      return;
    }

    std::vector< std::vector<const my_float_t*> > children_vectors(8);
    std::vector<const my_float_t*>::const_iterator pos_iter;
    for(pos_iter = pos_vec.begin(); pos_iter < pos_vec.end(); ++pos_iter){
      const my_float_t *pos = *pos_iter;
      size_t cube_idx = 0;
      if(*pos > cnode.centroid[0]) cube_idx += 1;
      ++pos;
      if(*pos > cnode.centroid[1]) cube_idx += 2;
      ++pos;
      if(*pos > cnode.centroid[2]) cube_idx += 4;

      if(*(cnode.children + cube_idx)) 
        (children_vectors.begin() + cube_idx)->push_back(*pos_iter);
    }

    octree_node_t<_Tp> **child_p = cnode.children;
    for(size_t i = 0; i < 8; ++i, ++child_p)
      if(*child_p) 
        bin_positions(**child_p, max_dist, children_vectors[i], bins);
  }

/*
  bool 
  point_in_node(const my_float_t *pt, const octree_node_t &cnode, 
                const my_float_t pt_tol)
  {
    my_float_t dist[3];
    for(size_t i = 0; i < 3; ++i){
      dist[i] = cnode.centroid[i] - pt[i];     
      dist[i] = (dist[i] < 0 ? -1.0 * dist[i] : dist[i]);
    }
 
    // If pt in box, return true
    if(dist[0] <= cnode.half_width && dist[1] <= cnode.half_width &&
       dist[2] <= cnode.half_width) return true;

    // If pt not in box + tol, return false
    if(dist[0] <= cnode.half_width + pt_tol && 
       dist[1] <= cnode.half_width + pt_tol &&
       dist[2] <= cnode.half_width + pt_tol) return false;
 
    // Project centroid to sphere with center at pt and radius pt_tol
    sphere S(pt, pt_tol);
    my_float_t pt_on_sphere[3];
    S.proj2surf(cnode.centroid, pt_on_sphere);

    // If projected centroid is in the box, return true
    // else return false 
    if(dist_squared(pt_on_sphere, pt) >= pt_tol * pt_tol) return false;
    return true; 
  }
*/                
 
  octree_node_t<_Tp> A_root;  //!< Root node of the tree
  uint A_max_numel; //!< Max number of elements in a cube at a level less than the max level
  my_float_t A_min_half_width;   //!< Smallest cube should be about this size

  static const std::string A_fname;  //!< Name of the source file
};


#endif
