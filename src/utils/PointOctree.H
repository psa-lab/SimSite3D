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
  
#ifndef POINT_OCTREE_HEADER_INCLUDED
#define POINT_OCTREE_HEADER_INCLUDED

#include <vector>
#include <queue>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <types.H>
#include <point.H>
#include <atom.H>

namespace SimSite3D{

// #define TRACE

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

//! Adaptive octree for the point_t datatype 
/*! It is reasonably efficient when compared with brute force, but a few more
 * checks and the number of returned points, etc would be reduced further.
 */
class PointOctree{

public:

  typedef enum{
    NO_INTERSECTION = 0,
    CENTER_IN_CUBE,
    POSSIBLE_INTERSECTION,
    CUBE_AND_SPHERE_INTERSECT
  }cube_sphere_intersect_t;

  PointOctree(atom_vci first, atom_vci last,
              const my_float_t min_half_width = 4.5, const uint max_numel = 1)
  {
    A_min_half_width = min_half_width;
    A_max_numel = max_numel;
    build(first, last);
  }

  ~PointOctree()
  {
    free_octree_node(&A_root);
  }

  void
  get_close_points(atom_vci a, const my_float_t max_dist, 
                   std::vector<atom_vci> *close_points)
  { 
    get_close_points(&A_root, a->pos, max_dist, close_points, 0);
  }

  void 
  get_close_points(const octree_node_t<atom_vci> *cnode, const my_float_t *pt,
                   const my_float_t pt_tol, 
                   std::vector<atom_vci> *close_points, int level)
  {
    if(cnode->data.size()){
      // Allow push_back to choose when to grow the vector and see how it goes
      //close_points->reserve(cnode->data.size() + close_points->size());
      std::vector<atom_vci>::const_iterator data_iter = cnode->data.begin();
      for( ; data_iter < cnode->data.end(); ++data_iter)
        close_points->push_back(*data_iter);
      return;
    }

    for(int i = 0; i < 8; ++i){
      if(cnode->children[i] == 0) continue;

      cube_sphere_intersect_t rv;
      rv = point_in_cube(cnode->children[i]->centroid, 
                         cnode->children[i]->half_width, pt, pt_tol);
      if(rv == CENTER_IN_CUBE ||
         (rv == POSSIBLE_INTERSECTION &&
          cube_and_sphere_intersect(cnode->children[i]->centroid,
                                    cnode->children[i]->half_width, 
                                    pt, pt_tol))){ 
        get_close_points(cnode->children[i], pt, pt_tol, close_points,
                         level + 1); 
      }
    }
  }

  void
  bin_trace()
  {
    print_bin_trace(&A_root, 0);
  }

private:

  void
  build(atom_vci first, atom_vci last)
  {
    // concept requirements
    __glibcxx_requires_valid_range(first, last);

    // Spin through the positions to find the min and max values
    my_float_t max_pt[3], min_pt[3];
    std::copy(first->pos, first->pos + 3, max_pt);
    std::copy(first->pos, first->pos + 3, min_pt);
    for(atom_vci p = first; p < last; ++p){
      for(uint i = 0; i < 3; ++i){
        if(max_pt[i] < p->pos[i]) max_pt[i] = p->pos[i];
        if(min_pt[i] > p->pos[i]) min_pt[i] = p->pos[i];
      }
    }

    // 
    // Use a square at the cost of more storage (larger tree) and possibly
    // a longer search time.  We do want a reasonable tolerance on the added
    // to each side so that points do not easily move outside
    //

    // Compute the centroid and add the tolerance to each coordinate of min 
    // and max to obtain the root side len.
    my_float_t max_side_len = 0;
    for(uint i = 0; i < 3; ++i){
      max_pt[i] += A_min_half_width;
      min_pt[i] -= A_min_half_width; 
      my_float_t len = max_pt[i] - min_pt[i];
      if(len > max_side_len) max_side_len = len;
    }

    //  
    // Should be an input variable at some point in time -- a Lys or Arg
    // sidechain can easily move more than 10.0 (A)
    //
    max_side_len += 40.0;

    // Tune size to be divisible by min half width otherwise we get 
    // a less tight set of boxes -- need to test on a number of different
    // point sets in order to get an idea of its impact
    my_float_t root_side_len = A_min_half_width;
    for( ; max_side_len > A_min_half_width; 
        max_side_len *= 0.5, root_side_len *= 2);

    for(uint i = 0; i < 3; ++i)
      A_root.centroid[i] = 0.5 * (max_pt[i] + min_pt[i]);
    A_root.half_width = 0.5 * root_side_len;

    // Copy data to root node
    A_root.data.reserve(last - first);
    for(atom_vci p = first; p < last; ++p) A_root.data.push_back(p);
    
    subdivide(&A_root);
  }

  void
  subdivide(octree_node_t<atom_vci>* cnode)
  {
    my_float_t h_width = 0.5 * cnode->half_width;
    if(cnode->data.size() <= A_max_numel || h_width < A_min_half_width) return;

    const my_float_t *centroid = cnode->centroid;
    std::vector<int> octant_count(8);
    std::fill(octant_count.begin(), octant_count.end(), 0);
    std::vector<int> partition(cnode->data.size());

    // Partition based on the position of each point
    std::vector<int>::iterator octant = partition.begin();
    std::vector<atom_vci>::const_iterator data_iter = cnode->data.begin();
    for( ; data_iter < cnode->data.end(); ++data_iter, ++octant){
      int idx = 0;
      //std::cout << **data_iter << "\n";
      for(int i = 2; i > -1; --i){
        idx <<= 1;
        if((*data_iter)->pos[i] >= centroid[i]) ++idx;
      }
      *octant = idx;
      ++octant_count[idx];
    }

    // Allocate the children for the octants with data and the childrens' 
    // data vectors
    for(int idx = 0; idx < 8; ++idx){
      if(octant_count[idx] == 0) continue;

      octree_node_t<atom_vci> *child = new octree_node_t<atom_vci>;
      child->parent = cnode;
      child->half_width = h_width;
      int dim_mask = 1;
      for(int j = 0; j < 3; ++j, dim_mask <<= 1){
        if(idx & dim_mask) child->centroid[j] = centroid[j] + h_width;
        else child->centroid[j] = centroid[j] - h_width;
      }
      child->data.reserve(octant_count[idx]);
      cnode->children[idx] = child;
    }

    // Partition the data based on the partition array
    octant = partition.begin();
    data_iter = cnode->data.begin();
    for( ; data_iter < cnode->data.end(); ++data_iter, ++octant){
       octree_node_t<atom_vci> *child = *(cnode->children + *octant);
       child->data.push_back(*data_iter);
    }
    cnode->data.clear();

    // Subdivide children
    for(int j = 0; j < 8; ++j)
      if(octant_count[j]) subdivide(cnode->children[j]);
#if 0

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
#endif
  }

#if 0
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
#endif

  //! Use a depth first search to bin the given positions -- of course,
  //! positions falling outside of the root node or in a an empty bin are
  //! discarded.
  /*! This function was not tested enough and it should also be timed to
   * check if it is really any faster than near_by_points for each point
   */
#if 0
  void
  bin_positions(const std::vector<const my_float_t*> &pos_vec, 
                const my_float_t max_dist,
                std::vector<bin_type> *bins)
  {
    bin_positions(A_root, max_dist, pos_vec, bins);
  }
#endif

  void
  free_octree_node(octree_node_t<atom_vci>* cnode)
  {
    bool leaf_node = false;
    if(cnode->data.size()){
      cnode->data.clear();
      leaf_node = true;
    }
    if(leaf_node) return;

    for(int i = 0; i < 8; ++i)
      if(cnode->children[i]){
        free_octree_node(cnode->children[i]);
        delete cnode->children[i];
        cnode->children[i] = 0;
      }
  }

#if 0
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
#endif

  inline cube_sphere_intersect_t
  point_in_cube(const my_float_t *centroid, const my_float_t half_width, 
                const my_float_t *pt, const my_float_t radius)
  {
    my_float_t dist[3];
    for(int i = 0; i < 3; ++i){
      dist[i] = centroid[i] - pt[i];
      dist[i] = (dist[i] < 0 ? -1.0 * dist[i] : dist[i]);
    }
  
    // If center in box, return true 
    if(dist[0] <= half_width && dist[1] <= half_width && dist[2] <= half_width)
      return CENTER_IN_CUBE;

    if(dist[0] <= half_width + radius && dist[1] <= half_width + radius &&
       dist[2] <= half_width + radius) return POSSIBLE_INTERSECTION;

    // If center not in box + tol, return false
    return NO_INTERSECTION;
  }

  cube_sphere_intersect_t 
  cube_and_sphere_intersect(const my_float_t *centroid, 
                            const my_float_t half_width,
                            const my_float_t *center, const my_float_t radius)
  {
  
    my_float_t sign[3];
    sign[0] = sign[1] = sign[2] = 0;

#ifdef TRACE
    std::cout << "\nCentroid: " << centroid[0] << " " << centroid[1] << " "
              << centroid[2] << "\n";
    std::cout << "half_width: " << half_width << "\n";
    std::cout << "center: " << center[0] << " " << center[1] << " "
              << center[2] << "\n";
    std::cout << "radius: " << radius << "\n";
#endif

    for(int i = 0; i < 3; ++i){
      my_float_t tmp = centroid[i] - half_width;
      if(tmp - radius <= center[i] && center[i] <= tmp) sign[i] = -1;
  
      tmp = centroid[i] + half_width;
      if(tmp <= center[i] && center[i] <= tmp + radius) sign[i] = 1;
    }
#ifdef TRACE
    std::cout << "Sign: " << sign[0] << " " << sign[1] << " " 
              << sign[2] << "\n";
#endif

    /* Center is in the rectangular solid defined by a face and ?? and is 
     * closer than radius to that face? */
    if((sign[0] && !sign[1] && !sign[2]) || (!sign[0] && sign[1] && !sign[2]) ||
       (!sign[0] && !sign[1] && sign[2])){
#ifdef TRACE
      std::cout << "Center of sphere is adjacent to a face\n";
#endif
      return CUBE_AND_SPHERE_INTERSECT;
    }

    /* Sphere contains a cube corner? */
    my_float_t dist[3];
    if(sign[0] && sign[1] && sign[2]){
      for(int i = 0; i < 3; ++i)
        dist[i] = centroid[i] + sign[i] * half_width - center[i];
#ifdef TRACE
      std::cout << "dist[]: " << dist[0] << " " << dist[1] << " "
                << dist[2] << "\n";
#endif
      if(radius*radius >= dist[0]*dist[0] + dist[1]*dist[1] + dist[2]*dist[2]){
#ifdef TRACE
        std::cout << "Sphere contains a cube corner\n";
#endif
        return CUBE_AND_SPHERE_INTERSECT;
      }else{
#ifdef TRACE
        std::cout << "Sphere does not contain a cube corner\n";
#endif
        return NO_INTERSECTION;
      }
    }

    if(!sign[0]){
      if(line_in_sphere(centroid, half_width, center, radius, 0, sign))
        return CUBE_AND_SPHERE_INTERSECT;
    }else if(!sign[1]){
      if(line_in_sphere(centroid, half_width, center, radius, 1, sign))
        return CUBE_AND_SPHERE_INTERSECT;
    }else if(!sign[2]){
      if(line_in_sphere(centroid, half_width, center, radius, 2, sign))
        return CUBE_AND_SPHERE_INTERSECT;
    }
  
#ifdef TRACE
    std::cout << "No intersection!\n";
#endif
    return NO_INTERSECTION;
  }

  void
  print_bin_trace(const octree_node_t<atom_vci> *cnode, int level)
  {
    if(cnode->data.size()){
      std::cout << "\nBin at level " << level << "\n";
      std::cout << "\nCentroid: " << cnode->centroid[0] << " " 
                << cnode->centroid[1] << " " << cnode->centroid[2] << "\n";
      std::cout << "\nHalf Width: " <<  cnode->half_width << "\n";
      std::cout << "\nPositions:\n";
  
      for(int i = 0; i < cnode->data.size(); ++i)
        std::cout << "\t " << cnode->data[i]->pos[0] << " "
                  << cnode->data[i]->pos[1] << " " 
                  << cnode->data[i]->pos[2] << "\n";
      return;
    }
  
    for(int i = 0; i < 8; ++i)
      if(cnode->children[i]) print_bin_trace(cnode->children[i], level + 1);
  }

  bool
  line_in_sphere(const my_float_t *centroid, const my_float_t half_width,
                 const my_float_t *center, const my_float_t radius, 
                 const int idx, const my_float_t *sign)
  {
    my_float_t pt[3];
    my_float_t RHS = radius*radius;
    my_float_t S1, S2;
  
    for(int i = 0; i < 3; ++i)
      if(i != idx){
        pt[i] = centroid[i] + sign[i]*half_width;
        RHS -= (pt[i] - center[i]) * (pt[i] - center[i]);
      }
#ifdef TRACE
    std::cout << "pt: " << pt[0] << " " << pt[1] << " " << pt[2] << "\n";
#endif
  
    if(RHS < 0) return false;
    RHS = std::sqrt(RHS);
#ifdef TRACE
    std::cout << "RHS: " << RHS << "\n";
#endif
  
    S1 = center[idx] - RHS;
    S2 = center[idx] + RHS;
    if((S1 < centroid[idx] - half_width && S2 < centroid[idx] - half_width) ||
       (S1 > centroid[idx] + half_width && S2 > centroid[idx] + half_width))
      return false;
    return true;
  }
 
  octree_node_t<atom_vci> A_root;  //!< Root node of the tree
  uint A_max_numel; //!< Max number of elements in a cube at a level less than the max level
  my_float_t A_min_half_width;   //!< Smallest cube should be about this size

  static const std::string A_fname;  //!< Name of the source file
};

} // end namesapce SimSite3D
#endif
