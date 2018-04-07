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
 * $Source: /psa/share/repository/pfizer_proj/src/basics/my_float_array.C,v $
 * $Revision: 1.1 $
 * $Author: vanvoor4 $
 * $Date: 2007-09-26 14:40:46 $
 * 
 * $Log: not supported by cvs2svn $
 *
 *
 */

#include <my_float_array.H>

using namespace SimSite3D;

my_float_array::my_float_array(uint init_num_pts, uint stride_in)
{
  _stride = stride_in;
  _reserved = _stride * init_num_pts;
  curr = 0;
  curr_end = 0;
  orig = 0;
  orig_end = 0;
  scratch = 0;

  if(init_num_pts > 0){
    curr = new my_float_t[_reserved];
    curr_end = curr;
  }
}

my_float_array::~my_float_array()
{
  _reserved = 0;
  if(orig) delete [] orig;
  orig = orig_end = 0;
  if(curr) delete [] curr;
  curr = curr_end = 0;
  if(scratch) delete [] scratch;
  scratch = 0;
}

const bool 
my_float_array::centroid_3D(my_float_t* C) const
{
  if(_stride != 3) return false;
  int p_max = curr_end - curr;

  std::fill(C, C + 3, 0);
  my_float_t* pt = curr;
  for(int p = 0; p < p_max; p += 3, pt += 3)
    for(uint i = 0; i < 3; ++i) C[i] += pt[i];
  uint n = size();
  for(uint i = 0; i < 3; ++i) C[i] /= n;
  return true;
}

void 
my_float_array::double_reserve()
{
  // Grab new mem, copy data from old to new and delete old
  int new_reserve = (_reserved > 0 ? 2*_reserved : 2*_stride);
  my_float_t* new_mem = new my_float_t[new_reserve];
  std::copy(curr, curr_end, new_mem);
  delete [] curr;
  curr = new_mem;
  new_mem = 0;
  curr_end = curr + _reserved;
  _reserved = new_reserve;
} 

void 
my_float_array::update_orig_pts()
{
  int n = curr_end - curr;
  int orig_n = orig_end - orig;
  if(orig_n >= n) return;
                                                                            
  my_float_t* new_mem = new my_float_t[n];
  if(orig){
    std::copy(orig, orig_end, new_mem);
    delete [] orig;
    orig = 0;
  }
  std::copy(curr + orig_n, curr_end, new_mem + orig_n);
  orig = new_mem;
  orig_end = new_mem + n;

  if(scratch) delete [] scratch;
  scratch = new my_float_t[n];
}

void
my_float_array::do_copy(const my_float_array& src)
{    
  _stride = src._stride;
  _reserved = src._reserved;
  curr = new my_float_t[_reserved];
  std::copy(src.curr, src.curr_end, curr);
  int N = src.curr_end - src.curr;
  curr_end = curr + N;
  if(src.orig && src.orig_end){
    orig = new my_float_t[N];
    std::copy(src.orig, src.orig_end, orig);
    orig_end = orig + N;
  }
  scratch = new my_float_t[N];
}

void
my_float_array::init()
{
  _stride = 3;
  _reserved = 0;
  orig = 0;
  orig_end = 0;
  curr = 0;
  curr_end = 0;
  scratch = 0;
}
