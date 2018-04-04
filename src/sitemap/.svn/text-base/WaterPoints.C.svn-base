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
 * $Source: /psa/share/repository/pfizer_proj/src/gen_points/WaterPoints.C,v $
 * $Revision: 1.1 $
 * $Author: vanvoor4 $
 * $Date: 2007-12-17 21:32:06 $
 *
 * $Log: not supported by cvs2svn $
 *
 *
 */

#include <WaterPoints.H>

using namespace SimSite3D;

const uint WaterPoints::num_pts = 20;
const my_float_t WaterPoints::relative_pos[] = 
  { 1.8618,  2.3524,  0.0000,
    0.9450,  2.8473,  0.0000,
    2.5541,  1.5737,  0.0000,
    1.7495,  2.2105,  1.0261,
    1.7495,  2.2105, -1.0261,
    1.8618, -2.3524,  0.0000,
    2.5541, -1.5737,  0.0000,
    0.9450, -2.8473,  0.0000,
    1.7495, -2.2105,  1.0261,
    1.7495, -2.2105, -1.0261,
   -1.7314,  0.0000,  2.4499,
   -1.6270, -1.0261,  2.3022,
   -1.6270,  1.0261,  2.3022,
   -0.7891,  0.0000,  2.8944,
   -2.4649,  0.0000,  1.7100,
   -1.7314,  0.0000, -2.4499,
   -1.6270, -1.0261, -2.3022,
   -1.6270,  1.0261, -2.3022,
   -2.4649,  0.0000, -1.7100,
   -0.7891,  0.0000, -2.8944
  };
/*
  { 3.0000,       0,       0,
    2.8191,  1.0261,       0,
    2.8191, -1.0261,       0,
    2.8191,       0,  1.0261,
    2.8191,       0, -1.0261,
   -0.6891, -2.9198,       0,
    0.3511, -2.9794,       0,
   -1.6462, -2.5080,       0,
   -0.6476, -2.7437,  1.0261,
   -0.6476, -2.7437, -1.0261,
   -1.0745,  1.3577,  2.4499,
   -1.8143,  0.6390,  2.3022,
   -0.2052,  1.9126,  2.3022,
   -0.4897,  0.6187,  2.8944,
   -1.5297,  1.9328,  1.7100,
   -1.0745,  1.3577, -2.4499,
   -1.8143,  0.6390, -2.3022,
   -0.2052,  1.9126, -2.3022,
   -1.5297,  1.9328, -1.7100,
   -0.4897,  0.6187, -2.8944
  };
*/

const my_float_t WaterPoints::relative_dir[] = 
  { 0.6206,  0.7841,  0.0000,
    0.3150,  0.9491,  0.0000,
    0.8514,  0.5246,  0.0000,
    0.5832,  0.7368,  0.3420,
    0.5832,  0.7368, -0.3420,
    0.6206, -0.7841,  0.0000,
    0.8514, -0.5246,  0.0000,
    0.3150, -0.9491,  0.0000,
    0.5832, -0.7368,  0.3420,
    0.5832, -0.7368, -0.3420,
   -0.5771,  0.0000,  0.8166,
   -0.5423, -0.3420,  0.7674,
   -0.5423,  0.3420,  0.7674,
   -0.2630,  0.0000,  0.9648,
   -0.8216,  0.0000,  0.5700,
   -0.5771,  0.0000, -0.8166,
   -0.5423, -0.3420, -0.7674,
   -0.5423,  0.3420, -0.7674,
   -0.8216,  0.0000, -0.5700,
   -0.2630,  0.0000, -0.9648
  };
/*
  { 1.0000,  0.0000,  0.0000,
    0.9397,  0.3420,  0.0000,
    0.9397, -0.3420,  0.0000,
    0.9397,  0.0000,  0.3420,
    0.9397,  0.0000, -0.3420,
   -0.2297, -0.9733,  0.0000,
    0.1170, -0.9931,  0.0000,
   -0.5487, -0.8360,  0.0000,
   -0.2159, -0.9146,  0.3420,
   -0.2159, -0.9146, -0.3420,
   -0.3582,  0.4526,  0.8166,
   -0.6048,  0.2130,  0.7674,
   -0.0684,  0.6375,  0.7674,
   -0.1632,  0.2062,  0.9648,
   -0.5099,  0.6443,  0.5700,
   -0.3582,  0.4526, -0.8166,
   -0.6048,  0.2130, -0.7674,
   -0.0684,  0.6375, -0.7674,
   -0.5099,  0.6443, -0.5700,
   -0.1632,  0.2062, -0.9648
  };
*/

WaterPoints::WaterPoints()
{
  // This is all hard coded for now -- we are looking for a rapid prototype
  // First 5 are for H1, 2nd 5 for H2, 3rd 5 for LP1, last 5 for LP2
  const my_float_t* pos = relative_pos;
  const my_float_t* dir = relative_dir;
  for(uint i = 0; i < num_pts; ++i, pos += 3, dir += 3){
    hbond_point_t tmp_pt;
    tmp_pt.act_type = (i < 10 ? ACCEPTOR : DONOR);
    std::copy(pos, pos + 3, tmp_pt.pos);
    std::copy(dir, dir + 3, tmp_pt.dir);
    water_pts.push_back(tmp_pt);
  }

  // Setup the "ideal" points
  pos = relative_pos;
  dir = relative_dir; 
  for(uint i = 0; i < 4; ++i, pos += 15){
    hbond_ideal_point_t tmp_pt;
    tmp_pt.pt_num = (i % 2 ? 2 : 1);
    tmp_pt.act_type = (i < 2 ? ACCEPTOR : DONOR);
    std::copy(pos, pos + 3, tmp_pt.pos);
    std::copy(dir, dir + 3, tmp_pt.dir);
    tmp_pt.fit_pts_beg = water_pts.begin() + 5*i;
    tmp_pt.fit_pts_end = water_pts.begin() + 5*(i+1);
    water_ideal_pts.push_back(tmp_pt);
  }
}
