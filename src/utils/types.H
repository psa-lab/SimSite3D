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
 * $Source: /psa/share/repository/pfizer_proj/src/basics/types.H,v $
 * $Revision: 1.13 $
 * $Author: vanvoor4 $ 
 * $Date: 2007-10-04 14:51:52 $
 * 
 * $Log: not supported by cvs2svn $
 * Revision 1.12  2007/09/20 17:36:33  vanvoor4
 * Removed underscore from define
 *
 * Revision 1.11  2007/08/21 16:42:09  vanvoor4
 * The interaction types got moved to gen_points an the point structure
 * got its own class.
 *
 * Revision 1.10  2007/02/07 16:02:19  vanvoor4
 * Added again as a quick fix solution.
 *
 * Revision 1.7  2006/11/17 21:36:50  vanvoor4
 * changed var for score_cutoff
 *
 * Revision 1.6  2006/11/16 20:20:31  vanvoor4
 * Changed the table buckets to vectors
 *
 * Revision 1.5  2006/11/15 20:09:26  vanvoor4
 * Added support for hbond direction info and my_float_t type.
 *
 * Revision 1.4  2006/08/25 15:08:02  vanvoor4
 * Added CVS Header
 *
 *
 */

#ifndef TYPES_H
#define TYPES_H

#include <limits>

typedef unsigned int uint;
typedef double my_float_t;
static const uint my_float_size = sizeof(my_float_t);
static const my_float_t my_float_max = std::numeric_limits<my_float_t>::max();
static const my_float_t my_float_min = std::numeric_limits<my_float_t>::min();
static const my_float_t my_float_eps = 
  std::numeric_limits<my_float_t>::epsilon();
static const my_float_t my_float_inf = 
  std::numeric_limits<my_float_t>::infinity();
static const my_float_t MIN_DOUBLE = 0.00000001;

#endif
