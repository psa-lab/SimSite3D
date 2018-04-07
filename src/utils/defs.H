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
 * $Source: /psa/share/repository/pfizer_proj/src/basics/defs.H,v $
 * $Revision: 1.10 $
 * $Author: vanvoor4 $ 
 * $Date: 2008-02-26 17:28:50 $
 * 
 * $Log: not supported by cvs2svn $
 * Revision 1.9  2007/12/17 21:07:59  vanvoor4
 * Comments and adding spacing to account for 4 char PDB metal strings
 *
 * Revision 1.8  2007/10/11 16:03:36  vanvoor4
 * Removed orbitals info from here to orbitals.H
 *
 * Revision 1.7  2007/08/21 15:31:33  vanvoor4
 * Added some enums and moved the triangle defines into the Sitemap
 * classes.
 *
 * Revision 1.6  2007/02/07 16:03:33  vanvoor4
 * It is my opinion that defs should be included as static constants in the
 * class where the are used rather than as global defs in a header file.
 *
 * Revision 1.5  2007/01/24 14:37:35  vanvoor4
 * updated to use enums.
 *
 * Revision 1.4  2007/01/03 16:57:09  vanvoor4
 * Changed defs for MOL2 atom orbitals to enum.
 *
 * Revision 1.3  2006/11/15 20:08:55  vanvoor4
 * Added support from type to string, and hbond direction info
 *
 * Revision 1.2  2006/08/25 15:07:45  vanvoor4
 * Added CVS Header
 *
 *
 */

#ifndef _DEFS_H
#define _DEFS_H

#include <string>
#include <vector>
#include <orbitals.H>


/* when matching ligand interaction centers onto template points, the
   matching for the hydrophobic points has to be less exact, h-bond
   point matchings are weight 1.0 during the DME, the least-squares-fit,
   and RMSD computations, this value is the relative weight for the
   hydrophobic points */
#define  HPHOB_MATCH_WEIGHT              0.3

/* maximal distance between atom centers when checking for hydro-compl */
#define  HYDRO_DIST                      4.5

/* maximal distance to look for polar ligand atoms for interface target
   atoms */
#define  HPHIL_MATCH_DIST                5.0

/* the range of donor acceptor distances in which a H-bond can be formed */
#define  MIN_HBOND_LENGTH_HYDROGEN       1.5
#define  MIN_HBOND_LENGTH                2.5
#define  MAX_HBOND_LENGTH                3.5

typedef enum{
  NONE=0,
  TRIANGLE_MATCH,
  DME,
  RMS_DEVI,
  VEC_DOT_PROD
}errorType;

#define  SMALL_FLOAT  1.0e-5
#define  SMALL_DOUBLE 1.0e-10

typedef enum atom_level_enum{
  UNKOWN_LEVEL = 0,
  ALPHA,
  BETA,
  GAMMA,
  DELTA,
  EPSILON,
  ZETA,
  ETA,
  LEVEL_ENUM_END
}atom_level_type;

#endif
