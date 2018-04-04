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
 * $Source: /psa/share/repository/pfizer_proj/src/gen_points/SitemapDefs.H,v $
 * $Revision: 1.3 $
 * $Author: vanvoor4 $
 * $Date: 2007-09-24 15:46:30 $
 * 
 * $Log: not supported by cvs2svn $
 * Revision 1.2  2007/08/21 18:42:16  vanvoor4
 * Added the working directory
 *
 * Revision 1.1  2007/08/09 17:39:31  vanvoor4
 * Initial checkin of enums and structs for the sitemaps dir
 *
 *
 * 
 * 
 */

#ifndef _SITEMAP_DEFS_HEADER_INCLUDED
#define _SITEMAP_DEFS_HEADER_INCLUDED

#include <basics.H>

namespace SimSite3D{

  //! Density at which to place hbonds
  typedef enum{
    HBOND_DENSITY_NOT_SET = 0,
    OPTIMUM_HBONDS,
    MIN_HBONDS,
    SPARSE_HBONDS,
    DENSE_HBONDS
  }hbond_method_t;

  typedef enum{
    HPHOB_METHOD_NOT_SET,
    THREE_PROTEIN_ATOMS, //!< If 3+ atoms are close ==> hphob pt.
    THREE_MORE_HPHOB, //!< If 3 more hphob than hphil are close 
    ATOM_CENTERS, //!< Use centers of Cs and Ses if they have no polar neighbor
    PSEUDO_SURFACE,
    BURIED_SURFACE_AREA  //!< not implemented yet.
  }hphob_method_t;

  typedef enum{
    SPHERE_SAMPLE_LEVEL_ZERO,
    SPHERE_SAMPLE_LEVEL_ONE,
    SPHERE_SAMPLE_LEVEL_TWO,
    SPHERE_SAMPLE_LEVEL_THREE,
  }sphere_sample_level_t;

  //! Container for args from command line, parameter files and env variables
  typedef struct{
    char* pts_fname;       //!< File name to store the sitemap
    char* protein_fname;   //!< Protein file to use to generate the sitemap
    char* ligand_fname;    //!< Ligand file to use to define sitemap volume
    char* ligs_dir;        //!< Default ligands directory
    char* proteins_dir;    //!< Default proteins directory
    char* work_dir;        //!< Scratch directory
    char* hbond_dens_str;  //!< string of the desired hbond density
    char* sphere_str;      //!< string of the center and radius sphere
    my_float_t cluster_diameter;  //!< Maximal cluster diameter for complete link
    my_float_t grid_spacing;   //!< The length of a side of a grid cube.
//    char* hphob_method_str;//!< Cstring holding the desired hphob method
    hbond_method_t hbond_method;  //!< Denotes how densely to represent hbonds
    hphob_method_t hphob_method;  //!< Denotes which hueristical method to use
    sphere_sample_level_t sphere_sample_level;
    int sphere_sample_level_int;
    bool normalize;        //!< True implies get mean and variance for sitemap
    verbose_level_t verbose;  //!< How much to write to stdout
  }sitemap_variables_t;

}

#endif
