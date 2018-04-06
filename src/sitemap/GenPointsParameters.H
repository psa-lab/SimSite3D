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

#ifndef GEN_POINTS_PARAMETERS_HEADER_INCLUDED
#define GEN_POINTS_PARAMETERS_HEADER_INCLUDED

#include <basics.H>
#include <BaseParameters.H>
#include <DiscreteSphere.H>

namespace SimSite3D{

  //! A data class (public members) used to hide the details, from the main
  //! application source, needed to get the parameters. 
  /*! Unlike the base class, the string parameters here are all cstrings to
   * support the popt library for parsing command line options
   */
  class GenPointsParameters : public BaseParameters{
  public:

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
      POPT_NULL_FLAG = 0,
      DISPLAY_VERSION,
      NO_NORMALIZATION,
      BASELINE_HPHOB_REP,
      PSEUDO_SURF_HPHOB_REP,
      MSMS_SURFACE,
      INCLUDE_WATERS, 
      INCLUDE_METALS,
      ALLOW_SMALL_SITEMAPS, 
      USE_BOUNDING_BOX_FOR_SITE_VOLUME
    }popt_flag_t;

    //! To be used by the SimSite3DPy search.parameters module
    GenPointsParameters();

    GenPointsParameters(const int argc, const char** argv);

    ~GenPointsParameters();

    //! Simple verification of command line parameters
    /*! Don't think this needs to be private
     * @return False implies unsafe to proceed, otherwise should be OK
     */
    status_t verify_params();

    std::string pts_fname;   //!< Path to sitemap points file
    std::string prot_fname;  //!< Path to the protein file used to build sitemap
    std::string lig_fname;   //!< Path to the ligand file used to build sitemap 
    std::string sphere_str;  //!< center and radius of sphere for bind vol
    std::string msms_binary;
    my_float_t cluster_diameter;
    my_float_t grid_spacing;
    hbond_method_t hbond_method; 
    hphob_method_t hphob_method;
    sphere_sample_level_t sphere_sample_level;
    bool normalize; //!< True imples normalize sitemap score against diverse set
    bool include_metals;  //!< True imples include metals in protein file as part of the protein, otherwise ignore all metal atoms
    std::vector<std::string> waters; //!< Water residues to include 
    verbose_level_t verbose_level;
    bool call_msms;  //!< True implies call $SIMSITE3D_INSTALL_DIR/bin/linux_msms to create the site's surface
    bool prune_to_lig; //!< True implies use the UnionOfBalls class rather than the RectangularSolid class for site map volume defined by ligand.
    my_float_t probe_radius; //!< Probe radius for MSMS
    my_float_t num_pts_per_area; //!< Number of points per (A)^2
    std::string user_provided_surf;  //!< Path to user provided msms_surf.vert file -- assumes a corresponding msms_surf.face file
    int min_res_per_chain; //!< Minimum number of residues per chain

  private:
    //! Initialize the class variables to standard nondefined values 
    //! (such as NULL).
    void init();
  
    //! Free up all C-style strings which are not NULL; calls init() at the end
    void free_cstrings();
 
    //! Given the command line arguments, parse them using the popt library
    /*! 
     * @param argc Number of command line arguments
     * @param argv Array of C-style strings holding the command line arguments
     * @return True if parsed successfully and may continue, else false
     */
    status_t get_opts(const int argc, const char** argv);
  
    char* prots_dir;       //!< Cmdline value for $SIMSITE3D_DBASE_PROTS
    char* ligs_dir;        //!< Cmdline value for $SIMSITE3D_DBASE_LIGS 
    char* hbond_dens_str;  //!< Cmdline value for hydrogen bond density
    char* A_scratch_dir;   //!< Cmdline value for $SIMSITE3D_SCRATCH_DIR
    char* A_proj_output;   //!< Cmdline value for $SIMSITE3D_PROJ_OUTPUT
    char* A_waters_str;    //!< Cmdline string for water residues to include
    static const std::string A_fname; //!< Name of source file
  };
}
#endif
