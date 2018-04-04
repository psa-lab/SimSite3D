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
 * $Source: /psa/share/repository/pfizer_proj/src/gen_points/HphobPoints.H,v $
 * $Revision: 1.10 $
 * $Author: vanvoor4 $
 * $Date: 2008-05-15 17:37:07 $
 * 
 * $Log: not supported by cvs2svn $
 * Revision 1.9  2007/12/17 21:25:57  vanvoor4
 * Changed enum for default sphere level
 *
 * Revision 1.8  2007/11/01 16:26:17  vanvoor4
 * Renamed header file and additional sphere levels
 *
 * Revision 1.7  2007/09/26 14:49:47  vanvoor4
 * Added support for hphob 3.0 (A) shell method.
 *
 * Revision 1.6  2007/08/29 20:27:43  vanvoor4
 * REmoved a redundant size() fcn
 *
 * Revision 1.5  2007/08/21 18:23:47  vanvoor4
 * Changes required to support the point_storage scheme
 *
 * Revision 1.4  2007/03/06 19:39:48  vanvoor4
 * Added support for changing the cluster radius.
 *
 * Revision 1.3  2007/02/07 15:31:35  vanvoor4
 * Reorganized to allow the testing any number of methods to specify
 * hydrophobicity.  Also had do deal with the unfortunate conflict of
 * point types from here and ../basics.
 *
 * Revision 1.2  2006/11/16 20:24:52  vanvoor4
 * Added CVS header
 *
 *
 */

#ifndef _HPHOB_POINTS_HEADER_INCLUDED
#define _HPHOB_POINTS_HEADER_INCLUDED

#include <basics.H>
#include <HbondPoints.H>  // To check for overlap with hbonds
#include <hphob_point_t.H>
#include <GenPointsParameters.H>
#include <hphob_triads.H>

namespace SimSite3D{


  typedef std::list<hphob_point_t>			list_of_hphobs;
  typedef point_storage<hphob_point_t>			hphob_point_vec;
  typedef hphob_point_vec::iterator			hphob_point_vi;
  typedef hphob_point_vec::const_iterator		hphob_point_vci;
  typedef point_storage<hphob_point_t, list_of_hphobs>	hphob_point_list;
  typedef hphob_point_list::iterator			hphob_point_li;
  typedef hphob_point_list::const_iterator		hphob_point_lci;

  //! Generate/Load hydrophobic points into a list
  class HphobPoints{
  public:

    //! Default constructor to read hphob points from a sitemap file set
    HphobPoints(std::istream& in, const uint num_lines, PDBBase& atoms, 
                PDBBase& points);

    //! Default constructor to generate hphob points
    HphobPoints(const GenPointsParameters::hphob_method_t method_in, 
                const double cluster_diameter_in = 3.5,
                const sphere_sample_level_t sample_level_in = 
                DISCRETE_SPHERE_LEVEL_ZERO);

    void init_positions_octree(const my_float_t max_dist, 
                               const my_float_t min_side_len,
                               const uint max_numel);
    
    //! Get a const iterator to the first hphob point
    hphob_point_lci 
    begin() const
    { return points.begin(); }

    //! Get a const iterator to one past the last hphob point
    hphob_point_lci 
    end() const
    { return points.end(); }

    //! Transform the hphob points given R and T
    void 
    transform(const my_float_t* R, const my_float_t* T)
    { points.transform(R, T); }

    //! Transform the hphob points by the inverse transform of the transform 
    //! given by R and T
    void 
    inverse_transform(const my_float_t* R, const my_float_t* T)
    { points.inverse_transform(R, T); }

    //! Revert the hphob points as read in or as created?
    void 
    revert()
    { points.revert(); }

    hphob_point_lci 
    closest_point(const my_float_t* pos, my_float_t* d) const
    { return points.closest_point(pos, d); }

    //! Compute the root mean squared deviation (RMSD) between the current
    //! and orignial positions of the hydrophobic points
    inline my_float_t
    compute_RMSD() const
    { return points.compute_RMSD(); }

    //! Number of hphob points
    /*! A helper function since list iterators do not support addition and
     * subtraction
     */
    size_t 
    size() const
    { return points.size(); }

    bool gen_points(const interact_atoms_vec& hphob_atoms,
                    const hphob_triad_vec& hphob_triads,
                    const std::vector<atom_vci>& rad_atoms,
                    BoundingVolume* site_vol);

    void cull_too_close_to_polar(hbond_fit_pt_vci hbonds_begin,
                                 hbond_fit_pt_vci hbonds_end);

    static bool atom_is_hphob(const residue_type res, const atom_type atom,
                              const hphob_triad_t** triad_ptr);

    //! Get transform to move current points to original position of
    //! the points
    inline void
    get_current_inverse_3D_transform(Quaternion *Q, my_float_t *T) const
    {
      points.get_current_inverse_3D_transform(Q, T);
    }

  private:

    bool grid_point_methods(const std::vector<atom_vci>& rad_atoms,
                            BoundingVolume* site_vol);

    bool pseudo_surface_points(const interact_atoms_vec& hphob_atoms,
                               const hphob_triad_vec& hphob_triads,
                               const std::vector<atom_vci>& rad_atoms,
                               BoundingVolume* site_vol);

    bool get_hphobic_atoms(const std::vector<atom_vci>& rad_atoms);

    void add_point(const my_float_t* pos, hphob_point_vec* pts_vec);
    void add_point(const my_float_t* pos, hphob_point_list* pts_list);

    static void init_hphob_nbrs_tbl();
    
    GenPointsParameters::hphob_method_t method;
    sphere_sample_level_t sphere_sample_level;
    double cluster_diameter;
    hphob_point_list points;

    static hphob_triad_map_t hphob_nbrs_tbl;

    static const my_float_t MIN_HYDRO_DIST; //!< closest a template point can be to center of any protein atom
    static const my_float_t MAX_HYDRO_DIST; //!< farthest a template point can be from the center of any protein atom
    static const my_float_t HPHIL_CUTOFF;  
    static const std::string A_fname; //!< Name of the source file
  };
}

#endif
