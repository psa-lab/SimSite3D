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

#ifndef _HBOND_POINTS_HEADER_INCLUDED
#define _HBOND_POINTS_HEADER_INCLUDED

#include <BoundingVolume.H>
#include <list>
#include <hbond_triads.H>
#include <PDBStructure.H>
#include <GenPointsParameters.H>
#include <interact_point_t.H>
#include <hbond_points.H>
#include <WaterPoints.H>

namespace SimSite3D{

  typedef struct {
    std::vector<my_float_t> alphas;
    std::vector<my_float_t> betas; 
    my_float_t ideal_alpha;
    my_float_t ideal_beta;
  }hbond_angles_t;
  typedef std::vector<hbond_angles_t> angles_vec;

  typedef struct {
    atom_type atom;
    const hbond_triad_t *triad;
    angles_vec *angles;
  }hbond_data_t;
  typedef std::map<residue_type, std::vector<hbond_data_t> > residue_table_t;

  //! Compute the hydrogen bonding template points based on the protein residues
  //! in the binding site.
  class HbondPoints{
  public:
    //! Build the hbonds_pts vector from a stored sitemap
    /*! This needs to be updated at some point so that the atoms and points
     * are correctly typed -- so that there is no possiblity of passing in the
     * wrong arguments
     */
    HbondPoints(std::istream& in, const uint num_lines,
                PDBBase& rad_atoms, PDBBase& points);

    //! Basic constructor
    /*! Sets a pointer to the bounding volume, loads the parameters for the
     * hbond geometry from parameter files and builds the parameters table
     * (keyed by the PDB residue name/type).
     *
     * @param density How coarse or finely to lay down template points
     * @param param_path Path to the SimSite3D parameters directory
     * @param vol_in Reference to the template bounding volume
     */
    HbondPoints(const GenPointsParameters::hbond_method_t hbond_density, 
                const std::string param_path, BoundingVolume* vol_in,
                const bool compute_volume = false);
  
    ~HbondPoints(){;}

    //! Given protein chains and the bounding volume compute the hbond
    //! template points
    /*! Get the rad and binding site atoms and compute the hbond positions
     * for the binding site residues.
     *
     * This method represents an initial attempt at getting water conserved
     * water molecules and required metals as part of the binding site.
     * If waters or metals are included, all waters or metals respectively
     * in the protein structure file and near the binding site will be 
     * considered as "part of the protein".
     *
     * Inclusion of water molecules requires correct protonation and 
     * orientation of all water molecules in the binding site.
     *
     * @param _begin Const iterator to the first protein chain
     * @param _end Const iterator to the last+1 protein chain
     * @param H2O_res_names The names of the water residues in the protein 
     *        structure to include as "part of the protein".
     * @param include_metals True implies include all metal atoms in the 
     *        the protein structure as "part of the protein".
     * @param bind_atoms Pointer to a vector to store the binding site 
     *        protein atoms
     */ 
    void find_binding_site_hbonds(const PDBStructure* prot,
                                  const std::vector<std::string>& H2O_res_names,
                                  const bool include_metals,
                                  const interact_atoms_vec& hbond_atoms,
                                  std::vector<atom_vci>* bind_atoms,
                                  std::vector<atom_vci>* rad_atoms);

    //! Get a constant iterator to the first hbond ideal point
    hbond_ideal_pt_vci 
    ideal_pts_beg() const 
    { return ideal_pts.begin(); } 

    //! Get a constant iterator to one past the last hbond ideal point
    hbond_ideal_pt_vci 
    ideal_pts_end() const 
    { return ideal_pts.end(); } 

    //! Get a constant iterator to the first hbond fit point
    hbond_fit_pt_vci 
    fit_pts_beg() const 
    { return fit_pts.begin(); } 

    //! Get a constant iterator to one past the last hbond fit point
    hbond_fit_pt_vci 
    fit_pts_end() const 
    { return fit_pts.end(); } 

    const size_t
    fit_pts_size() const
    { return fit_pts.size(); }

    //! Transform the positions of the hbond points and rotate their associated
    //! normal vectors (directions)
    void 
    transform(const my_float_t* R, const my_float_t* T)
    { 
      fit_pts.transform(R, T); 
      ideal_pts.transform(R, T);
    }

    //! Transform the positions of the hbond points and rotate their associated
    //! normal vectors (directions) by the inverse transform 
    void 
    inverse_transform(const my_float_t* R, const my_float_t* T)
    { 
      fit_pts.inverse_transform(R, T); 
      ideal_pts.inverse_transform(R, T);
    }

    //! Revert the positions of the hbond points and rotate their associated
    //! normal vectors (directions) to the positions and normals first stored.
    void 
    revert()
    { 
      fit_pts.revert(); 
      ideal_pts.revert(); 
    }

    //! Get transform to move current fit points to original position of
    //! the fit points
    inline void
    get_current_inverse_3D_transform(Quaternion *Q, my_float_t *T) const
    {
      fit_pts.get_current_inverse_3D_transform(Q, T);
    }

    //! Compute the root mean squared deviation (RMSD) between the current
    //! and orignial positions of the fit points
    inline my_float_t
    compute_RMSD() const
    { return fit_pts.compute_RMSD(); }

    hbond_fit_pt_vci
    closest_fit_pt(const my_float_t* pos, my_float_t* d) const
    { return fit_pts.closest_point(pos, d); }

    hbond_ideal_pt_vci
    closest_ideal_pt(const my_float_t* pos, my_float_t* d) const
    { return ideal_pts.closest_point(pos, d); }

    void
    close_fit_pts(const my_float_t *pos, const my_float_t radius,
                  hbond_fit_pt_vec::float_const_iter_map *pts_map) const
    { fit_pts.close_points(pos, radius, pts_map); }

    bool
    fail() 
    { return A_fail; }

    static const atom_vci NO_CARBON_NEIGHBOR;
  private:

    //! Get points for the specified waters near the binding site
    /*! Check if water atoms are near the binding site based on the location of
     * the oxygen atoms of H2O -- assume that if a chainID is specified with
     * the water residue #, the chainID always preceeds the residue #.  
     * Similarly if an insertion code (iCode) is specified with a water 
     * residue #, the iCode always follows the residue #
     */
    void get_nearby_waters(const PDBStructure* prot,
                           const std::vector<std::string>& H2O_res,
                           std::vector<atom_vci>* nearby_waters,
                           std::vector<atom_vci>* bind_atoms,
                           std::vector<atom_vci>* rad_atoms);

    void init_iterators(hbond_ideal_pt_vec* ideal_pts_p,
                        hbond_fit_pt_vec* fit_pts_p,
                        std::vector<int> num_fit_pts);
  
    //! Set down the potential hbonding template points for the given atom
    /*! Using the given atom and residue, set down the hbond template points
     * for that atom and remove those template points with heavy overlap of
     * of the protein
     *
     * @param chain Const iterator to the current chain
     * @param residue Const iterator to the current residue
     * @param B_atom Const iterator to the N or O atom requiring template 
     *        points
     * @param rad_atoms Const vector specifying the rad atoms
     * @param ideal_pts_p Pointer to the vector holding the ideal hbond points
     * @param fit_pts_p Pointer to the vector holding the hbond points for 
     *        fitting
     * @return False if we could not lay down template points for the given 
     *   N or O (not including mainchain Ns and Os)
     */
    bool get_atom_hbonds(const chain_const_iter chain, 
                         const residue_vci residue,
                         const atom_vci B_atom, 
                         const std::vector<atom_vci>& rad_atoms,
                         hbond_ideal_pt_vec* ideal_pts_p,
                         hbond_fit_pt_vec* fit_pts_p,
                         std::vector<int>* num_fit_pts);

    //! The positions of the hbonding template points are computed based on 
    //! the hbonding atom and two of its neighbors.
    /*! Given the hbonding atom and 2 of its neighbors, compute the template
     * points by transforming the angles and distances from the model
     * to the structure file coordinates of the given atoms.
     * 
     * @param atom The N or O atom requiring hbond template points
     * @param hbond_info Struct with the atom name, favored angles, 
     *    interaction distance and template point type
     * @param positions Array holding the positions of neighbor atom A, 
     *        current N or O and neighbor atom C
     * @param rad_atoms Vector denoting the rad atoms
     * @param ideal_pts_p Pointer to the vector holding the ideal hbond points
     * @param fit_pts_p Pointer to the vector holding the hbond points for 
     *        fitting
     * @return Currently returns true 
     */
    bool calc_positions(const atom_vci atom, const atom_vci carbon_nbr,
                        const hbond_data_t& hbond_info, 
                        const my_float_t **positions,
                        const std::vector<atom_vci>& rad_atoms,
                        hbond_ideal_pt_vec* ideal_pts_p,
                        hbond_fit_pt_vec* fit_pts_p,
                        std::vector<int>* num_fit_pts);

    bool calc_volume(const atom_vci atom, const atom_vci carbon_nbr,
                     const hbond_data_t& hbond_info, 
                     const my_float_t **positions,
                     const std::vector<atom_vci>& rad_atoms);

    void compute_offset(const my_float_t bond_len, const my_float_t alpha,
                        const my_float_t beta, const atom_vci atom,
                        my_float_t *offset);

    //! Add template points for "bound" water molecules that have the potential
    //! to interact with molecules in the sitemap volume
    /*! 
     * @param water_O Iterator to the water oxygen atom that requires template 
     *        points
     * @param rad_atoms Vector holding the rad atoms
     * @param water_hbond_pts Class to get the water hbond points
     * @param ideal_pts_p Pointer to the vector holding the ideal hbond points
     * @param fit_pts_p Pointer to the vector holding the hbond points for 
     *        fitting
     * @param num_fit_pts Pointer to the vector holding the number of fit points
     *        for a given object (metal atom, lone pair or hydrogen atom).
     */
    void add_water_points(const atom_vci water_O, 
                          const std::vector<atom_vci>& rad_atoms,
                          WaterPoints* water_hbond_pts,
                          hbond_ideal_pt_vec* ideal_pts_p,
                          hbond_fit_pt_vec* fit_pts_p,
                          std::vector<int>* num_fit_pts);

    //! Add template points for "included" metal atoms that have the potential
    //! to interact with molecules in the sitemap volume
    /*! This method is similar to that of SLIDE with a couple of exceptions:
     * 1) To stay consistent with the SimSite3D handling of points, no template
     *    points can be closer than 2.5 (A) to the center of a metal atom 
     *    (except for the template points arising from a given metal atom).
     * 2) The distribution of points on the sphere is somewhat different and
     *    arises from a systematic discretization of a sphere
     *
     * @param metal Iterator to the metal atom that needs template points
     * @param rad_atoms Vector holding the rad atoms
     * @param ideal_pts_p Pointer to the vector holding the ideal hbond points
     * @param fit_pts_p Pointer to the vector holding the hbond points for 
     *        fitting
     * @param num_fit_pts Pointer to the vector holding the number of fit points
     *        for a given object (metal atom, lone pair or hydrogen atom).
     */
    void add_metal_points(const atom_vci metal, 
                          const std::vector<atom_vci>& rad_atoms,
                          hbond_ideal_pt_vec* ideal_pts_p,
                          hbond_fit_pt_vec* fit_pts_p,
                          std::vector<int>* num_fit_pts);

    //! Wrapper function to setup temporaries and transform a given point
    /*! This functions main purpose is to reduce the number of lines in 
     * calc_positions
     * 
     * @param bond_len Distance to put template points from interacting heavy 
     *        atom's center
     * @param alpha In-plane angle
     * @param beta Out-of-plane angle
     * @param R_inv The rotation matrix used to rotate the point -- named R 
     *        inverse because of the construction of this matrix was by first
     *        finding R (in calc_positions).
     * @param pt The point to rotate
     */
    void transform_point(const my_float_t bond_len, const my_float_t alpha,
                         const my_float_t beta, const my_float_t* R_inv,
                         my_float_t *pt);


    //! Given the hbond triad get the postions of the atoms nbr_A and nbr_C.
    /*!
     * @param residue Const iterator to the current residue
     * @param triad Const reference to the current atom triplet
     * @param positions Pointer to array to store the positions for nbr_A and
     *        nbr_C
     * @return True if both neighbors were found (in the given residue), else
     *         return false.
     */
    bool sidechain_nhbr_pos(residue_vci residue, 
                            const hbond_triad_t& triad,
                            const my_float_t **positions, atom_vci* carbon_nbr);
  
    typedef bool (*add_angles_t)(std::map<residue_type, angles_vec>* angles,
                                 const std::vector<std::string>& tokens);

    //! Read in the ideal and favored hbond angles for the binding residues
    /*!
     * @param hbond_density Sampling density for hbond angles
     * @param param_path Path to the SimSite3D parameters directory
     * @return True if both files were found and had no parsing errors, 
     *         otherwise false
     */
    bool load_hbond_angles(GenPointsParameters::hbond_method_t hbond_density, 
                           std::string param_path);

    static bool add_ideal_angles(std::map<residue_type, angles_vec>* angles,
                                 const std::vector<std::string>& tokens);

    static bool add_favored_angles(std::map<residue_type, angles_vec>* angles,
                                   const std::vector<std::string>& tokens);

    //! Read in a file which holds the possible angles for favored hbonds
    //! based on the interacting atom and 2 of its neighbors
    /*! Two types of angle files should be read.  One is the ideal angles 
     * file that is expected to have only the angles for which the DHA and
     * HAP (hydrogen-acceptor-preacceptor) angles are optimum (typically
     * angle DHA = Pi).  Only one angle pair per H or lone pair is expected -- 
     * specifying more than one angle pair per H or lone pair will yield 
     * undetermined results.
     *
     * The second angles file should be one which holds one or more angle pairs
     * for a given lone pair or H.  The angle pairs represent a discrete 
     * sampling of the favorable angles allowed/used by a given lone pair
     * or H.
     * @param fname Path to the hbond angles parameter file
     * @param add_angles Pointer to the function to use to load the angles.
     *        The current options are add_ideal_angles() or 
     *        add_favored_angles().  
     * @return True if the file exists and no parsing errors occured, otherwise
     *         false
     */
    bool read_angles_file(std::string fname, add_angles_t add_angles);

    //! Build a table keyed by residue type that holds the residues and 
    //! corresponding atoms which have been parameterized to make hbonds
    void build_residue_table();

    //! Add hbond fit point to the fit points vector if it is inside the 
    //! binding site volume and does not have significant overlap with the 
    //! protein
    /*!
     * @param atom Iterator to the atom that could make the interaction
     * @param pt_pos Position of the template point
     * @param pt_act_type Interaction type of the template point
     * @param rad_atoms Atoms near the binding site
     * @param fit_pts_p Pointer to the vector holding the hbond fit points
     * @return True if point was added to the fit points vector, else false
     */
    bool add_fit_point(const atom_vci atom, const my_float_t* pt_pos, 
                       const interactionType pt_act_type, 
                       const std::vector<atom_vci>& rad_atoms, 
                       hbond_fit_pt_vec* fit_pts_p);

    //! Add hbond ideal point to the ideal points vector
    /*!
     * @param atom Iterator to the atom that could make the interaction
     * @param pt_pos Position of the template point
     * @param pt_act_type Interaction type of the template point
     * @param pt_num Ideal point number
     * @param num_pts Number of fit points corresponding to the ideal point
     * @param ideal_pts_p Pointer to the vector holding the hbond fit points
     * @param num_fit_pts Vector holding the number of fit points for each
     *        ideal point
     */
    void add_ideal_point(const atom_vci atom, const atom_vci carbon_nbr,
                         const my_float_t* pt_pos, 
                         const interactionType pt_act_type, 
                         const uint pt_num, const uint num_pts,
                         hbond_ideal_pt_vec* ideal_pts_p,
                         std::vector<int>* num_fit_pts);

    bool load_carbon_nbrs(const std::string rad_fname);
  
    BoundingVolume* bounding_volume;//!< Pointer to the template bounding volume
    hbond_fit_pt_vec fit_pts;
    hbond_ideal_pt_vec ideal_pts;
    std::map<residue_type, angles_vec> A_O_angles; //!< O hbond angles table -- keyed by PDB residue
    std::map<residue_type, angles_vec> A_N_angles; //!< N hbond angles table -- keyed by PDB residue
    residue_table_t hbond_residue_table; //!< Table holding the Hbond geometry -- keyed by PDB residue -- for a given residue and atom

    bool A_fail;
    bool A_compute_volume;
    my_float_array A_hbond_vol;

    static const my_float_t MAXRADDIST;  //!< So called rad file distance
    static const my_float_t MAXBINDDIST; //!< So called bind file distance
    static const my_float_t MIN_VDW_DIST;//!< Min distance for template pt to prot atom 
    static const my_float_t MIN_VDW_DIST_2; //!< Square of min vdw dist
    static const my_float_t MIN_VDW_DIST_METAL ;//!< Min distance for template pt to metal atom 
    static const my_float_t MIN_VDW_DIST_METAL_2; //!< Square of min vdw dist to metals
    static const std::string _fname;     //!< "HbondPoints.C"

  };
}

#endif
