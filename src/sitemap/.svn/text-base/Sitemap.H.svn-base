/******************************************************************************
 * Copyright (c) 2006,2007, Michigan State University (MSU) Board of Trustees.
 *   All rights reserved.
 *
 * This file is part of the ASCbase Software project.
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

#ifndef SITEMAP_HEADER_INCLUDED
#define SITEMAP_HEADER_INCLUDED

#include <SitemapAtomsFile.H>
#include <SitemapPointsFile.H>
#include <GenPointsParameters.H>
#include <interact_point_t.H>
#include <HbondSurfaces.H>
#include <SimpleTrimesh.H>

namespace ASCbase{

  typedef dir_point_storage<interact_point_t>		interact_pts_vec;
  typedef interact_pts_vec::const_iterator		interact_pts_vci;

  //! Both generates and loads sitemaps.  Kinda glues together the different
  //! parts in an ad hoc manner.  This works, but the interface is quite 
  //! clunky.
  class Sitemap{
  public:

    //! Read/load a sitemap
    Sitemap(const std::string path, const std::string struct_id, 
            const BaseParameters& args, const bool normalize = true,
            const bool load_hbond_surfaces = false,
            const bool hydrophobic_query = false,
            const verbose_level_t verbosity = VERBOSE_SILENT);

    //! First step in creating a sitemap
    Sitemap(const GenPointsParameters& usr_args);

    //! Typical cleanup of dynamic mem
    virtual ~Sitemap();

    //! Check if class is viable
    /*!
     * @return False if no fatal errors have occured, else true
     */
    bool fail() const { return A_fail; }
 
    //! Write the three sitemap files
    bool write_files(const std::string struct_id, const std::string path);

    //! Determine if the site map has enough points given the current
    //! implementation of the site maps, alignment method, and alignment 
    //! techniques scoring
    bool has_enough_points();

    //! Get the filename of the PDB structure used to generate the sitemap
    std::string
    atoms_file_name() const
    { return A_atoms_fname; }

    //! Get the filename of the mol2 ligand used to generate the sitemap
    std::string
    ligand_file_name() const
    { 

      // Check if a ligand was specified at creation time
      if(A_lig_fname == "NONE") return "NONE";

      // Check A_site_path first, since the ligand could be in the same
      // directory as the site map
      std::string my_lig = A_site_path + "/" + A_lig_fname;
      if(!normal_file_exists(my_lig, false)){
        my_lig = A_dbase_ligs + "/" + A_lig_fname;

        if(!normal_file_exists(my_lig, false)){
          std::cerr << "Unable to find the file: " << my_lig << "\n";
          return "NOT_FOUND";
        }
      }
      return my_lig;
    }   

    const SitemapAtomsFile&
    interacting_atoms() const
    { return *A_site_atoms; } 

    SitemapAtomsFile&
    interacting_atoms()
    { return *A_site_atoms; } 

    //const SitemapAtomsFile&
    const PDBStructure & 
    bind_site_atoms() const
    { return *A_bind_site_res; }

    //SitemapAtomsFile&
    PDBStructure&
    bind_site_atoms()
    { return *A_bind_site_res; }

    const SitemapPointsFile&
    sitemap_points() const
    { return *points; } 

    interact_pts_vci 
    interact_pts_beg() const
    { return interact_points.begin(); }

    interact_pts_vci 
    interact_pts_end() const
    { return interact_points.end(); }

    void
    transform_interact_pts(my_float_t* R, my_float_t* T)
    { interact_points.transform(R,T); }

    void
    revert_interact_pts()
    { interact_points.revert(); }

    void
    interact_pts_centroid_3D(my_float_t* C)
    { interact_points.centroid_3D(C); }

    //! Get a reference to the sitemap hbond points
    HbondPoints&
    hbond_points()
    { return *hbond_pts; }

    //! Get a const reference to the sitemap hbond points
    const HbondPoints&
    hbond_points() const
    { return *hbond_pts; }

    //! Get a reference to the sitemap hphob points
    HphobPoints&
    hphob_points()
    { return *hphob_pts; }

    //! Get a const reference to the sitemap hphob points
    const HphobPoints&
    hphob_points() const
    { return *hphob_pts; }

    size_t
    fit_points_size() const
    { 
      return (hbond_pts->fit_pts_end() - hbond_pts->fit_pts_beg()) 
        + hphob_pts->size();
    }

    const my_float_t*
    get_norm_moments() const
    { return A_norm_moments; }

    void
    transform_bind_site_res_pos(const my_float_t* R, const my_float_t* T)
    { A_bind_site_res->transform(R, T); }

    void
    revert_bind_site_res_pos()
    { A_bind_site_res->revert(); }

    const BoundingVolume&
    bounding_volume_handle() const
    { return *A_site_vol; }

    const BoundingVolume&
    site_volume_estimate_handle() const
    { return *A_site_vol_est; }

    inline void
    revert()
    {
      hbond_pts->revert();
      hphob_pts->revert();
      A_bind_site_res->revert();
    }

    //! Compute the root mean squared deviation (RMSD) between the current
    //! and orignial positions of the hydrogen bonding and hydrophobic points
    inline my_float_t
    compute_RMSD() const
    { 
      // This function was built for ease of use and not speed
      my_float_t hbond_sum = 0.0;
      if(hbond_pts->fit_pts_size() > 0) hbond_sum = hbond_pts->compute_RMSD();
      my_float_t hphob_sum = 0.0;
      if(hphob_pts->size() > 0) hphob_sum = hphob_pts->compute_RMSD();
      hbond_sum *= hbond_sum * hbond_pts->fit_pts_size();
      hphob_sum *= hphob_sum * hphob_pts->size();
      size_t N = hbond_pts->fit_pts_size() + hphob_pts->size();
      if(N == 0) return my_float_max;
      return std::sqrt((hbond_sum + hphob_sum) / N);
    }

    //! Transform the objects that comprise the model of the binding site
    inline void
    transform(const my_float_t* R, const my_float_t* T)
    {
      hbond_pts->transform(R, T);
      hphob_pts->transform(R, T);
      A_bind_site_res->transform(R, T);
    }

    //! Transform the objects that comprise the model of the binding site using
    //! the inverse of the transform defined by R & T
    inline void
    inverse_transform(const my_float_t* R, const my_float_t* T)
    {
      hbond_pts->inverse_transform(R, T);
      hphob_pts->inverse_transform(R, T);
      A_bind_site_res->inverse_transform(R, T);
    }

    HbondSurfaces<hbond_surface_t>&
    hbond_surfaces() const
    { return *A_hbond_surfaces; }

    std::string
    site_path() const
    { return A_site_path; }

    std::string
    get_struct_id() const
    { return A_struct_id; }


    //! Compute the maximum value for each feature and store in the given 
    //! vector
    void
    compute_max_feature_values(std::vector<my_float_t> *max_feat_vals) const;

    //! Compute the fraction of the volume of this site's volume estimate
    //! that is inside the other site's volume estimate
    my_float_t
    compute_site_vol_intersect(const Sitemap& other) const;

  protected:
    void set_fail_flag(const bool val) { A_fail = val; }

    bool load_volume(std::string path, std::string struct_id, 
                     std::string dbase_sites, std::string vol_tag,
                     BoundingVolume **vol, geometry::sphere_t *S = 0);

    HbondPoints *hbond_pts;      //!< Hydrophillic portion of sitemap
    HphobPoints *hphob_pts;      //!< Hydrophobic portion of sitemap
    interact_pts_vec interact_points;
    BoundingVolume *A_site_vol_est;   //!< Upper estimate of the site volume

  private:
    //! Initialize class pointer variables to zero (NULL)
    void init();

    bool
    get_user_surf_patches(const std::string &lig_fname,
                          const std::string &scratch_dir,
                          const std::string &sphere_str,
                          const std::string &user_vert_fname);

    bool 
    get_MSMS_surf_patches(const std::string &lig_fname,
                          const std::string &scratch_dir,
                          const std::string &sphere_str,
                          const std::string &install_dir,
                          const bool include_metals, 
                          const std::string &msms_binary,
                          const std::vector<std::string> &waters);

    //! Get the normalization stats for the given sitemap file
    bool get_norm_stats(const std::string csv_fname, 
                        const BaseParameters& args,
                        bool hydrophobic_query = false);

    //! Setup the sitemap volume 
    /*!
     * @param ligand_fname Path to the ligand to define the binding site volume
     * @param grid_spacing The rectuangular grid sampling distance
     * @param sphere_str Sphere as center and radius in a string
     * @param use_union_of_balls If a ligand is given and use_union_of_balls
     *   is true, take the union of balls (1 centered at each ligand atom) as
     *   the volume estimate and not the ligand bounding box.
     * @return Pointer to the volume object or NULL if an error occured
     */
    BoundingVolume* setup_volume(const std::string ligand_fname,
                                 const my_float_t grid_spacing,
                                 const std::string sphere_str,
                                 const bool use_union_of_balls);

    void get_prot_interact_atoms(const PDBStructure* prot,
                                 BoundingVolume* site_vol,
                                 std::vector<atom_vci>* bind_atoms,
                                 std::vector<atom_vci>* rad_atoms, 
                                 interact_atoms_vec* hbond_atoms,
                                 interact_atoms_vec* hphob_atoms,
                                 hphob_triad_vec* hphob_triads,
                                 const int min_chain_size);

    //! Calculate the hydrogen bonding template points for the protein in the 
    //! binding volume
    /*!
     * @param hbond_density Density to discretely sample hydrogen bonding volume
     * @param params_dir Directory holding the hbond point sampling files
     * @param site_vol Pointer to the sitemap volume
     * @param include_waters True implies include all water molecules in 
     *        the protein structure as "part of the protein".
     * @param include_metals True implies include all metal atoms in the 
     *        the protein structure as "part of the protein".
     * @param rad_atoms Pointer to vector to hold rad atoms.
     */
    void
    find_binding_site_hbonds(const GenPointsParameters::hbond_method_t hbond_density,
                             const std::string params_dir,
                             const std::vector<std::string>& waters,
                             BoundingVolume* site_vol,
                             const interact_atoms_vec& hbond_atoms,
                             std::vector<atom_vci>* bind_atoms,
                             std::vector<atom_vci>* rad_atoms);

    //! Calculate the template points for all metal atoms in the binding site
    bool find_binding_site_metals();

    //! Calculate the hydrophobic bonding template points for the protein in
    //! the binding volume
    /*!
     * @param hphob_method which hydrophobic model/method to use
     * @param cluster_diameter The maximum cluster diameter for complete
     *        link clustering of hphob points.
     * @param rad_atoms Reference to the rad atoms (returned from hbonds method)
     * @return !fail()
     */
    bool find_binding_site_hphobs(const GenPointsParameters::hphob_method_t hphob_method,
                                  const my_float_t cluster_diameter,
                                  const interact_atoms_vec& hphob_atoms,
                                  const hphob_triad_vec& hphob_triads,
                                  const std::vector<atom_vci>& bind_atoms, 
                                  BoundingVolume* site_vol,
                                  const sphere_sample_level_t sample_level =
                                  DISCRETE_SPHERE_LEVEL_TWO);

    //! Write out the largest value possible for each feature 
    void write_feature_scales(std::ostream &out);

    void clean_up_string(std::string& S);

    void build_points_header(const GenPointsParameters& args);

    bool gen_lig_box(const std::string fname, my_float_t* corners);

    void get_site_volume_estimate(hbond_fit_pt_vci hbond_pts_beg,
                                  hbond_fit_pt_vci hbond_pts_end, 
                                  hphob_point_lci hphob_points_beg,
                                  hphob_point_lci hphob_points_end, size_t npts,
                                  my_float_t vol_est_tol);

    bool load_ligand(const std::string &fname, CoordFile** lig);
  
    bool
    get_prot_atom_close_to_ligand(const std::string ligand_fname, 
                                  std::vector<atom_vci>& bind_atoms, 
                                  atom_vci* close_atom);


    verbose_level_t verbose;
    BoundingVolume *A_site_vol;
    SitemapAtomsFile *A_site_atoms; //!< Atoms from residues near the binding site
    //SitemapAtomsFile *A_bind_site_res;
    PDBStructure *A_bind_site_res;
    SitemapPointsFile *points;
    PDBStructure *prot_atoms;
    atom_map_t interact_atoms;
    my_float_t A_norm_moments[2];
    my_float_t A_hphob_norm_moments[2];
    HbondSurfaces<hbond_surface_t>* A_hbond_surfaces;   
    SimpleTrimesh *old_binding_site_mesh;

    geometry::sphere_t A_site_vol_est_as_sphere; //!< I don't have the time to compute the intersection of the different bounding volume classes -- use this instead

    std::string A_site_path;     //!< Directory holding the sitemap
    std::string A_struct_id;     //!< Structure ID of the site map
    std::string A_dbase_ligs;    //!< Directory holding dbase ligands
    std::string A_lig_fname;     //!< Name of the ligand file
    std::string A_atoms_fname;     //!< Name of the sitemap atoms file
    bool A_include_metals;
    bool A_MSMS_surf;
    bool A_user_surf;
    my_float_t A_probe_radius; //!< Probe radius for MSMS
    my_float_t A_num_pts_per_area; //!< Number of points per (A)^2
    int A_min_res_per_chain;

    bool A_fail;
    static const int MIN_NUMBER_POINTS__POLAR_SITE;
    static const int MIN_NUMBER_POINTS__HPHOB_SITE;
    static const my_float_t TOLERANCE;
    static const my_float_t MAXRADDIST;  //!< Max rad file dist?
    static const my_float_t MAXBINDDIST;  //!< Max interaction with binding site dist
    static const std::string A_fname; //!< Name of source file
  };

}

#endif
