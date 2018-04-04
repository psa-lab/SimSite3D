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

/*  
 * $Source: /psa/share/repository/pfizer_proj/src/search/ModelSitemap.H,v $
 * $Revision: 1.4 $
 * $Author: vanvoor4 $
 * $Date: 2009-01-12 21:22:53 $
 * 
 * $Log: not supported by cvs2svn $
 * Revision 1.3  2007/11/01 16:40:27  vanvoor4
 * Changed due to the need to possibly pass the parameters class
 * to normalize_sitemaps
 *
 * Revision 1.2  2007/09/26 14:27:15  vanvoor4
 * Added initial support for testing octree
 *
 * Revision 1.1  2007/08/21 19:27:19  vanvoor4
 * initial checkin
 *
 *  
 *  
 *  
 */

#ifndef MODEL_SITEMAP_HEADER_FILE_INCLUDED
#define MODEL_SITEMAP_HEADER_FILE_INCLUDED

#include <Sitemap.H>
#include <ModelHbondSurfaces.H>
#include <TransformableTrimesh.H>

namespace ASCbase{

  //! A triangle has 3 labeled vertices and 3 edges
  typedef struct{
    my_float_t edge_len[3];   //!< Lengths of the triangle edges
    interact_pts_vci vertices[3];      //!< Iterators to the vertices
  }triangle_t;
    
  // This is a bit old -- I wouldn't use a list any more.  I wonder what the
  // current penalty is for using a list.  It might be minimal since this is
  // not an inner processing loop
  typedef triangle_t* triangle_ptr;
  typedef std::vector<std::list<triangle_t>::iterator> idx_tbl_lvl_short;
  typedef idx_tbl_lvl_short::const_iterator bucket_vci;
  typedef std::vector<idx_tbl_lvl_short> idx_tbl_lvl_long;
  typedef std::vector<idx_tbl_lvl_long> idx_tbl_lvl_perimeter;
  typedef std::vector<idx_tbl_lvl_perimeter> idx_tbl_lvl_class;
  typedef std::vector<idx_tbl_lvl_class> index_table_t;

  //! Create and provide an interface to the index table for the model (query) 
  //! sitemap
  class ModelSitemap : public Sitemap{
  public:
    ModelSitemap(const std::string path, const std::string struct_id,
                 const BaseParameters& args, 
                 const my_float_t max_corr_surf_pt_dist,
                 const bool normalize = true,
                 const bool load_hbond_surfaces = false,
                 const bool hydrophobic_query = false);

    virtual
    ~ModelSitemap();

    //! Given c,p,l,s get the iterators for the corresponding bucket
    /*! Given the indexes into the index table, set the pointers to the
     * begin and end of the bucket.
     *
     * @param c hash_class index 
     * @param p perimeter level index 
     * @param l longest edge length index
     * @param s shortest edge length index
     * @param bucket_begin Pointer to the beginning of the bucket
     * @param bucket_end Pointer to the end of the bucket
     * @return true if bucket is nonempty
     */
    virtual bool 
    get_bucket_iters(const uint c, const my_float_t perimeter,
                     const my_float_t max_len, const my_float_t min_len,
                     bucket_vci* bucket_begin, bucket_vci* bucket_end);

    bool report_table_stats(std::ostream& out);

    void init_positions_octrees(const my_float_t max_dist, 
                                const my_float_t min_side_len,
                                const uint max_numel);

    inline void 
    transform(const my_float_t* R, const my_float_t* T,
              const bool transform_mesh = true,
              const bool transform_hbond_surfaces = true)
    {
      if(A_model_hbond_surfaces && transform_hbond_surfaces) 
        A_model_hbond_surfaces->transform(R, T);
      if(A_binding_site_mesh && transform_mesh) 
        A_binding_site_mesh->transform(R, T);
      Sitemap::transform(R, T);
    }

    inline void 
    inverse_transform(const my_float_t* R, const my_float_t* T,
                      const bool transform_mesh = true, 
                      const bool transform_hbond_surfaces = true)
    {
      if(A_model_hbond_surfaces && transform_hbond_surfaces) 
        A_model_hbond_surfaces->inverse_transform(R, T);
      if(A_binding_site_mesh && transform_mesh) 
        A_binding_site_mesh->inverse_transform(R, T);
      Sitemap::inverse_transform(R, T);
    }

    inline void 
    revert(const bool transform_mesh = true, 
           const bool transform_hbond_surfaces = true)
    {
      if(A_model_hbond_surfaces && transform_hbond_surfaces) 
        A_model_hbond_surfaces->revert();
      if(A_binding_site_mesh && transform_mesh) A_binding_site_mesh->revert();
      Sitemap::revert();
    }

    //! Get transform to move current hbond fit points to original position of
    //! the hbond fit points
    inline void
    get_current_inverse_3D_transform(Quaternion *Q, my_float_t *T) const
    {

      // If ICP was used (only reason to currently call this function), 
      // we should be able to use surface points -- once things are allowed
      // some flexibility we will need to choose something like the fixed
      // alpha carbons to compute the rigid body transformation
      if(A_binding_site_mesh)
        A_binding_site_mesh->get_current_inverse_3D_transform(Q, T);
      else{
        // Sites need not have hbond points
        if(hbond_pts->fit_pts_size()) 
          hbond_pts->get_current_inverse_3D_transform(Q, T);
        else
          hphob_pts->get_current_inverse_3D_transform(Q, T);
      }
    }

    const ModelHbondSurfaces*
    model_hbond_surfaces() const
    { return A_model_hbond_surfaces; }

    //! A constant reference to the binding site mesh surface
    const geometry::TransformableTrimesh&
    binding_site_mesh_handle() const
    { return *A_binding_site_mesh; }

    //! A reference to the binding site mesh surface (not constant)
    /*! Use this function with care -- at the present it should only be
     * used when initializing a SurfDepsJoints instance
     */
    geometry::TransformableTrimesh&
    mutable_binding_site_mesh_handle()
    { return *A_binding_site_mesh; }

    ModelHbondSurfaces*
    mutable_model_hbond_surfaces()
    { return A_model_hbond_surfaces; }

    bool
    write_msms_caps(std::string &ofname_pref)
    {
      if(A_model_hbond_surfaces) 
        return A_model_hbond_surfaces->write_msms_caps(ofname_pref);
      std::cerr << "Cannot write msms caps if surfaces do not exist\n";
      return false;
    }

    //! Compute the maximum value for each feature and store in the given 
    //! vector
    void
    compute_max_feature_vals(std::vector<my_float_t> *max_feat_vals,
                             const my_float_t max_corr_surf_pt_dist) const;

    //! Get a const ref to the vector holding the max feature values
    const std::vector<my_float_t>&
    get_max_feature_vals() const
    { return A_max_feature_vals; }



    static const int pmax;  //!< Maximum number of buckets for the perimeter
    static const int lmax;  //!< Maximum number of buckets for the longest side
    static const int smax;  //!< Maximum number of buckets for the shortest side
    static const my_float_t MIN_PERIMETER;
    static const my_float_t MAX_PERIMETER;
    static const my_float_t MIN_LONGEST_SIDE;
    static const my_float_t MAX_LONGEST_SIDE;
    static const my_float_t MIN_SHORTEST_SIDE;
    static const my_float_t MAX_SHORTEST_SIDE;
    static const uint BUCKETS_PER_A_PERIMETER;
    static const uint BUCKETS_PER_A_LENGTH;
    static const uint NUMBER_OF_TYPE_HASH_CLASSES;
    static const std::string _fname;   //!< Source file name

  private:
    //! Frees any allocated objects & calls init() afterwards
    void clear();
  
    //! Initializes pointers to 0 and other variables to 0 or another useful 
    //! value
    void init();

    void create_index_table();

    void compute_combinations(std::vector<uint>* arr_p, int index, int length);

    // List is used because we have interators pointing to triangles and
    // are adding to the list at the same time.  And at no point do we 
    // look up items in the list by value, but rather, through the index table.
    std::list<triangle_t> triangles;
    index_table_t index_table;
    interact_pts_vci A_pts_begin;  //!< Iterator to the first split point (orig positions)
    ModelHbondSurfaces* A_model_hbond_surfaces;
    geometry::TransformableTrimesh *A_binding_site_mesh;
    std::vector<my_float_t> A_max_feature_vals;

  };

}

#endif
