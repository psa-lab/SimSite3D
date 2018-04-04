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

#ifndef DATABASE_SITEMAP_HEADER_FILE_INCLUDED
#define DATABASE_SITEMAP_HEADER_FILE_INCLUDED

#include <Sitemap.H>
#include <TransformableTrimesh.H>
#include <ImmovableTrimesh.H>

namespace SimSite3D{

//! Create and provide an interface to the database sitemaps.
/*! We need this class because currently the database sitemaps require their
 *  surfaces stored using the ImmovableTrimesh class
 */
class DbaseSitemap : public Sitemap{
public:
  DbaseSitemap(const std::string path, const std::string struct_id,
               const BaseParameters& args, 
               const my_float_t max_corr_surf_pt_dist_in=1.5,
               const bool check_all_triangles = false, 
               const bool load_hbond_surfaces = false)
    : Sitemap(path, struct_id, args, false, load_hbond_surfaces, false,
              VERBOSE_SILENT)
  { 
    init();
    A_max_corr_surf_pt_dist = max_corr_surf_pt_dist_in;
    if(args.load_surf_files && 
       !load_surface_mesh(struct_id, check_all_triangles)){
      set_fail_flag(true);
      return;
    }
  }

  ~DbaseSitemap()
  { clear(); } 

  const geometry::SimpleTrimeshTwo&
  binding_site_mesh_handle() const
  { return *A_binding_site_mesh; }

  const my_float_t
  max_corr_surf_pt_dist() const
  { return A_max_corr_surf_pt_dist; }

private:

  bool
  load_surface_mesh(const std::string struct_id, const bool check_all_triangles)
  {
    std::string surf_fname = site_path() + "/" + struct_id + "_surf";
    if(check_all_triangles){
      std::cout << "Allocating TransformableTrimesh for: " << struct_id << "\n";
      A_binding_site_mesh = new geometry::TransformableTrimesh(surf_fname);
    }else{
      std::cout << "Allocating ImmovableTrimesh for: " << struct_id << "\n";
       A_binding_site_mesh = 
         new geometry::ImmovableTrimesh(surf_fname, A_max_corr_surf_pt_dist);
    }

    if(A_binding_site_mesh->fail()){
      std::cerr << "Unable to load the surface mesh" << std::endl;
      return false;
    }
    return true;
  }

  //! Frees any allocated objects & calls init() afterwards
  void 
  clear()
  { 
    if(A_binding_site_mesh) delete A_binding_site_mesh;
    init(); 
  }

  //! Initializes pointers to 0 and other variables to 0 or another useful 
  //! value
  void 
  init()
  { 
    A_binding_site_mesh = 0; 
    A_max_corr_surf_pt_dist = 1.5;
  }

  my_float_t A_max_corr_surf_pt_dist;
  geometry::SimpleTrimeshTwo *A_binding_site_mesh;
};

} // end namespace SimSite3D

#endif
