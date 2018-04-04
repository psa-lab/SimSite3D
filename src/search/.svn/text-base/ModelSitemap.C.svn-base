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
 * $Source: /psa/share/repository/pfizer_proj/src/search/ModelSitemap.C,v $
 * $Revision: 1.6 $
 * $Author: vanvoor4 $
 * $Date: 2009-01-12 21:22:53 $
 * 
 * $Log: not supported by cvs2svn $
 * Revision 1.5  2008/07/29 14:54:58  vanvoor4
 * Added support for binding site surface (mesh)
 *
 * Revision 1.4  2008/05/13 16:35:55  vanvoor4
 * Don't build the octree for now
 *
 * Revision 1.3  2007/11/01 16:40:54  vanvoor4
 * Changed due to the fact we may need to pass the parameters class to
 * normalize_sitemaps()
 *
 * Revision 1.2  2007/09/26 14:27:20  vanvoor4
 * Added initial support for testing octree
 *
 * Revision 1.1  2007/08/21 19:27:26  vanvoor4
 * initial checkin
 *
 *  
 *  
 *  
 */

#include <ModelSitemap.H>
#include <hash_class.H>
#include <iomanip>

using namespace SimSite3D;

const std::string ModelSitemap::_fname = "ModelSitemap.C";
const int ModelSitemap::pmax = static_cast<int>(ceil(BUCKETS_PER_A_PERIMETER * (MAX_PERIMETER - MIN_PERIMETER)));
const int ModelSitemap::lmax = static_cast<int>(ceil(BUCKETS_PER_A_LENGTH * (MAX_LONGEST_SIDE - MIN_LONGEST_SIDE)));
const int ModelSitemap::smax = static_cast<int>(ceil(BUCKETS_PER_A_LENGTH * (MAX_SHORTEST_SIDE - MIN_SHORTEST_SIDE)));

#if 1
const my_float_t ModelSitemap::MIN_PERIMETER = 9.0;
const my_float_t ModelSitemap::MAX_PERIMETER = 13.0;
const my_float_t ModelSitemap::MIN_LONGEST_SIDE = 3.5;
const my_float_t ModelSitemap::MAX_LONGEST_SIDE = 4.5;
const my_float_t ModelSitemap::MIN_SHORTEST_SIDE = 1.8;
const my_float_t ModelSitemap::MAX_SHORTEST_SIDE = 3.5;
#else
const my_float_t ModelSitemap::MIN_PERIMETER = 9.0;
const my_float_t ModelSitemap::MAX_PERIMETER = 16.0;
const my_float_t ModelSitemap::MIN_LONGEST_SIDE = 4.0;
const my_float_t ModelSitemap::MAX_LONGEST_SIDE = 7.0;
const my_float_t ModelSitemap::MIN_SHORTEST_SIDE = 1.8;
const my_float_t ModelSitemap::MAX_SHORTEST_SIDE = 4.0;
#endif

const uint ModelSitemap::BUCKETS_PER_A_PERIMETER = 2;
const uint ModelSitemap::BUCKETS_PER_A_LENGTH = 4;
const uint ModelSitemap::NUMBER_OF_TYPE_HASH_CLASSES = 10;

// Note -- we pass false to the base class for load_hbond_surfaces since 
// the class ModelHbondSurfaces takes care of everything for the query site
ModelSitemap::ModelSitemap(const std::string path, const std::string struct_id,
                           const BaseParameters& args, 
                           const my_float_t max_corr_surf_pt_dist,
                           const bool normalize, const bool load_hbond_surfaces,
                           const bool hydrophobic_query)
 : Sitemap(path, struct_id, args, normalize, false, hydrophobic_query, 
           VERBOSE_SILENT)
{
  init();
  if(fail()) return;

  if(!load_volume(path, struct_id, args.dbase_sites, "<site_vol_est>", 
                  &A_site_vol_est)){
    std::string msg = "Unable to get the query sitemap volume estimate\n";
    err_msg(_fname, "cstr", msg);
    set_fail_flag(true);
    return;
  }

  A_pts_begin = interact_pts_beg();
  create_index_table(); 
  if(!triangles.size()){
    set_fail_flag(true);
    warn(_fname, "cstr()", std::string("Query sitemap (") + struct_id +
         ") points are too sparse");
    report_table_stats(std::cout);
    return;
  }

  if(args.load_surf_files){
    std::cout << "Allocating TransformableTrimesh for: " << struct_id << "\n";
    std::string surf_fname = site_path() + "/" + struct_id + "_surf";
    A_binding_site_mesh = new geometry::TransformableTrimesh(surf_fname);
    if(A_binding_site_mesh->fail()){
      std::cerr << "Unable to load the surface mesh" << std::endl;
      set_fail_flag(true);
      return;
    }
  }else
    std::cout << "You have chosen to not load the molecular surfaces" 
              << std::endl;

  if(hydrophobic_query) std::cout << "Query site map is hydrophobic\n";

  if(load_hbond_surfaces){
        // Verify that path to sitemap exists
    std::cout << "Loading model hbond surfaces for: " << struct_id << "\n";
    std::string my_struct = site_path() + "/" + struct_id;
    std::string my_site = my_struct + "_surf_caps.csv";
    if(!normal_file_exists(my_site, false)){
      std::cerr << "Unable to find the file: " << my_site << "\n"
                << "Skipping the binding site of " << struct_id << "\n";
      set_fail_flag(true);
      return;
    }
    A_model_hbond_surfaces =
      new ModelHbondSurfaces(my_site, Sitemap::bind_site_atoms(), VERBOSE_SILENT);
  }

  compute_max_feature_vals(&A_max_feature_vals, max_corr_surf_pt_dist);
  std::cout << "Inside ModelSitemap::ModelSitemap, num points is: " 
            << sitemap_points().num_atoms() << "\n";

}

ModelSitemap::~ModelSitemap()
{
  clear();
}

void
ModelSitemap::clear()
{
  if(A_model_hbond_surfaces) delete A_model_hbond_surfaces;
  if(A_binding_site_mesh) delete A_binding_site_mesh;
  init();
}

void
ModelSitemap::init()
{
  A_model_hbond_surfaces = 0;
  A_binding_site_mesh = 0;
}

bool 
ModelSitemap::get_bucket_iters(const uint c, const my_float_t perimeter,
                               const my_float_t max_len,
                               const my_float_t min_len,
                               bucket_vci* bucket_begin, bucket_vci* bucket_end)
{
  int p = -1;
  int l = -1;
  int s = -1;
  if(perimeter >= MIN_PERIMETER && perimeter < MAX_PERIMETER)
    p = static_cast<uint>(BUCKETS_PER_A_PERIMETER *(perimeter - MIN_PERIMETER));
  if(max_len >= MIN_LONGEST_SIDE && max_len < MAX_LONGEST_SIDE)
    l = static_cast<uint>(BUCKETS_PER_A_LENGTH * (max_len - MIN_LONGEST_SIDE));
  if(min_len >= MIN_SHORTEST_SIDE && min_len < MAX_SHORTEST_SIDE)
    s = static_cast<uint>(BUCKETS_PER_A_LENGTH * (min_len - MIN_SHORTEST_SIDE));
  if(p == -1 || l == -1 || s == -1) return false;

  idx_tbl_lvl_short& bucket = index_table[c][p][l][s];
  if(bucket.empty()) return false;
  *bucket_begin = bucket.begin();
  *bucket_end = bucket.end();
  return true; 
}

void
ModelSitemap::create_index_table()
{
  // init
  index_table.resize(NUMBER_OF_TYPE_HASH_CLASSES);
  for(uint c = 0; c < NUMBER_OF_TYPE_HASH_CLASSES; ++c){
    index_table[c].resize(pmax);

    for(int p = 0; p < pmax; ++p){
      index_table[c][p].resize(lmax);

      for(int l = 0; l < lmax; ++l) index_table[c][p][l].resize(smax);
    }
  }
 
  // create
  std::vector<uint> tmp_array(interact_pts_end() - A_pts_begin);
  for(uint i = 0; i < tmp_array.size(); ++i) tmp_array[i] = i;
  for(uint i = 2; i < tmp_array.size(); i++){
    tmp_array[2] = tmp_array[i];
    compute_combinations(&tmp_array, 2, 3);
  }
}

void
ModelSitemap::compute_combinations(std::vector<uint>* arr_p,
                                   int index, int length)
{
  std::vector<uint>& arr = *arr_p;

  // Get itertors, edge lengths and perimeter for the given vertices
  triangle_t triangle; 
  for(uint i = 0; i < 3; ++i) triangle.vertices[i] = A_pts_begin + arr[i];
  my_float_t perimeter = 0;
  for(uint i = 0; i < 3; ++i){
    triangle.edge_len[i] = 
      dist(triangle.vertices[i]->pos, triangle.vertices[(i+1) % 3]->pos);
    perimeter += triangle.edge_len[i];
  }

  // Find min and max edge lengths -- eqilateral triangle may cause undefined
  // results?
  my_float_t min_len = triangle.edge_len[0];
  my_float_t max_len = triangle.edge_len[0];
  for(uint i = 1; i < 3; ++i){
    if(triangle.edge_len[i] > max_len) max_len = triangle.edge_len[i];
    else if(triangle.edge_len[i] < min_len) min_len = triangle.edge_len[i];
  }

  // Find the triangles bucket based on its perimeter and edge lengths
  int perimeter_class = -1;
  if(perimeter >= MIN_PERIMETER && perimeter <= MAX_PERIMETER)
    perimeter_class = 
      static_cast<int>(BUCKETS_PER_A_PERIMETER * (perimeter - MIN_PERIMETER));
  int longest_class = -1;
  if(max_len >= MIN_LONGEST_SIDE && max_len <= MAX_LONGEST_SIDE)
    longest_class = 
      static_cast<int>(BUCKETS_PER_A_LENGTH * (max_len - MIN_LONGEST_SIDE));
  int shortest_class = -1;
  if(min_len >= MIN_SHORTEST_SIDE && min_len <= MAX_SHORTEST_SIDE)
    shortest_class = 
      static_cast<int>(BUCKETS_PER_A_LENGTH * (min_len - MIN_SHORTEST_SIDE));
    
  // Find the buckets which the triangle fits in -- including one above and
  // below the edge lengths and perimeter classes (if they exist).
  if(perimeter_class != -1 && longest_class != -1 && shortest_class != -1){
    bool first = true;
    int c = 0;
    for(int i = 0; i < 3; ++i) c += triangle.vertices[i]->act_type;
    c = hash_class[c];
    int p_start = (0 > perimeter_class - 1 ? 0 : perimeter_class - 1);
    int p_stop = (pmax < perimeter_class + 2 ? pmax : perimeter_class + 2);
    for(int p = p_start; p < p_stop; p++){
      int l_start = (0 > longest_class - 1 ? 0 : longest_class - 1);
      int l_stop = (lmax < longest_class + 2 ? lmax : longest_class + 2);
      for(int l = l_start; l < l_stop; l++){
        int s_start = (0 > shortest_class - 1 ? 0 : shortest_class - 1);
        int s_stop = (smax < shortest_class + 2 ? smax : shortest_class + 2);
        for(int s = s_start; s < s_stop; s++){
          if(first){
            triangles.push_front(triangle);
            first = false;
          }
  
          // Push back an iterator to the stored triangle
          index_table[c][p][l][s].push_back(triangles.begin());
        }
      }
    }
  }

  if(index > 0){
    int orig = arr[index-1]; /* store original entry at position index-1 */
    while(arr[index-1] < arr[index] - 1) {
      /* increase value at position index - 1 as long as it is
         still smaller than the value at position index */
      arr[index-1]++;
      /* recursively increase all elements in the array at
         positions [0..index-2] */
      compute_combinations(arr_p, index - 1, length);
    }
    /* restore original value at position index-1 */
    arr[index-1] = orig;
  }
}

bool
ModelSitemap::report_table_stats(std::ostream& out)
{
  uint num_tri_ptrs = 0;
  uint empty = 0;
  uint count = 0;
  uint max = 0;
  for(uint c = 0; c < 10; c++ )
    for(int p = 0; p < pmax; p++)
      for(int l = 0; l < lmax; l++)
        for(int s = 0; s < smax; s++)
          if(index_table[c][p][l][s].size() > 0 ){
            num_tri_ptrs += index_table[c][p][l][s].size();
            count++;
            if(index_table[c][p][l][s].size() > max )
              max = index_table[c][p][l][s].size();
          }
          else empty++;

  out << "Multi-level indexing stats\n";
  out << std::setw(10) << triangles.size() << " template triangles\n";
  out << std::setw(10) << num_tri_ptrs << " triangle pointers\n"; 
  out << std::setw(10) << empty << " empty buckets\n";
  out << std::setw(10) << count << " non-empty-buckets\n";

  if(num_tri_ptrs != 0)
    out << std::setw(10) << num_tri_ptrs / count << " triangles per bucket "
        << "on average\n" 
        << std::setw(10) << max << " maximal number of triangles\n";
  else{
    std::cerr << "No triangles identified for the query (model) sitemap.\n"
              << "Program exiting\n";
    return false;
  }
  return true;
}

void
ModelSitemap::compute_max_feature_vals(std::vector<my_float_t> *max_feat_vals,
                                       const my_float_t max_corr_surf_pt_dist)
const
{
  Sitemap::compute_max_feature_values(max_feat_vals);
  if(A_binding_site_mesh){
    max_feat_vals->push_back(A_binding_site_mesh->number_of_vertices());
    max_feat_vals->push_back(A_binding_site_mesh->number_of_faces());
  
    max_feat_vals->push_back(max_corr_surf_pt_dist);
    max_feat_vals->push_back(max_corr_surf_pt_dist);
    max_feat_vals->push_back(A_binding_site_mesh->get_total_SA());
  
    // Currently max value for norm penalty is 4.0
    max_feat_vals->push_back(4.0);
  }else{
    for(int i = 0; i < 5; ++i) max_feat_vals->push_back(-1.0);
  }

  if(A_model_hbond_surfaces){
    int num_pts = A_model_hbond_surfaces->num_A_D_points();
    if(num_pts < 1) num_pts = 1; 
    for(int i = 0; i < 6; ++i) max_feat_vals->push_back(num_pts);

    num_pts = A_model_hbond_surfaces->num_polar_points();
    if(num_pts < 1) num_pts = 1; 
    for(int i = 0; i < 6; ++i) max_feat_vals->push_back(num_pts);

    // still missing the max value for term 27 (complementary surface area 
    // of the caps)
    max_feat_vals->push_back(-1.0);
  }
}

