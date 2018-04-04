

#ifndef SIMPLE_TRIMESH_HEADER_INCLUDED
#define SIMPLE_TRIMESH_HEADER_INCLUDED

#include <map>
#include <types.H>
#include <BoundingVolume.H>
#include <SimpleOctree.H>


namespace ASCbase{

  //! Simple surface class built for speed of computing the shortest point to
  //! surface distance.
  class SimpleTrimesh{
  public:
    typedef struct{
      const my_float_t *vertices[3];
      const my_float_t *normals[3];
      my_float_t area;
    }triangle_t;

    typedef struct{
      my_float_t model_pt[3];
      my_float_t target_pt[3];
    }correspond_t;
    typedef std::multimap<my_float_t, correspond_t> correspond_map;
    typedef correspond_map::const_iterator correspond_mmci;

    SimpleTrimesh(const std::string fname, bool check_vertex_order = false);

    SimpleTrimesh(SimpleTrimesh& src, BoundingVolume &bvol);

    SimpleTrimesh(SimpleTrimesh& src, const my_float_t* center,
                  const my_float_t radius);
 
    ~SimpleTrimesh();

    void 
    clip(const my_float_t* center, const my_float_t radius);

    void 
    intersect(const BoundingVolume *bvol);

    bool 
    write(const std::string ofname) const ;

    my_float_t
    get_total_SA() const;

    bool old_compare(const SimpleTrimesh& target, const my_float_t max_pt_dist, 
                     size_t *num_points, my_float_t* RMSE, 
                     size_t *num_faces = 0, my_float_t *area = 0) const;

    void compare(const SimpleTrimesh& target, const my_float_t max_pt_dist, 
                 size_t *num_points, my_float_t* RMSE, correspond_map *cmap = 0,
                 size_t *num_faces = 0, my_float_t *area = 0) const;

    void 
    transform(const my_float_t* R, const my_float_t* T);

    void
    inverse_transform(const my_float_t* R, const my_float_t* T);

    void 
    revert();

    const my_float_t*
    vertices_begin() const
    { return A_vertices; }

    std::vector<triangle_t>::const_iterator
    faces_begin() const
    { return A_triangles.begin(); }

    std::vector<triangle_t>::const_iterator
    faces_end() const
    { return A_triangles.end(); }

    size_t
    number_of_vertices() const
    { return A_num_verts; }

    size_t
    number_of_faces() const
    { return A_triangles.size(); }

    bool
    fail() const
    { return A_fail; }

    // Called from target sitemaps
    void 
    build_vert_face_map();

    //! To be called from the query sitemap
    void
    build_vert_octree();

    typedef std::vector<triangle_t> face_vec;
    typedef face_vec::const_iterator face_vci;

  private:
    typedef struct{
      const my_float_t* old_n;
      const my_float_t* v;
      const my_float_t* n;
    }locations;

    typedef struct{
      my_float_t corr_pt[3];
      my_float_t distance;
    }pt_plus_dist;
    
    typedef std::map<const my_float_t*, pt_plus_dist> corr_map;
  
    void init();

    bool read_vert_file(const std::string fname);

    bool read_face_file(const std::string fname, bool check_vertex_order);

    size_t intersect(BoundingVolume &bvol, size_t *new_num_verts, 
                     my_float_t **new_vertices, my_float_t **new_normals, 
                     std::vector<triangle_t> *new_faces);

    size_t clip(const my_float_t* center, const my_float_t radius,
                size_t *new_num_verts, my_float_t **new_vertices, 
		my_float_t **new_normals, std::vector<triangle_t> *new_faces);

    void copy_kept(std::map<const my_float_t*, locations> &V_keep,
                   my_float_t **new_vertices, my_float_t **new_normals, 
                   std::vector<triangle_t> *new_faces);

    void compute_areas();

    size_t A_num_verts;
    my_float_t *A_vertices;
    my_float_t *A_normals;
    my_float_t *A_original_vertices;
    my_float_t *A_original_normals;

    typedef std::map<const my_float_t*, std::vector<face_vci> > vert_face_map;
    typedef vert_face_map::iterator vert_face_mi;
    typedef vert_face_map::const_iterator vert_face_mci;

    face_vec A_triangles;
    vert_face_map A_map_vert_to_faces;
    point_vec_t A_vert_vector;
    SimpleOctree<point_vci> *A_vert_octree;
    bool A_fail;
    static const std::string A_fname;
  };
}
#endif
