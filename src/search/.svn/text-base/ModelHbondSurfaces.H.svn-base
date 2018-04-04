
#ifndef MODEL_HBOND_SURFACES_HEADER_INCLUDED
#define MODEL_HBOND_SURFACES_HEADER_INCLUDED

#include <SimpleTrimesh.H>
#include <HbondSurfaces.H>
#include <TriMeshSphere.H>

namespace ASCbase{

//! A simple class to model a hydrogen bond spherical cap by a triangle mesh
//! (surface for potential complementary hydrogen bonding atom to sit)
/*!
 *
 * I am currently in a bind and to get this working quickly, I am 
 * abusing the format to handle metals as well.  In the case of metals
 * we do NOT have neighbor atom or an initial cap cutting plane (metals are
 * modeled as not having a prefered direction).
 */
class model_hbond_surf_t : public hbond_surface_t{
public:

  model_hbond_surf_t(const std::string& data_line, PDBBase& prot_atoms,
                     const uint level=2) ;
  
  model_hbond_surf_t(const model_hbond_surf_t& src)
    : hbond_surface_t(src), A_mesh_surf(src.A_mesh_surf)
  { 
    init();
    do_copy(src); 
  }

  model_hbond_surf_t& operator=(const model_hbond_surf_t& src)
  {
    if(this == &src) return *this;

    hbond_surface_t::operator=(src);
    A_mesh_surf = src.A_mesh_surf;
    do_copy(src);
    return *this;
  }

  virtual ~model_hbond_surf_t()
  { 
    if(A_dists) delete [](A_dists);
    if(A_closest_pts) delete [](A_closest_pts);
    init(); 
  }

  inline void
  transform(const my_float_t* R, const my_float_t* T)
  { A_mesh_surf.transform(R, T); }

  inline void
  inverse_transform(const my_float_t* R, const my_float_t* T)
  { A_mesh_surf.inverse_transform(R, T); }

  inline void
  revert()
  { A_mesh_surf.revert(); }

  //! Compute the best value of the 6 features for each point & the estimate
  //! of the complimentary surface area
  /*!
   * The features for each point are:
   *   <ol>
   *     <li>Positive value if point was matched</li>
   *     <li>Best dot product for pts within tolerace</li>
   *     <li>Linear weighting of closest point within tolerance</li>
   *     <li>Squared weighting of closest point within tolerance</li>
   *     <li>Dot product times linear weighting of closest point within 
   *         tolerance</li>
   *     <li>Dot product times squared weighting of closest point within 
   *         tolerance</li>
   *   </ol>
   */
  void
  compute_best_terms(std::vector<hbond_surface_t>::iterator surfs_begin,
                     std::vector<hbond_surface_t>::iterator surfs_end,
                     my_float_t *AA_DD_terms, my_float_t *doneptor_terms,
                     my_float_t *surface_area, const my_float_t dist_tol=1.0) 
                     const;

  const uint
  num_terms() const
  { return A_num_terms; }

  const uint
  num_verts() const
  { return A_mesh_surf.num_verts(); }

  //! Get a constant iterator to the first vertex in the vector
  const geometry::vertex_vci
  verts_begin() const
  { return A_mesh_surf.verts_begin(); }

  //! Get a constant iterator to the end of the vertex vector
  const geometry::vertex_vci
  verts_end() const
  { return A_mesh_surf.verts_end(); }

  //! Get a constant pointer to the first point in the closest points array
  const my_float_t*
  closest_pts_begin() const
  { return A_closest_pts; }

  //! Get a constant pointer to the end of the closest points array
  const my_float_t*
  closest_pts_end() const
  { return A_closest_pts +  3 * num_verts(); }

  //! Get a constant pointer to the first distance in the distances array
  const my_float_t*
  closest_pts_dists_begin() const
  { return A_dists; }

  //! Get a constant pointer to the end of the distances array
  const my_float_t*
  closest_pts_dists_end() const
  { return A_dists + num_verts(); }

  void
  write(std::ostream &out, const interactionType point_type, 
        const char delim='|') const
  { A_mesh_surf.write(out, point_type, delim); }

  inline void
  write_msms_headers(std::ostream &vert_out, std::ostream &face_out) const
  { A_mesh_surf.write_msms_headers(vert_out, face_out); }

  inline void
  write_msms_cap(std::ostream &vert_out, std::ostream &face_out) const
  { A_mesh_surf.write_msms_cap(vert_out, face_out, act_type); }

private:
  bool
  adjust_terms(const my_float_t dot_prod, const my_float_t distance,
               const my_float_t dist_tol, my_float_t* row) const
  {
    // Compute the best value for each term
    row[0] = 1;
    row[1] = (row[1] < dot_prod ? dot_prod : row[1]);
        
    my_float_t linear_w = (dist_tol - distance) / dist_tol;
    if(row[2] < linear_w){
      row[2] = linear_w;
      row[3] = linear_w * linear_w;
    }

/*
    std::cout << "linear W: " << linear_w << std::endl;
    std::cout << "dot_prod: " << dot_prod << std::endl;
    std::cout << "row[4]: " << row[4] << std::endl;
*/
    my_float_t tmp = linear_w * dot_prod;
    if(row[4] < tmp){
      row[4] = tmp;
      row[5] = tmp * linear_w;
      return true;
    }
    return false;
  }

  void
  do_copy(const model_hbond_surf_t& src)
  {
    if(src.A_closest_pts){
      if(A_closest_pts) delete [](A_closest_pts);
      size_t N = 3 * src.A_mesh_surf.num_verts();
      A_closest_pts = new my_float_t[N];
      std::copy(src.A_closest_pts, src.A_closest_pts + N, A_closest_pts);
    }else A_closest_pts = 0;

    if(src.A_dists){
      if(A_dists) delete [](A_dists);
      size_t N = src.A_mesh_surf.num_verts();
      A_dists = new my_float_t[N];
      std::copy(src.A_dists, src.A_dists + N, A_dists);
    }else A_dists = 0;
  }

  inline void
  init()
  {
    A_closest_pts = 0;
    A_dists = 0; 
  }

  geometry::TriMeshSphere A_mesh_surf;
  my_float_t *A_closest_pts;  //! Current closest point for each cap point
  my_float_t *A_dists;        //! Current distances for corresponding points
  static const uint A_num_terms;
};

class ModelHbondSurfaces : public HbondSurfaces<model_hbond_surf_t>{

public:
  typedef HbondSurfaces<model_hbond_surf_t>::surfaces_vci m_surf_vci;

  ModelHbondSurfaces(const std::string fname, PDBBase& rad_atoms,
                     const verbose_level_t verbosity = VERBOSE_SILENT);

  virtual ~ModelHbondSurfaces()
  { ; }

  bool
  write_msms_caps(std::string &ofname_pref)
  {

    std::ofstream face_file, vert_file;
    if(!open_ofstream(vert_file, ofname_pref + ".vert")) return false;
    if(!open_ofstream(face_file, ofname_pref + ".face")) return false;

    m_surf_vci hb_surf = surf_caps_begin();
    hb_surf->write_msms_headers(vert_file, face_file);
    for( ; hb_surf != surf_caps_end(); ++hb_surf){
      // This is bad, but need to be fast for now -- later we should
      // just renumber the verts
      hb_surf->write_msms_cap(vert_file, face_file);
      vert_file << "\n<+++++ END +++++>";
      face_file << "\n<+++++ END +++++>";
    }
    return true;
  }

  //! Get the number of acceptor and donor points
  const size_t
  num_A_D_points() const
  { return A_num_A_D_points; }

  //! Get the number of polar points (acceptor, donor, and doneptor)
  const size_t
  num_polar_points() const
  { return A_num_polar_points; }
  
private:
  void
  init()
  {
    A_num_A_D_points = 0;
    A_num_polar_points = 0;
  }

  size_t A_num_A_D_points;
  size_t A_num_polar_points;

};


} // end namespace ASCbase
#endif // endif header file included
