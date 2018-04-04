#ifndef HBOND_SURFACES_HEADER_FILE_INCLUDED
#define HBOND_SURFACES_HEADER_FILE_INCLUDED

#include <sphere.H>
#include <hbond_points.H>
#include <BaseParameters.H>
#include <HbondBase.H>

namespace ASCbase{

//! A simple class to model a hydrogen bond spherical cap (surface for potential
//! complementary hydrogen bonding atom to sit)
/*!
 *
 * I am currently in a bind and to get this working quickly, I am 
 * abusing the format to handle metals as well.  In the case of metals
 * we do NOT have neighbor atom or an initial cap cutting plane (metals are
 * modeled as not having a prefered direction).  However, we do desire to keep
 * the method transformation invariant; to do this we put the closest atom to
 * the metal as the C_nbr_atom and the 2nd closest as the 2nd nbr atom 
 * (this is confusing but quick & dirty).
 */
class hbond_surface_t : public hbond_ideal_pt_base{
public: 
  //! Constructor to create the representation
  hbond_surface_t(const atom_vci hbond_atom, const atom_vci C_nbr_atom,
                  const atom_vci second_nbr_atom, 
                  BoundingVolume &site_vol, const int cap_number,
                  const bool include_metals = false,
                  const alloc_t a = ALLOC_POSITION)
    : hbond_ideal_pt_base(hbond_atom, C_nbr_atom, second_nbr_atom, site_vol,
                          cap_number, include_metals, a),
      A_surf(hbond_atom->pos, SURF_SPHERE_RAD)
  {
  }

  //! Constructor for an existing surface (typically read from file)
  hbond_surface_t(const std::string& data_line, PDBBase& prot_atoms);

  hbond_surface_t(const hbond_surface_t& src)
    : hbond_ideal_pt_base(src)
  { do_copy(src); }

  const hbond_surface_t&
  operator=(const hbond_surface_t& src)
  {
    if(this != &src){
      do_copy(src);
      hbond_ideal_pt_base::operator= (src);
    }
    return *this;
  }

  virtual ~hbond_surface_t()
  { ; }

  bool
  closest_point(const my_float_t* pt, my_float_t* closest_pt, 
                const my_float_t tol = 1.0);

  const my_float_t*
  ideal_dir() const
  { return A_cut_plane.normal(); }

  inline void
  transform(const my_float_t* R, const my_float_t* T)
  { std::cerr << "Unimplemented function in hbond_surface_t\n"; }

  inline void
  inverse_transform(const my_float_t* R, const my_float_t* T)
  { std::cerr << "Unimplemented function in hbond_surface_t\n"; }

  inline void
  revert()
  { std::cerr << "Unimplemented function in hbond_surface_t\n"; }


private:
  void 
  do_copy(const hbond_surface_t& src)
  {
    A_surf = src.A_surf;
    A_cut_plane = src.A_cut_plane;
    A_cap_plane_rad = src.A_cap_plane_rad;
    A_circles = std::vector<geometry::iCircle>(src.A_circles.begin(), 
                                               src.A_circles.end());
  }

  geometry::sphere_t A_surf;
  geometry::plane_t A_cut_plane;
  my_float_t A_cap_plane_rad; //!< Radius of cap restricted to cutting plane
  std::vector<geometry::iCircle> A_circles;
  static const my_float_t SURF_SPHERE_RAD;  //!< Radius of sphere used to define the cap
};

typedef std::vector<hbond_surface_t>::iterator hbond_surf_vi;

// Concentrate on the searching for now
template<class surf_class>
class HbondSurfaces : public HbondBase{
public:
  typedef surf_class						value_type;
  typedef typename std::vector<surf_class>      		surfaces_vec;
  typedef typename std::vector<surf_class>::iterator		surfaces_vi;
  typedef typename std::vector<surf_class>::const_iterator	surfaces_vci;

  HbondSurfaces(const std::string fname, PDBBase& rad_atoms,
                const verbose_level_t verbosity = VERBOSE_SILENT)
  {
    std::ifstream surf_caps_file;
    if(!open_ifstream(surf_caps_file, fname)) return;

    for(std::string line; std::getline(surf_caps_file, line); ){
      if(line.length() == 0 or line[0] == '#') continue;

      surf_class tmp_surf_cap(line, rad_atoms);
      A_surf_caps.push_back(tmp_surf_cap);
    }
  }

  virtual ~HbondSurfaces()
  { ; }

  //! Get a constant iterator to the first surface cap
  surfaces_vi
  surf_caps_begin()
  { return A_surf_caps.begin(); }

  //! Get a constant iterator to one past the last surface cap
  surfaces_vi
  surf_caps_end()
  { return A_surf_caps.end(); }

  //! Get a constant iterator to the first surface cap
  surfaces_vci
  surf_caps_begin() const 
  { return A_surf_caps.begin(); }

  //! Get a constant iterator to one past the last surface cap
  surfaces_vci
  surf_caps_end() const
  { return A_surf_caps.end(); }

  inline void
  transform(const my_float_t* R, const my_float_t* T)
  {
    for(surfaces_vi V = A_surf_caps.begin(); V < A_surf_caps.end(); ++V)
      V->transform(R, T);
  }

  inline void
  inverse_transform(const my_float_t* R, const my_float_t* T)
  {
    for(surfaces_vi V = A_surf_caps.begin(); V < A_surf_caps.end(); ++V)
      V->inverse_transform(R, T);
  }

  inline void
  revert()
  {
    for(surfaces_vi V = A_surf_caps.begin(); V < A_surf_caps.end(); ++V)
      V->revert();
  }
  
private:
  surfaces_vec A_surf_caps;  //!< Caps used to model hbond volumes
  //static const std::string A_fname;
};


} // end namespace ASCbase

#endif
