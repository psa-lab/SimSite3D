#ifndef FACE_ATTRIB_HEADER_FILE_INCLUDED
#define FACE_ATTRIB_HEADER_FILE_INCLUDED

#include <dir_point.H>
#include <dir_point_storage.H>

namespace SimSite3D{
namespace geometry{

//! A class to hold attributes associated with a vertex
class VertAttrib : public dir_point_t{
public:
  typedef dir_point_storage<VertAttrib>::iterator vi; //!< Iterator into a dir_point_storage holding vertices
  typedef dir_point_storage<VertAttrib>::const_iterator vci; //!< Const iterator into a dir_point_storage holding vertices

  //! Default constructor
  VertAttrib(alloc_t a = ALLOC_POSITION) : dir_point_t(a)
  { ; }

  //! Copy constructor
  VertAttrib(const VertAttrib& V) : dir_point_t(V)
  {
    do_copy(V);
  }

  //! Assignment operator
  const VertAttrib&
  operator=(const VertAttrib& V)
  {
    if(&V != this){
      do_copy(V);
      dir_point_t::operator= (V);
    }
    return *this;
  }

  virtual ~VertAttrib()
  {
  }


  static dir_point_storage<VertAttrib> NULL_VERTICES_STORAGE;
  static const vi NULL_VI;
  static const vci NULL_VCI;

private:
  inline void
  do_copy(const VertAttrib& V)
  {
  }
};
 
} // end namespace geometry
} // end namespace SimSite3D

#endif
