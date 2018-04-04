#include <boost/python.hpp>
#include <HbondSurfaces.H>

using namespace boost::python;
using SimSite3D::HbondSurfaces;
using SimSite3D::hbond_surface_t;

BOOST_PYTHON_MODULE(_HbondSurfaces)
{
  class_< HbondSurfaces<hbond_surface_t> >("HbondSurfaces", 
    init<const std::string, SimSite3D::PDBBase&>())
  ;
}
