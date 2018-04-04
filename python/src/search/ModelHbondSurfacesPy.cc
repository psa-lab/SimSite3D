#include <boost/python.hpp>
#include <ModelHbondSurfaces.H>

using namespace boost::python;
using SimSite3D::ModelHbondSurfaces;
using SimSite3D::HbondSurfaces;
using SimSite3D::model_hbond_surf_t;



BOOST_PYTHON_MODULE(_ModelHbondSurfaces)
{
  class_< model_hbond_surf_t >
    ("model_hbond_surf_t", init<const std::string&, SimSite3D::PDBBase&>())
  ;





  class_< ModelHbondSurfaces,
          bases< SimSite3D::HbondSurfaces<model_hbond_surf_t> > >
    ("ModelHbondSurfaces", init<const std::string, SimSite3D::PDBBase&>())
  ;
}
