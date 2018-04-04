#include <boost/python.hpp>
#include <ModelHbondSurfaces.H>

using namespace boost::python;
using ASCbase::ModelHbondSurfaces;
using ASCbase::HbondSurfaces;
using ASCbase::model_hbond_surf_t;



BOOST_PYTHON_MODULE(_ModelHbondSurfaces)
{
  class_< model_hbond_surf_t >
    ("model_hbond_surf_t", init<const std::string&, ASCbase::PDBBase&>())
  ;





  class_< ModelHbondSurfaces,
          bases< ASCbase::HbondSurfaces<model_hbond_surf_t> > >
    ("ModelHbondSurfaces", init<const std::string, ASCbase::PDBBase&>())
  ;
}
