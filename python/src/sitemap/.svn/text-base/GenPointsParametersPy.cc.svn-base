#include <boost/python.hpp>
#include <GenPointsParameters.H>

using namespace boost::python;
using SimSite3D::GenPointsParameters;

class verify_wrapper{
public:
  static void
  verify_params(GenPointsParameters &P)
  {
    P.verify_params();
  }
};

BOOST_PYTHON_MODULE(_parameters)
{
  class_<GenPointsParameters, bases<SimSite3D::BaseParameters> >
    ("parameters", init< > ())
    .def("verify_params", &verify_wrapper::verify_params)
    .def_readwrite("pts_fname", &GenPointsParameters::pts_fname)
    .def_readwrite("prot_fname", &GenPointsParameters::prot_fname)
    .def_readwrite("lig_fname", &GenPointsParameters::lig_fname)
    .def_readwrite("sphere_str", &GenPointsParameters::sphere_str)
    .def_readwrite("normalize", &GenPointsParameters::normalize)
    .def_readwrite("include_metals", &GenPointsParameters::include_metals)
    .def_readwrite("call_msms", &GenPointsParameters::call_msms)
    // All other variables will take more thought -- trying to be fast here
  ;
}
