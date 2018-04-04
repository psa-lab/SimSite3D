#include <boost/python.hpp>
#include "../utils/stream_basicsPy.hh"
#include <ProtLigScoreParameters.H>

using namespace boost::python;
using SimSite3D::ProtLigScoreParameters;

BOOST_PYTHON_MODULE(_parameters)
{
  class_<ProtLigScoreParameters, bases<SimSite3D::BaseParameters> >
    ("parameters", init< > ())
    .def_readwrite("prot_fname", &ProtLigScoreParameters::prot_fname)
    .def_readwrite("lig_fname", &ProtLigScoreParameters::lig_fname)
    .def_readwrite("lig_list_fname", &ProtLigScoreParameters::lig_list_fname)
    .def_readwrite("build_interact_tbl", 
                   &ProtLigScoreParameters::build_interact_tbl)
    .def_readwrite("print_interactions", 
                   &ProtLigScoreParameters::print_interactions)
  ;
}
