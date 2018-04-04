#include <boost/python.hpp>
#include <MatchTriangles.H>

using namespace boost::python;
using SimSite3D::MatchTriangles;
using SimSite3D::ModelSitemap;

BOOST_PYTHON_MODULE(_MatchTriangles)
{
  class_<MatchTriangles>("MatchTriangles", 
		         init<const my_float_t, const my_float_t, const bool>())
    .def("align", &MatchTriangles::align<SimSite3D::rigid_align_t>)
    ;
}
