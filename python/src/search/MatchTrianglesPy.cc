#include <boost/python.hpp>
#include <MatchTriangles.H>

using namespace boost::python;
using ASCbase::MatchTriangles;
using ASCbase::ModelSitemap;

BOOST_PYTHON_MODULE(_MatchTriangles)
{
  class_<MatchTriangles>("MatchTriangles", 
		         init<const my_float_t, const my_float_t, const bool>())
    .def("align", &MatchTriangles::align<ASCbase::rigid_align_t>)
    ;
}
