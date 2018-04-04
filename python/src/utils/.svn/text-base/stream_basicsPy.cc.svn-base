#include <boost/python.hpp>
#include <stream_basics.H>
#include "stream_basicsPy.hh"

using namespace boost::python;

// Hijacked from Py++
BOOST_PYTHON_MODULE(_stream_basics)
{
  class_<stl_ofstream, boost::noncopyable >( "stl_ofstream", no_init )
    .def( init< >() )
    .def("open", &stl_ofstream::open)
    .def("close", &stl_ofstream::close)
    .def("write", &stl_ofstream::write)
  ;
}
