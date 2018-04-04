
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <types.H>

namespace bp = boost::python;

BOOST_PYTHON_MODULE(_stl_containers)
{
    { //::std::vector< double >
        typedef bp::class_< std::vector< double > > vector_less__double__greater__exposer_t;
        vector_less__double__greater__exposer_t vector_less__double__greater__exposer = vector_less__double__greater__exposer_t( "vector_less__double__greater_" );
        bp::scope vector_less__double__greater__scope( vector_less__double__greater__exposer );
        vector_less__double__greater__exposer.def( bp::vector_indexing_suite< ::std::vector< double >, true >() );
    }

    { //::std::vector< bool >
        typedef bp::class_< std::vector< bool > > vector_less__bool__greater__exposer_t;
        vector_less__bool__greater__exposer_t vector_less__bool__greater__exposer = vector_less__bool__greater__exposer_t( "vector_less__bool__greater_" );
        bp::scope vector_less__bool__greater__scope( vector_less__bool__greater__exposer );
        vector_less__bool__greater__exposer.def( bp::vector_indexing_suite< ::std::vector< bool >, true >() );
    }

#if 0
  class_<std::vector<bool> >("vec_bool")
    .def(vector_indexing_suite<std::vector<bool> >())
  ;

  class_<std::vector<my_float_t> >("vec_float")
    .def(vector_indexing_suite<std::vector<my_float_t> >())
  ;
#endif
}	
