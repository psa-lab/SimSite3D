#include <boost/python.hpp>
#include <DbaseSitemap.H>

using namespace boost::python;
using ASCbase::DbaseSitemap;

BOOST_PYTHON_MODULE(_DbaseSitemap)
{
  class_<DbaseSitemap, bases<ASCbase::Sitemap> >("DbaseSitemap", 
		       init<const std::string, const std::string, 
		            const ASCbase::BaseParameters&, const my_float_t, 
                            const bool>())
    .def(init<const std::string, const std::string,
              const ASCbase::BaseParameters&, const my_float_t, const bool, 
              const bool>())
  ;
}
