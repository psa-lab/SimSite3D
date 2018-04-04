#include <boost/python.hpp>
#include <DbaseSitemap.H>

using namespace boost::python;
using SimSite3D::DbaseSitemap;

BOOST_PYTHON_MODULE(_DbaseSitemap)
{
  class_<DbaseSitemap, bases<SimSite3D::Sitemap> >("DbaseSitemap", 
		       init<const std::string, const std::string, 
		            const SimSite3D::BaseParameters&, const my_float_t, 
                            const bool>())
    .def(init<const std::string, const std::string,
              const SimSite3D::BaseParameters&, const my_float_t, const bool, 
              const bool>())
  ;
}
