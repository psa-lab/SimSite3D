#include <boost/python.hpp>
#include <ModelSitemap.H>

using namespace boost::python;
using ASCbase::ModelSitemap;

BOOST_PYTHON_MODULE(_ModelSitemap)
{
  class_<ModelSitemap, bases<ASCbase::Sitemap> >("ModelSitemap", 
		       init<const std::string, const std::string, 
		            const ASCbase::BaseParameters&, const bool>())
    .def(init<const std::string, const std::string,
              const ASCbase::BaseParameters&, const bool, const bool, 
              const bool>())
    .def("get_bucket_iters", &ModelSitemap::get_bucket_iters)
  ;
}
