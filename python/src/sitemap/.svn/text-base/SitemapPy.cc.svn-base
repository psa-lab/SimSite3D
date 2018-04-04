#include <boost/python.hpp>
#include <Sitemap.H>

using namespace boost::python;
using ASCbase::Sitemap;
using ASCbase::HbondPoints;
using ASCbase::GenPointsParameters;

// From Py++
typedef HbondPoints const & ( Sitemap::*hbond_points_function_type )(  ) const;

BOOST_PYTHON_MODULE(_Sitemap)
{
  class_<Sitemap>("sitemap", 
		  init<const std::string, const std::string, 
		       const ASCbase::BaseParameters&, const bool >())
    .def(init<const GenPointsParameters&>())
    .def(init<const std::string, const std::string, 
         const ASCbase::BaseParameters&, const bool, const bool, const bool>())
    // Something screwy with the conversion of this function -- 
    // the Python function is returning True when Sitemap::A_fail is false
    //.def("fail", &Sitemap::fail)
    //.def("hbond_points", hbond_points_function_type(&Sitemap::hbond_points),
    //     return_value_policy<copy_const_reference>())
    // No idea if this will bite me later
    .def("atoms_file_name", &Sitemap::atoms_file_name)
    .def("hbond_points", hbond_points_function_type(&Sitemap::hbond_points),
        return_value_policy<reference_existing_object>())
    .def("write_files", &Sitemap::write_files)
  ;
}
