#include <boost/python.hpp>
#include <boost/python/iterator.hpp>
#include <HbondPoints.H>

using namespace boost::python;
using SimSite3D::HbondPoints;
using SimSite3D::PDBBase;
using SimSite3D::hbond_ideal_pt_vci;



// Declared as no_init for now -- this class is included so that I can use a
// pointer to HbondPoints and get the ideal points
BOOST_PYTHON_MODULE(_HbondPoints)
{
  //class_<HbondPoints, boost::noncopyable>("hbond_points", no_init)
  class_<HbondPoints>("hbond_points", no_init)
    .def( init<std::istream&, const uint, PDBBase&, PDBBase&>() )
    .def("ideal_pts", 
         range<return_value_policy<reference_existing_object> > 
           (&HbondPoints::ideal_pts_beg, &HbondPoints::ideal_pts_end))
    .def("fit_pts",
         range<return_value_policy<reference_existing_object> >
           (&HbondPoints::fit_pts_beg, &HbondPoints::fit_pts_end))
    .def("transform", &HbondPoints::transform)
    .def("revert", &HbondPoints::revert)
  ;
}
