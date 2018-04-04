#include <boost/python.hpp>
#include <atom.H>
#include "__array_1.pypp.hpp"

using namespace boost::python;
using SimSite3D::point_t;
using SimSite3D::atom_t;

// Machinery to get at the position stored in point_t -- not foolproof since
// it doesn't support python type iteration but seems to work using numpy.array
// and has no problems using a loop "for i in 3: do_something(pos[i]) "
struct point_t_wrapper : point_t, wrapper<point_t>{
  point_t_wrapper(point_t const & arg ) 
  : point_t( arg ) , wrapper< point_t >(){
    // copy constructor
  }

  point_t_wrapper()
  : point_t() , wrapper< point_t >(){
    // null constructor
  }

  static pyplusplus::containers::static_sized::array_1_t< double, 3>
  pyplusplus_pos_wrapper( point_t & inst ){
    return pyplusplus::containers::static_sized::array_1_t<double, 3>(inst.pos);
  }
};


// Note one must be careful here since the point_t and additional machinery
// has a nonstandard allocation scheme so that the positions of atoms, etc
// can be stored in a contiguous hunk of memory to allow for rapid operations
// (i.e. transformation of coordinates).
BOOST_PYTHON_MODULE(_atom)
{
  { // point_t
    typedef class_<point_t_wrapper> point_t_exposer_t;
    point_t_exposer_t point_t_exposer = point_t_exposer_t("point_t");
    scope point_t_scope(point_t_exposer);
    // What do these lines do?
    //point_t_exposer.def( boost.python.self != boost.python.self );
    //point_t_exposer.def( boost.python.self == boost.python.self );
    pyplusplus::containers::static_sized::register_array_1< double, 3 >
      ( "__array_1_double_3" );
    { //SimSite3D::point_t::pos [variable], type=my_float_t[3]
      typedef pyplusplus::containers::static_sized::array_1_t< double, 3> 
        ( *array_wrapper_creator )( point_t & );
      point_t_exposer.add_property("pos", 
        make_function(array_wrapper_creator(&point_t_wrapper::pyplusplus_pos_wrapper), 
        with_custodian_and_ward_postcall< 0, 1 >()) );
    }
  }


  // UPdate this with enums at some point
  class_<atom_t, bases<point_t> >("atom", init<>())
    .def_readwrite("atom_num", &atom_t::atom_num)
    .def_readwrite("res_num", &atom_t::res_num)
    .def_readwrite("chainID", &atom_t::chainID)
    .def_readwrite("altLoc", &atom_t::altLoc)
    .def_readwrite("iCode", &atom_t::iCode)
    .def_readwrite("occupancy", &atom_t::occupancy)
    .def_readwrite("tempFactor", &atom_t::tempFactor)
    .def_readwrite("vdw_radius", &atom_t::vdw_radius)
    .def_readwrite("name_str", &atom_t::name_str)
    .def_readwrite("res_str", &atom_t::res_str)
    .def_readwrite("is_hetero", &atom_t::is_hetero)
  ;
}
