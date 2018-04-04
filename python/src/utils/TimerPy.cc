#include <boost/python.hpp>
#include <Timer.H>
using namespace boost::python;
using ASCbase::Timer;

object
wrap(Timer& T)
{
  double r = -1.0;
  double v = -1.0;
  double p = -1.0;
  T.get(&r, &v, &p);
  return object(make_tuple(r, v, p));
}

BOOST_PYTHON_MODULE(_system_timers)
{
  class_<Timer, boost::noncopyable>("system_timers", init<>())
    .def("start", &Timer::start)
//    .def("get", &Timer::get_real_timer)
    .def("get", &wrap)
    .def("fail", &Timer::fail)
    .def("started", &Timer::started)
  ;

}
