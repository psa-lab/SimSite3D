#include <VertAttrib.H>

using namespace ASCbase;
using namespace ASCbase::geometry;

dir_point_storage<VertAttrib> VertAttrib::NULL_VERTICES_STORAGE;
const VertAttrib::vi VertAttrib::NULL_VI = 
  VertAttrib::NULL_VERTICES_STORAGE.end();
const VertAttrib::vci VertAttrib::NULL_VCI = 
  VertAttrib::NULL_VERTICES_STORAGE.end();
