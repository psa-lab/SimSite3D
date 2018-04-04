#include <boost/python.hpp>
#include <HbondSurfacesScore.H>
#include "../utils/stream_basicsPy.hh"

using namespace boost::python;
using ASCbase::HbondSurfacesScore;
using ASCbase::DbaseSitemap;
using ASCbase::ModelSitemap;
using ASCbase::SearchParameters;

// template this
class out_wrapper{
public:
  static void 
  write_header(HbondSurfacesScore& S, stl_ofstream& out)
  {
    S.write_score_header(out.get());
  }

  static void
  score_aligns(HbondSurfacesScore& S, ASCbase::rigid_align_vec& aligns, 
               DbaseSitemap* search, stl_ofstream& out)
  {
    S.score_alignments(aligns, search, out.get()); 
  }
};

BOOST_PYTHON_MODULE(_HbondSurfacesScore)
{
  class_<HbondSurfacesScore>("HbondSurfacesScore", 
                             init<ModelSitemap*, const SearchParameters&>())
    //.def("write_score_header", &HbondSurfacesScore::write_score_header)
    .def("write_score_header", &out_wrapper::write_header)
    //.def("score_alignments", &HbondSurfacesScore::score_alignments)
    .def("score_alignments", &out_wrapper::score_aligns)
  ;
}
