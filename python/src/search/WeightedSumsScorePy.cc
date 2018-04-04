#include <boost/python.hpp>
#include <WeightedSumsScore.H>
#include <point_and_surf_score.H>
#include <NoTier1Score.H>
#include <ScoreRigidAlignments.H>
#include "../utils/stream_basicsPy.hh"

using namespace boost::python;
using ASCbase::rigid_align_t;
using ASCbase::ScoreRigidAlignments;
using ASCbase::WeightedSumsScore;
using ASCbase::point_and_surf_score;
using ASCbase::NoTier1Score;
using ASCbase::DbaseSitemap;
using ASCbase::ModelSitemap;
using ASCbase::SearchParameters;

// Wrapper class to handle passing in an std::ofstream object
template < class tier1_SF, class tier2_SF, typename align_T >
class out_wrapper{
public:
  static void 
  write_header(ScoreRigidAlignments< tier1_SF, tier2_SF, align_T>& S, 
               stl_ofstream& out)
  {
    S.write_score_header(out.get());
  }

  static void
  score_aligns(ScoreRigidAlignments< tier1_SF, tier2_SF, align_T >& S, 
               std::vector<align_T>& aligns, 
               DbaseSitemap* search, stl_ofstream& out)
  {
    S.score_alignments(aligns, search, out.get()); 
  }
};

// Test this method first.  If it works as expected we an rename this file
// and put in the other scoring functions or if we really want to 
// (which I don't at this time) take the effort to understand (if it is even
// possible) to specify the template paramters from python code (I recall
// this being impossible as the boost::python interface requires explicit
// instances of the C++ classes for all desired combinations of template 
// parameters).

// Weighted sums or single tiered scoring using site map points
BOOST_PYTHON_MODULE(_WeightedSumsScore)
{
  class_< ScoreRigidAlignments< NoTier1Score, WeightedSumsScore, 
                                rigid_align_t > >
    ("WeightedSumsScore", init<ModelSitemap*, const SearchParameters&>())

    .def("write_score_header", 
         &out_wrapper< NoTier1Score, WeightedSumsScore, rigid_align_t >
           ::write_header)
    .def("score_alignments", 
         &out_wrapper< NoTier1Score, WeightedSumsScore , rigid_align_t>
           ::score_aligns)
  ;
}
