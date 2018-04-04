#include <boost/python.hpp>
#include "../utils/stream_basicsPy.hh"
#include <SearchParameters.H>

using namespace boost::python;
using ASCbase::SearchParameters;

class params_out_wrap{
public:
  static void
  report(SearchParameters &P, stl_ofstream& out)
  {
    P.report_parameters(out.get());
  }
};
#if 0
// template this
template<T>
class out_wrapper_t{
public:
  static void
  write_header(T& S, stl_ofstream& out)
  {
    S.write_score_header(out.get());
  }

  static void
  score_aligns(T& S, ASCbase::rigid_align_vec& aligns,
               Sitemap* search, stl_ofstream& out)
  {
    S.score_alignments(aligns, search, out.get());
  }
};
#endif


BOOST_PYTHON_MODULE(_parameters)
{
  class_<SearchParameters, bases<ASCbase::BaseParameters> >
    ("parameters", init< > ())
    .def("report", &params_out_wrap::report)
    // We won't be using this variable from the python interface
    .def_readwrite("ext_score_method", &SearchParameters::ext_score_method)

    .def_readwrite("model_file_name", &SearchParameters::model_file_name)

    // The next three variables are not needed by the python interface once
    // we have finished time testings
    .def_readwrite("ofname", &SearchParameters::ofname)
    .def_readwrite("dbase_file_name", &SearchParameters::dbase_file_name)
    .def_readwrite("score_str", &SearchParameters::score_str)

    .def_readwrite("db_index_fname", &SearchParameters::db_index_fname)
    .def_readwrite("dbstart", &SearchParameters::db_start)
    .def_readwrite("dbstop", &SearchParameters::db_stop)
    .def_readwrite("normalize", &SearchParameters::normalize)
    .def_readwrite("num_scores_to_keep", &SearchParameters::num_scores_to_keep)
    .def_readwrite("score_cutoff", &SearchParameters::score_cutoff)
    .def_readwrite("min_num_atoms", &SearchParameters::min_num_atoms)
    .def_readwrite("dmetol", &SearchParameters::dmetol)
    .def_readwrite("lsetol", &SearchParameters::lsetol)
    .def_readwrite("ligand_rmsd", &SearchParameters::ligand_rmsd)
    .def_readwrite("sitemap_rmsd", &SearchParameters::sitemap_rmsd)
    .def_readwrite("write_ligands", &SearchParameters::write_ligands)
    .def_readwrite("align_to_query", &SearchParameters::align_to_query)
    .def_readwrite("num_rand_aligns", &SearchParameters::num_rand_aligns)
    .def_readwrite("time_process", &SearchParameters::time_process)
    .def_readwrite("allow_hphob_triangles", 
		   &SearchParameters::allow_hphob_triangles)
    .def_readwrite("do_internal_prot_lig_score", 
                   &SearchParameters::do_internal_prot_lig_score)
    .def_readwrite("fine_tune_all_tier2_alignments",
                   &SearchParameters::fine_tune_tier2_alignments)
    .def_readwrite("fine_tune_best_tier2_alignment",
                   &SearchParameters::fine_tune_best_tier2_alignment)
    .def_readwrite("save_rigid_scores",
                   &SearchParameters::save_rigid_scores)
    .def_readwrite("max_corr_surf_pt_dist",
                   &SearchParameters::max_corr_surf_pt_dist)
  ;
}
