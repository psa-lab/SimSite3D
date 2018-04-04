#include <boost/python.hpp>
#include <mol2File.H>
using namespace boost::python;
using ASCbase::mol2File;

BOOST_PYTHON_MODULE(_mol2File)
{
  class_<mol2File>("mol2", init<const std::string>())
    .def("write", &mol2File::write)
  ;
#if 0
  class_<ASCbase::ASCbaseSearchParameters, bases<ASCbase::BaseParameters> >
    ("parameters", init< > ())
    .def_readwrite("ext_score_method", 
		   &ASCbase::ASCbaseSearchParameters::ext_score_method)
    .def_readwrite("model_file_name", 
		   &ASCbase::ASCbaseSearchParameters::model_file_name)
    .def_readwrite("dbase_file_name", 
		   &ASCbase::ASCbaseSearchParameters::dbase_file_name)
    .def_readwrite("ofname", &ASCbase::ASCbaseSearchParameters::ofname)
    .def_readwrite("score_str", &ASCbase::ASCbaseSearchParameters::score_str)
    .def_readwrite("normalize", &ASCbase::ASCbaseSearchParameters::normalize)
    .def_readwrite("num_scores_to_keep", 
		   &ASCbase::ASCbaseSearchParameters::num_scores_to_keep)
    .def_readwrite("score_cutoff", 
		   &ASCbase::ASCbaseSearchParameters::score_cutoff)
    .def_readwrite("min_num_atoms", 
		   &ASCbase::ASCbaseSearchParameters::min_num_atoms)
    .def_readwrite("dmetol", &ASCbase::ASCbaseSearchParameters::dmetol)
    .def_readwrite("lsetol", &ASCbase::ASCbaseSearchParameters::lsetol)
    .def_readwrite("ligand_rmsd", 
		   &ASCbase::ASCbaseSearchParameters::ligand_rmsd)
    .def_readwrite("sitemap_rmsd", 
		   &ASCbase::ASCbaseSearchParameters::sitemap_rmsd)
    .def_readwrite("write_ligands", 
		   &ASCbase::ASCbaseSearchParameters::write_ligands)
    .def_readwrite("align_to_query", 
		   &ASCbase::ASCbaseSearchParameters::align_to_query)
    .def_readwrite("num_rand_aligns", 
		   &ASCbase::ASCbaseSearchParameters::num_rand_aligns)
    .def_readwrite("time_process", 
		   &ASCbase::ASCbaseSearchParameters::time_process)
    .def_readwrite("allow_hphob_triangles", 
		   &ASCbase::ASCbaseSearchParameters::allow_hphob_triangles)
  ;
#endif
}
