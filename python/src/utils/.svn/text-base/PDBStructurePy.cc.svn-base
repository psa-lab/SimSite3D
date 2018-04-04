#include <boost/python.hpp>
#include <PDBStructure.H>
using namespace boost::python;
using SimSite3D::PDBStructure;

BOOST_PYTHON_MODULE(_PDBStructure)
{
  class_<PDBStructure>("PDBStructure", init<const std::string>())
  ;
#if 0
  class_<SimSite3D::SimSite3DSearchParameters, bases<SimSite3D::BaseParameters> >
    ("parameters", init< > ())
    .def_readwrite("ext_score_method", 
		   &SimSite3D::SimSite3DSearchParameters::ext_score_method)
    .def_readwrite("model_file_name", 
		   &SimSite3D::SimSite3DSearchParameters::model_file_name)
    .def_readwrite("dbase_file_name", 
		   &SimSite3D::SimSite3DSearchParameters::dbase_file_name)
    .def_readwrite("ofname", &SimSite3D::SimSite3DSearchParameters::ofname)
    .def_readwrite("score_str", &SimSite3D::SimSite3DSearchParameters::score_str)
    .def_readwrite("normalize", &SimSite3D::SimSite3DSearchParameters::normalize)
    .def_readwrite("num_scores_to_keep", 
		   &SimSite3D::SimSite3DSearchParameters::num_scores_to_keep)
    .def_readwrite("score_cutoff", 
		   &SimSite3D::SimSite3DSearchParameters::score_cutoff)
    .def_readwrite("min_num_atoms", 
		   &SimSite3D::SimSite3DSearchParameters::min_num_atoms)
    .def_readwrite("dmetol", &SimSite3D::SimSite3DSearchParameters::dmetol)
    .def_readwrite("lsetol", &SimSite3D::SimSite3DSearchParameters::lsetol)
    .def_readwrite("ligand_rmsd", 
		   &SimSite3D::SimSite3DSearchParameters::ligand_rmsd)
    .def_readwrite("sitemap_rmsd", 
		   &SimSite3D::SimSite3DSearchParameters::sitemap_rmsd)
    .def_readwrite("write_ligands", 
		   &SimSite3D::SimSite3DSearchParameters::write_ligands)
    .def_readwrite("align_to_query", 
		   &SimSite3D::SimSite3DSearchParameters::align_to_query)
    .def_readwrite("num_rand_aligns", 
		   &SimSite3D::SimSite3DSearchParameters::num_rand_aligns)
    .def_readwrite("time_process", 
		   &SimSite3D::SimSite3DSearchParameters::time_process)
    .def_readwrite("allow_hphob_triangles", 
		   &SimSite3D::SimSite3DSearchParameters::allow_hphob_triangles)
  ;
#endif
}
