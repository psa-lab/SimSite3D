#include <boost/python.hpp>
#include <BaseParameters.H>
using namespace boost::python;

BOOST_PYTHON_MODULE(_base_parameters)
{
  class_<SimSite3D::BaseParameters>("base_parameters", init< >())
    .def_readonly("install_dir", &SimSite3D::BaseParameters::install_dir)
    // While testing keep this next line
    .def_readwrite("load_surf_files", &SimSite3D::BaseParameters::load_surf_files)
    .def_readwrite("dbase_sites", &SimSite3D::BaseParameters::dbase_sites)
    .def_readwrite("dbase_ligs", &SimSite3D::BaseParameters::dbase_ligs)
    .def_readwrite("dbase_prots", &SimSite3D::BaseParameters::dbase_prots)
    .def_readwrite("diverse_sites", &SimSite3D::BaseParameters::diverse_sites)
    .def_readwrite("diverse_ligs", &SimSite3D::BaseParameters::diverse_ligs)
    .def_readwrite("proj_output", &SimSite3D::BaseParameters::proj_output)
    .def_readwrite("scratch_dir", &SimSite3D::BaseParameters::scratch_dir)
    .def_readwrite("require_min_npts", 
                   &SimSite3D::BaseParameters::require_min_npts)
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
