#include <boost/python.hpp>
#include <boost/python/iterator.hpp>
#include <ProtLigScore.H>

using namespace boost::python;
using SimSite3D::ProtLigScore;
using SimSite3D::PDBStructure;
using SimSite3D::mol2File;


BOOST_PYTHON_MODULE(_ProtLigScore)
{
  // Note: I changed the default value of the third parameter to true, I hope
  // the default value does not get changed back to false -- this parameter 
  // could go away since it introduces relatively little extra storage
  //
  // Use return_by_value since under the present design, speed is not an
  // issue when using this class and it saves messing around with 
  // to_python for C++ std::string
  // NOTE: we are using the default parameters for ProtLigScore -- save the
  // interactions and compute the "charge sums" on the ligand.  For that
  // reason DO NOT compute the charge sums on the ligand before calling
  // this constructor or the charges will be incorrect/inconsistent.
  class_<ProtLigScore>("prot_lig_score",
                       init<PDBStructure &, mol2File &>())
    .def("lig_act_strings", 
         range<return_value_policy<return_by_value> > 
           (&ProtLigScore::lig_act_strings_begin, 
            &ProtLigScore::lig_act_strings_end))
  ;
}

