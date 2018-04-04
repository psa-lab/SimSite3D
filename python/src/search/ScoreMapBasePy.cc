// This file has been generated by Py++.
// Modified by Jeff Van Voorst

#include <boost/python.hpp>
#include "__array_1.pypp.hpp"
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <ScoreMapBase.H>

namespace bp = boost::python;

static void
rigid_vec_append(std::vector< SimSite3D::rigid_align_t > & x,
                 SimSite3D::rigid_align_t const& v)
{
  x.push_back(v);
}

static void
rigid_vec_clear(std::vector< SimSite3D::rigid_align_t > & x)
{
  x.clear();
}

static int
rigid_vec_size(std::vector< SimSite3D::rigid_align_t > & x)
{
  return x.size();
} 

struct rigid_align_t_wrapper : SimSite3D::rigid_align_t, bp::wrapper< SimSite3D::rigid_align_t > {

    rigid_align_t_wrapper(SimSite3D::rigid_align_t const & arg )
    : SimSite3D::rigid_align_t( arg )
      , bp::wrapper< SimSite3D::rigid_align_t >(){
        // copy constructor
        
    }

    rigid_align_t_wrapper()
    : SimSite3D::rigid_align_t()
      , bp::wrapper< SimSite3D::rigid_align_t >(){
        // null constructor
        
    }

    static pyplusplus::containers::static_sized::array_1_t< double, 9>
    pyplusplus_R_wrapper( ::SimSite3D::rigid_align_t & inst ){
        return pyplusplus::containers::static_sized::array_1_t< double, 9>( inst.R );
    }

    static pyplusplus::containers::static_sized::array_1_t< double, 3>
    pyplusplus_T_wrapper( ::SimSite3D::rigid_align_t & inst ){
        return pyplusplus::containers::static_sized::array_1_t< double, 3>( inst.T );
    }

#if 0
    static pyplusplus::containers::static_sized::array_1_t< double, 3>
    pyplusplus_Q_T_wrapper( ::SimSite3D::rigid_align_t & inst ){
        return pyplusplus::containers::static_sized::array_1_t< double, 3>( inst.Q_T );
    }

    static pyplusplus::containers::static_sized::array_1_t< double, 3>
    pyplusplus_Trans_C_wrapper( ::SimSite3D::rigid_align_t & inst ){
        return pyplusplus::containers::static_sized::array_1_t< double, 3>( inst.Trans_C );
    }

    static pyplusplus::containers::static_sized::array_1_t< unsigned int, 3>
    pyplusplus_pts_idx_wrapper( ::SimSite3D::rigid_align_t & inst ){
        return pyplusplus::containers::static_sized::array_1_t< unsigned int, 3>( inst.pts_idx );
    }

    static pyplusplus::containers::static_sized::array_1_t< double, 3>
    pyplusplus_tri_params_wrapper( ::SimSite3D::rigid_align_t & inst ){
        return pyplusplus::containers::static_sized::array_1_t< double, 3>( inst.tri_params );
    }
#endif

};

BOOST_PYTHON_MODULE(_ScoreMapBase){
#if 0
// moved to basics/_stl_containers.cc
    { //::std::vector< double >
        typedef bp::class_< std::vector< double > > vector_less__double__greater__exposer_t;
        vector_less__double__greater__exposer_t vector_less__double__greater__exposer = vector_less__double__greater__exposer_t( "vector_less__double__greater_" );
        bp::scope vector_less__double__greater__scope( vector_less__double__greater__exposer );
        vector_less__double__greater__exposer.def( bp::vector_indexing_suite< ::std::vector< double >, true >() );
    }

    { //::std::vector< bool >
        typedef bp::class_< std::vector< bool > > vector_less__bool__greater__exposer_t;
        vector_less__bool__greater__exposer_t vector_less__bool__greater__exposer = vector_less__bool__greater__exposer_t( "vector_less__bool__greater_" );
        bp::scope vector_less__bool__greater__scope( vector_less__bool__greater__exposer );
        vector_less__bool__greater__exposer.def( bp::vector_indexing_suite< ::std::vector< bool >, true >() );
    }
#endif
    { //::std::vector< rigid_align_t >
        typedef bp::class_< std::vector< SimSite3D::rigid_align_t > > vector_less__rigid_align_t__greater__exposer_t;
        vector_less__rigid_align_t__greater__exposer_t vector_less__rigid_align_t__greater__exposer = vector_less__rigid_align_t__greater__exposer_t( "vector_less__rigid_align_t__greater_" );
        bp::scope vector_less__rigid_align_t__greater__scope( vector_less__rigid_align_t__greater__exposer );
        vector_less__rigid_align_t__greater__exposer.def( bp::vector_indexing_suite< ::std::vector< SimSite3D::rigid_align_t >, true >() );
        // Added by Jeff Van Voorst
        vector_less__rigid_align_t__greater__exposer.def("append", &rigid_vec_append, bp::with_custodian_and_ward<1,2>());
        vector_less__rigid_align_t__greater__exposer.def("clear", &rigid_vec_clear);
        vector_less__rigid_align_t__greater__exposer.def("size", &rigid_vec_size);

    }

    { //::SimSite3D::rigid_align_t
        typedef bp::class_< rigid_align_t_wrapper > rigid_align_t_exposer_t;
        rigid_align_t_exposer_t rigid_align_t_exposer = rigid_align_t_exposer_t( "rigid_align_t" );
        bp::scope rigid_align_t_scope( rigid_align_t_exposer );
        rigid_align_t_exposer.def( bp::self != bp::self );
        rigid_align_t_exposer.def( bp::self == bp::self );
        pyplusplus::containers::static_sized::register_array_1< double, 3 >( "__array_1_double_3" );
#if 0
        { //SimSite3D::rigid_align_t::Q_T [variable], type=my_float_t[3]
        
            typedef pyplusplus::containers::static_sized::array_1_t< double, 3> ( *array_wrapper_creator )( ::SimSite3D::rigid_align_t & );
            
            rigid_align_t_exposer.add_property( "Q_T"
                , bp::make_function( array_wrapper_creator(&rigid_align_t_wrapper::pyplusplus_Q_T_wrapper)
                                    , bp::with_custodian_and_ward_postcall< 0, 1 >() ) );
        }
#endif
        pyplusplus::containers::static_sized::register_array_1< double, 9 >( "__array_1_double_9" );
        { //SimSite3D::rigid_align_t::R [variable], type=my_float_t[9]
        
            typedef pyplusplus::containers::static_sized::array_1_t< double, 9> ( *array_wrapper_creator )( ::SimSite3D::rigid_align_t & );
            
            rigid_align_t_exposer.add_property( "R"
                , bp::make_function( array_wrapper_creator(&rigid_align_t_wrapper::pyplusplus_R_wrapper)
                                    , bp::with_custodian_and_ward_postcall< 0, 1 >() ) );
        }
        { //SimSite3D::rigid_align_t::T [variable], type=my_float_t[3]
        
            typedef pyplusplus::containers::static_sized::array_1_t< double, 3> ( *array_wrapper_creator )( ::SimSite3D::rigid_align_t & );
            
            rigid_align_t_exposer.add_property( "T"
                , bp::make_function( array_wrapper_creator(&rigid_align_t_wrapper::pyplusplus_T_wrapper)
                                    , bp::with_custodian_and_ward_postcall< 0, 1 >() ) );
        }
#if 0
        { //SimSite3D::rigid_align_t::Trans_C [variable], type=my_float_t[3]
        
            typedef pyplusplus::containers::static_sized::array_1_t< double, 3> ( *array_wrapper_creator )( ::SimSite3D::rigid_align_t & );
            
            rigid_align_t_exposer.add_property( "Trans_C"
                , bp::make_function( array_wrapper_creator(&rigid_align_t_wrapper::pyplusplus_Trans_C_wrapper)
                                    , bp::with_custodian_and_ward_postcall< 0, 1 >() ) );
        }
#endif
        rigid_align_t_exposer.def_readwrite( "ext_scores", &SimSite3D::rigid_align_t::ext_scores );
        rigid_align_t_exposer.def_readwrite( "frag_atoms_flags", &SimSite3D::rigid_align_t::frag_atoms_flags );
        rigid_align_t_exposer.def_readwrite( "match_print", &SimSite3D::rigid_align_t::match_print );
        pyplusplus::containers::static_sized::register_array_1< unsigned int, 3 >( "__array_1_unsigned_int_3" );
// Hmm, maybe Gentoo and Enterprise Linux have defined size_t to be 
// different types -- lets just grab unsigned int since we don't necessarily
// need the possible additional space of size_t
#if 0
        { //SimSite3D::rigid_align_t::pts_idx [variable], type=size_t[3]
        
            typedef pyplusplus::containers::static_sized::array_1_t< unsigned int, 3> ( *array_wrapper_creator )( ::SimSite3D::rigid_align_t & );
            
            rigid_align_t_exposer.add_property( "pts_idx"
                , bp::make_function( array_wrapper_creator(&rigid_align_t_wrapper::pyplusplus_pts_idx_wrapper)
                                    , bp::with_custodian_and_ward_postcall< 0, 1 >() ) );
        }
        rigid_align_t_exposer.def_readwrite( "score", &SimSite3D::rigid_align_t::score );
        rigid_align_t_exposer.def_readwrite( "terms", &SimSite3D::rigid_align_t::terms );
        { //SimSite3D::rigid_align_t::tri_params [variable], type=my_float_t[3]
        
            typedef pyplusplus::containers::static_sized::array_1_t< double, 3> ( *array_wrapper_creator )( ::SimSite3D::rigid_align_t & );
            
            rigid_align_t_exposer.add_property( "tri_params"
                , bp::make_function( array_wrapper_creator(&rigid_align_t_wrapper::pyplusplus_tri_params_wrapper)
                                    , bp::with_custodian_and_ward_postcall< 0, 1 >() ) );
        }
#endif
	// Added by Jeff
	rigid_align_t_exposer.def_readwrite("frag_file", &SimSite3D::rigid_align_t::frag_file );
        rigid_align_t_exposer.def_readwrite( "hb_caps_score", &SimSite3D::rigid_align_t::hb_caps_score );
        rigid_align_t_exposer.def_readwrite( "hb_caps_match_print", &SimSite3D::rigid_align_t::hb_caps_match_print );
#if 0
        rigid_align_t_exposer.def_readwrite("affiscore",  &SimSite3D::rigid_align_t::affiscore);
        rigid_align_t_exposer.def_readwrite("orientscore",  &SimSite3D::rigid_align_t::orientscore);
        rigid_align_t_exposer.def_readwrite("ligand_efficiency",  &SimSite3D::rigid_align_t::ligand_efficiency);
#endif
    }
}
