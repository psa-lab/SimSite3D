/******************************************************************************
 * Copyright (c) 2006,2007, Michigan State University (MSU) Board of Trustees.
 *   All rights reserved.
 *
 * This file is part of the ASCbase Software project.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 * Authors: Jeffrey Van Voorst, vanvoor4@msu.edu
 *          Leslie Kuhn, Ph.D., KuhnL@msu.edu 
 *****************************************************************************/

/*
 * $Source: /psa/share/repository/pfizer_proj/src/basics/mol2_atoms.H,v $
 * $Revision: 1.5 $
 * $Author: vanvoor4 $
 * $Date: 2008-02-26 19:09:24 $
 * 
 * $Log: not supported by cvs2svn $
 * Revision 1.4  2007/10/12 17:12:11  vanvoor4
 * Added S.o as an alias for S.O.  Also, given the Tripos Mol2
 * file format pdf, S.O and S.O2 are the given atom types.
 *
 * Revision 1.3  2007/10/11 17:22:30  vanvoor4
 * S.o2 is more likely to be encountered than S.O2 not sure why
 * S.O2 was here in first place
 *
 * Revision 1.2  2007/10/11 16:03:04  vanvoor4
 * Added ASCbase namespace and changed AM to AMIDE
 *
 * Revision 1.1  2007/08/21 15:38:14  vanvoor4
 * Initial checkin
 *
 *
 * 
 * 
 */

#ifndef MOL2_ATOMS_HEADER_INCLUDED
#define MOL2_ATOMS_HEADER_INCLUDED

#include <AtomTypes.H>

namespace ASCbase{

  //!<Sybil bond types
  typedef enum mol2_bond_enum{
    UNKNOWN_BOND = 0,
    SINGLE_BOND,
    DOUBLE_BOND,
    TRIPLE_BOND,
    AMIDE_BOND,
    AROMATIC_BOND,
    DUMMY_BOND,
    NOT_CONNECTED
  }bond_type;
  
  typedef struct{
    std::string bond_type_str;
    bond_type type;
  }mol2_bond_info_t;
  
  const uint mol2_bond_table_size = 8;
  const mol2_bond_info_t mol2_bond_table[] = {
    { "un", UNKNOWN_BOND  },
    { "1" , SINGLE_BOND   },
    { "2" , DOUBLE_BOND   },
    { "3" , TRIPLE_BOND   },
    { "am", AMIDE_BOND    },
    { "ar", AROMATIC_BOND },
    { "du", DUMMY_BOND    },
    { "nc", NOT_CONNECTED }
  };
  
  typedef struct{
    std::string atom_type_str;
    atom_type atom;
    orbit_type orbit;
    my_float_t vdw_radii;
  }mol2_atom_info_t;
  
  //! MOL2 file van der Waals Radii Definitions 
  /*! The radii are the same as used in SLIDE 3.0.2
   */
  const uint mol2_atom_table_size = 79;
  const mol2_atom_info_t mol2_atom_table[] = {
    { "Al"   , AL, DEFAULT_ORBIT, 1.500 },
    { "Ag"   , AG, DEFAULT_ORBIT, 1.550 },
    { "As"   , AS, DEFAULT_ORBIT, 0.830 },
    { "Au"   , AU, DEFAULT_ORBIT, 1.900 },
    { "B"    , B,  DEFAULT_ORBIT, 1.700 },
    { "Ba"   , BA, DEFAULT_ORBIT, 3.100 },
    { "Be"   , BE, DEFAULT_ORBIT, 1.900 },
    { "Bi"   , BI, DEFAULT_ORBIT, 1.800 },
    { "Br"   , BR, DEFAULT_ORBIT, 1.950 },
    { "C.1"  , C , SP1, 1.850 },
    { "C.2"  , C , SP2, 1.850 },
    { "C.3"  , C , SP3, 1.800 },
    { "C.cat", C , CAT, 1.800 },
    { "C.ar" , C , AR , 1.800 },
    { "C"    , C , DEFAULT_ORBIT, 1.850 },
    { "Ca"   , CALCIUM, DEFAULT_ORBIT, 1.600 },
    { "Cd"   , CADMIUM, DEFAULT_ORBIT, 1.750 }, // vdw could be mixed up with CD
    { "Ce"   , CE, DEFAULT_ORBIT, 1.860 },
    { "Cl"   , CL, DEFAULT_ORBIT, 2.030 },
    { "Co"   , CO, DEFAULT_ORBIT, 1.130 },
    { "Cr"   , CR, DEFAULT_ORBIT, 1.130 },
    { "Cs"   , CS, DEFAULT_ORBIT, 3.010 },
    { "Cu"   , CU, DEFAULT_ORBIT, 1.150 },
    { "Er"   , ER, DEFAULT_ORBIT, 1.590 },
    { "F"    , F , DEFAULT_ORBIT, 1.550 },
    { "Fe"   , FE, DEFAULT_ORBIT, 1.950 },
    { "Ga"   , GA, DEFAULT_ORBIT, 1.550 },
    { "Ge"   , GE, DEFAULT_ORBIT, 2.720 },
    { "H"    , H , DEFAULT_ORBIT, 1.000 },
    { "Hg"   , HG, DEFAULT_ORBIT, 2.000 },
    { "I"    , I , DEFAULT_ORBIT, 2.350 },
    { "In"   , IN, DEFAULT_ORBIT, 1.900 },
    { "K"    , K , DEFAULT_ORBIT, 2.390 },
    { "La"   , LA, DEFAULT_ORBIT, 1.830 },
    { "Li"   , LI, DEFAULT_ORBIT, 1.220 },
    { "Mg"   , MG, DEFAULT_ORBIT, 1.500 },
    { "Mn"   , MN, DEFAULT_ORBIT, 1.190 },
    { "Mo"   , MO, DEFAULT_ORBIT, 1.750 },
    { "N.1"  , N , SP1, 1.750 },
    { "N.2"  , N , SP2, 1.750 },
    { "N.3"  , N , SP3, 1.850 },
    { "N.4"  , N , SP4, 1.850 },
    { "N.pl3", N , PL3, 1.750 },
    { "N.ar" , N , AR , 1.750 },
    { "N.am" , N , AMIDE, 1.750 },
    { "N"    , N , DEFAULT_ORBIT, 1.750 },
    { "Na"   , NA, DEFAULT_ORBIT, 1.600 },
    { "Ni"   , NI, DEFAULT_ORBIT, 1.240 },
    { "O.3"  , O , SP3, 1.650 },
    { "O.2"  , O , SP2, 1.600 },
    { "O.co2", O , CO2, 1.600 },
    { "O"    , O , DEFAULT_ORBIT, 1.600 },
    { "Os"   , OS, DEFAULT_ORBIT, 1.580 },
    { "P.3"  , P , SP3, 2.100 },  // vdw radius may be off
    { "P"    , P , DEFAULT_ORBIT, 2.100 },
    { "Pd"   , PD, DEFAULT_ORBIT, 1.440 },
    { "Rb"   , RB, DEFAULT_ORBIT, 2.650 },
    { "Re"   , RE, DEFAULT_ORBIT, 1.300 },
    { "Rh"   , RH, DEFAULT_ORBIT, 1.220 },
    { "Ru"   , RU, DEFAULT_ORBIT, 1.200 },
    { "S.3"  , S , SP3, 2.000 },
    { "S.2"  , S , SP2, 2.000 },
    { "S.o"  , S , O_ORBITAL  , 2.000 },
    { "S.O"  , S , O_ORBITAL  , 2.000 },
    { "S.o2" , S , O2 , 2.000 },
    { "S.O2" , S , O2 , 2.000 },
    { "S"    , S , DEFAULT_ORBIT, 2.000 },
    { "Si"   , SI, DEFAULT_ORBIT, 2.200 },
    { "Sn"   , SN, DEFAULT_ORBIT, 2.300 },
    { "Sr"   , SR, DEFAULT_ORBIT, 2.700 },
    { "Te"   , TE, DEFAULT_ORBIT, 2.000 },
    { "Ti"   , TI, DEFAULT_ORBIT, 1.950 },
    { "Tl"   , TL, DEFAULT_ORBIT, 1.710 },
    { "U"    , U , DEFAULT_ORBIT, 1.750 },
    { "V"    , V , DEFAULT_ORBIT, 1.060 },
    { "Y_"   , Y_, DEFAULT_ORBIT, 1.610 },
    { "Zn"   , ZN, DEFAULT_ORBIT, 1.150 },
    { "Zr"   , ZR, DEFAULT_ORBIT, 1.420 },
    { "Un"   , UNKNOWN_ATOM, UNKNOWN_ORBIT, 2.000 }
  };
}  
#endif
