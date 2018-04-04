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
 * $Source: /psa/share/repository/pfizer_proj/src/basics/AtomTypes.H,v $
 * $Revision: 1.1 $
 * $Author: vanvoor4 $
 * $Date: 2008-05-13 17:06:55 $
 * 
 * $Log: not supported by cvs2svn $
 * 
 */

#ifndef ATOM_TYPES_HEADER_INCLUDED
#define ATOM_TYPES_HEADER_INCLUDED

#include <orbitals.H>
#include <types.H>

namespace ASCbase{

  //! The majority of atom types
  typedef enum atom_type_enum{
    NULL_ATOM_TYPE = 0,
    UNKNOWN_ATOM,
    CA,          //!< alpha carbon
    CALCIUM,     //!< calcium is written out to reduce conflicts with alpha C
    CB,          //!< beta carbon
    C,           //!< carbon
    O,           //!< oxygen
    N,           //!< nitrogen
    CADMIUM,     //!< cadmium is written out to reduce conflicts with delta C
    CD,          //!< delta carbon
    CD1,         //!< first delta carbon
    CD2,         //!< second delta carbon
    CE,          //!< epsilon carbdon
    CE1,         //!< first epsilon carbon
    CE2,         //!< second epsilon carbon
    CE3,         //!< third epsilon carbon
    CG,          //!< gamma carbond
    CG1,         //!< first gamma carbond
    CG2,         //!< second gamma carbond
    CH2,
    CH3,
    CZ,
    CZ2,
    CZ3,
    ND1,
    ND2,
    NE,
    NE1,
    NE2,
    NH1,
    NH2,
    NZ,
    OD1,
    OD2,
    OE,    //!< For PCA
    OE1,
    OE2,
    OG,
    OG1,
    OH,
    OXT,
    SD,
    SG,
    P,       //!< Phopshorus
    S,
    H,
    D,
    SI,
    CO,
    AG,
    SB, 
    SE,
    SN,
    SR, 
    TE,
    AU,
    NI,
    IN,
    B,
    BA,
    BE,
    BI,
    BR,
    TL,
    OS,
    CL,
    HG,
    LA,
    K,
    FE,
    MO,
    I,
    AS,
    U,
    F,
    RE,
    RU,
    ER,      //!< Erbium
    LI,
    MN,
    PD,
    CU,
    AL,
    MG,
    NA,
    CR,
    CS,
    TI,
    ZN,     //!< Zinc
    RB,
    V,
    RH,
    GE,
    ZR,
    GA,
    Y_,
    AD1,      //!< First ambigious atom (is it O or N for ASN?)
    AD2,      //!< Second ambigious atom (is it O or N for ASN?)
    AE1,      //!< First ambigious atom (is it O or N for GLN?)
    AE2,      //!< Second ambigious atom (is it O or N for GLN?)
    HETATM,
    DU,       //!< dummy atom 
    ANY       //!< Used mainly for rules which do not depend on atom name/type
  }atom_type;

  //! Atom level
  typedef enum atom_level_enum{
    UNKOWN_LEVEL = 0,
    ALPHA,
    BETA,
    GAMMA,
    DELTA,
    EPSILON,
    ZETA,
    ETA,
    LEVEL_ENUM_END
  }atom_level_type;

  //! Silly little struct to help convert from string to atom_type
  typedef struct{
    atom_type atom;
    std::string name;
  }atom_conv_type;

  typedef enum{
    WARN_UNKNOWN_ATOMS,
    SILENTLY_IGNORE_UNKNOWN_ATOMS
  }atoms_warn_type;

  //! a PDB atom for a unknown residue, a PDB hetatm, a PDB H atom, or a mol2 
  //! atom
  /*! The Van der Waals radii of a number of metals are the recommended values
   * given in table 12 of S. S. Batsanov, "Van der Waals radii of elements from
   * the data of structural inorganic chemistry", Russian Chemical Bulletin,
   * V. 44, 1995.
   */
  typedef struct{
    atom_type atom;             //!< ASCbase atom identifier (enum)
    orbit_type orbit;           //!< ASCbase atom orbital identifier (enum)
    my_float_t vdw_radius;      //!< Van der Waals radius
  }atom_info_t;

  //! A data class to look up the "hard coded" data for a particular atom based
  //! on its name and orbital.
  class AtomTypes{
  public:
    //! Get the the atom's Van der Waals radius based on the atom's name and
    //! orbital
    /*!
     * @param atom Name type of the atom
     * @param orbit Orbital type of the atom
     * @return constant pointer to the atom's atom_info_t
     */
    static const atom_info_t*    
    lookup_vdw_rad(const atom_type atom, const orbit_type orbit = DEFAULT_ORBIT)
    {
       if(vdw_table.empty()) build_vdw_table();
       std::multimap<atom_type, const atom_info_t*>::const_iterator tbl_iter;
       tbl_iter = vdw_table.find(atom);
       if(tbl_iter == vdw_table.end()) return 0; 

       for(;tbl_iter != vdw_table.end() && tbl_iter->first == atom; ++tbl_iter){
         if(tbl_iter->second->orbit == orbit) return tbl_iter->second;
       }
       return 0;
    }

  private:
    //! Setup the map used to search for the atom's structure
    static void build_vdw_table();

    static std::multimap<atom_type, const atom_info_t*> vdw_table;
    static const uint A_atoms_array_size;
    static const atom_info_t A_atoms_array[];
  };
}
#endif
