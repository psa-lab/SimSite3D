/******************************************************************************
 * Copyright (c) 2006,2007, Michigan State University (MSU) Board of Trustees.
 *   All rights reserved.
 *
 * This file is part of the SimSite3D Software project.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 * Authors: Jeffrey Van Voorst, jeff.vanvoorst@gmail.com
 *          Leslie Kuhn, Ph.D., KuhnL@msu.edu 
 *****************************************************************************/
  
/*
 * $Source: /psa/share/repository/pfizer_proj/src/basics/atom.H,v $
 * $Revision: 1.5 $
 * $Author: vanvoor4 $
 * $Date: 2008/04/17 18:34:51 $
 * 
 * $Log: atom.H,v $
 * Revision 1.5  2008/04/17 18:34:51  vanvoor4
 * Added some statics to the atom class for NULL iterator
 * types
 *
 * Revision 1.4  2008/02/26 19:10:59  vanvoor4
 * Changed some comments.
 *
 * Revision 1.3  2007/12/17 21:07:15  vanvoor4
 * Added a filed for hetero versus atom
 *
 * Revision 1.2  2007/08/21 16:38:21  vanvoor4
 * Moved the checking functions such as is_water(), etc to this class
 *
 * Revision 1.1  2007/08/21 15:30:01  vanvoor4
 * Initial checkin
 *
 *
 * 
 * 
 */ 

#ifndef ATOM_HEADER_INCLUDED
#define ATOM_HEADER_INCLUDED

#include <PDB_residues.H>
#include <point_storage.H>

namespace SimSite3D{


//! Generic atom type which uses the point_t class to store its position
class atom_t : public point_t{
public:
  atom_t(alloc_t a = ALLOC_POSITION) : point_t(a)
  {
    atom_num = 0;
    res_num = 0;
    name = NULL_ATOM_TYPE;
    res = NULL_RESIDUE;
    chainID = ' ';
    altLoc = ' ';
    iCode = ' ';
    occupancy = 0.00;
    tempFactor = 0.0;
    orbit = UNKNOWN_ORBIT;
    vdw_radius = 1000.0;
    is_hetero = false;
    act_type = UNKNOWN_INTERACTION;
    name_str = "";
    res_str = "";
    charge = 0.0;
    orig_charge = 0.0;
    subst_id = 0;
    hydro = 82;
  }

  // Must be a faster way to copy, but copying should not play a sigificant 
  // role in the speed for SimSite3D
  atom_t(const atom_t& a) : point_t(a)
  {
    do_copy(a);
  }

  const atom_t& 
  operator=(const atom_t& a)
  {
    if(&a != this){
      do_copy(a);
      point_t::operator= (a);
    }
    return *this; 
  }

  ~atom_t() { ; }

  //! Allows sorting of an atom_vec_t by atom number
  static bool 
  cmp(const atom_t &A, const atom_t &B)
  { return A.atom_num < B.atom_num; }

  //! Allows flexibility in how waters are designated
  const bool
  is_water() const
  {
    if(res == HOH) return true;
    else return false;
  } 

  //! Allows flexibility in how metals are designated
  const bool
  is_metal() const
  { return ( (res == PDB_METAL) ? true : false); }

  const bool
  is_nitrogen() const
  {
     if(name == N || name == ND1 || name == ND2 || name == NE || name == NE1
       || name == NE2 || name == NH1 || name == NH2 || name == NZ) 
       return true;
     return false;
  } 

  const bool
  is_oxygen() const
  {
    if(name == O || name == OD1 || name == OD2 || name == OE1 || name == OE2
       || name == OG || name == OG1 || name == OH) return true;
    return false;
  } 

  const bool
  is_hydrophobic() const
  {
    if(name == C || name == CA || name == CB || name == CD || name == CD1
       || name == CD2 || name == CE || name == CE1 || name == CE2 
       || name == CE3 || name == CG || name == CG1 || name == CG2 
       || name == CH2 || name == CH3 || name == CZ || name == CZ2 
       || name == CZ3 || name == S) return true;
    return false; 
  }

  const bool
  is_hbonder() const
  {
    if(res != PRO && is_nitrogen()) return true;
    else if(is_oxygen()) return true;
    return false;
  }

  friend std::ostream& operator<<(std::ostream &out, const atom_t &a)
  {
    // Cheap hack for mol2 versus PDB
    if(a.res_str.length() == 0){
      if(a.name_str.length()) out << a.name_str;
      else{
        if(a.chainID != ' ') out << "(" << a.chainID << ")";
        out << PDB_residues::residue_to_string(a.res) << a.res_num;
        if(a.iCode != ' ') out << a.iCode;
        out << " " << PDB_residues::atom_to_string(a.name);
      }
    }else{
      if(a.chainID != ' ') out << "(" << a.chainID << ")";
      out << a.res_str << a.res_num;
      if(a.iCode != ' ') out << a.iCode;
      out << " " << a.name_str;
    }
    return out;
  }

  uint atom_num;       //!< Atom (serial) number
  uint res_num;        //!< Residue number
  atom_type name;      //!< Atom type
  residue_type res;    //!< Residue type
  char chainID;        //!< Chain id
  char altLoc;         //!< Alterate location code
  char iCode;          //!< Insertion code
  my_float_t occupancy; //!< Occupancy value
  my_float_t tempFactor; //!< Temperature factor
  orbit_type orbit;      //!< Atom hybridization/orbital
  my_float_t vdw_radius; //!< Atom's Van der Waals radius
  std::string name_str;  //!< For PDB HETATMs', PDB Hs' and mol2 atoms' names
  std::string res_str;   //!< For PDB HETATMs
  bool is_hetero;        //!< True implies PDB HETATM, otherwise false 
  interactionType act_type; //!< Donor, acceptor, doneptor, hphob, or nothing
  my_float_t charge;     //!< Protein charges + SimSite3D summed atomic charges on ligands
  my_float_t orig_charge;  //!< Hold the ligand atomic charges as read in

  int subst_id; //!< mol2 junk -- ignored for now
  std::string subst_name; //!< mol2 junk -- ignored for now
  int hydro;    //!< Hydrophobicity value

  static const atom_t NULL_ATOM;  //!< the "standard" null atom
  static const point_storage<atom_t> NULL_ATOM_VECTOR;  //!< the "standard" null atom vector
  static const point_storage<atom_t>::const_iterator NULL_ATOM_VCI;  //!< the "standard" null atom vector constant iterator

private:
  void
  do_copy(const atom_t& a)
  {
    atom_num = a.atom_num;
    res_num = a.res_num;
    name = a.name;
    res = a.res;
    chainID = a.chainID;
    altLoc = a.altLoc;
    iCode = a.iCode;
    occupancy = a.occupancy;
    tempFactor = a.tempFactor;
    orbit = a.orbit;
    vdw_radius = a.vdw_radius;
    name_str = a.name_str;
    res_str = a.res_str;
    is_hetero = a.is_hetero;
    act_type = a.act_type;
    charge = a.charge;
    orig_charge = a.orig_charge;
    subst_id = a.subst_id;
    subst_name = a.subst_name;
    hydro = a.hydro;
  };

};

typedef point_storage<atom_t> atom_vec_t;
typedef atom_vec_t::iterator atom_vi;
typedef atom_vec_t::const_iterator atom_vci;
typedef std::map<uint, atom_vci> atom_map_t;


} // end of namespace SimSite3D
#endif
