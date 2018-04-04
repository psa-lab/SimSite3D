/******************************************************************************
 * Copyright (c) 2007, Michigan State University (MSU) Board of Trustees.
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
 * $Source: /psa/share/repository/pfizer_proj/src/basics/PDB_metals.H,v $
 * $Revision: 1.4 $
 * $Author: vanvoor4 $
 * $Date: 2008-02-26 18:19:23 $ 
 * 
 * $Log: not supported by cvs2svn $
 * Revision 1.3  2008/01/04 20:54:24  vanvoor4
 * Added some checking enums
 *
 * Revision 1.2  2007/12/20 01:57:40  vanvoor4
 * Removed class prefix
 *
 * Revision 1.1  2007/12/17 21:10:15  vanvoor4
 * Initial checkin
 * 
 * 
 */

#ifndef PDB_METALS_HEADER_INCLUDED
#define PDB_METALS_HEADER_INCLUDED

#include <AtomTypes.H>
#include <point.H>

namespace ASCbase{

  typedef enum{
    WARN_UKNOWN_METALS,
    SILENTLY_IGNORE_UKNOWN_METALS
  }metal_warn_type;

  //! Information about supported PDB metal types
  typedef struct{
    std::string pdb_metal_str;    //!< PDB 4 character atom name string
    atom_type metal_name;         //!< Internal atom name 
    interactionType act_type;     //!< "Type" of metal interaction
    my_float_t interact_pt_rad;   //!< Template point distance from metal center
    my_float_t min_act_dist; //!< Min distance between metal and acceptor atom
    my_float_t max_act_dist; //!< Max distance between metal and acceptor atom
  }pdb_metal_info_t;

  //! Simple class to facilitate lookup of metals using 4 char PDB metal names
  class PDB_metals{
  public:
  
    //! Look up metal info for the given 4 char PDB metal name
    /*! NOTE:  because of the small number of elements in the table (at this 
     * point 9 metal types) and the relatively low occurance of metals in PDB
     * files, just use a linear scan.
     *
     * @param metal_name 4 char PDB metal name
     * @param warn_if_unknown Should we use warn() if a metal is not in the 
     *        table?
     * @return Const ptr to metal info struct, otherwise 0
     */
    static const pdb_metal_info_t* 
    lookup(const std::string metal_name, const std::string metal_res,
           metal_warn_type warn_if_unknown = WARN_UKNOWN_METALS);

    //! Look up metal info for the atom_type
    /*! NOTE:  because of the small number of elements in the table (at this 
     * point 9 metal types) and the relatively low occurance of metals in PDB
     * files, just use a linear scan.
     *
     * @param metal enum atom_type for the desired metal
     * @param warn_if_unknown Should we use warn() if a metal is not in the 
     *        table?
     * @return Const ptr to metal info struct, otherwise 0
     */
    static const pdb_metal_info_t* 
    lookup(const atom_type metal, 
           metal_warn_type warn_if_unknown = WARN_UKNOWN_METALS);

    static const my_float_t MAX_METAL_1_SQUARED_LENGTH;
    static const my_float_t MAX_METAL_2_SQUARED_LENGTH;
    static const my_float_t MIN_METAL_1_SQUARED_LENGTH;
    static const my_float_t MIN_METAL_2_SQUARED_LENGTH;

  private:
    static void
    set_metal_res_names();

    static std::map<atom_type, std::vector<std::string> > A_metal_res_names;
    static const uint A_metal_array_size;   //!< Number of els in metal array
    static const pdb_metal_info_t A_metal_array[]; //!< Array holding metal info
    static const std::string A_fname;   //!< Name of the source file
  };
}  
#endif
