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
 * $Source: /psa/share/repository/pfizer_proj/src/search/ExternalScoringFunction.H,v $
 * $Revision: 1.4 $
 * $Author: vanvoor4 $
 * $Date: 2007-11-01 16:41:35 $
 * 
 * $Log: not supported by cvs2svn $
 * Revision 1.3  2007/08/21 19:21:30  vanvoor4
 * changed license header
 *
 * Revision 1.2  2007/06/13 17:21:46  vanvoor4
 * Added support for pre and post options with respect to ligand and
 * protein specifications on scoring function command line
 *
 * Revision 1.1  2007/06/07 20:53:45  vanvoor4
 * initial checkin
 *
 * 
 * 
 */

#ifndef _EXTERNAL_SCORING_FUNCTION_HEADER_FILE
#define _EXTERNAL_SCORING_FUNCTION_HEADER_FILE

#include <types.H>
#include <string>
#include <vector>

namespace ASCbase{

  typedef enum{
    LIGAND_FIRST,
    PROTEIN_FIRST,
  }prot_lig_order_t;

  //! Base class: wrapper and parser for external scoring functions
  class ExternalScoringFunction{
  public:
    //! Constructor
    /*!
     * @param cmdline Command line to run the chosen external protein-ligand
     *  scoring function with the variables $PROTEIN and $LIGAND denoting
     *  where to stuff the values for the proteins and ligands to score.
     */
    ExternalScoringFunction(std::string cmdline);

    //! Virtual destruction
    virtual ~ExternalScoringFunction() {;}

    //! Determine if calling score will actually do something 
    bool fail() const { return a_fail; }

    //! Score protein-ligand pair using chosen external scoring method
    /*! Score using chosen external scoring method -- force derived classes
     * to override with their version.
     * @param prot_path String holding the path of the protein PDB file
     * @param lig_path String holding the path of the ligand mol2 file
     * @param scores Pointer to vector holding the scores of the pair
     * @return True if prot_path and lig_path exist and are not directories,
     *         otherwise false
     */
    virtual bool score(std::string prot_path, std::string lig_path, 
                       std::vector<my_float_t> *scores) = 0;

    //! Get the name(s) of the external scoring function(s) used
    /*! To be overloaded by the derived classes
     * @param names Pointer to vector of strings holding the names of the 
     *        scoring function(s) in the same order as the scores stored in the 
     *        ext_scores vector in rigid_align_t
     */
    virtual void sf_names(std::vector<std::string> *names) = 0;

  protected:
    //! Run the chosen external scoring method
    /*! Run the chosen external scoring method 
     * @param prot_path String holding the path of the protein PDB file
     * @param lig_path String holding the path of the ligand mol2 file
     * @return True if prot_path and lig_path exist and are not directories,
     *         otherwise false
     */
    bool run(std::string prot_path, std::string lig_path);

    bool a_fail;              //!< True means unable to run or failure occured
  private:
    std::vector<std::string> cmdline_toks;
    std::vector<std::string>::const_iterator prot_str_iter;
    std::vector<std::string>::const_iterator lig_str_iter;

  };
}
#endif
