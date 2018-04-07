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
 * $Source: /psa/share/repository/pfizer_proj/src/search/DrugScoreInterface.H,v $
 * $Revision: 1.4 $
 * $Author: vanvoor4 $
 * $Date: 2007-11-01 16:39:47 $
 * 
 * $Log: not supported by cvs2svn $
 * Revision 1.3  2007/08/21 19:20:44  vanvoor4
 * Changed license header
 *
 * Revision 1.2  2007/06/13 17:22:39  vanvoor4
 * Added support for pre and post options with respect to ligand and
 * protein specifications on scoring function command line
 *
 * Revision 1.1  2007/06/07 20:53:16  vanvoor4
 * initial checkin
 *
 * 
 * 
 */

#ifndef _DRUGSCORE_INTERFACE_HEADER_FILE
#define _DRUGSCORE_INTERFACE_HEADER_FILE

#include <ExternalScoringFunction.H>

namespace SimSite3D{

  //! Base class: wrapper and parser for external scoring functions
  class DrugScoreInterface : public ExternalScoringFunction{
  public:
    //!  
    /*!
     */
    DrugScoreInterface(std::string cmdline);

    //! Do nothing destructor
    ~DrugScoreInterface() {;}

    //! Wrapper to ExternalScoringFunction::run and parse output score file
    /*! Currently is operating under the assumption that PAIR_10 is enough.
     * A change in options would cause the DrugScore output files to have a
     * different prefix.  Assume the .scr file holds the score only
     * 
     * @param prot_path String holding the path of the protein PDB file
     * @param lig_path String holding the path of the ligand mol2 file
     * @oaram scores Pointer to vector to hold the score from DrugScore
     * @return True if DrugScore seemed to work, otherwise false
     */
    bool score(std::string prot_path, std::string lig_path, 
               std::vector<my_float_t> *scores);

    //! Store the string "DrugScore" in names
    /*! 
     * @param names Pointer to the vector to write "DrugScore" as the first and
     *        only scoring function
     */
    void sf_names(std::vector<std::string> *names);

  private:
    static const std::string _fname;
  };
}
#endif
