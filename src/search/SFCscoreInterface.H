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
 * $Source: /psa/share/repository/pfizer_proj/src/search/SFCscoreInterface.H,v $
 * $Revision: 1.4 $
 * $Author: vanvoor4 $
 * $Date: 2007-11-01 16:39:36 $
 * 
 * $Log: not supported by cvs2svn $
 * Revision 1.3  2007/08/21 19:30:53  vanvoor4
 * changed license header
 *
 * Revision 1.2  2007/06/13 17:18:46  vanvoor4
 * Added support for .spf file to come after the protein and ligands.
 *
 * Revision 1.1  2007/06/07 20:54:05  vanvoor4
 * initial checkin
 *
 * 
 * 
 */

#ifndef _SFCSCORE_INTERFACE_HEADER_FILE
#define _SFCSCORE_INTERFACE_HEADER_FILE

#include <ExternalScoringFunction.H>

namespace SimSite3D{

  //! Wrapper and parser for the external scoring functions in SFCscore
  /*! This is a hack to get SFCscore to work with SimSite3D.  At present, running
   * with a specified weight file is not supported (by SimSite3D; i.e. this class)
   */
  class SFCscoreInterface : public ExternalScoringFunction{
  public:
    //! Constructor
    /*!
     * @param cmdline String holding the SFCscore command line
     */
    SFCscoreInterface(std::string cmdline);

    //! Do nothing destructor
    ~SFCscoreInterface() {;}

    //! Wrapper to ExternalScoringFunction::run and parse output score file
    /*! 
     * 
     * @param prot_path String holding the path of the protein PDB file
     * @param lig_path String holding the path of the ligand mol2 file
     * @oaram scores Pointer to vector to hold the score(s) from SFCscore
     * @return True if SFCscore seemed to work, otherwise false
     */
    bool score(std::string prot_path, std::string lig_path, 
               std::vector<my_float_t> *scores);

    //! Get a list of the scoring functions calculated by SFCscore
    /*! Get the scoring functions listed in the .spf file.
     *
     * @param name Pointer to vector of strings to hold the names of the 
     *        scoring functions listed in the .spf file
     */
    void sf_names(std::vector<std::string> *names);

  private:
    static const std::string _fname; //!< "SFCscoreInterface.C"
    std::vector<std::string> sfc_names; //!< Vector holding names of SFC scoring functions requested in the spf file

    //! Get the names of the scoring functions in the .spf file
    /*! Given the .spf file as passed to the constructor, parse the file
     * looking for the line starting with "SCORE" (no preceeding spaces).
     * The listed scoring functions will be used to determine the files to 
     * parse and list in the score listing (report).
     *
     * @param spf_fname Path to the .spf file--typically specified in the 
     *        options field of
     *        the SimSite3DSoftParams/external_scoring_functions.txt file
     * @return True if the .spf file was opened for reading, otherwise false
     */
    bool get_scoring_functions(std::string spf_fname);
  };
}
#endif
