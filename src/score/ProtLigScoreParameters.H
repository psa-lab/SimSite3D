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
 * Authors: Jeffrey Van Voorst, vanvoor4@msu.edu
 *          Leslie Kuhn, Ph.D., KuhnL@msu.edu 
 *****************************************************************************/
  
/*
 * $Source: /psa/share/repository/pfizer_proj/src/gen_points/GenPointsParameters.H,v $
 * $Revision: 1.6 $
 * $Author: vanvoor4 $
 * $Date: 2009-01-12 21:07:26 $
 * 
 * $Log: not supported by cvs2svn $
 * Revision 1.5  2008/07/29 16:03:01  vanvoor4
 * Changed status values from popt, etc to be an enum rather than boolean
 *
 * Revision 1.4  2008/07/29 14:46:56  vanvoor4
 * Changed the default parameters to use the hydrophobic spheres as the
 * default hydrophobic representation.
 *
 * Revision 1.3  2007/12/17 21:29:16  vanvoor4
 * Added support for waters on command line
 *
 * Revision 1.2  2007/11/14 16:23:03  vanvoor4
 * Name change of some variables and removed the conf file option from the
 * command line.
 *
 * Revision 1.1  2007/11/01 16:09:07  vanvoor4
 * Initial Checkin
 * 
 *
 */

#ifndef PROT_LIG_SCORE_PARAMETERS_HEADER_INCLUDED
#define PROT_LIG_SCORE_PARAMETERS_HEADER_INCLUDED

#include <basics.H>
#include <BaseParameters.H>

namespace SimSite3D{

//! A data class (public members) used to hide the details, from the main
//! application source, needed to get the parameters. 
class ProtLigScoreParameters : public BaseParameters{
public:

  typedef enum{
    PRINT_VERSION = 1,
    PRINT_HELP,
    BUILD_INTERACTIONS_TABLE,
    PRINT_INTERACTIONS
  }popt_args_t;

  //! To be used by the SimSite3DPy score.parameters module
  ProtLigScoreParameters();

  //! Set the score parameters using the values in argv
  ProtLigScoreParameters(const int argc, const char** argv);

  ~ProtLigScoreParameters();

  std::string prot_fname;
  std::string lig_fname;
  std::string lig_list_fname;
  bool build_interact_tbl;
  bool print_interactions;

private:
  //! Initialize the class variables to standard nondefined values 
  //! (such as NULL).
  void init();

  //! Given the command line arguments, parse them using the popt library
  /*! 
   * @param argc Number of command line arguments
   * @param argv Array of C-style strings holding the command line arguments
   * @return True if parsed successfully and may continue, else false
   */
  status_t get_opts(const int argc, const char** argv);

  //! Simple verification of command line parameters
  /*!
   * @return False implies unsafe to proceed, otherwise should be OK
   */
  status_t verify_parameters();
  
  static const std::string A_fname; //!< Name of source file
};
} // end namespace SimSite3D
#endif
