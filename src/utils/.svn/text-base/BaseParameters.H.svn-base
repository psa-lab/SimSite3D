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
 * $Source: /psa/share/repository/pfizer_proj/src/basics/BaseParameters.H,v $
 * $Revision: 1.3 $
 * $Author: vanvoor4 $
 * $Date: 2008-07-29 16:02:19 $
 * 
 * $Log: not supported by cvs2svn $
 * Revision 1.2  2007/12/17 21:00:13  vanvoor4
 * Added some comments
 *
 * Revision 1.1  2007/11/01 16:02:18  vanvoor4
 * Initial checkin
 * 
 *
 */

#ifndef BASE_PARAMETERS_HEADER_INCLUDED
#define BASE_PARAMETERS_HEADER_INCLUDED

#include <string>
#include <map>

namespace ASCbase{

  //! A method to load database, diverse sitemaps, etc directories from 
  //! environment variables and parameter files
  /*! The order is to first read the $ASCBASE_INSTALL_DIR environment variable.
   * Second, parse $ASCBASE_INSTALL_DIR/ASCbaseSoftParams/ascbase_software.conf
   * for system wide default values.  Thirdly, check for any updated parameters
   * in the environment.  Finally, if load_conf_file is called by a derived 
   * class, the parsed values override all previously stored values.
   */
  class BaseParameters{
  public:
    typedef enum{
      FATAL_ERROR,
      INVALID_PARAMETER,
      INITIALIZING, 
      DISPLAY_HELP_ONLY,
      READY
    }status_t;

    BaseParameters();

    ~BaseParameters()
    { ; }

    static bool get_env_var(const std::string var, std::string *val);

    status_t
    status() const
    { return A_status; }

    std::string dbase_sites;    //!< Directory holding the database sitemaps
    std::string dbase_ligs;     //!< Directory holding the database ligands
    std::string dbase_prots;    //!< Directory holding the database proteins
    std::string diverse_sites;  //!< Directory holding the diverse sitemaps
    std::string diverse_ligs;//!< Directory holding ligands for diverse sitemaps
    std::string proj_output; //!< Directory to store results 
    std::string scratch_dir; //!< Directory where ASCbase can create temp files

    // for now we keep this around till we fold the information into 
    // header files.
    std::string install_dir; 

    // While testing keep this here
    bool load_surf_files;

    // needed in this class so that Sitemap can see the variable
    bool require_min_npts;

  protected:
    bool load_conf_file(std::string conf_fname);

    void
    print_version(std::ostream &out, std::string prog_name);

    status_t A_status;  //!< Parameter status

  private:

    void load_environment();

    void init_str_to_var_map();


    std::map<std::string, std::string*> A_str_to_var;
    static const std::string A_fname;  //!< Name of source file
  };
}
#endif
