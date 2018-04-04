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
 * $Source: /psa/share/repository/pfizer_proj/src/gen_points/SitemapPointsFile.H,v $
 * $Revision: 1.3 $
 * $Author: vanvoor4 $
 * $Date: 2009-01-12 21:12:10 $
 * 
 * $Log: not supported by cvs2svn $
 * Revision 1.2  2007/08/29 20:23:14  vanvoor4
 * Changed handling of hbond point types
 *
 * Revision 1.1  2007/08/21 18:43:01  vanvoor4
 * initial checkin
 *
 *
 * 
 * 
 */

#ifndef _SITEMAP_POINTS_FILE_HEADER_INCLUDED
#define _SITEMAP_POINTS_FILE_HEADER_INCLUDED

#include <PDBBase.H>
#include <HbondPoints.H>
#include <HphobPoints.H>

namespace ASCbase{

  //! Creation and storage of sitemap points files
  class SitemapPointsFile : public PDBBase {
  public:
    //! Cstr to create a sitemap points file
    SitemapPointsFile();

    //! Cstr to load an existing sitemap points file
    SitemapPointsFile(const std::string filename, 
                      const verbose_level_t verbosity = VERBOSE_SILENT);

    ~SitemapPointsFile();

    //! Add hydrogen bonding points
    bool add_points(const hbond_fit_pt_vci _begin, const hbond_fit_pt_vci _end);

    //! Add hydrophobic points 
    bool add_points(const hphob_point_lci _begin, const hphob_point_lci _end);

    const uint num_polar_points() const 
    { 
      return counts[0] + counts[1] + counts[2]; 
    }

    const uint num_hphobic_points() const 
    { 
      return counts[3];
    }
  private:
    //! Shamelessly assumes header comes before data
//    bool read_header(std::ifstream &fin);

    std::vector<uint> counts;  
    static const std::string _fname;  //!< Source file name
  };
}

#endif
