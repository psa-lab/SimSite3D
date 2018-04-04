/******************************************************************************
 * Copyright (c) 2006-2011, Michigan State University (MSU) Board of Trustees.
 * This file is part of the SimSite3D software project.
 *
 * Authors: Jeffrey Van Voorst, vanvoor4@msu.edu
 *          Leslie Kuhn, Ph.D., KuhnL@msu.edu 
 *
 *  SimSite3D is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  SimSite3D is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc.,
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *  You may also visit http://www.gnu.org/licenses/gpl-2.0.html.
 * 
 * $Date: 2011-02-25 10:12:03 -0500 (Fri, 25 Feb 2011) $: Date of last commit
 * $Author: vanvoor4 $: Author of last commit
 * $Rev: 1620 $: svn revision of last commit
 *
 * svn file: $Id: search_sitemaps.C 1620 2011-02-25 15:12:03Z vanvoor4 $
 * file location: $URL: file:///psa/share/repository/SimSite3D/branches/surfaces-branch/src/search/search_sitemaps.C $
 *****************************************************************************/
#include <Search.H>

using namespace ASCbase;

int main(const int argc, const char **argv)
{
  std::cout << "\n" << argv[0] << " (" << PACKAGE_NAME << ") " 
            << PACKAGE_VERSION << "\n\n";

  // Do not return -1 here since system, fork, etc return -1 on failure, and we
  // wish to distinguish between system and program failure
  SearchParameters my_params(argc, argv);
  BaseParameters::status_t status = my_params.status();
  if(status == BaseParameters::DISPLAY_HELP_ONLY) return 0;
  else if(status != BaseParameters::READY){
    std::cerr << "\n" << argv[0]
              << " *FAILED* \n\tCould not initialize parameters\n";
    return 1;
  }
  Search my_search(&my_params);
  if(my_search.fail()){
    std::cerr << "\n" << argv[0]
              << " *FAILED* \n\tCould not initialize search\n";
    return 1;
  }

  if(!my_search.run()) return 1;
  return 0; 
}
