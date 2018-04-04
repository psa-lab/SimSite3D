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
 * $Date: 2011-02-25 10:00:51 -0500 (Fri, 25 Feb 2011) $: Date of last commit
 * $Author: vanvoor4 $: Author of last commit
 * $Rev: 1618 $: svn revision of last commit
 *
 * svn file: $Id: gen_points.C 1618 2011-02-25 15:00:51Z vanvoor4 $
 * file location: $URL: file:///psa/share/repository/SimSite3D/branches/surfaces-branch/src/sitemap/gen_points.C $
 *****************************************************************************/
#include <GenPointsParameters.H>
#include <Sitemap.H>
#include <normalize_sitemap.H>

int main(const int argc, const char** argv)
{
  std::cout << argv[0] << " (" << PACKAGE_NAME << ") " << PACKAGE_VERSION
            << "\n\n";

  // Do not return -1 here since system, fork, etc return -1 on failure, and we
  // wish to distinguish between system and program failure
  ASCbase::GenPointsParameters my_params(argc, argv);
  ASCbase::BaseParameters::status_t status = my_params.status();
  if(status == ASCbase::BaseParameters::DISPLAY_HELP_ONLY) return 0;
  else if(status != ASCbase::BaseParameters::READY){
    std::cerr << "\n" << argv[0] 
              << " *FAILED* \n\tCould not initialize parameters\n";
    return 1;
  }
  ASCbase::Sitemap my_site(my_params);
  if(my_site.fail()){
    std::cerr << "\n" << argv[0] 
              << " *FAILED* \n\tCould not initialize sitemap\n";
    return 1;
  }

  std::string path;
  std::string struct_id;
  get_path_and_struct_id(my_params.pts_fname, &path, &struct_id);

  if(my_params.require_min_npts && !my_site.has_enough_points()){ 
    std::cerr << "\nUnable to generate site map for " << struct_id << "\n";
    return 1;
  }

  my_site.write_files(struct_id, path);
  std::cout << "\nFinished generating site map for " << struct_id << "\n";

  int rv = 0;
  if(my_params.normalize){
    std::string xml_fname = path + "/" + struct_id + "_s.csv";
    std::cout << "file name is: " << xml_fname << "\n";
    rv = normalize_sitemap(xml_fname, my_params);
  }

  if(path.length() <= 2)
    std::cout << "  The site map files may be found in your current "
              << "directory\n\n";
  else std::cout << "  The site map files may be found at " << path << "\n\n";


  return rv;
}
