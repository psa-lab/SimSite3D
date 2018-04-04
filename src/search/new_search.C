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
 * $Date: 2012-02-22 00:16:10 -0500 (Wed, 22 Feb 2012) $: Date of last commit
 * $Author: vanvoor4 $: Author of last commit
 * $Rev: 1632 $: svn revision of last commit
 *
 * svn file: $Id: new_search.C 1632 2012-02-22 05:16:10Z vanvoor4 $
 * file location: $URL: file:///psa/share/repository/SimSite3D/branches/surfaces-branch/src/search/new_search.C $
 *****************************************************************************/
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <openssl/evp.h>
#include <sstream>
#include <Search.H>

using namespace SimSite3D;

bool 
compute_hash(const SearchParameters &params, std::string hash_val)
{
  std::string path;
  std::string struct_id;
  get_path_and_struct_id(params->model_file_name, &path, &struct_id);

  EVP_MD_CTX ctx;
  EVP_MD_CTX_init(&ctx);
  if(EVP_DigestInit_ex(&ctx, EVP_sha1, NULL) == 0) return false;

  // loop over the file types "_s.csv", "_s.pdb", "_surf.face", "_surf.vert"
  // for each one read it in, and do an EVP_DigestUpdate(&ctx, *, size);      
  // maybe we just ....


  unsigned char md_value[EVP_MAX_MD_SIZE];
  int md_len;
  if(EVP_DigestFinal_ex(&ctx, md_value, &md_len) == 0) return false;
  EVP_MD_CTX_cleanup(&ctx);

// hash each sitemap file except for the "new" sitemap normalization data
// file.
// we allso need to hash 
// params.score_str
// params.max_tier1_aligns
// params.min_num_atoms
// params.dmetol
// params.lsetol
// params.allow_hphob_triangles
// params.use_hbond_surfaces
// params.fine_tune_tier2_alignments
// params.fine_tune_best_tier2_alignment
// params.scale_terms
// params.do_IK
// params.max_corr_surf_pt_dist
// params.fine_tune_ratio
// params.IK_type
}

// in the child we need to change the values for:
// params.require_min_npts to false
// params.ext_score_method to nothing
// params.ofname to a temp file
// params.db_index_fname to nothing
// params.db_start to nothing
// params.db_stop to nothing
// params.normalize to false
// params.num_scores_to_keep to 1
// params.score_cutoff to 1000.0
// params.min_num_atoms to 0
// params.align_to_query to false
// params.num_rand_aligns to zero
// ok need to move timing functions  ......
// params.add_struct_id_field to false
// params.do_internal_prot_lig_score to false
// params.check_all_triangles to false
// params.save_rigid_scores to false
// params.alignments_fname to default value
// params.dbase_sites to right value
// params.dbase_ligs to right value
// params.dbase_prots to right value
// params.proj_output to right value
// params. time_process is false

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

  if(my_params.normalize == true){
    
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
