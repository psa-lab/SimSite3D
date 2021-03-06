/******************************************************************************
 * Copyright (c) 2006-2011, Michigan State University (MSU) Board of Trustees.
 * This file is part of the SimSite3D software project.
 *
 * Authors: Jeffrey Van Voorst, jeff.vanvoorst@gmail.com
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
 * $Date: 2011-02-25 10:01:22 -0500 (Fri, 25 Feb 2011) $: Date of last commit
 * $Author: vanvoor4 $: Author of last commit
 * $Rev: 1619 $: svn revision of last commit
 *
 * svn file: $Id: prot_lig_score.C 1619 2011-02-25 15:01:22Z vanvoor4 $
 * file location: $URL: file:///psa/share/repository/SimSite3D/branches/surfaces-branch/src/score/prot_lig_score.C $
 *****************************************************************************/
#include <ProtLigActTable.H>
#include <ProtLigScoreParameters.H>

using namespace SimSite3D;

typedef struct{
  std::string lig_id;
  my_float_t affi;
  my_float_t affieff;
  my_float_t orient;
  std::vector<size_t> atom_idz;
  //std::vector<atom_vci> prot_polar_atoms;
  //std::vector<residue_vci> prot_hphob_SC;
  //std::vector<atom_vci> prot_atoms;
}score_t;

// -- NOTE::::::::
// the ligand iterators will be invalidated anyhow -- we should only keep
// protein atom info -- is protein atom indexes faster than protein atom
// iterators?

// idea -- save only the polar atoms & residues hit for each ligand
// sort the vectors 
// save the uniques in a map --
// use residue->CB (with lookup map) to save hphob residues.

int main(const int argc, const char **argv)
{
  std::cout << "\n" << argv[0] << " (" << PACKAGE_NAME << ") " 
            << PACKAGE_VERSION << "\n\n";

  // Do not return -1 here since system, fork, etc return -1 on failure, and we
  // wish to distinguish between system and program failure
  ProtLigScoreParameters my_params(argc, argv);
  BaseParameters::status_t status = my_params.status();
  if(status == BaseParameters::DISPLAY_HELP_ONLY) return 0;
  else if(status != BaseParameters::READY){
    std::cerr << "\n" << argv[0]
              << " *FAILED* \n\tCould not initialize parameters\n";
    return 1;
  }
  
  PDBStructure prot(my_params.prot_fname);
  if(prot.fail()){
    std::cerr << " *FAILURE* \n\tCould not load the protein structure\n";
    return 1;
  }
  if(my_params.lig_fname.length()){
    mol2File lig(my_params.lig_fname);
    if(lig.fail()){
      std::cerr << " *FAILURE* \n\tCould not load the ligand structure\n";
      return 1;
    }
    
    lig.calc_charge_sums();
    ProtLigScore my_score(prot, lig, my_params.print_interactions);
    my_score.score_header_lines(std::cout);
    my_score.report_scores(std::cout);
    if(my_params.print_interactions){
      std::cout << "\nLigand -- protein interactions\n";
      my_score.write_lig_acts();
    }
  }else{
    std::ifstream lig_list_file;
    open_ifstream(lig_list_file, my_params.lig_list_fname);
    if(lig_list_file.fail()) return 1;

    ProtLigActTable table(prot);
    ProtLigScore::score_header_lines(std::cout);
    for(std::string line; std::getline(lig_list_file, line); ){
      mol2File lig(line);
      if(lig.fail()){
        std::cerr << "Could not load the ligand file\n"
                  << "\t" << line << "\n\tSkipping . . .\n";
        continue;
      }

      lig.calc_charge_sums();
      ProtLigScore my_score(prot, lig);
      my_score.report_scores(std::cout);
      if(my_params.build_interact_tbl) table.add_row(my_score);
    }
    if(my_params.build_interact_tbl) table.print_table(std::cout);
  }

  return 0; 
}
