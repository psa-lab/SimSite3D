
#include <sstream>
#include <ProtLigScore.H>

using namespace SimSite3D;

const my_float_t ProtLigScore::MIN_BURY_SQUARED_DIST = 1.9 * 1.9;
const my_float_t ProtLigScore::MAX_BURY_SQUARED_DIST = 4.5 * 4.5;
const my_float_t ProtLigScore::HYDRO_SQUARED_DIST = 4.5 * 4.5;
const my_float_t ProtLigScore::MAX_INTERACT_SQUARED_DIST = 4.5 * 4.5;
const my_float_t ProtLigScore::AFFI_64_WEIGHTS[] = 
  {-5.218,-0.239,-0.083, 0.134,0.000,0.000,0.000 };
const my_float_t ProtLigScore::ORIENT_99_WEIGHTS[] = 
  {-3.917,-0.061,-0.023,-0.189,0.000,0.000,0.000 };


ProtLigScore::ProtLigScore(PDBStructure &prot, mol2File &lig,
                           bool save_prot_lig_acts, bool compute_charge_sums)
{
  A_prot = &prot;
  A_lig = &lig;
  A_hbond_geometry.set_protein(A_prot);
  A_raw_terms.num_metal_hbonds = 0;
  A_raw_terms.num_salt_bridges = 0;
  A_raw_terms.num_hbonds = 0;
  A_raw_terms.num_unsat_polar = 0;
  A_raw_terms.num_unsat_charge = 0;
  A_raw_terms.total_target_hydro = 0;
  A_raw_terms.contact_hphob_hphob = 0;
  A_raw_terms.total_hphob_hphob_comp = 0;
  A_num_lig_carbons = 0;
  A_num_exposed_lig_carbons = 0;
  A_affiscore = my_float_max;
  A_orientscore = my_float_max;
  A_affiefficientscore = my_float_max;

  if(compute_charge_sums) lig.calc_charge_sums();

  // Ignore hydrogen atoms --- interactions are only associated with heavy atoms
  atom_vci lig_atom = A_lig->atoms_begin();
  for( ; lig_atom < A_lig->atoms_end(); ++lig_atom)
    if(lig_atom->name != H)
      A_ligand_flag.insert(atom_act_pair(lig_atom, INITIAL));

  // I need to get a much better handle on PDB chains versus small molecules
  // and peptides
  residue_vci res = A_prot->residues_begin();
  for(res = A_prot->residues_begin(); res < A_prot->residues_end(); ++res){
    for(atom_vci atom = res->atoms_begin; atom < res->atoms_end; ++atom)
      if(atom->name != H)
        A_target_flag.insert(atom_act_pair(atom, INITIAL));
  }
  std::vector<atom_vci>::const_iterator metal_iter = A_prot->metals_begin();
  for( ; metal_iter < A_prot->metals_end(); ++metal_iter)
    A_target_flag.insert(atom_act_pair(*metal_iter, INITIAL));

  A_keep_prot_lig_acts = save_prot_lig_acts;
  if(A_keep_prot_lig_acts){
    A_hphob_prot_lig_SC.resize(lig.num_atoms());
    A_lig_exposed_Cz.resize(lig.num_atoms());
    std::fill(A_lig_exposed_Cz.begin(), A_lig_exposed_Cz.end(), false);
  }

  compute_terms();
  compute_scores();

#if 0
  print_hbonds(std::cout);

  std::cout << "# metal hbonds      : " << A_raw_terms.num_metal_hbonds << "\n";
  std::cout << "# salt bridges      : " << A_raw_terms.num_salt_bridges << "\n";
  std::cout << "# hydrogen bonds    : " << A_raw_terms.num_hbonds << "\n";
  std::cout << "# usat polar        : " << A_raw_terms.num_unsat_polar << "\n";
  std::cout << "# usat charge       : " << A_raw_terms.num_unsat_charge << "\n";
  std::cout << "# total_target_hydro: " << A_raw_terms.total_target_hydro << "\n";
  std::cout << "# contact_hphob_hphob: " << A_raw_terms.contact_hphob_hphob << "\n";
  std::cout << "# total_hphob_hphob_comp: " << A_raw_terms.total_hphob_hphob_comp << "\n";

#endif
  //print_atom_act_classes(&A_ligand_flag);
  //print_atom_act_classes(&A_target_flag);


  //A_prot = 0;
  //A_lig = 0;
}
#if 0
    NUMBER_OF_METAL_HBONDS = 0, -- needs to be done
    NUMBER_OF_SALT_BRIDGES, -- needs metal SBs
    NUMBER_OF_HBONDS, -- done
    NUMBER_OF_UNSAT_POLAR, -- done
    NUMBER_OF_UNSAT_CHARGE, -- done
    TOTAL_TARGET_HYDRO, -- done
    CONTACT_HPHOB_HPHOB, -- done
    TOTAL_HPHOB_HPHOB_COMP -- done
#endif

ProtLigScore::~ProtLigScore()
{

}

void
ProtLigScore::compute_scores()
{
  A_orientscore = ORIENT_99_WEIGHTS[0] +
                  ORIENT_99_WEIGHTS[1] * A_raw_terms.total_target_hydro +
                  ORIENT_99_WEIGHTS[2] * A_raw_terms.contact_hphob_hphob +
                  ORIENT_99_WEIGHTS[3] *(A_raw_terms.num_metal_hbonds + 
                                         A_raw_terms.num_salt_bridges + 
                                         A_raw_terms.num_hbonds);
  A_affi_terms.resize(4);
  A_affi_terms[0] = AFFI_64_WEIGHTS[0];
  A_affi_terms[1] = AFFI_64_WEIGHTS[1] * A_raw_terms.total_hphob_hphob_comp;
  A_affi_terms[2] = AFFI_64_WEIGHTS[2] * (A_raw_terms.num_metal_hbonds +
                                          A_raw_terms.num_salt_bridges +
                                          A_raw_terms.num_hbonds);
  A_affi_terms[3] = AFFI_64_WEIGHTS[3] * (A_raw_terms.num_unsat_polar +
                                          A_raw_terms.num_unsat_charge);
  A_affiefficientscore = A_affi_terms[1] + A_affi_terms[2] + A_affi_terms[3];
  A_affiscore = A_affi_terms[0] + A_affiefficientscore;
  A_affiefficientscore /= A_ligand_flag.size();
}

void
ProtLigScore::score_header_lines(std::ostream &out, bool classic)
{
  out << "Score fields\n" 
      << " 1) Ligand file name\n"
      << " 2) OrientScore (orientation score)\n"
      << " 3) Ligand efficency (AffiScore / N)\n"
      << " 4) AffiScore (affinity score) ~ Kcals/mol\n"
      << " 5) Weighted contribution of hydrophobic complementarity\n"
      << " 6) Weighted contribution of protein-ligand salt bridges, "
      << "hydrogen bonds, and\n    ligand-metal interactions\n"
      << " 7) Weighted contribution of protein and ligand unsatisfied polar "
      << "and \n    charged atoms\n"
      << " 8) Total hydrophobic complementarity\n"
      << " 9) Count of protein-ligand hydrogen bonds\n"
      << "10) Count of protein-ligand salt bridges\n"
      << "11) Count of metal-ligand interactions\n"
      << "12) Count of interfacial protein and ligand unsatisfied polar atoms\n"
      << "13) Count of interfacial protein and ligand unsatisfied charged "
      << "atoms\n"
      << "14) Percentage of buried ligand carbon atoms\n";
  if(classic)
    out << "15) No bump checking was performed\n"
        << "16) Total overlap was not checked\n";
}

void
ProtLigScore::classic_report(std::ostream &out)
{
  out.precision(3);
  // Uncomment for fixed float representation
  //out.setf(std::ios_base::fixed, std::ios_base::floatfield);
  out << A_lig->name() << ": O: " << A_orientscore
      << " A/N: " << A_affiefficientscore << " A: " << A_affiscore << " [" 
      << A_affi_terms[1] << " " << A_affi_terms[2] << " " 
      << A_affi_terms[3] << "] " << A_raw_terms.contact_hphob_hphob << " "
      << A_raw_terms.num_hbonds << " " << A_raw_terms.num_salt_bridges << " "
      << A_raw_terms.num_metal_hbonds  << " "
      << A_raw_terms.num_unsat_polar << " "
      << A_raw_terms.num_unsat_charge << " "
      << (1.0 - A_num_exposed_lig_carbons / 
         static_cast<my_float_t>(A_num_lig_carbons))*100.0 
      // We don't do bumping -- no overlap or number of bumps to output
      << "% 0 0.0" << "\n";
}

void
ProtLigScore::report_scores(std::ostream &out, char delim)
{
  out.precision(3);
  // Uncomment for fixed float representation
  //out.setf(std::ios_base::fixed, std::ios_base::floatfield);
  out << A_lig->name() << delim << A_orientscore << delim
      << A_affiefficientscore << delim << A_affiscore << delim
      << A_affi_terms[1] << delim << A_affi_terms[2] << delim
      << A_affi_terms[3] << delim << A_raw_terms.contact_hphob_hphob << delim
      << A_raw_terms.num_hbonds << delim 
      << A_raw_terms.num_salt_bridges << delim
      << A_raw_terms.num_metal_hbonds << delim
      << A_raw_terms.num_unsat_polar << delim
      << A_raw_terms.num_unsat_charge << delim
      << (1.0 - A_num_exposed_lig_carbons / 
         static_cast<my_float_t>(A_num_lig_carbons))*100.0 
      // We don't do bumping -- no overlap or number of bumps to output
      << "%" << delim << "\n";
}

void
ProtLigScore::gen_lig_act_strings(const char field_delim,
                                  const char prot_delim)
{
  A_lig_act_strings.clear();
  A_lig_act_strings.resize(A_lig->num_atoms());

  // Exposed ligand carbon atoms
  for(size_t i = 0; i < A_lig_act_strings.size(); ++i)
    if(A_lig_exposed_Cz[i]) A_lig_act_strings[i] += std::string("EXPOSED");

  // Prot-lig salt bridges
  for(size_t i = 0; i < A_prot_lig_SB.size(); ++i){
    atom_vci &atom = A_prot_lig_SB[i].A;
    std::vector<std::string>::iterator my_str;
    my_str = A_lig_act_strings.begin() +
      (A_prot_lig_SB[i].B - A_lig->atoms_begin());
    std::ostringstream ostr;
    if(my_str->size()) ostr << prot_delim;
    if(atom->chainID != ' ') ostr << "(" << atom->chainID << ")";
    ostr << PDB_residues::residue_to_string(atom->res) << atom->res_num;
    if(atom->iCode != ' ') ostr << atom->iCode;
    ostr << " " << PDB_residues::atom_to_string(atom->name);

    std::cout << "SB: " << ostr.str() << std::endl;
    *my_str += ostr.str();
  }

  // Prot-lig hbonds
  for(size_t i = 0; i < A_prot_lig_hbonds.size(); ++i){
    atom_vci &atom = A_prot_lig_hbonds[i].A;
    std::vector<std::string>::iterator my_str;
    my_str = A_lig_act_strings.begin() +
      (A_prot_lig_hbonds[i].B - A_lig->atoms_begin());
    std::ostringstream ostr;
    if(my_str->size()) ostr << prot_delim;
    if(atom->chainID != ' ') ostr << "(" << atom->chainID << ")";
    ostr << PDB_residues::residue_to_string(atom->res) << atom->res_num;
    if(atom->iCode != ' ') ostr << atom->iCode;
    ostr << " " << PDB_residues::atom_to_string(atom->name);

//    std::cout << "Hbond: " << ostr.str() << std::endl;
    *my_str += ostr.str();
  }

  // Ligand hydrophobic interactions with protein hydrophobic side chain atoms
  std::vector<std::string>::iterator my_str = A_lig_act_strings.begin();
  for(size_t i = 0; i < A_hphob_prot_lig_SC.size(); ++i, ++my_str){
    // Get unique side chains for each lig atom
    std::vector<residue_vci>::const_iterator unq_res_end;
    unq_res_end = std::unique(A_hphob_prot_lig_SC[i].begin(),
                              A_hphob_prot_lig_SC[i].end());
    std::vector<residue_vci>::const_iterator res_iter;
    res_iter = A_hphob_prot_lig_SC[i].begin();
    for( ; res_iter < unq_res_end; ++res_iter){
      const residue_vci &res = *res_iter;
      std::ostringstream ostr;
      if(my_str->size()) ostr << prot_delim;
      if(res->chainID != ' ') ostr << "(" << res->chainID << ")";
      ostr << PDB_residues::residue_to_string(res->name) << res->number;
      if(res->icode != ' ') ostr << res->icode;

      *my_str += ostr.str();
    }
  }
}

void
ProtLigScore::write_lig_acts(std::ostream &out, char field_delim,
                             char prot_delim)
{
  if(!A_lig_act_strings.size()) gen_lig_act_strings(field_delim, prot_delim);

  // Write it out as a row
  //   Write column headings
  atom_vci lig_atom = A_lig->atoms_begin();
  std::vector<std::string>::const_iterator prot_act = A_lig_act_strings.begin();
  for( ; lig_atom < A_lig->atoms_end(); ++lig_atom, ++prot_act)
    if(lig_atom->name != H) out << lig_atom->name_str << field_delim;
  out << "\n";
  //   Write binary values
  lig_atom = A_lig->atoms_begin();
  prot_act = A_lig_act_strings.begin();
  for( ; lig_atom < A_lig->atoms_end(); ++lig_atom, ++prot_act)
    if(lig_atom->name != H){
      if(prot_act->length() && *prot_act != "EXPOSED") 
        out << "1" << field_delim;
      else out << "0" << field_delim;
    }
  out << "\n\n";
 
  // Write it out as columns
  lig_atom = A_lig->atoms_begin();
  prot_act = A_lig_act_strings.begin();
  for( ; lig_atom < A_lig->atoms_end(); ++lig_atom, ++prot_act)
    if(lig_atom->name != H)
      out << lig_atom->name_str << field_delim 
          << *prot_act << field_delim << "\n";
}

void
ProtLigScore::compute_terms()
{
  atom_vci lig_atom = A_lig->atoms_begin();
  for( ; lig_atom < A_lig->atoms_end(); ++lig_atom){
    if(lig_atom->name == H) continue;
 
    int sum_hphob_hphob = 0;
    int num_hphob_hphob_contact_1_lig_atm = 0;
    int num_lig_nbrs = 0;
    if(lig_atom->name == C) ++A_num_lig_carbons;

    metal_interactions(lig_atom, &num_lig_nbrs);
    prot_lig_atom_interactions(lig_atom, &num_lig_nbrs, &sum_hphob_hphob,
                               &num_hphob_hphob_contact_1_lig_atm);

    if(A_ligand_flag[lig_atom] != INITIAL){
      if(lig_atom->hydro < 0 && sum_hphob_hphob != 0)
        sum_hphob_hphob += lig_atom->hydro;
      if(sum_hphob_hphob != 0 && num_hphob_hphob_contact_1_lig_atm != 0){
        A_raw_terms.total_hphob_hphob_comp += 
          1.0 - 
          std::fabs(lig_atom->hydro - 
                    (sum_hphob_hphob - lig_atom->hydro) / 
                    (static_cast<my_float_t>(num_hphob_hphob_contact_1_lig_atm))
                   )/635.0;
      }
    }

    // If this ligand atom does not have any protein neighbor atoms, it is 
    // exposed
    if(lig_atom->name == C && num_lig_nbrs == 0){
      ++A_num_exposed_lig_carbons;
      if(A_keep_prot_lig_acts)
        *(A_lig_exposed_Cz.begin() + (lig_atom - A_lig->atoms_begin())) = true;
    }
  }

  flag_intra_ligand_hbonds();
  flag_intra_protein_hbonds();
  sum_target_hydrophobicity();

  int num_prot_unsat_charge, num_prot_unsat_polar;
  int num_lig_unsat_charge, num_lig_unsat_polar;
  flag_unsat_polar_atoms(&A_target_flag, &num_prot_unsat_charge, 
                         &num_prot_unsat_polar);
  flag_unsat_polar_atoms(&A_ligand_flag, &num_lig_unsat_charge, 
                         &num_lig_unsat_polar);

  A_raw_terms.num_unsat_charge = num_prot_unsat_charge + num_lig_unsat_charge;
  A_raw_terms.num_unsat_polar = num_prot_unsat_polar + num_lig_unsat_polar; 
  A_raw_terms.num_salt_bridges = A_prot_lig_SB.size();
  A_raw_terms.num_hbonds = A_prot_lig_hbonds.size();
}

void
ProtLigScore::metal_interactions(atom_vci lig_atom, int *num_lig_nbrs)
{
  atom_vci prev_prot_atom = atom_t::NULL_ATOM_VCI; 
  std::vector<atom_vci>::const_iterator metal_iter = A_prot->metals_begin();
  for( ; metal_iter < A_prot->metals_end(); ++metal_iter){
    atom_vci metal = *metal_iter;

    // use octree here after errors are gone
    my_float_t sq_dist = dist_squared(lig_atom->pos, metal->pos);
    if(MAX_INTERACT_SQUARED_DIST < sq_dist) continue;

    // Search for direct-interfacial atoms
    if(HYDRO_SQUARED_DIST < sq_dist) continue;
    ++(*num_lig_nbrs);
    if(A_ligand_flag[lig_atom] == INITIAL) 
      A_ligand_flag[lig_atom] = INTERFACIAL_ATOM;
    if(A_target_flag[metal] == INITIAL) 
      A_target_flag[metal] = INTERFACIAL_ATOM;

    // Does this make sense?  I have seen residues duplicated, or portions
    // of residues duplicated -- i.e. two possible positions are not
    // on consecutive atom lines.
    //
    // we don't want duplicated interaction between the same pair of atoms.
    // Some side chains in protein may have more than one orientations
    if(prev_prot_atom != atom_t::NULL_ATOM_VCI &&
       metal->res == prev_prot_atom->res &&
       metal->name == prev_prot_atom->name &&
       metal->altLoc != prev_prot_atom->altLoc) continue;
    prev_prot_atom = metal;

    if(A_hbond_geometry.metal_hbond(*lig_atom, *metal, sq_dist)){
      A_ligand_flag[lig_atom] = METAL_DIRECT_HBOND;
      A_target_flag[metal] = METAL_DIRECT_HBOND;
      ++(A_raw_terms.num_metal_hbonds);
    }else if(A_hbond_geometry.metal_salt_bridge(*lig_atom, *metal, sq_dist)){
      A_ligand_flag[lig_atom] = SALT_BRIDGE;
      A_target_flag[metal] = SALT_BRIDGE;
      ++(A_raw_terms.num_salt_bridges);
    }
  }
}

void
ProtLigScore::prot_lig_atom_interactions(atom_vci lig_atom, int *num_lig_nbrs,
                                         int *sum_hphob_hphob,
                                         int *num_hphob_hphob_contact_1_lig_atm)
{
  atom_vci prev_prot_atom = atom_t::NULL_ATOM_VCI; 
  residue_vci res = A_prot->residues_begin();
  for(res = A_prot->residues_begin(); res < A_prot->residues_end(); ++res){
    atom_vci prot_atom = res->atoms_begin;
    for( ; prot_atom < res->atoms_end; ++prot_atom){
      if(prot_atom->name == H) continue;

      // use octree here after errors are gone
      my_float_t sq_dist = dist_squared(lig_atom->pos, prot_atom->pos);
      if(MAX_INTERACT_SQUARED_DIST < sq_dist) continue;

      // Search for direct-interfacial atoms
      if(HYDRO_SQUARED_DIST < sq_dist) continue;
      ++(*num_lig_nbrs);
      if(A_ligand_flag[lig_atom] == INITIAL) 
        A_ligand_flag[lig_atom] = INTERFACIAL_ATOM;
      if(A_target_flag[prot_atom] == INITIAL) 
        A_target_flag[prot_atom] = INTERFACIAL_ATOM;

      // we don't want duplicated interaction between the same pair of atoms.
      // Some side chains in protein may have more than one orientations
      if(prev_prot_atom != atom_t::NULL_ATOM_VCI &&
         prot_atom->res == prev_prot_atom->res &&
         prot_atom->name == prev_prot_atom->name &&
         prot_atom->altLoc != prev_prot_atom->altLoc) continue;
      prev_prot_atom = prot_atom;
    
      // Check for polar interaction
      if((lig_atom->act_type != NOTHING || lig_atom->charge != 0.0) &&
         ((prot_atom->act_type != NOTHING && prot_atom->act_type != HYDROPHOB)
         || prot_atom->charge != 0.0))
        polar_interaction(prot_atom, lig_atom, sq_dist);

      // Check for hydrophobic interaction and/or contacts
      if(lig_atom->hydro < 0 && prot_atom->hydro < 0){
        ++A_raw_terms.contact_hphob_hphob; 
        *sum_hphob_hphob += prot_atom->hydro;
        ++(*num_hphob_hphob_contact_1_lig_atm);
        // Keep it like this for now since it is how we did it for the
        // Thrombin screening paper
        if(prot_atom->act_type == HYDROPHOB){
          A_hphob_prot_SC[res] = true;
          if(A_keep_prot_lig_acts)
            (A_hphob_prot_lig_SC.begin() + 
             (lig_atom - A_lig->atoms_begin()))->push_back(res);
        }
      }
    }
  }
}

void
ProtLigScore::polar_interaction(atom_vci prot_atom, atom_vci lig_atom,
                                my_float_t sq_dist)
{
  my_float_t cos_theta, cos_delta;
  if(A_hbond_geometry.salt_bridge(*prot_atom, *lig_atom, sq_dist)){
    A_ligand_flag[lig_atom] = SALT_BRIDGE;
    A_target_flag[prot_atom] = SALT_BRIDGE;
    salt_bridge_t sb;
    sb.A = prot_atom;
    sb.B = lig_atom;
    sb.sq_dist = sq_dist;
    A_prot_lig_SB.push_back(sb);
#if 0
  std::cout << "\n" << PDB_residues::residue_to_string(prot_atom->res) << " "
            << prot_atom->res_num << " " 
            << PDB_residues::atom_to_string(prot_atom->name) 
            << " | " << lig_atom->name_str  << "|"
            << std::sqrt(sq_dist) << "| ** SB ** \n";
#endif
  }else if(A_hbond_geometry.prot_lig_hbond(prot_atom, lig_atom, *A_prot, *A_lig,
                                           sq_dist, &cos_theta, &cos_delta)){
    A_ligand_flag[lig_atom] = DIRECT_HBOND;
    A_target_flag[prot_atom] = DIRECT_HBOND;
    hbond_t hb;
    hb.A = prot_atom;
    hb.B = lig_atom;
    hb.sq_dist = sq_dist;
    hb.cos_DHA_angle = cos_theta;
    hb.cos_preacc_angle = cos_delta;
    A_prot_lig_hbonds.push_back(hb);
  }
}

void
ProtLigScore::flag_intra_ligand_hbonds()
{
  my_float_t cos_DHA, cos_H__A_AA;

  atom_act_mi lig_atom_act_iter = A_ligand_flag.begin();;
  for( ; lig_atom_act_iter != A_ligand_flag.end(); ++lig_atom_act_iter){
    atom_vci cur_atom = lig_atom_act_iter->first;
    atom_act_class_t &cur_act_class = lig_atom_act_iter->second;
    if(cur_act_class != INTERFACIAL_ATOM || cur_atom->act_type == NOTHING || 
       cur_atom->act_type == HYDROPHOB) continue;

    atom_vci lig_atom = A_lig->atoms_begin();
    for( ; lig_atom < A_lig->atoms_end(); ++lig_atom){
      if(lig_atom->name == H || cur_atom == lig_atom) continue;

      my_float_t sq_dist = dist_squared(cur_atom->pos, lig_atom->pos);
      if(A_hbond_geometry.salt_bridge(*cur_atom, *lig_atom, sq_dist)){
        cur_act_class = INTRA_LIGAND_SALT_BRIDGE;

        atom_act_mi lig_atom_iter = A_ligand_flag.find(lig_atom);
        if(lig_atom_iter->second == INTERFACIAL_ATOM)
          lig_atom_iter->second = INTRA_LIGAND_SALT_BRIDGE;
      }else if(A_hbond_geometry.intra_lig_hbond(cur_atom, lig_atom, sq_dist,
                                                *A_lig, &cos_DHA, 
                                                &cos_H__A_AA)){
        cur_act_class = INTRA_LIGAND_HBOND;
        atom_act_mi lig_atom_iter = A_ligand_flag.find(lig_atom);
        if(lig_atom_iter->second == INTERFACIAL_ATOM)
          lig_atom_iter->second = INTRA_LIGAND_HBOND;
      }
    }
  }
}

void
ProtLigScore::flag_intra_protein_hbonds()
{
  my_float_t cos_DHA, cos_H__A_AA;

  atom_act_mi cur_atom_act_iter = A_target_flag.begin();;
  for( ; cur_atom_act_iter != A_target_flag.end(); ++cur_atom_act_iter){
    atom_vci cur_atom = cur_atom_act_iter->first;
    atom_act_class_t &cur_act_class = cur_atom_act_iter->second;
    if(cur_act_class != INTERFACIAL_ATOM || cur_atom->act_type == NOTHING || 
       cur_atom->act_type == HYDROPHOB) continue;
 
    residue_vci res = A_prot->residues_begin();
    for(res = A_prot->residues_begin(); res < A_prot->residues_end(); ++res){
      atom_vci prev_prot_atom = atom_t::NULL_ATOM_VCI; 
      atom_vci prot_atom = res->atoms_begin;
      for( ; prot_atom < res->atoms_end; ++prot_atom){
        if(prot_atom->name == H) continue;

        my_float_t sq_dist = dist_squared(cur_atom->pos, prot_atom->pos);
        if(A_hbond_geometry.salt_bridge(*cur_atom, *prot_atom, sq_dist)){
          cur_act_class = INTRA_TARGET_SALT_BRIDGE;

          atom_act_mi prot_atom_iter = A_target_flag.find(prot_atom);
          if(prot_atom_iter->second == INTERFACIAL_ATOM)
            prot_atom_iter->second = INTRA_TARGET_SALT_BRIDGE;
        }else if(A_hbond_geometry.intra_prot_hbond(cur_atom, prot_atom, sq_dist,
                                                  *A_prot, &cos_DHA, 
                                                  &cos_H__A_AA)){
          cur_act_class = INTRA_TARGET_HBOND;
          atom_act_mi prot_atom_iter = A_target_flag.find(prot_atom);
// this condition does not exist in other code
// It is required to omit it in order to have both codes give the same scores
// Unfortunately, this means that it is possible that protein atoms making
// intraprotein hbonds and are more than 4.5 (A) away from
// ligand heavy atoms still contribute to the total_target_hydro term if the
// other protein atom participating in the hbond is less than 4.5 (A) from at
// least 1 ligand atom.
          //if(prot_atom_iter->second == INTERFACIAL_ATOM)
            prot_atom_iter->second = INTRA_TARGET_HBOND;
        }
      }
    }
  }
}

void
ProtLigScore::sum_target_hydrophobicity()
{
  atom_act_mci atom_act_iter = A_target_flag.begin();
  for( ; atom_act_iter != A_target_flag.end(); ++atom_act_iter)
    if(atom_act_iter->second != INITIAL)
      A_raw_terms.total_target_hydro += (400.0 - atom_act_iter->first->hydro);
  A_raw_terms.total_target_hydro *= 0.001;
}

void
ProtLigScore::flag_unsat_polar_atoms(atom_act_map *my_map, 
                                     int *num_unsat_charge,
                                     int *num_unsat_polar)
{
  *num_unsat_charge = 0;
  *num_unsat_polar = 0;
  for(atom_act_mi iter = my_map->begin(); iter != my_map->end(); ++iter){
    atom_vci atom = iter->first;
    atom_act_class_t &act_class = iter->second;
    if(act_class != INTERFACIAL_ATOM || atom->act_type == NOTHING ||
       atom->act_type == HYDROPHOB || atom->act_type == UNKNOWN_INTERACTION) 
      continue;

    if(atom->charge != 0){
      act_class = UNSAT_CHARGE;
      ++(*num_unsat_charge);
    }else{
      act_class = UNSAT_POLAR;
      ++(*num_unsat_polar);
    }
  }
}

void
ProtLigScore::print_hbonds(std::ostream &out) 
{
  std::cout << "protein - ligand hbonds\n";
  for(size_t i = 0; i < A_prot_lig_hbonds.size(); ++i)
  {
    hbond_t &hb = A_prot_lig_hbonds[i];
    out << PDB_residues::residue_to_string(hb.A->res) 
        << hb.A->res_num << " " << PDB_residues::atom_to_string(hb.A->name)
        << " --- " << hb.B->name_str << "()     " 
        << std::sqrt(hb.sq_dist) << " "
        << std::acos(hb.cos_DHA_angle) * 180.0 / M_PI << " "
        << std::acos(hb.cos_preacc_angle)  * 180.0 / M_PI<< "\n";
  }
}

void
ProtLigScore::print_atom_act_classes(atom_act_map *my_map)
{
  for(atom_act_mi iter = my_map->begin(); iter != my_map->end(); ++iter){
    atom_vci atom = iter->first;
    const atom_act_class_t &act_class = iter->second;
    if(act_class == INITIAL) continue;

    if(atom->name_str.length()) std::cout << atom->name_str;
    else
      std::cout << PDB_residues::residue_to_string(atom->res) << " " 
                << atom->res_num 
                << " " << PDB_residues::atom_to_string(atom->name);

    switch(act_class){
    case INITIAL:
      std::cout << "  INITIAL\n";
      break;
    case INTERFACIAL_ATOM:
      std::cout << "  INTERFACIAL_ATOM\n";
      break;
    case METAL_DIRECT_HBOND:
      std::cout << "  METAL_DIRECT_HBOND\n";
      break;
    case SALT_BRIDGE:
      std::cout << "  SALT_BRIDGE\n";
      break;
    case DIRECT_HBOND:
      std::cout << "  DIRECT_HBOND\n";
      break;
    case UNSAT_CHARGE:
      std::cout << "  UNSAT_CHARGE\n";
      break;
    case UNSAT_POLAR:
      std::cout << "  UNSAT_POLAR\n";
      break;
    case INTRA_LIGAND_HBOND:
      std::cout << "  INTRA_LIGAND_HBOND\n";
      break;
    case INTRA_LIGAND_SALT_BRIDGE:
      std::cout << "  INTRA_LIGAND_SALT_BRIDGE\n";
      break;
    case INTRA_TARGET_HBOND:
      std::cout << "  INTRA_TARGET_HBOND\n";
      break;
    case INTRA_TARGET_SALT_BRIDGE:
      std::cout << "  INTRA_TARGET_SALT_BRIDGE\n";
      break;
    default:
      std::cout << "  *****\n";
      break;
    }

  }
}
