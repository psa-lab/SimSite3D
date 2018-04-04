
#include <HbondGeometry.H>
#include <cstdio>

#define TRACE 0

using namespace ASCbase;

const my_float_t HbondGeometry::MIN_SALT_BRIDGE_SQUARED_LENGTH = 2.5 * 2.5; 
const my_float_t HbondGeometry::MAX_SALT_BRIDGE_SQUARED_LENGTH = 4.5 * 4.5;
const my_float_t HbondGeometry::MIN_HBOND_SQUARED_LENGTH = 2.5 * 2.5; 
const my_float_t HbondGeometry::MAX_HBOND_SQUARED_LENGTH = 3.5 * 3.5;

const my_float_t HbondGeometry::MIN_DHA_ANGLE = 120.0; // a bit steep
const my_float_t HbondGeometry::MIN_H_A_AA_ANGLE = 90.0; // a bit steep for sp3,
// but maybe the occurance of less than 90.0 degrees for sp3 is limited anyhow 
// (at least if we were considering D..A-AA instead of H..A-AA.

const my_float_t HbondGeometry::COSINE_MIN_DHA_ANGLE = 
  std::cos(deg2rad(HbondGeometry::MIN_DHA_ANGLE));
const my_float_t HbondGeometry::COSINE_MIN_H_A_AA_ANGLE = 
  std::cos(deg2rad(HbondGeometry::MIN_H_A_AA_ANGLE));

const my_float_t HbondGeometry::A_SC_H_positions[] = {
  -0.463070309734338, 0.908826654672135, 0.0,  // Arg NH1
  -0.463070309734338, 0.908826654672135, 0.0,  // Arg NH2
  -0.431070626975513, 0.924433942777383, 0.0,  // Asn ND2
  -0.431070626975513, 0.924433942777383, 0.0,  // Gln NE2
  -0.414871375937316, 0.931816366795453, 0.0,  // Lys NZ
  -0.280676836533827, 0.918052565724514, 0.0,  // Ser OG
  -0.264611861584319, 0.922811228100786, 0.0,  // Thr OG1
  -0.296656314599950, 0.913014255643347, 0.0   // Tyr OH
};

const my_float_t HbondGeometry::A_lig_H_positions[] = {
  -0.337144927826108, 0.952067906003100, 0.0, // SP3/SP4 Nitrogen TETRAHEDRAL
  -0.320454584864420, 0.904935831448491, 0.0, // SP3 Oxygen TETRAHEDRAL
  -0.337144927826108, 0.0, 0.952067906003100  // N.sp3 TETRAHEDRAL_2
};

const my_float_t HbondGeometry::A_main_chain_NH_position[] = {
  -0.414871375937316, 0.931816366795453, 0
};

const size_t HbondGeometry::A_SC_donors_array_size = 12;
const HbondGeometry::sidechain_donor_t HbondGeometry::A_SC_donors_array[] = {
  { ARG, { CD, NE, CZ}, NULL },
  { ARG, { CZ, NH1, NH2}, &A_SC_H_positions[0] },
  { ARG, { CZ, NH2, NH1}, &A_SC_H_positions[3] },
  { ASN, { CG, ND2, OD1}, &A_SC_H_positions[6] },
  { GLN, { CD, NE2, OE1}, &A_SC_H_positions[9] },
  { HIS, { CG, ND1, CE1}, NULL },
  { HIS, {CE1, NE2, CD2}, NULL },
  { LYS, { CE, NZ, NULL_ATOM_TYPE}, &A_SC_H_positions[12] },
  { SER, { CB, OG, NULL_ATOM_TYPE}, &A_SC_H_positions[15] },
  { THR, { CB, OG1, NULL_ATOM_TYPE}, &A_SC_H_positions[18] },
  { TRP, {CD1, NE1, CE2}, NULL },
  { TYR, { CZ, OH, NULL_ATOM_TYPE}, &A_SC_H_positions[21] }
};
HbondGeometry::SC_tbl_res_lvl HbondGeometry::A_SC_donors;

const size_t HbondGeometry::A_SC_preacc_array_size = 13;
const HbondGeometry::sidechain_preacc_t HbondGeometry::A_SC_preacc_array[] = {
  { ASN, OD1,  CG },
  { ASP, OD1,  CG },
  { ASP, OD2,  CG },
  { CYS,  SG,  CB },
  { GLN, OE1,  CD },
  { GLU, OE1,  CD },
  { GLU, OE2,  CD },
  { HIS, ND1,  CG },  // Do we want to check CE1 as well?
  { HIS, NE2, CD2 },  // Do we want to check CE1 as well?
  { MET,  SD,  CG },  // Do we want to check CE as well?
  { SER,  OG,  CB },
  { THR, OG1,  CB },
  { TYR,  OH,  CZ }
};
HbondGeometry::preacc_tbl_res_lvl HbondGeometry::A_SC_preaccs;

bool
HbondGeometry::intra_prot_hbond(atom_vci A, atom_vci B, my_float_t sq_dist,
                                PDBStructure &prot, my_float_t *cos_theta, 
                                my_float_t *cos_delta)
{
  *cos_theta = *cos_delta = 1.0;
  my_float_t H_pos[3];

  if(sq_dist < MIN_HBOND_SQUARED_LENGTH || MAX_HBOND_SQUARED_LENGTH < sq_dist ||
     !is_hbond_interaction(A->act_type, B->act_type)) 
    return false;

  atom_vci acceptor = A;
  atom_vci donor = B;
  if((A->act_type == DONOR || A->act_type == DONEPTOR) && B->act_type != DONOR){
   donor = A;
   acceptor = B;
  }

  bool rv = false;
  if(protein_donor_cos_theta(donor, acceptor, prot, H_pos, cos_theta) &&
     prot_acceptor_cos_delta(acceptor, prot, H_pos, cos_delta)) rv = true;

#if TRACE
  std::cout << PDB_residues::residue_to_string(A->res) << " "
            << A->res_num << " " 
            << PDB_residues::atom_to_string(A->name) 
            << " | " << PDB_residues::residue_to_string(B->res) << " "
            << B->res_num << " " 
            << PDB_residues::atom_to_string(B->name)  << "|"
            << std::sqrt(sq_dist) << " " << std::acos(*cos_theta) * 180.0/ M_PI
            << " " << std::acos(*cos_delta) * 180.0/M_PI << "\n\n";

#endif
  return rv;
}

bool
HbondGeometry::intra_lig_hbond(atom_vci A, atom_vci B, my_float_t sq_dist,
                               mol2File &ligand, my_float_t *cos_theta,
                               my_float_t *cos_delta)
{
  *cos_theta = *cos_delta = 1.0;
  my_float_t H_pos[3];

  if(sq_dist < MIN_HBOND_SQUARED_LENGTH || MAX_HBOND_SQUARED_LENGTH < sq_dist ||
     !is_hbond_interaction(A->act_type, B->act_type)) 
    return false;

  atom_vci acceptor = A;
  atom_vci donor = B;
  if((A->act_type == DONOR || A->act_type == DONEPTOR) && B->act_type != DONOR){
   donor = A;
   acceptor = B;
  }

  if(ligand_donor_cos_theta(donor, acceptor, ligand, H_pos, cos_theta) &&
     lig_acceptor_cos_delta(acceptor, ligand, H_pos, cos_delta)) return true;
  return false;
} 

bool
HbondGeometry::prot_lig_hbond(atom_vci prot_atom, atom_vci lig_atom,
                              PDBStructure &prot, mol2File &ligand,
                              my_float_t sq_dist, 
                              my_float_t *cos_theta, my_float_t *cos_delta)
{
  *cos_theta = *cos_delta = 1.0;
  my_float_t H_pos[3];

  if(sq_dist < MIN_HBOND_SQUARED_LENGTH || MAX_HBOND_SQUARED_LENGTH < sq_dist ||
     !is_hbond_interaction(prot_atom->act_type, lig_atom->act_type)) 
    return false;

  bool rv = false;
  // Protein atom is the proton donor
  if(lig_atom->act_type == ACCEPTOR ||
     (lig_atom->act_type == DONEPTOR && 
      (prot_atom->act_type == DONOR || prot_atom->res == PDB_METAL))){
 
    if(protein_donor_cos_theta(prot_atom, lig_atom, prot, H_pos, cos_theta) &&
       lig_acceptor_cos_delta(lig_atom, ligand, H_pos, cos_delta)){
      rv =  true;
    }

  // Ligand atom is the proton donor
  }else{
    if(ligand_donor_cos_theta(lig_atom, prot_atom, ligand, H_pos, cos_theta) &&
       prot_acceptor_cos_delta(prot_atom, prot, H_pos, cos_delta)){
      rv = true;
      }
  }
#if TRACE
  std::cout << PDB_residues::residue_to_string(prot_atom->res) << " "
            << prot_atom->res_num << " " 
            << PDB_residues::atom_to_string(prot_atom->name) 
            << " | " << lig_atom->name_str  << "|"
            << std::sqrt(sq_dist) << " " << std::acos(*cos_theta) * 180.0/ M_PI
            << " " << std::acos(*cos_delta) * 180.0/M_PI;
  if(rv) std::cout << " -- is hbond\n\n";
  else std::cout << " -- NO\n\n";
#endif
  return rv;
}

bool
HbondGeometry::protein_donor_cos_theta(atom_vci prot_atom, atom_vci acceptor,
                                       PDBStructure &prot, my_float_t *H_pos,
                                       my_float_t *cos_theta)
{
  *cos_theta = 1.0;
  std::fill(H_pos, H_pos + 3, 0.0);

  // assume we have hbond -- should check before entering 
  residue_vci res = prot.get_residue(prot_atom);
  if(res == PDBStructure::NULL_RESIDUE_VCI) return false;

  const my_float_t *local_H_pos = 0;
  std::vector<atom_vci> atoms(3, atom_t::NULL_ATOM_VCI);
  const sidechain_donor_t *SC_donor = NULL;
  atoms[1] = prot_atom;

  // Main chain nitrogen
  if(prot_atom->name == N){
    atoms[2] = res->get_atom(CA);
    residue_vci prev_res = (res == prot.residues_begin() ? res : res - 1);
    atom_vci C_atom = prev_res->get_atom(C);

    // N-terminus or chain break -- allow H to rotate to point to acceptor
    //
    // Residue number seems to be unreliable
    //if(res == prot.residues_begin() || res->number - 1 != prev_res->number){
    // N-C bond length is 1.32 (A) -- allow up to 1.5 (A) to be nice
    if(res == prot.residues_begin() || atom_t::NULL_ATOM_VCI == C_atom || 
       2.25 < dist_squared(C_atom->pos, prot_atom->pos)){
      atoms[0] = atoms[2];
      atoms[2] = atom_t::NULL_ATOM_VCI;
      local_H_pos = A_main_chain_NH_position;

    // Current residue is bonded to the previous residue
    }else{
      atoms[0] = C_atom;
      // If this residue does not have an alpha-carbon (Litian He 03/17/2004)
      // use the C & O from previous residue
      // Note: this sounds dubious because without alpha carbons many tools
      // will choke and numerous other issues will result as well from missing
      // alpha carbons (Jeff, Jan 20, 2009).
      if(atoms[2] == atom_t::NULL_ATOM_VCI) atoms[2] = prev_res->get_atom(O);
    }
  }else if((SC_donor = get_sidechain_donor(prot_atom->res, prot_atom->name))){
    local_H_pos = SC_donor->local_H_pos;
    atoms[0] = res->get_atom(SC_donor->atoms[0]);
    // Planar nitrogen
    if(SC_donor->atoms[2] != NULL_ATOM_TYPE)
      atoms[2] = res->get_atom(SC_donor->atoms[2]);
    // Sp3 atoms
    else atoms[2] = atom_t::NULL_ATOM_VCI;
  }else{
    std::cerr << "unsupported protein side chain donor atom\n";
    std::cerr << PDB_residues::residue_to_string(prot_atom->res)
              << prot_atom->res_num
              << " " << PDB_residues::atom_to_string(prot_atom->name) << "\n";
    return false;
  }

  // If there is 1 hydrogen atom and it lies in the plane defined by the
  // three heavy atoms, we can bisect the angle atoms[0] - atoms[1] - atoms[2]
  // and project the position of the hydrogen using the vector of bisection
  // on the oblique angle side at a distance length away from atoms[1].
  //
  // Seems I had forgotten that the two bonds (vectors) need not have the
  // same length.  
  if(!local_H_pos){
    // All should be nitrogens -- length is 1.02 for all N-H bonds in proteins
    my_float_t length = 1.02;
    my_float_t cos_angle = 1.0;
    cosine_angle(atoms[0]->pos, atoms[1]->pos, atoms[2]->pos, &cos_angle);
    
    // We might be able to reduce the computation later using trig identities
    my_float_t NH_angle = M_PI - std::acos(cos_angle) * 0.5;
    my_float_t my_H_pos[3]; 
    my_H_pos[0] = length * std::cos(NH_angle);
    my_H_pos[1] = length * std::sin(NH_angle);
    my_H_pos[2] = 0.0;

    // Note that there is only 1 hydrogen position in this case
    compute_H_position(atoms[0], atoms[1], atoms[2], acceptor, my_H_pos, 
                       FLIP_NBR_Y, H_pos);

  }else
    // Note here we either have sp2, sp3 or pl3 types -- i.e. choose the better
    // position for the acceptor or point the hydrogen to the acceptor
    compute_H_position(atoms[0], atoms[1], atoms[2], acceptor, local_H_pos, 
                       FLIP_ACC_Y, H_pos);

  cosine_angle(acceptor->pos, H_pos, prot_atom->pos, cos_theta);        
#if TRACE
  std::cout << "H pos: " << H_pos[0] << " " << H_pos[1] << " " << H_pos[2] 
            << "\n";
#endif
  if(*cos_theta > COSINE_MIN_DHA_ANGLE) return false;
  return true;
}

bool
HbondGeometry::ligand_donor_cos_theta(atom_vci lig_atom, atom_vci acceptor, 
                                      mol2File &ligand, my_float_t *H_pos,
                                      my_float_t *cos_theta)
{
  // Set things to arbitrary values to help with debugging
  *cos_theta = 1.0; 
  std::fill(H_pos, H_pos + 3, 0.0);

  const std::vector<atom_vci> &lig_nbrs = ligand.get_nbrs(lig_atom);

  // Count number of heavy atoms bound to lig_atom and keep any of the heavy
  // atom neighbors (if there is at least 1)
  int num_H_nbrs = 0;
  std::vector<atom_vci> hvy_nbrs;
  atom_vci H_nbr = atom_t::NULL_ATOM_VCI;
  std::vector<atom_vci>::const_iterator nbr_iter;
  for(nbr_iter = lig_nbrs.begin(); nbr_iter < lig_nbrs.end(); ++nbr_iter){
    if((*nbr_iter)->name == H){ 
      H_nbr = *nbr_iter;
      ++num_H_nbrs;
    }else hvy_nbrs.push_back(*nbr_iter);
  }

  // Kinda keep given hydrogen positions for ar, am, pl3 & sp2. 
  if(lig_atom->name == N && 
     (lig_atom->orbit == AR || lig_atom->orbit == AMIDE || 
      lig_atom->orbit == PL3 || lig_atom->orbit == SP2)){

    // Allow the H(s) to rotate if the donor is bonded to only 1 heavy 
    // neighbor atom -- we can use the given H positions if we use only the
    // second part of the if statement.
    if(hvy_nbrs.size() == 1 && num_H_nbrs > 0){
      my_float_t local_H_pos[3];
      compute_local_H_pos(hvy_nbrs.front(), lig_atom, H_nbr, local_H_pos);
      compute_H_position(hvy_nbrs.front(), lig_atom, atom_t::NULL_ATOM_VCI, 
                         acceptor, local_H_pos, FLIP_ACC_Y, H_pos);
      cosine_angle(lig_atom->pos, H_pos, acceptor->pos, cos_theta);        

    // Use the H(s) as given in the ligand
    }else{
      // Compute the "best" angle -- (we seek as close to 180 degrees as
      // possible) the closer cosine is to -1.0 the better the angle.
      for(nbr_iter = lig_nbrs.begin(); nbr_iter < lig_nbrs.end(); ++nbr_iter)
        if((*nbr_iter)->name == H){
          my_float_t tmp;
          cosine_angle(lig_atom->pos, (*nbr_iter)->pos, acceptor->pos, &tmp);
          if(tmp < *cos_theta){ 
            *cos_theta = tmp;
            std::copy((*nbr_iter)->pos, (*nbr_iter)->pos + 3, H_pos);
          }
        }
    }
  }else{
    
    /////////////////////////
    // NOTE:
    // Need to actually move the hydrogens if other scoring depends on it or we
    // wish to write out the ligands with move hydrogen positions
    /////////////////////////
   
    if(hvy_nbrs.size() >= 4){ 
      std::cerr << "The hydrogen bond donor atom ("
                << ligand.atom_type_to_str(lig_atom) 
                << ") -- has more than 3 heavy neighbors; \n"
                << "  it cannot make a hydrogen bond\n";
      return false ;

    // TETRAHEDRAL_3_NEIGHBORS; 
    // Hydrogen atom position is fixed if N.sp4
    }else if(hvy_nbrs.size() == 3){
      if(lig_atom->orbit != SP4){
        std::cerr << "The hydrogen bond donor atom ("
                  << ligand.atom_type_to_str(lig_atom) 
                  << ") -- has 3 heavy neighbors; \n"
                  << "  it cannot make a hydrogen bond\n";
        return false ;
      }

#if 0
      // Use the given hydrogen position -- as with the previous case
      // when computing the centroid will not correspond to a vector/position
      // that is as close to tetrahedral as possible since we may have 
      // different bond lengths depending on the heavy neighbors.
      
      // Compute the centroid of the 3 hvy neighbors
      my_float_t centroid[3];
      std::fill(centroid, centroid + 3, 0.0);
      for(nbr_iter == lig_nbrs.begin(); nbr_iter < lig_nbrs.end(); ++nbr_iter)
        if((*nbr_iter)->name != H)
          my_axpy(3, 1.0, (*nbr_iter)->pos, 1, centroid, 1);
      for(size_t i = 0; i < 3; ++i) centroid[i] /= 3.0;

      // and push out
      // the hydrogen along the vector from the centroid to the donor atom
      // at a distance of about 1.0 (A).
      my_float_t H_dir[3];
      unit_vector(H_dir, lig_atom->pos, centroid);
#endif      
      cosine_angle(lig_atom->pos, H_nbr->pos, acceptor->pos, cos_theta); 
      std::copy(H_nbr->pos, H_nbr->pos + 3, H_pos);

    // TETRAHEDRAL_2_NEIGHBORS; 
    }else if(hvy_nbrs.size() == 2){
      // Special case -- need to place down the 1 or 2 tetrahedral hydrogen 
      // atoms they are fixed in this case -- but we are in the local XZ plane 
      // and not the local XY plane
      // We don't particularly care at this point whether there are 1 or 2
      // hydrogen atoms (we assume the hydrogen will move to make an hbond if
      // possible -- could be quite incorrect, but wth).  Also, if we want to 
      // write out the ligands with "moved" hydrogen atoms we will need to pay
      // closer attention
      compute_H_position(hvy_nbrs.front(), lig_atom, hvy_nbrs.back(),
                         acceptor, &A_lig_H_positions[6], FLIP_Z, H_pos);
      cosine_angle(lig_atom->pos, H_pos, acceptor->pos, cos_theta);        

    // TETRAHEDRAL
    }else if(hvy_nbrs.size() == 1){
      if(lig_atom->name == N && (lig_atom->orbit == SP3 || 
                                 lig_atom->orbit == SP4))
        compute_H_position(hvy_nbrs.front(), lig_atom, atom_t::NULL_ATOM_VCI,
                           acceptor, &A_lig_H_positions[0], FLIP_ACC_Y, H_pos);
      else if(lig_atom->name == O && lig_atom->orbit == SP3)
        compute_H_position(hvy_nbrs.front(), lig_atom, atom_t::NULL_ATOM_VCI,
                           acceptor, &A_lig_H_positions[3], FLIP_ACC_Y, H_pos);
      else{
        std::cerr << "Unexpected hydrogen bond donor atom (" 
                  << ligand.atom_type_to_str(lig_atom)
                  << ")-- cannot compute hydrogen bonding geometry\n";
        return false;
      }
      cosine_angle(lig_atom->pos, H_pos, acceptor->pos, cos_theta);        

    // In keeping with the current line of thinking, if the atom does not have
    // any heavy neighbors we allow it to rotate as much as it likes and allow
    // the theta angle to be 180 degrees.
    }else{
      *cos_theta = -1.0;
      my_float_t DH_dir[3];
      unit_vector(DH_dir, acceptor->pos, lig_atom->pos);
      std::copy(lig_atom->pos, lig_atom->pos + 3, H_pos);
      if(lig_atom->name == O) my_axpy(3, 0.96, DH_dir, 1, H_pos, 1);
      else if(lig_atom->name == N) my_axpy(3, 1.02, DH_dir, 1, H_pos, 1);
      else{
        std::cerr << "Unexpected hydrogen bond donor atom (" 
                  << ligand.atom_type_to_str(lig_atom)
                  << ")-- cannot compute hydrogen bonding geometry\n";
        return false;
      }
    }
  }

#if TRACE
  std::cout << "H pos: " << H_pos[0] << " " << H_pos[1] << " " << H_pos[2] 
            << "\n";
#endif
  if(*cos_theta > COSINE_MIN_DHA_ANGLE) return false;
  return true;
}

bool 
HbondGeometry::prot_acceptor_cos_delta(atom_vci prot_acceptor, 
                                       PDBStructure &prot, my_float_t *H_pos,
                                       my_float_t *cos_delta)
{
  // NOTE: do NOT modify H_pos

  *cos_delta = 1.0;
  residue_vci res = prot.get_residue(prot_acceptor);
  const sidechain_preacc_t *SC_info = NULL;
  atom_vci AA;
  if(prot_acceptor->name == O || prot_acceptor->name == OXT) 
    AA = res->get_atom(C);
  else{
    SC_info = get_sidechain_preacc(prot_acceptor->res, prot_acceptor->name);
    if(SC_info) AA = res->get_atom(SC_info->preacceptor);
    else{
      std::cerr << "Could find sidechain preacceptor for "
                << prot_acceptor->res << prot_acceptor->res_num << " "
                << prot_acceptor->name << "\n"
                << "Ignoring preacceptor angle conditions for this hbond\n";
      *cos_delta = -1.0;
       return true;
    }
  }
  cosine_angle(AA->pos, prot_acceptor->pos, H_pos, cos_delta); 
  if(*cos_delta <= COSINE_MIN_H_A_AA_ANGLE) return true;
  return false;
}

bool
HbondGeometry::lig_acceptor_cos_delta(atom_vci lig_acceptor, mol2File &ligand,
                                      my_float_t *H_pos, my_float_t *cos_delta)
{
  // NOTE: do NOT modify H_pos

  // Preacceptor check -- allow any neighboring ligand heavy atom to be a 
  // preacceptor -- will need to search for the best if we want to use
  // the Proflex energy function or similar idea
  *cos_delta = 1.0;
  const std::vector<atom_vci> &lig_nbrs = ligand.get_nbrs(lig_acceptor);
  std::vector<atom_vci>::const_iterator nbr_iter;
  for(nbr_iter = lig_nbrs.begin(); nbr_iter < lig_nbrs.end(); ++nbr_iter)
    if((*nbr_iter)->name != H){
      cosine_angle((*nbr_iter)->pos, lig_acceptor->pos, H_pos, cos_delta);
      if(*cos_delta <= COSINE_MIN_H_A_AA_ANGLE) return true;
    }
  return false;
}


void
HbondGeometry::compute_H_position(atom_vci hvy_nbr_A, atom_vci donor, 
                                  atom_vci hvy_nbr_B, atom_vci acceptor, 
                                  const my_float_t *local_H_pos, 
                                  const local_coord_flip coord_flip, 
                                  my_float_t *new_H_pos)
{
  // 1) Compute translation of donor to origin & the following rotation to 
  // get acceptor and hvy_nbr in the XY plane with hvy_nbr on the 
  // positive X-axis
#if 0
  std::cout << "hvy nbr pos: " << hvy_nbr->pos[0] << " "
    << hvy_nbr->pos[1] << " " << hvy_nbr->pos[2] << "\n";
  std::cout << "donor pos: " << donor->pos[0] << " "
    << donor->pos[1] << " " << donor->pos[2] << "\n";
  std::cout << "acceptor pos: " << acceptor->pos[0] << " "
    << acceptor->pos[1] << " " << acceptor->pos[2] << "\n";
#endif

  my_float_t UU[3], VV[3], AA[3];
  for(uint i = 0; i < 3; i++){
    UU[i] = hvy_nbr_A->pos[i] - donor->pos[i];
    if(hvy_nbr_B != atom_t::NULL_ATOM_VCI)
      VV[i] = hvy_nbr_B->pos[i] - donor->pos[i];
    else VV[i] = acceptor->pos[i] - donor->pos[i];
    AA[i] = acceptor->pos[i] - donor->pos[i];
  }
  my_float_t R[9];
  // Premultiply R^t by global coord to get to local coordinate system
  // Premulitply R by local coord to get to global coordinates
  get_local_orientation(UU, VV, R);
#if 0
  std::cout << "\nUU:" << UU[0] << ", " << UU[1] << ", " << UU[2] << "\n";
  std::cout << "VV:" << VV[0] << ", " << VV[1] << ", " << VV[2] << "\n";
  std::cout << "R: " << R[0] << ", " <<  R[1] << ", " << R[2] << ",\n";
  std::cout << "   " << R[3] << ", " <<  R[4] << ", " << R[5] << ",\n";
  std::cout << "   " << R[6] << ", " <<  R[7] << ", " << R[8] << "\n";
  std::cout << "Local H: " << local_H_pos[0] << ", " << local_H_pos[1] << ", "
            << local_H_pos[2] << "\n";
#endif

  // 2) The donor - C* bond lies on the X-axis.  Check which half plane contains
  // the heavy neighbor atom or acceptor and move the local position to the 
  // global H position.
  // Assumption: local H position is always in the 2nd quadrant of the XY plane
  my_float_t tmp[3];
  std::copy(local_H_pos, local_H_pos + 3, tmp);
  my_float_t val = 0.0;
  // Flip away from hvy_nbr_B
  if(coord_flip == FLIP_NBR_Y){
    for(size_t i = 0; i < 3; ++i) val += VV[i] * R[3 + i];
    if(val > 0) tmp[1] *= -1.0;
  // Flip towards Acceptor
  }else if(coord_flip == FLIP_ACC_Y){
    for(size_t i = 0; i < 3; ++i) val += AA[i] * R[3 + i];
    if(val < 0) tmp[1] *= -1.0;
  // Flip towards acceptor
  }else if(coord_flip == FLIP_Z){
    for(size_t i = 0; i < 3; ++i) val += AA[i] * R[6 + i];
    if(val < 0) tmp[2] *= -1.0;
  }
  std::copy(donor->pos, donor->pos + 3, new_H_pos);
  my_gemm(1, 3, 3, 1.0, tmp, 3, R, 3, new_H_pos, 3, 1.0);
}  

bool
HbondGeometry::compute_local_H_pos(atom_vci hvy_nbr, atom_vci lig_atom, 
                                   atom_vci H_atom, my_float_t *local_H_pos)
{
  // 1) Get angle between H-N-C* where C* is whatever heavy atom
  //    is bound to N.
  my_float_t cos_angle;
  my_float_t bond_len = 
    cosine_angle(H_atom->pos, lig_atom->pos, hvy_nbr->pos, &cos_angle);
  if(cos_angle >= 0.0){
    std::cerr << "H-N-C* angle is 90 degrees or less\n";
    return false;
  }

  // 2) Calculate the H position(s) in the local plane 
  // (2nd quadrant of XY plane)
  //std::cout << "computeing local H position\n";
  //local_H_pos[0] = -1.0 * cos_angle * bond_len;
  local_H_pos[0] = cos_angle * bond_len;
  local_H_pos[1] =  
    std::sqrt(bond_len*bond_len - local_H_pos[0]*local_H_pos[0]);
  local_H_pos[2] = 0;  // in the XY plane
  return true;
}

my_float_t
HbondGeometry::cosine_angle(my_float_t *A_pos, my_float_t *B_pos, 
                            my_float_t *C_pos, my_float_t *cos_angle)
{
  my_float_t BA_vec[3], BC_vec[3]; 
  my_float_t len_BA = unit_vector(BA_vec, A_pos, B_pos);
  unit_vector(BC_vec, C_pos, B_pos);
  *cos_angle = dot(BA_vec, BC_vec);
  return len_BA;
}

void
HbondGeometry::build_SC_donor_tables()
{
  if(A_SC_donors.size()) return;
  for(size_t i = 0; i < A_SC_donors_array_size; ++i)
    A_SC_donors[A_SC_donors_array[i].residue][A_SC_donors_array[i].atoms[1]] =
      A_SC_donors_array + i;
}

void
HbondGeometry::build_SC_preacc_tables()
{
  if(A_SC_preaccs.size()) return;
  for(size_t i = 0; i < A_SC_preacc_array_size; ++i)
    A_SC_preaccs[A_SC_preacc_array[i].residue][A_SC_preacc_array[i].acceptor] =
      A_SC_preacc_array + i;
}
