
#include <ProtLigActTable.H>

using namespace ASCbase;

void
ProtLigActTable::build_hphob_SC_to_atom_map(PDBStructure &prot)
{
  if(A_hphob_SC_to_atom.size()) return;

  typedef std::pair<residue_vci, atom_vci> SC_atom_pair_t;
  residue_vci res = prot.residues_begin();
  for(res = prot.residues_begin(); res < prot.residues_end(); ++res){
    if(res->name == GLY || res->name == SER) continue;

    A_hphob_SC_to_atom.insert(SC_atom_pair_t(res, res->get_atom(CB)));
  }
}
