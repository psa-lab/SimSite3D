import os
import tempfile
import re
from sys import stderr
from numpy import *
from ASCbasePy import *

# The openeye modules (wrappers to C++ libraries) are stored somewhere on
# the PSA system.  When writing this script they resided at
# /soft/linux64/openeye/wrappers/v1.7.2.4/python.  Pleae append their 
# current location to your $PYTHONPATH variable
from openeye.oechem import *
from openeye.oequacpac import *

def run(hetgrp, mol2_fname="", use_lig_expo=True):
  if(use_lig_expo):
    return __run(hetgrp, mol2_fname=mol2_fname)
  else:
    return __run_no_lig_expo(hetgrp, mol2_fname=mol2_fname)
################################################################################

################################################################################
def __run_no_lig_expo(hetgrp, mol2_fname=""):
  """
  Bond order is likely bogus -- this function exists for volumetric calcs 
  and should only be used if one does not have internet access to Ligand
  Expo
  """
  mol_charge = "/home/vanvoor4/Desktop/quacpac_test/openeye/bin/molcharge"

  # Write temporary pdb file
  (pdb_fd, pdb_out_fname) = tempfile.mkstemp(suffix=".pdb")
  pdb_file = os.fdopen(pdb_fd, "w+")

  for a in hetgrp.atoms:
    if(not a.is_pdb_hydrogen()):
      print >> pdb_file, a
  pdb_file.close()

  # molcharge it -- no am1bcc since we don't know bond order anyhow
  mol_charge_cmd = "%s -in %s -out %s >& /dev/null " \
    % (mol_charge, pdb_out_fname, mol2_fname)
  os.system(mol_charge_cmd)
  os.unlink(pdb_out_fname)
################################################################################

################################################################################
def __run(hetgrp, mol2_fname, charge_type="MMFF94", check=True):
  """
  Note: this function uses OpenEye Python wrappers to oechem & Quacpac 
  """

  # Get the ideal cif from the RCSB Ligand Expo
  lig_id = hetgrp.atoms[0].resName.upper().strip()
  if(not os.path.exists("%s.cif" % (lig_id))):
    os.system("wget http://ligand-expo.rcsb.org/reports/%s/%s/%s.cif" % \
              (lig_id[0], lig_id, lig_id))

  # Load the ideal cif and use it to convert the pdb hetgrp to sdf format
  my_cif = utils.pdb_cif.molecule("%s.cif" % (lig_id))
  if(my_cif.fail): 
    print "Unable to load:", lig_id
    return False
  my_chem_comp = my_cif.chemical_component
  (pdb_atom_names, my_sdf) = my_chem_comp.pdb_hetgrp_to_sdf(hetgrp)

  # Check if an error occured
  if(len(pdb_atom_names) == 0):
    # If PDB_atom_names is [], then we had a key mismatch
    return False
    

  # Read the SDF string into OE data structure
  oe_ims = oemolistream()
  oe_ims.SetFormat(OEFormat_MDL)
  oe_ims.openstring("%s\n" % (my_sdf))
  oe_mol_graph = OEGraphMol()
  OEReadMolecule(oe_ims, oe_mol_graph)

  # Assign partial charges and protons using the molcharge library
  if(charge_type.upper() == "AM1BCC"):
    charge_type = OECharges_AM1BCC
    print "am1bcc charges will be used"
  else:
    charge_type = OECharges_MMFF94
    print "MMFF94 charges will be used"
  noHydrogen = False
  debug = False
  OEAssignPartialCharges(oe_mol_graph, charge_type, noHydrogen, debug)

  # Assign the PDB names to the heavy atoms, and arbitrary naming to hydrogen
  # atoms
  aitr = oe_mol_graph.GetAtoms()
  hydrogen_count = 1
  hvy_count = 1
  while(aitr.IsValid()):
    a = aitr.Target()
    if(a.IsHydrogen()):
      aitr.Target().SetName("H%02d" % (hydrogen_count))
      hydrogen_count += 1
    else:
      aitr.Target().SetName(pdb_atom_names[hvy_count - 1])
      hvy_count += 1
    aitr.Next()

  # Use low level writer so that the atom names assigned in this script are
  # used
  ofs = oemolostream(mol2_fname)
  OEWriteMol2File(ofs, oe_mol_graph)

  if(check):
    my_mol2 = utils.mol2.molecule(mol2_fname)
    mol2_atoms = {}
    for atom in my_mol2.atoms:
      if(not atom.name == "H" and not atom.name == "D"): 
        mol2_atoms[atom.name_str] = array(atom.position)

    atom_names = []
    for atom in hetgrp.atoms:
      if(atom.is_pdb_hydrogen()): continue

      # check for alt loc
      id = atom.name.strip() 
      if(id in atom_names): continue
      atom_names.append(id)

      if(not id in mol2_atoms):
        print >> stderr, """
ERROR: Could not find atom with name %s in the mol2 file
""" % (id)
        return False

      pdb_pos = array(atom.position)
      my_diff = pdb_pos - mol2_atoms[id] 
      sq_dist = dot(my_diff, my_diff)
      # Be generous with distance tolerance -- should be zero
      if(sq_dist > 0.01):
        print >> stderr, """
ERROR: After conversion, the atom with the name %s in the mol2 file has
a significantly different position (%f).
""" % (id, sqrt(sq_dist))
        return False
    
  return True 
  
################################################################################

################################################################################
def __old_run(hetgrp, mol2_fname=""):
  """
  Bond order should be ok
  """
  babel = "/soft/linux64/openeye/examples/oechem-utilities/babel2"
  # test new version of molchareg
  #mol_charge = "/soft/linux/bin/molcharge"
  mol_charge = "/home/vanvoor4/Desktop/quacpac_test/openeye/bin/molcharge"

  # Get the ideal sdf from the RCSB Ligand Expo
  lig_id = hetgrp.atoms[0].resName.upper()
  if(not os.path.exists("%s_ideal.sdf" % (lig_id))):
    os.system("wget http://ligand-expo.rcsb.org/reports/%s/%s/%s_ideal.sdf" % \
              (lig_id[0], lig_id, lig_id))

  ideal_sdf = utils.sdf.molecule("%s_ideal.sdf" % (lig_id))
  mol_out = utils.sdf.molecule()
  unk_atom = False
  for a,b in zip(hetgrp.atoms, ideal_sdf.atoms):
# More nuisance PDB nonstandardness -- 1q0n at this time has the following
# lines
#HETATM 1364 HN61 PH2   181      -1.622   9.734  10.333  1.00  8.74            H
#HETATM 1365 HN62 PH2   181      -1.896   8.790  11.718  1.00  8.74            H
#HETATM 1367 H111 PH2   181      -7.074   6.066   5.108  1.00  9.55            H
#HETATM 1368 H112 PH2   181      -5.921   5.006   5.081  1.00  9.55            H
    if(b.element == "H" or a.is_pdb_hydrogen() or
       re.match("H[ONG][0-9]{2}", a.name) or re.match("H[0-9]{3}", a.name)):
      continue

    atom = utils.sdf.atom()
    if(a.name.startswith(" C") and b.element == "C"): atom.element = "C"
    elif(a.name.startswith(" N") and b.element == "N"): atom.element = "N"
    elif(a.name.startswith(" O") and b.element == "O"): atom.element = "O"
    elif(a.name.startswith(" P") and b.element == "P"): atom.element = "P"
    elif(a.name.startswith(" S") and b.element == "S"): atom.element = "S"
    elif(a.name.startswith("CL") and b.element == "CL"): atom.element = "CL"
    else:
      unk_atom = True
      print >> stderr, "Unknown atom type\nPDB: %s\nSDF: %s\n\tSkipping %s . . ." % (a, b, mol2_fname)
      break
 
    atom.position = a.position[:]
    atom.other = b.other
    mol_out.atoms.append(atom)

  if(unk_atom): return False

  # write out all bonds with both atom numbers <= mol_out.num_atoms
  num_atoms = len(mol_out.atoms)
  for bond in ideal_sdf.bonds:
    if(bond.atom1 <= num_atoms and bond.atom2 <= num_atoms):
      mol_out.bonds.append(bond)

  mol_out.build_header(ideal_sdf.header_lines[0].strip(), 
                       ["  ", "PDB mol with bond order from Ligand Expo " + \
                        "ideal sdf and implicit hydrogens"])
  (sdf_fd, sdf_out_fname) = tempfile.mkstemp(suffix=".sdf")
  sdf_file = os.fdopen(sdf_fd, "w+")

  print >> sdf_file, mol_out
  sdf_file.close()
 
  babel_cmd = "%s -h -i sdf %s -o mol2 %s" \
    % (babel, sdf_out_fname, mol2_fname)
  #mol_charge_cmd = "%s -in %s -out %s " \
  mol_charge_cmd = "%s -am1bcc -in %s -out %s" \
    % (mol_charge, sdf_out_fname, mol2_fname)
  os.system(mol_charge_cmd)
  os.unlink(sdf_out_fname)

  # Need to "rename" the mol2 atoms heavy atoms to have teh same name as in
  # the pdb file
  mol2_mol = utils.mol2.molecule(mol2_fname)
  bad_order = False
  for a,b in zip(hetgrp.atoms, mol2_mol.atoms):
    if(b.name == "H"): continue

    # braindead check for different ordering of atoms
    for x,y in zip(a.position, b.position):
      if(x-y < -0.01 or x-y > 0.01):
        print "atoms are not in the same order or have differing positions"
        os.unlink(mol2_fname)
        bad_order = True
        break;

    b.name_str = a.name.strip()
 
  if(not bad_order):  
    mol2_file = open(mol2_fname, "w+") 
    print >> mol2_file, mol2_mol
    return True
  else: 
    return False


if(__name__ == "__main__"):

  for f in os.listdir("."):
    if(not f.endswith("_l.pdb")): continue
  
    print "now processing", f[:-4]
    # Load the pdb hetgroup as bound in the xtal -- with hydrogen atoms removed
    pdb_lig = utils.pdb.hetgroup(f)
    run(pdb_lig, mol2_fname = f[:-3] + "mol2")
    print 
