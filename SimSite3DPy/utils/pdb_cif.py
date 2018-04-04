from sys import stderr
from datetime import datetime
import sys,os
from numpy import *
import sdf, pdb

################################################################################

################################################################################
class molecule:

  def __init__(self, fname):
    """
    Assumption: the first line of the CIF file contains either data_XXXX or
    data_YYY where XXXX is a valid 4 char pdb structure code or YYY is a valid
    3 char pdb HET code (think Ligand Expo)
    """
    self.fail = False

    try:
      infile = open(fname, "r")
    except IOError, (errno, strerror):
      print >> stderr, "Unable to open the file", fname
      print >> stderr, "error(%s): %s" % (errno, strerror)
      print >> stderr, "	in (%s)" % (self.__module__)
      self.fail = True
      return

    # Check the file type -- read it from the first line in the CIF file
    line = infile.next()
    toks = line.rstrip(" \n").split("_")
    if(not toks[0] == "data"): 
      print >> stderr, """Unknown file type:
Expected data_XXX or data_YYYY where XXX is a valid Ligand Expo identifier or 
YYYY is a valid PDB structure identifier.  Received (%s)
""" % (line.rstrip("\n"))
      self.fail = True
      return


    if(len(toks[1]) <= 3):
      print "Assuming this is a Ligand Expo CIF file"
      self.chemical_component = chem_comp(infile, toks[1])
      if(self.chemical_component.fail):
        self.fail = True
        return
    elif(len(toks[1]) == 4):
      print "Assuming this is an RCSB PDB CIF file"
      print "Reading of mmCIF structure files is unimplemented\n"
      self.fail = True
      return
    else:
      print >> stderr, """Unknown file type:
Expected data_XXX or data_YYYY where XXX is a valid Ligand Expo identifier or 
YYYY is a valid PDB structure identifier.  Received (%s)
""" % (line.rstrip("\n"))
      self.fail = True
      return
    
################################################################################

################################################################################
class atom:

  def __init__(self, name, pos, element):
    self.name = name
    self.pos = array(pos)
    self.element = element

class bond:
  def __init__(self, first_atom, second_atom, order):
    self.first_atom = first_atom
    self.second_atom = second_atom
    self.order = order
 
  def __repr__(self):
    return ",".join([self.first_atom, self.second_atom, self.order])

################################################################################

################################################################################
class chem_comp:
  ##############################################################################

  ##############################################################################
  def __init__(self, cif_file, mol_name=""):
    self.fail = False
    line = cif_file.next()
    if(not line.startswith("#")):
      print >> stderr, "Expected a '#' line\n"
      self.fail = True
      return

    self.mol_name = mol_name
    self.comp = {}
    self.comp_atoms = (0,)
    self.comp_bonds = (0,)

    for line in cif_file:
      # Read the _chem_comp section
      if(line.startswith("_chem_comp.")):
        (line_type, field, string) = self.__data_line(line)
        self.comp[field] = string
        self.read_chem_comp_section(cif_file, field)
      elif(line.startswith("loop_")):
        (section_id, fields, data) = self.__loop_section(cif_file)
        if(section_id == "_chem_comp_atom"): self.comp_atoms = (fields, data)
        elif(section_id == "_chem_comp_bond"): self.comp_bonds = (fields, data)
      if(self.fail): return

    if(self.comp["pdbx_ideal_coordinates_missing_flag"] == "N"):
      self.init_data_structs()
    elif(self.comp["pdbx_ideal_coordinates_missing_flag"] == "Y"
       and self.comp["pdbx_model_coordinates_missing_flag"] == "N"):
      self.init_data_structs(use_ideal_coords=False)
      print "Ideal coordinates are missing -- using model coordinates"
    else:
      print "In", mol_name
      print "Ideal AND model coordinates are missing -- EPIC FAIL"
      self.fail = True
      return
    self.fail = False
  ##############################################################################

  ##############################################################################
  def read_chem_comp_section(self, cif_file, prev_field):
    for line in cif_file:
      if(line.startswith("#")): return

      (line_type, field, string) = self.__data_line(line)
      if(not line_type == "" and not line_type == "_chem_comp"):
        print >> stderr, "Error parsing file; expected _chem_comp item"
        self.fail = True
        return

      if(field == ""):
        self.comp[prev_field] += string
      else:
        self.comp[field] = string
        prev_field = field
  ##############################################################################

  ##############################################################################
  def init_data_structs(self, use_ideal_coords=True):
    """
    Use the strings in self.comp_atoms to get the atom names and coordinates
    Use the strings in self.comp_bonds to get the bond orders
      
    It may be that the ideal coordinates are missing, in that case we 
    need to try the model coordinates
   
    Each loop section has 1 field heading per line with possible contination
    to the next line.  After that there are space delimited rows till we hit
    the end.
    """

    # Setup the atom information
    atom_fields = self.comp_atoms[0]
    pos_idx = -1
    atom_id_idx = -1 
    element_idx = -1
    for i, field in zip(range(len(atom_fields)), atom_fields):
      if(field == "atom_id"): atom_id_idx = i
      if(use_ideal_coords and field == "pdbx_model_Cartn_x_ideal"): pos_idx = i
      if(not use_ideal_coords and field == "model_Cartn_x"): pos_idx = i
      if(field == "type_symbol"): element_idx = i

    atom_data = self.comp_atoms[1]
    self.atoms = []
    for a in atom_data:
      #print " ".join(a)
      pos = array([ float(a[pos_idx + i]) for i in range(3) ])
      self.atoms.append(atom(a[atom_id_idx].strip("\""), pos, a[element_idx]))

    # Setup the bond information
    bond_fields = self.comp_bonds[0]
    first_atom_idx = -1
    second_atom_idx = -1
    bond_order_idx = -1
    for i, field in zip(range(len(bond_fields)), bond_fields):
      if(field == "atom_id_1"): first_atom_idx = i
      if(field == "atom_id_2"): second_atom_idx = i
      if(field == "value_order"): bond_order_idx = i

    bond_data = self.comp_bonds[1]
    self.bonds = [ bond(b[first_atom_idx].strip("\""), 
                        b[second_atom_idx].strip("\""), 
                        b[bond_order_idx]) for b in bond_data ]

  ##############################################################################

  ##############################################################################
  def __data_line(self, line):

    line = line.rstrip(" \n")

    # continuation line
    if(line.startswith(";")): return ("", "", line.lstrip(";"))
     
    dot_pos = line.find(".")
    line_type = ""
    if(dot_pos > -1):
      line_type = line[:dot_pos]
      line = line[dot_pos+1:]
    toks = line.split()
    if(len(toks) == 1): return (line_type, toks[0], "")
    return (line_type, toks[0], toks[1])
  ##############################################################################

  ##############################################################################
  def __loop_section(self, cif_file):
    fields = []
    data = []
    prev_section_id = "" 
    
    for line in cif_file:
      line = line.rstrip(" \n")
      if(line.startswith("#")): return (prev_section_id, fields, data)

      if(line.startswith("_")):
        dot_pos = line.find(".")
        section_id = line[:dot_pos]
        if(prev_section_id == ""): prev_section_id = section_id
        elif(not prev_section_id == section_id):
          print >> stderr, "Error parsing file; expected a %s line, but received a %s line" % (prev_section_id, section_id)
          return(prev_section_id, fields, data)

        fields.append(line[dot_pos+1:])

      # Assume this is a data line
      else:
        if(line.startswith(";")):
          toks = line[1:].split()
          if(len(data) and len(toks)): data[-1].extend(toks)
          elif(len(toks)): data.append(toks)
        else: data.append(line.split())
  ##############################################################################

  ##############################################################################
  def pdb_hetgrp_to_sdf(self, hetgrp):
    """
    Convert the chem_comp atom and bond data from cif format to sdf and use
    the positions of the hetgrp
 
    Note: we use SDF since we want to have explicit bond order, but use 
    molcharge assign the orbitals and partial charges
    """

    N = len(self.atoms)
    atoms_idz = dict([ (a.name, i) for i,a in zip(range(N), self.atoms)])

    # Map idz from full ligand to ligand fragment and set sdf atoms
    my_sdf = sdf.molecule()
    N = len(hetgrp.atoms)
    kept_idz = {}
    kept_names = []
    for i, a in zip(range(N), hetgrp.atoms):
      atom_name = a.name.strip()

      # simple check to avert alt loc problems -- we generally assume
      # the first location is the best (could be an incorrect assumption,
      # but we want it automated).
      if(atom_name in kept_names): continue
      
      # annoying check, but we want to get things running 
      if(not atom_name in atoms_idz): 
        print >> stderr, "The SDF file does not have the atom:",  atom_name
        print >> stderr, "Cannot process the file\n"
        return ([], my_sdf)

#      print i, atom_name
      kept_idz[atoms_idz[atom_name]] = i
      kept_names.append(atom_name)
      
      sdf_atom = sdf.atom()
      sdf_atom.other = ["0", "0", "0", "0", "0"]
      # Take the element from the corresponding cif atom
      sdf_atom.element = self.atoms[atoms_idz[atom_name]].element 
      # Take the atom names and coordinates from the given hetgrp
      sdf_atom.position = a.position
      my_sdf.atoms.append(sdf_atom)

#    print kept_names
    # Check if the het group in the PDB file is missing any heavy atoms
    for cif_atom in self.atoms:
      if(not cif_atom.name.startswith("H") and 
         not cif_atom.name.startswith("D")):
        if(not cif_atom.name in kept_names):
          print >> stderr, "The given hetgrp does not have the atom:", 
          print >> stderr, cif_atom.name, "Cannot process the file\n"
          return ([], my_sdf)
 
    # Setup the bonds for the sdf molecule
    for bond in self.bonds:
      atom1 = atoms_idz[bond.first_atom]
      atom2 = atoms_idz[bond.second_atom]
      if(not atom1 in kept_idz or not atom2 in kept_idz): continue

      sdf_bond = sdf.bond()
      if(bond.order == "SING"): sdf_bond.type = 1
      elif(bond.order == "DOUB"): sdf_bond.type = 2
      elif(bond.order == "TRIP"): sdf_bond.type = 3
      else: 
        print >> stderr, "Unknown bond order: %s" % (bond.order)
        sdf_bond.type = -1
      sdf_bond.other = ["0", "0", "0"]
      # Adjust atom indices as we use 0 indexed in the code, but SDF 
      # uses 1 indexed
      sdf_bond.atom1 = kept_idz[atom1] + 1
      sdf_bond.atom2 = kept_idz[atom2] + 1
      my_sdf.bonds.append(sdf_bond)

    now = datetime.now()    
    timestamp = "Created by pdb_model_to_sdf at %s on %s"  % \
      (str(now).split(".")[0].split()[1], now.strftime("%A, %B %d, %Y"))

    my_sdf.build_header(self.mol_name, [timestamp, ""])
    return (kept_names, my_sdf)

  ##############################################################################

  ##############################################################################
  def to_sdf(self):
    """
    Convert the chem_comp atom and bond data from cif format to sdf using
    the coordinates in the cif file
 
    Note: the current purpose of this method is PDB cif -> sdf for FragGen
    """

    my_sdf = sdf.molecule()
    N = len(self.atoms)
    atoms_idz = dict([ (a.name, i) for i,a in zip(range(1, N+1), self.atoms)])

    for cif_atom in self.atoms:
      sdf_atom = sdf.atom()
      sdf_atom.other = ["0", "0", "0", "0", "0"]
      sdf_atom.element = cif_atom.element
      sdf_atom.position = cif_atom.pos[:]
      my_sdf.atoms.append(sdf_atom)
      
    for cif_bond in self.bonds:
      sdf_bond = sdf.bond()
      if(cif_bond.order == "SING"): sdf_bond.type = 1
      elif(cif_bond.order == "DOUB"): sdf_bond.type = 2
      elif(cif_bond.order == "TRIP"): sdf_bond.type = 3
      else: 
        print >> stderr, "Unknown bond order: %s" % (cif_bond.order)
        sdf_bond.type = -1
      sdf_bond.other = ["0", "0", "0"]
      sdf_bond.atom1 = atoms_idz[cif_bond.first_atom]
      sdf_bond.atom2 = atoms_idz[cif_bond.second_atom]
      my_sdf.bonds.append(sdf_bond)

    now = datetime.now()    
    timestamp = "Created by pdb_cif.to_sdf at %s on %s"  % \
      (str(now).split(".")[0].split()[1], now.strftime("%A, %B %d, %Y"))
    my_sdf.build_header(self.mol_name, [timestamp, ""])

    return my_sdf
################################################################################

################################################################################

if(__name__ == "__main__"):
  # The openeye modules (wrappers to C++ libraries) are stored somewhere on
  # the PSA system.  When writing this script they resided at
  # /soft/linux64/openeye/wrappers/v1.7.2.4/python.  Pleae append their 
  # current location to your $PYTHONPATH variable
  from openeye.oechem import *
  from openeye.oequacpac import *
  import tempfile
 
  my_prot = pdb.residues("/psa/pdb/pdb1qpg.ent")
  MAP = None
  for hetgrp in my_prot.hetgroups:
    if(hetgrp.hetID == "MAP"): MAP = hetgrp

  my_mol = molecule("MAP.cif")
  if(my_mol.fail): sys.exit(-1)
  my_chem_comp = my_mol.chemical_component
 
  # This is the ideal PDBX to sdf using only the heavy atoms in the hetgrp
  (pdb_atom_names, my_sdf) = my_chem_comp.pdb_hetgrp_to_sdf(MAP)

  # Read the SDF string into OE data structure
  oe_ims = oemolistream()
  oe_ims.SetFormat(OEFormat_MDL)
  oe_ims.openstring("%s\n" % (my_sdf))
  oe_mol_graph = OEGraphMol()
  OEReadMolecule(oe_ims, oe_mol_graph) 

  # Assign partial charges using the molcharge library
# For now we don't need the accuracy of AM1BCC and it is very slow to calculate
#    charge_type = OECharges_AM1BCC
  charge_type = OECharges_MMFF94
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
  ofs = oemolostream("test.mol2")
  OEWriteMol2File(ofs, oe_mol_graph)
