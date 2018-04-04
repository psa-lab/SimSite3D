
from sys import stderr
from numpy import sqrt, dot
import re
import os

# For now this makes things much more complicated than they need be
#class hetnames:
#
#  def __init__(self):
#    self.C = re.compile(" C")
#    self.N = re.compile(" N")
#    self.P = re.compile(" P")
#    self.O = re.compile(" O")
#    self.S = re.compile(" S")

class metal_table_entry:

  def __init__(self, pdb_metal, metal_name, act_type, pt_rad, min_act_dist,
               max_act_dist, vdw_radius):
    self.pdb_metal_name = pdb_metal
    self.element = metal_name
    self.act_type = act_type
    self.pt_rad = pt_rad
    self.min_act_rad = min_act_dist
    self.max_act_rad = max_act_dist
    self.radius = vdw_radius

class metal_lookup:

  def __init__(self):
    self.table = {}
    # must be a better way to have a data file read at run time
    fname = os.environ["ASCBASE_INSTALL_DIR"]
    fname += "/ASCbasePy/utils/pdb_metals.csv"
    infile = open(fname, "r")

    for line in infile:
      if(line.startswith("#")): continue
      self.__add_line(line)

  def __add_line(self, line):
    toks = [ t.strip() for t in line.split(",") ]
    pdb_metal_name = toks[0].strip("\"")
    self.table[pdb_metal_name] = \
      metal_table_entry(pdb_metal_name, toks[1], toks[2], float(toks[3]),
                        float(toks[4]), float(toks[5]), -1.0) 

  def add_type_and_orbit_to_metal(self, atom):
    if(not atom.name in self.table):
      print >> stderr, "Warning: an unsupported metal atom (%s)" % (atom.name),
      print >> stderr, "was found"
      print >> stderr, "Skipping the following atom line:\n%s\n" % (atom)
      return False

    tmp = self.table[atom.name]
    atom.radius = tmp.radius
    #atom.orbital = tmp.orbital
    atom.interact_type = tmp.act_type
    self.min_act_rad = tmp.min_act_rad
    self.max_act_rad = tmp.max_act_rad
    return True


class residue_table_entry:

  def __init__(self, info_tuple):
    self.atom, self.residue, self.radius, self.orbital, self.letter, \
      self.interact_type, self.carbon_nbr = info_tuple

class residue_lookup:

  def __init__(self):
    self.table = {}
    self.table[" O  "] = {}

    # Waters are handled as a special case -- other hetatoms are not 
    # handled at all at this point in time (nor are metals)
    self.table[" O  "]["HOH"] = \
      residue_table_entry(("O", "HOH", 1.36, "SP3", "ALPHA", "DONEPTOR", None))
    self.table[" O  "]["OH2"] = \
      residue_table_entry(("O", "HOH", 1.36, "SP3", "ALPHA", "DONEPTOR", None))
  
    # must be a better way to have a data file read at run time
    infile = \
      open(os.environ["ASCBASE_INSTALL_DIR"] + "/ASCbasePy/utils/pdb_residues.csv", "r")
#"/home/vanvoor4/code/python/trunk/coord_files/pdb_residues.csv", "r") 

    main_chain = {}
    # MAIN_CHAIN stuff
    for i in range(5):
      line = infile.next()
      self.__add_line(line, main_chain)
    self.__add_main_chain_atoms(main_chain, "GLY")
      
    prev_res = "GLY"
    for line in infile:
      curr_res = self.__add_line(line, self.table)
      if(curr_res != prev_res):
        prev_res = curr_res
        self.__add_main_chain_atoms(main_chain, curr_res)

  def __add_main_chain_atoms(self, main_chain, res):
    for a in main_chain:
      if(not a in self.table): self.table[a] = {}
      if(res == "PRO" and a == " N  "):
        tmp = main_chain[a]["MAIN_CHAIN"]
        self.table[" N  "]["PRO"] = \
          residue_table_entry(("N", "PRO", tmp.radius, tmp.orbital, \
                               tmp.letter, "NOTHING", tmp.carbon_nbr ))
      else: self.table[a][res] = main_chain[a]["MAIN_CHAIN"]
  
  def __add_line(self, line, the_dict): 
    toks = [ t.strip() for t in line.split(",") ]
    a_name = toks[0].strip("\"")
    r_name = toks[1].strip("\"")
    
    if(not a_name in the_dict):
      the_dict[a_name] = {}
    if(not len(toks[8])):
      toks[8] = None
    the_dict[a_name][r_name] = \
      residue_table_entry((toks[2], toks[3], float(toks[4]), toks[5], 
                           toks[6], toks[7], toks[8]))
    return r_name

  def add_type_and_orbit_to_atom(self, atom):
    if(not atom.name in self.table or \
       not atom.resName in self.table[atom.name]):
      print >> stderr, "Warning: a nonstandard residue (%s)" % (atom.resName),
      print >> stderr, "and atom (%s) pair was found" % (atom.name)
      print >> stderr, "Skipping the following atom line:\n%s\n" % (atom)
      return False

    tmp = self.table[atom.name][atom.resName]
    atom.radius = tmp.radius
    atom.orbital = tmp.orbital
    atom.interact_type = tmp.interact_type
    atom.carbon_nbr = tmp.carbon_nbr 
    return True

class atom:
  
  def __init__(self, line):
    self.record_name = line[0:6]
    self.serial = int(line[6:11])
    self.name = line[12:16]
    self.altLoc = line[16]
    self.resName = line[17:20]
    self.chainID = line[21]
    self.resSeq = int(line[22:26])
    self.iCode = line[26] 
    self.position = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
    self.occupancy = float(line[54:60])
    self.tempFactor = float(line[60:66])
  
  def __cmp__(self, other):
    if(self.serial < other.serial):
      return -1
    elif(self.serial > other.serial):
      return 1
    else:
      return 0

  def __repr__(self):
    return "%s%5d %s%s%s %s%4d%s   %8.3f%8.3f%8.3f%6.2f%6.2f" % \
      (self.record_name, self.serial, self.name, self.altLoc, self.resName,
       self.chainID, self.resSeq, self.iCode, self.position[0], 
       self.position[1], self.position[2], self.occupancy, self.tempFactor)

  def is_altloc(self, other):
    if(not self.name == other.name): return False
    if(not self.resSeq == other.resSeq): return False
    if(not self.chainID == other.chainID): return False
    if(not self.iCode == other.iCode): return False

    # Assume we have the same atom, check altLoc
    if(not self.altLoc == other.altLoc): return True

    # Return False -- these 2 atoms are the same, and are not alternate 
    # locations
    return False

  def is_pdb_hydrogen(self):
    c13 = self.name[0]
    c14 = self.name[1].upper()
    if((c13 == ' ' or c13.isdigit()) and (c14 == 'H' or c14 == 'D')):
      return True
    # Some hydrogens are H followed by a character & 2 integers (H[A-Z][0-9]{2})
    # Others are H followed by 3 integers (H[0-9]{3})
    if(re.match("H[A-Z][0-9]{2}", self.name) or \
       re.match("H[0-9]{3}", self.name) or \
      # this could be dangerous, but the hydrogen atom naming seems to
      # be upto whoever deposited the HET entries in Ligand Expo
       (self.name[0] == 'H' and not self.name[1] == " " and \
        not self.name[2] == " " and not self.name[3] == " ")): return True
    return False
################################################################################

################################################################################
class seqres_entries:
 
  def __init__(self):
    self.seqres = {}

  def add_line(self, line):
    chainID = line[11]
    num_res = int(line[13:17])
    
    res_names = line.rstrip(" \n")[19:].split()
    if(not chainID in self.seqres): self.seqres[chainID] = []
    self.seqres[chainID].extend(res_names)
    if(len(self.seqres) > num_res):
      print "somehow we have more residues names than is specified by SEQRES"
     
################################################################################

################################################################################
class residue:

  def __init__(self):
    self.atoms = []
    self.alt_atoms_locs = {}
    self.iCode = ""
    self.resSeq = -1
    self.resName = ""
    self.chainID = ""
    
  def __len__ (self):
    return len(self.atoms) 

  def __iter__(self):
    for atom in self.atoms:
      yield atom
    return

  def __repr__(self):
    str = ""
    for atom in self.atoms:
      if(len(str)): str += "\n"
      str += "%s" % (atom)
    return str

  def append(self, atom):
    if(len(self.atoms) == 0):
      self.iCode = atom.iCode
      self.resSeq = atom.resSeq
      self.resName = atom.resName
      self.chainID = atom.chainID
    elif(self.iCode != atom.iCode or self.resSeq != atom.resSeq or \
         self.chainID != atom.chainID):
      print >> stderr, "Attempted append an atom to the incorrect residue"

    # take care of the alternate location, by choosing the first location
    # We probably should be using both locations ...
    if(not atom.altLoc == " "):
      for a in self.atoms:
        if(a.is_altloc(atom)):
          if(not atom.altLoc in self.alt_atoms_locs): 
            self.alt_atoms_locs[atom.altLoc] = []
          self.alt_atoms_locs[atom.altLoc].append(atom) 
          return 

    self.atoms.append(atom)

  # Bad name since this is a check to see if atom should be appended to the
  # residue rather than is the atom in the residue class
  def __contains(self, atom):
    if(self.resSeq == atom.resSeq and self.iCode == atom.iCode and \
       self.chainID == atom.chainID):
      return True
    else:
      return False

  def is_water(self):
    if(self.resName == "HOH"):
      return True
    return False

  def __cmp__(self, other):
    rv = self.resSeq - other.resSeq
    if(rv):
      return rv
 
    if(self.iCode != other.iCode):
      if(self.iCode < other.iCode):
        return -1
      else:
        return 1
    return 0

class residues:

  def __init__(self, fname):
    pdb_residue_tbl = residue_lookup()
    pdb_metals = metal_lookup()

    self.residues = []
    self.metals = []
    self.waters = [] 
    self.hetgroups = hetgroups()
    self.fail = False
    self.seqres = seqres_entries()
    
    try:
      pdb_file = open(fname, "r")
    except IOError, (errno, strerror):
      print >> stderr, "Unable to open the file", fname
      print >> stderr, "error(%s): %s" % (errno, strerror)
      print >> stderr, "in %s" % (self.__module__)
      self.fail = True
      return

    first = True
    current_residue = residue()
    for line in pdb_file:
      if(line.startswith("HET   ")): self.hetgroups.add_HET_record(line)
      elif(line.startswith("SEQRES")): self.seqres.add_line(line)
      elif(line.startswith("ATOM  ") or line.startswith("HETATM")): break

    # At the moment it seems like the safest thing to do is reset the file
    # to the beginning since ftell & fseek work with the current location in
    # the file.  However, we don't know how far to back up since we don't know
    # how much of the file is buffered.
    pdb_file.seek(0, 0)

    for line in pdb_file:
      if(line.startswith("MODEL ")):
        model_num = int(line[6:].strip(" \t\n"))
        if(model_num == 1):  
          print "Warning!\n\tFound multiple models, using the first model only" 
        else: break

      if(line.startswith("HETATM")):
        tmp_atom = atom(line)
        # Ignores hydrogens
        if(tmp_atom.is_pdb_hydrogen()): continue

        if(tmp_atom.resName == "HOH"): self.waters.append(tmp_atom)
        # We want to distinguish between metals "by themselves" and metals
        # that are covalently a part of ligands
        elif(tmp_atom.name in pdb_metals.table and \
             tmp_atom.resName.strip() == tmp_atom.name.strip()):
          pdb_metals.add_type_and_orbit_to_metal(tmp_atom)
          self.metals.append(tmp_atom)
        else:
          self.hetgroups.add_HETATOM(line)
#          print "Warning: an unsupported hetatom (%s) was found" % \
#            (tmp_atom.name)
#          print "Skipping the following atom line:\n%s\n" % (tmp_atom)
        continue

      elif(not line.startswith("ATOM  ")): # and not line.startswith("HETATM")):
        continue

      tmp_atom = atom(line)
      # Ignores hydrogens
      if(tmp_atom.is_pdb_hydrogen()):
        continue

      # the residue table only has entries for the 20 standard residues
      if(not pdb_residue_tbl.add_type_and_orbit_to_atom(tmp_atom)):
        print >> stderr, "Warning: an unsupported ATOM (%s) was found" % \
          (tmp_atom.name)
        print >> stderr, "Skipping the following atom line:\n%s\n" % (tmp_atom)
        continue

      if(first):
        current_residue.append(tmp_atom)
        first = False
      else:
        if(not current_residue._residue__contains(tmp_atom)):
          if(current_residue.is_water()):
            self.waters.append(current_residue)
          else: 
            self.residues.append(current_residue)
          current_residue = residue()
        current_residue.append(tmp_atom)
      
    # Add the last residue
    if(current_residue.is_water()):
      self.waters.append(current_residue)
    else: 
      self.residues.append(current_residue)
################################################################################
    
################################################################################
  def write_atoms_in_canonical_order(self, ofile):
    """
    One of the many assuptions of the IK code is that the atoms with in
    a residue are in canonical form.  This makes a number of items much
    easier to work with.  Unfortunately, the MD files from Prof. Bob Cukier
    have N -> CA -> side chain -> C -> O

    At the present we use implicit hydrogen atoms -- ignore all hydrogen atoms
    just for kicks since I don't want to handle the "old" hydrogen
    naming convention, the "new" PDB convention and Prof. Bob Cukier's 
    hydrogen naming convention
    """
    
    main_chain_atoms = [" N  ", " CA ", " C  ", " O  "]
    my_atoms = { "ALA": main_chain_atoms[:], "ARG": main_chain_atoms[:],
                 "ASN": main_chain_atoms[:], "ASP": main_chain_atoms[:],
                 "CYS": main_chain_atoms[:], "GLN": main_chain_atoms[:],
                 "GLU": main_chain_atoms[:], "GLY": main_chain_atoms[:],
                 "HIS": main_chain_atoms[:], "ILE": main_chain_atoms[:],
                 "LEU": main_chain_atoms[:], "LYS": main_chain_atoms[:],
                 "MET": main_chain_atoms[:], "PHE": main_chain_atoms[:],
                 "PRO": main_chain_atoms[:], "SER": main_chain_atoms[:],
                 "THR": main_chain_atoms[:], "TRP": main_chain_atoms[:],
                 "TYR": main_chain_atoms[:], "VAL": main_chain_atoms[:] }
    my_atoms["ALA"].append(" CB ")
    my_atoms["ARG"].extend([" CB ", " CG ", " CD ", " NE ", " CZ ", " NH1", " NH2"])
    my_atoms["ASN"].extend([" CB ", " CG ", " OD1", " ND2"])
    my_atoms["ASP"].extend([" CB ", " CG ", " OD1", " OD2"])
    my_atoms["CYS"].extend([" CB ", " SG "])
    my_atoms["GLN"].extend([" CB ", " CG ", " CD ", " OE1", " NE2"])
    my_atoms["GLU"].extend([" CB ", " CG ", " CD ", " OE1", " OE2"])
#    my_atoms["GLY"]
    my_atoms["HIS"].extend([" CB ", " CG ", " ND1", " CD2", " CE1", " NE2"])
    my_atoms["ILE"].extend([" CB ", " CG1", " CG2", " CD1"])
    my_atoms["LEU"].extend([" CB ", " CG ", " CD1", " CD2"])
    my_atoms["LYS"].extend([" CB ", " CG ", " CD ", " CE ", " NZ "])
    my_atoms["MET"].extend([" CB ", " CG ", " SD ", " CE "])
    my_atoms["PHE"].extend([" CB ", " CG ", " CD1", " CD2", " CE1", " CE2", " CZ "])
    my_atoms["PRO"].extend([" CB ", " CG ", " CD "])
    my_atoms["SER"].extend([" CB ", " OG "])
    my_atoms["THR"].extend([" CB ", " OG1", " CG2"])
    my_atoms["TRP"].extend([" CB ", " CG ", " CD1", " CD2", " NE1", " CE2", " CE3", " CZ2", " CZ3", " CH2"])
    my_atoms["TYR"].extend([" CB ", " CG ", " CD1", " CD2", " CE1", " CE2", " CZ ", " OH "])
    my_atoms["VAL"].extend([" CB ", " CG1", " CG2"])

    for res in self.residues:
      if(res.resName in my_atoms):
        for atom_name in my_atoms[res.resName]:
	  for a in res:
	    if(a.name == atom_name): print >> ofile, a
      else:
        for a in res:
	  print >> ofile, a

    for m in self.metals: print >> ofile, m
    for H2O in self.waters: print >> ofile, H2O

#    self.residues = []
#    self.metals = []
#    self.waters = [] 


################################################################################
    
################################################################################
  def __len__ (self):
    return len(self.residues) 
################################################################################
    
################################################################################
  
  # Use a generator since I am not familiar enough with using __iter__
  # and next to produce "different iterators"
  def __iter__(self):
    for r in self.residues:
      yield r
    return

################################################################################
    
################################################################################
  def __getitem__(self, index):
    return self.residues[index]

################################################################################
    
################################################################################
class HET_record:
  """
  Data class
  Holds information from a PDB HET record (line)
  """

  def __init__(self, line):
     """
     Given a "HET   " record, split the line into its fields
     """
     if(len(line) < 31):
       print >> stderr, "Bad HET line\n\t%s" % (line)
       return

     self.hetID = line[7:10]
     self.chainID = line[12]
     self.seqNum = int(line[13:17])
     self.iCode = line[17]
     self.numHetAtoms = int(line[20:25])
     self.text = line[30:].strip("\n")
     self.atoms = []

  def het_line(self):
    my_str = "HET    %3s  %s%4d%s   %5d     %s" % \
      (self.hetID, self.chainID, self.seqNum, self.iCode, self.numHetAtoms,
       self.text.rstrip())
    return my_str

  def __eq__(self, other):
    if(not self.hetID == other.hetID): return False
    if(not self.chainID == other.chainID): return False
    if(not self.seqNum == other.seqNum): return False
    if(not self.iCode == other.iCode): return False
    if(not self.numHetAtoms == other.numHetAtoms): return False
    return True

  def __repr__(self):
    my_str = "HET    %3s  %s%4d%s   %5d     %s" % \
      (self.hetID, self.chainID, self.seqNum, self.iCode, self.numHetAtoms,
       self.text.rstrip())
    for a in self.atoms: my_str += "\n%s" % (a.__repr__())
    return my_str

class hetgroups:
  """
  Simple class to handle the matching between HET entries and the
  corresponding HETATM records
  """

  def __init__(self):
    self.HET_records = {}
    
  def add_HET_record(self, line):
    """
    Given a "HET   " record, store the information in an dictionary with 
    columns 8-18 as the index
    """ 
    key = line[7:18]
    if(key in self.HET_records):
      print >> stderr, "HET records are not unique;",
      print >> stderr, "some HET groups could be missed"
      return 

    self.HET_records[key] = HET_record(line)

  # At the present we are more concerned with getting the prototype running
  # and will not be storing other HET information such as name, and synonyms

  def add_HETATOM(self, line):

    hetatm = atom(line)
    idx = "%s  %s%4d%s" % (hetatm.resName, hetatm.chainID, hetatm.resSeq, 
                           hetatm.iCode)
    if(not idx in self.HET_records):
#      AAz = [ "ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", 
#              "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", 
#              "TYR", "VAL"]
#      if(hetatm.resName in AAz):
#        print >> stderr, """
#For some reason we got a standard residue in pdb.add_HETATM()
#Skipping . . .
#"""
#        return False

      print >> stderr, """
Could not find a corresponding HET record for the following HETATM line:
\t%s
Making up our own HET entry
Maybe it corresponds to a MODRES record or it is one of the standard
amino or nucleic acids?""" % (line.strip("\n"))

      self.add_HET_record("HET    %3s  %s%4d%s      1     Added by ASCbasePy" % \
                          (hetatm.resName, hetatm.chainID, hetatm.resSeq, 
                           hetatm.iCode))
    self.HET_records[idx].atoms.append(hetatm) 
    return True

  def __iter__(self):
    for idx, hetgrp in self.HET_records.iteritems():
      yield hetgrp
    return

class hetgroup:
  """
  Simple class to load all heavy atoms from a PDB file
  """

  def __init__(self, fname=""):
    if(len(fname)): self.__load_from_file(fname)


  def __load_from_file(self, fname):
    self.num_atoms = 0
    self.atoms = []

    try:
      pdb_file = open(fname, "r")
    except IOError, (errno, strerror):
      print >> stderr, "Unable to open the file", fname
      print >> stderr, "error(%s): %s" % (errno, strerror)
      print >> stderr, "in %s" % (self.__module__)
      return

    for line in pdb_file:
      if(not line.startswith("HETATM") and not line.startswith("ATOM  ")):
        continue

      tmp_atom = atom(line)
      # Ignore hydrogens
      if(tmp_atom.is_pdb_hydrogen()): continue

      self.atoms.append(tmp_atom)
    self.num_atoms = len(self.atoms)

class my_pair:
  
  def __init__(self, a=None, b=None):
    self.a = a
    self.b = b

# returns (RMSD, A \ A&B, B \ A&B)
def atom_serial_rmsd(A, B):
  AorB = {}

  for A_res in A:
    for atom in A_res:
      AorB[atom.serial] = my_pair(a=atom)

  for B_res in B:
    for atom in B_res:
      if(atom.serial in AorB):
        AorB[atom.serial].b = atom
      else:
        AorB[atom.serial] = my_pair(b=atom)

  rmsd = 0.0
  A_only = []
  B_only = []
  N = 0
  for serial, atom_pair in AorB.iteritems():
    a, b = atom_pair.a, atom_pair.b
    if(not a is None and not b is None):
      rmsd += sum([ (x-y)*(x-y) for x,y in zip(a.position, b.position) ])
      N += 1
    elif(not atom_pair.a is None):
      A_only.append(a)
    elif(not atom_pair.b is None):
      B_only.append(b)

  rmsd = sqrt(rmsd/float(N))

  return (rmsd, A_only, B_only)

# Very simple C_alpha RMSD -- feed residues for which you wish to compute
# C_alpha rmsd and omit the rest
def Calpha_rmsd(A,B):

  rmsd = 0.0
  N = 0
  for (A_res, B_res) in zip(A,B):
    A_CA = None
    B_CA = None

    for a in A_res:
      if(a.name == " CA "):
        A_CA = a
        break

    for b in B_res:
      if(b.name == " CA "):
        B_CA = b
        break

    if(not A_CA is None and not B_CA is None):
      rmsd += sum([ (x-y)*(x-y) for x,y in zip(A_CA.position, B_CA.position) ])
      N += 1

  return sqrt(rmsd/float(N))

# Very simple main chain RMSD -- feed residues for which you wish to compute
# main chain rmsd and omit the rest
def main_chain_rmsd(A,B):
  A_atoms = { " N  ":[], " CA ":[], " C  ":[], " O  ":[] }
  B_atoms = { " N  ":[], " CA ":[], " C  ":[], " O  ":[] }
  rmsd = 0.0
  N = 0

  for (A_res,B_res) in zip(A,B):

    for a in A_res:
      if(a.name in A_atoms):
        A_atoms[a.name] = a.position

    for b in B_res:
      if(b.name in B_atoms):
        B_atoms[b.name] = b.position

    for name in A_atoms:
      rmsd += sum([ (x-y)*(x-y) for x,y in zip(A_atoms[name], B_atoms[name]) ])
      N += 1

  return sqrt(rmsd/float(N))

# Very simple heavy atom RMSD -- feed residues for which you wish to compute
# the heavy atom rmsd and omit the rest
def hvy_atom_rmsd(A,B):
  rmsd = 0.0
  N = 0

  #self.resName

  for (A_res,B_res) in zip(A,B):
    if(not A_res.resName == B_res.resName):
      print "Cannot compute heavy atom rmsd if the sequences are different\n"
      return 1E+37

    A_pos = dict([ (a.name, a.position) for a in A_res])
    B_pos = dict([ (b.name, b.position) for b in B_res])

    positions = []
    for a in A_res:
      if(a.is_pdb_hydrogen()):
        continue

#      positions.append((A_pos[a], B_pos[a]))
      rmsd += sum([ (x-y)*(x-y) for x,y in zip(A_pos[a.name], B_pos[a.name]) ])
      N += 1

  return sqrt(rmsd/float(N))

###############################################################################
def get_prot_polar_atoms(residues, min_pt, max_pt, tol=3.5):
  """
  Get the polar protein atoms with interactions falling inside the box 
  defined by min_pt - tol and max_pt + tol
  """
  p_min_pt = min_pt[:]
  p_max_pt = max_pt[:]
  for i in range(3):
    p_min_pt[i] -= tol
    p_max_pt[i] += tol
    
  polar_atoms = [] 
  for res in residues:
    for atom in res.atoms:
      if(atom.interact_type == "ACCEPTOR" or atom.interact_type == "DONOR" or \
         atom.interact_type == "DONEPTOR"):

        keep_atom = True
        for i in range(3):
          if(atom.position[i] < p_min_pt[i] or p_max_pt[i] < atom.position[i]):
            keep_atom = False

        if(keep_atom): polar_atoms.append(atom)
  return polar_atoms
###############################################################################

###############################################################################
def get_hbond_atom_nbrs(residues, hbond_atom, C_name, other_name):
  """
  Get the carbon nbr and other atom with respect to hbond_atom
  """
  # Not sure about the "python" way to do this, but ...
  for i in range(len(residues)):
    if(residues[i]._residue__contains(hbond_atom)):
      if(hbond_atom.name == " O  "):
        #print "carbonyl O"
        for tmp_atom in residues[i]:
          if(tmp_atom.name.strip() == C_name): C_nbr = tmp_atom
    
        if(i+1 < len(residues)):
          for tmp_atom in residues[i+1]:
            if(tmp_atom.name.strip() == other_name): other_nbr = tmp_atom

           # distance check here -- to ensure residues are consecutive
      elif(hbond_atom.name == " N  "):
        #print "amide N"
        #print my_ideal_pts[0].atoms
        if(i > 0):
          for tmp_atom in residues[i-1]:
            if(tmp_atom.name.strip() == C_name): C_nbr = tmp_atom
            elif(tmp_atom.name.strip() == other_name): other_nbr = tmp_atom

           # distance check here -- to ensure residues are consecutive
      else:
        #print "other atom"
        for tmp_atom in residues[i]:
          if(tmp_atom.name.strip() == C_name): C_nbr = tmp_atom
          elif(tmp_atom.name.strip() == other_name): other_nbr = tmp_atom
      break
  return (C_nbr, other_nbr)
###############################################################################

###############################################################################
def renumber(fname):
  
  infile = open(fname, "r")

  lines = []
  cnt = 0
  for line in infile:
    if(line.startswith("ATOM  ") or line.startswith("HETATM")):
      cnt += 1
      line = "%s%5d%s" % (line[:6], cnt, line[11:].rstrip("\n"))
    lines.append(line)

  if(cnt > 99999):
    print "Too many atoms to renumber"
  else:
    infile.close()
    out_file = open(fname, "w+")
    for line in lines: print >> out_file, line
###############################################################################

###############################################################################
def transform(in_fname, out_fname, R, T):
  """
  Assumptions 
    Assumes R is for premultiplication by coordinates and 
    R and T are 3x3 and 3x1 numpy arrays 
    apply transformation to all ATOM and HETATM records
    transformation is x' = xR + T
  """
  
  infile = open(in_fname, "r")
  out = open(out_fname, "w+") 

  for line in infile:
    if(line.startswith("ATOM  ") or line.startswith("HETATM")):
      my_atom = atom(line)
      my_atom.position = dot(my_atom.position, R) + T
      print >> out, my_atom
    else:
      print >> out, line.rstrip("\r\n")
