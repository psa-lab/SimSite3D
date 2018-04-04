class molecule:

  def __init__(self, fname=None):

    self.num_atoms = 0
    self.num_bonds = 0
    self._fail = False
    if(fname == None): 
      self.atoms = []
      self.bonds = []
      return

    try:
      mol2_file = file(fname, "r")
    except IOError, (errno, strerror):
      print "Unable to open the file", fname
      print "error(%s): %s" % (errno, strerror)
      self._fail = True
      return

    line = mol2_file.readline()
    while line:
      if(line.startswith("@<TRIPOS>MOLECULE")):
        self.__mol_section(mol2_file)
      elif(line.startswith("@<TRIPOS>ATOM")):
        self.__atom_section(mol2_file)
      elif(line.startswith("@<TRIPOS>BOND")):
        self.__bond_section(mol2_file)

      line = mol2_file.readline()
 
  def fail(self):
    return self._fail

  def __mol_section(self, mol2_file):
    self.mol_name = mol2_file.readline().strip()
    line = mol2_file.readline()
    toks = line.split()
    self.num_atoms = int(toks[0])
    self.num_bonds = int(toks[1])
       
  def __atom_section(self, mol2_file):
    self.atoms = []
    for i in range(self.num_atoms):
      self.atoms.append(atom(mol2_file.readline()))

  def __bond_section(self, mol2_file):
    self.bonds = []
    for i in range(self.num_bonds):
      self.bonds.append(bond(mol2_file.readline()))

  def __repr__(self, fixed_format=False):
    mol_str = """@<TRIPOS>MOLECULE
%s
%5d %5d     0     0     0
SMALL
USER_CHARGES

@<TRIPOS>ATOM""" % (self.mol_name, len(self.atoms), len(self.bonds))
    for a in self.atoms: mol_str += "\n" + str(a)
    mol_str += "\n@<TRIPOS>BOND"
    for b in self.bonds: mol_str += "\n" + str(b)
    return mol_str

class atom:

  def __init__(self,line=None):
    if(line == None):
      return

    toks = line.split()
    self.serial = int(toks[0])
    self.name_str = toks[1]
    self.position = [ float(s) for s in toks[2:5] ]

    # work around for fixed format Pfizer mol2 files
    if(len(toks[5]) > 5):
      toks.append(toks[5][5:])

    tmp = toks[5].split(".")
    self.name = tmp[0]
    if(len(tmp) > 1): self.orbital = tmp[1].rstrip()
    else: self.orbital = ""
    self.rest_of_line = " ".join(toks[6:])
    
  def __cmp__(self, other):
    if(self.serial < other.serial):
      return -1
    elif(self.serial > other.serial):
      return 1
    else:
      return 0

  def __repr__(self, fixed_format=False):
    if(self.orbital): atom_type = "%s.%s" % (self.name, self.orbital)
    else: atom_type = self.name
    
    if(fixed_format):
      return "%7d %4s %8.4f %8.4f %8.4f %-5s%s" % \
        (self.serial, self.name_str, self.position[0], self.position[1],
         self.position[2], atom_type, self.rest_of_line)
    else:
      return "%7d %4s %8.4f %8.4f %8.4f %-5s %s" % \
        (self.serial, self.name_str, self.position[0], self.position[1],
         self.position[2], atom_type, self.rest_of_line)

  def is_hydrogen(self):
    if(self.name == H or self.name == D):
      return True
    return False

class bond:

  def __init__(self, line):
    toks = line.split()
    self.number = int(toks[0])
    self.atom1 = int(toks[1])
    self.atom2 = int(toks[2])
    self.type = toks[3]
    self.status_bits = ""
    if(len(toks) > 3):
      self.status_bits = " ".join(toks[4:])

  def __repr__(self):
    return "%6d %4d %4d %s %s" % \
      (self.number, self.atom1, self.atom2, self.type, self.status_bits)
