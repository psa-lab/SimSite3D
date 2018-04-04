class sdf_molecules:

  def __init__(self, fname=None):
    
    self.mols = []

    try:
      sdf_file = file(fname, "r")
    except IOError, (errno, strerror):
      print "Unable to open the file", fname
      print "error(%s): %s" % (errno, strerror)
      return
  
    # It is likely that there is a more elegant way to do this, I don't 
    # have the time to search for it
    while 1:
      try: my_mol = molecule(sdf_file=sdf_file)
      except StopIteration:
        return

      if(len(my_mol.header_lines)):
        self.mols.append(my_mol)

  def sort_by_score(self, score_type="Tanimoto"):
    """
    We want to be able to sort in descending order by score
    The problem is FragGen can write out different scores and
    they are labeled by the type -- this means we must have as input
    the name of the score field we wish to sort
    """

    tups = [ (float(my_mol.features[score_type][0].strip(" \n")), my_mol) for my_mol in self.mols ]
    tups.sort(cmp=lambda x,y: cmp(x[0],y[0]), reverse=True)
    self.mols = [ x[1] for x in tups ]

class molecule:

  def __init__(self, fname=None, sdf_file=None):

    self.num_atoms = 0
    self.num_bonds = 0
    self.atoms = []
    self.bonds = []
    self.header_lines = []
    self.features = {}

    if(sdf_file == None):
      if(fname == None): return
      try:
        sdf_file = file(fname, "r")
      except IOError, (errno, strerror):
        print "Unable to open the file", fname
        print "error(%s): %s" % (errno, strerror)
        return

    try: self.__header_section(sdf_file)
    except StopIteration:
      self.header_lines = []
      raise StopIteration
      return
  
#    if(len(self.header_lines) < 4): 

    self.__atom_section(sdf_file)
    self.__bond_section(sdf_file)
    self.__feature_section(sdf_file)

  def build_header(self, mol_name, comments):
    self.num_atoms = len(self.atoms)
    self.num_bonds = len(self.bonds)
    self.header_lines = ["", "", "", ""]
    self.header_lines[0] = mol_name
    self.header_lines[1] = comments[0]
    self.header_lines[2] = comments[1]
    self.header_lines[3] = "%3d%3d  0  0  0  0  0  0  0  0  0" % \
      (self.num_atoms, self.num_bonds)

  def __header_section(self, sdf_file):
    for i in range(4):
      try: line = sdf_file.next()
      except StopIteration: raise StopIteration

      self.header_lines.append(line.rstrip())

    self.num_atoms = int(line[0:3])
    self.num_bonds = int(line[3:6])
       
  def __atom_section(self, sdf_file):
    for i in range(self.num_atoms):
      self.atoms.append(atom(sdf_file.next()))

  def __bond_section(self, sdf_file):
    for i in range(self.num_bonds):
      self.bonds.append(bond(sdf_file.next()))
 
  def __feature_section(self, sdf_file):
    tag = ""
    for line in sdf_file:
      if(line.startswith("> <")): 
        tag = line[3: line.rfind(">")]
        self.features[tag] = []
      elif(line.startswith("$$$$")): break
      elif(line.startswith("M  END")): continue
      else: self.features[tag].append(line.rstrip("\n"))

  def __repr__(self):
    lines = self.header_lines[:]
    for a in self.atoms:
      lines.append("%s" % (a))
    for b in self.bonds:
      lines.append("%s" % (b))

    lines.append("M END")
    
    for tag, data in self.features.iteritems():
      lines.append("> <%s>" % (tag))
      lines.extend(data)

    lines.append("$$$$")
    return "\n".join(lines)

class atom:

  def __init__(self, line=None):
    if(line == None):
      self.position = []
      self.element = ""
      self.other = []
      return

    toks = line.split()
    self.position = [ float(s) for s in toks[0:3] ]
    self.element = toks[3]
    self.other = toks[4:]
    
  def __repr__(self):
    return " %9.4f %9.4f %9.4f %-2s  %s" % \
      (self.position[0], self.position[1], self.position[2], self.element,
       "  ".join(self.other))

  def is_hydrogen(self):
    if(self.name == H or self.name == D):
      return True
    return False

class bond:

  def __init__(self, line=None):
    if(line == None):
      self.atom1 = -1
      self.atom2 = -1
      self.type = ""
      self.other = []
      return

    toks = line.split()
    self.atom1 = int(toks[0])
    self.atom2 = int(toks[1])

    # Crude way to check if RCSB sdf files have anything other than numeric
    # bond types
    self.type = int(toks[2])
    self.other = toks[3:]

  def __repr__(self):
    return "%3d%3d  %d  %s" % \
      (self.atom1, self.atom2, self.type, "  ".join(self.other))
