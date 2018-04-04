import gzip 
from sys import stderr

class info:

  def __init__(self, header):
    self.headers = [header]
################################################################################

################################################################################
class row:

  def __init__(self, entry_line):
    toks = entry_line.split("\t")
    
    self.id = toks[0]
    self.header = toks[1].title()
    self.date = toks[2]
    self.compound = toks[3].title()
    self.source = toks[4].capitalize()
    self.authors = toks[5]
    if(toks[6] == "NOT"): self.resolution = -1.0
    else: self.resolution = float(toks[6])
    self.exp_method = toks[7].rstrip("\n")

    self.fastas = {}
    self.name = ""
  
  def add_fasta(self, fasta_seq, chainID, name):
    self.fastas[chainID] = fasta_seq
    self.name = name.title()

  def __str__(self):
    if(self.resolution > 0.0): res_str = "%.2f" % (self.resolution)
    else: res_str = "NOT"
    toks = [self.id, self.name, self.header, self.date, self.compound, 
            self.source, self.authors, res_str, self.exp_method]
    my_str = " ".join(toks)
    for chain in self.fastas:
      my_str += "\n>%s_%s\n%s" % (self.id, chain, self.fastas[chain])
    return my_str
################################################################################

################################################################################
class table:
  
  def __init__(self, entries_fname, fastas_fname):
    self.records = {}
    self.load_entries(entries_fname)
    self.load_fasta_sequences(fastas_fname)

  def load_entries(self, fname):
    try:
      if(fname.endswith(".gz")): infile = gzip.open(fname, "r")
      else: infile = open(fname, "r")
    except IOError, (errno, strerror):
      print "Unable to open the file: %s" % (fname)
      print "Error(%s): %s" % (errno, strerror)
#      print "in %s" % (self.__module__)
      return {}

    # Skip the first 2 lines
    for i in range(2): infile.next()

    for line in infile:
      self.records[line[0:4].lower()] = row(line) 

  def load_fasta_sequences(self, fname):
    try:
      if(fname.endswith(".gz")): infile = gzip.open(fname, "r")
      else: infile = open(fname, "r")
    except IOError, (errno, strerror):
      print "Unable to open the file: %s" % (fname)
      print "Error(%s): %s" % (errno, strerror)
#      print "in %s" % (self.__module__)
      return {}

    for desc_line in infile:
      fasta = infile.next().strip("\n") 

      pdb_id = desc_line[1:5].lower()
      chain = desc_line[6]
      name = desc_line[desc_line.rfind("length:") + 7:].lstrip("0123456789")
      name = name.strip(" \n")
      self.records[pdb_id].add_fasta(fasta, chain, name)

################################################################################

################################################################################
  def get_FASTAs(self, pdb_ids, chainIDs=[]):
    """
    Return the sequences in FASTA format for all the unique chains associated
    with each pdb_id.  Or if chainIDs has a length, return the sequence for
    each (pdb_id, chainID) pair
    """
    seqs = {}
 
    if(not len(chainIDs)): 
      for pdb_code in pdb_ids:
        if(not pdb_code.lower() in self.records):
          print >> stderr, "Warning: Unable to find entry for:", pdb_code
          continue

        pdb_rec = self.records[pdb_code.lower()]
  
        for chain, fasta in pdb_rec.fastas.iteritems():
          my_header = ">%s:%s|PDBID|CHAIN|SEQUENCE" % (pdb_code.upper(), chain)
          if(fasta in seqs): seqs[fasta].headers.append(my_header)
          else: seqs[fasta] = info(my_header)
    else:
      for pdb_code, chainID in zip(pdb_ids, chainIDs):
        if(not pdb_code.lower() in self.records):
          print >> stderr, "Warning: Unable to find entry for:", pdb_code
          continue

        pdb_rec = self.records[pdb_code.lower()]
  
        for chain, fasta in pdb_rec.fastas.iteritems():
          if(chain == chainID):
            my_header = ">%s:%s|PDBID|CHAIN|SEQUENCE" % \
              (pdb_code.upper(), chain)
            if(fasta in seqs): seqs[fasta].headers.append(my_header)
            else: seqs[fasta] = info(my_header)
            break

    return seqs

################################################################################

################################################################################
if(__name__ == "__main__"):
  
  entries_fname = "entries.idx" 
  fastas_fname = "pdb_seqres.txt.gz"
  t = table(entries_fname, fastas_fname)

  print t.records["1ua2"]
  
