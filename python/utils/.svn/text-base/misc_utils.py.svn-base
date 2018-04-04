from ASCbasePy.utils import Quaternion
from numpy import zeros, eye, isnan
from sys import stderr

class trans_data:

  def __init__(self):
    self.T = zeros((3,))
    self.R = eye(3)

################################################################################

################################################################################
def load_lig_atoms_to_align(fname):
  """
  Load the names of the ligand heavy atoms to align.
  Assumptions:
    1 atom name per line
  """
  try:
    infile = open(fname, "r")
  except IOError, (errno, strerror):
    print >> stderr, "Unable to open the file", fname
    print >> stderr, "error(%s): %s" % (errno, strerror)
    return []

  atom_names = []
  for line in infile:
    if(not line.startswith("#") and len(line)):
      atom_names.append(line.strip("\n"))

  return atom_names
################################################################################

################################################################################
def load_align_data(fname):
  from numpy import array

  try:
    infile = open(fname, "r")
  except IOError, (errno, strerror):
    print "Unable to open the file", fname
    print "error(%s): %s" % (errno, strerror)
    return {}

  xforms = {} 
  for line in infile:
    if(line.startswith("#") or line == "\n"): continue

    toks = line.split("|")
    #(id, rot, trans, rmsd) = (toks[0], toks[1], toks[2], toks[3])
    my_trans = trans_data()
    rot = array([ float(s) for s in toks[1].split() ])
    if(rot.shape[0] == 4):
      q = Quaternion(Q=rot)
      my_trans.R = q.get_ortho_rot_mat()
    elif(rot.shape[0] == 9): my_trans.R = rot.reshape((3,3)) 
    else:  
      print >> stderr, "Incorrect size for rotation"
      continue

    my_trans.T = array([ float(s) for s in toks[2].split() ]) 
    xforms[toks[0]] = my_trans 
   
  return xforms 

################################################################################

################################################################################
def load_saved_hetgrps(pkl_fname):
  """
  Load a pickled dictionary of hetgroups
  """
  from os import path
  import cPickle

  if(path.exists(pkl_fname)):
    try:
      pkl_in = file(pkl_fname, "rb")
    except IOError, (errno, strerror):
      print "Unable to open the file", pkl_fname
      print "error(%s): %s" % (errno, strerror)
      print "in %s" % (self.__module__)
      return {}

    hetgrps = cPickle.load(pkl_in)
    pkl_in.close()
  else:
    print  """
Unable to load the het groups file 
  (%s) 
because it doesn't exist
""" % (pkl_fname)
    return {}

  return hetgrps
################################################################################

################################################################################
def load_uniq_csv(fname, max_above_50_cnt=4, max_seq_id=90.0):
  """
  Assumptions:
    The records in the file are sorted from most distinct to least distinct
    in terms of sequence identity
  """
  try:
    csv_file = file(fname, "r")
  except IOError, (errno, strerror):
    print "Unable to open the file", fname
    print "error(%s): %s" % (errno, strerror)
    return {}

  for line in csv_file:
    if(line.startswith("PDB ID")): break

  above_50_count = 0
  records = []
  chains = []
  for line in csv_file:
    line = line.strip("\n")
    if(len(line) == 0 or line.startswith("#")): continue
    toks = line.split("|")
    if(len(toks) < 3):
      print "Line is too short: %s", line
      continue

    # Rows may have been added for the SCOP/CATH classifications and we do not
    # want to include the "extra" rows
    if(len(toks[0]) == 0 and len(toks[1]) == 0): continue

    if(len(toks[2]) == 0 or float(toks[2]) <= 50.0):
      chains.append(":".join(toks[0:2]))
      records.append(toks)
    elif(float(toks[2]) < max_seq_id):
      chains.append(":".join(toks[0:2]))
      records.append(toks)
      above_50_count += 1
      if(above_50_count >= max_above_50_cnt): break

  return (chains, records)
################################################################################

################################################################################
def load_FASTAs(FASTA_dir):
  import os,gzip
  
  if(not os.path.isdir(FASTA_dir)):
    print "Cannot find the directory: %s" % (FASTA_dir)
    return {}

  seqs = {}
  seq = ""
  header = ""
  first_header = True
  for fname in os.listdir(FASTA_dir):
    if(fname.find("fasta") == -1): continue
    if(fname == "fasta.tst"): continue

    my_fname = "%s/%s" % (FASTA_dir, fname)
    print  my_fname
    if(fname.endswith(".gz")): infile = gzip.open(my_fname)
    else: infile = file(my_fname)

    # PDB FASTAs have a header followed by one or more sequence lines --
    # this is because of the 80 character per line limit (a good thing)
    for line in infile:
      if(line.startswith(">")):
        if(first_header): first_header = False
        else:
          if(seq in seqs): seqs[seq].headers.append(header)
          else: seqs[seq] = info(header)
        header = line.strip("\n")
        seq = ""
      else:
        seq += line.strip("\n")

  if(seq in seqs): seqs[seq].headers.append(header)
  else: seqs[seq] = info(header)
  return seqs
################################################################################

################################################################################
def load_oocalc_csv(csv_fname):
  """
  Assumptions: CSV file is created by copying and pasting PDB Ligand Expo
  search results table into oocalc.  The CSV file is delimited by pipes ('|') 
  and has a break in the line due to a return between the crystallographic 
  article and the authors.  There are no blank lines and the first line contains
  the column headings.
  """

  try:
    csv_file = file(csv_fname, "r")
  except IOError, (errno, strerror):
    print "Unable to open the file: %s" % (csv_fname)
    print "Error(%s): %s" % (errno, strerror)
    print "in %s" % (self.__module__)
    return {}

  headings = csv_file.next().split("|")
  lig_data = {}
  for line in csv_file:
    toks = line.split("|")
    lig_data[toks[0].upper()] = toks[1:]

  return lig_data
################################################################################

################################################################################
def load_pdb_file(pdb_id, save_dir=".", local_pdb_dir="/psa/pdb",
                  pdb_ftp_url = "ftp://ftp.wwpdb.org/pub/pdb/data/" + \
                  "structures/divided/pdb"):
  """
  Load the pdb file ("pdbXXXX.ent") from local_pdb_dir; if not found, look for 
  it in save_dir; if not in save_dir, get it from RCSB
  """
  from os import path
  from os import system
  from ASCbasePy.utils import pdb

  # Load the corresponding PDB file 
  prot_fname = "%s/pdb%s.ent" % (local_pdb_dir, pdb_id.lower())
  if(not path.exists(prot_fname) and
     not path.exists("%s/pdb%s.ent" % (save_dir, pdb_id.lower()))):

    prot_url = "%s/%s/pdb%s.ent.gz" % \
      (pdb_ftp_url, pdb_id.lower()[1:3], pdb_id.lower())
    system("wget " + prot_url)
    system("gunzip pdb%s.ent.gz" % (pdb_id.lower()))
    system("mv pdb%s.ent %s" % (pdb_id.lower(), save_dir))
    prot_fname = "%s/pdb%s.ent" % (save_dir, pdb_id.lower())
  elif(not path.exists(prot_fname)):
    prot_fname = "%s/pdb%s.ent" % (save_dir, pdb_id.lower())
  prot = pdb.residues(prot_fname)
  if(prot.fail):
    print "Could not open the protein file %s" % (prot_fname)
    return ("",None)

  return (prot_fname, prot)
################################################################################

################################################################################
def load_clustal_pairwise_seq(fname):
  """
  This works for clustal.out from EBI web interface
  If using the installed clustal, we must save at least
 "
Sequence format is
Sequence 1: 
 .
 .
 .
Sequence N:

Aligning...

Sequences (1:2) Aligned. Score: 1 
 .
 .
 .
Sequences (N-1:N) Aligned. Score: 100
"
  """
  try:
    infile = file(fname, "r")
  except IOError, (errno, strerror):
    print "Unable to open the file: %s" % (fname)
    print "Error(%s): %s" % (errno, strerror)
    print "in %s" % (self.__module__)
    return {}

  # Get line before the first sequence line
  for line in infile:
    if(line.startswith("Sequence format is")): break
  
  # Get sequences -- assume no spaces
  ids = []
  for line in infile:
    toks = line.split()
    if(len(toks) < 2): break
    ids.append(line.split()[2].split("|")[0])
  
  # Between the ids and the pairwise sequence lines there are some additional
  # lines
  for line in infile:
    if(line.startswith("Aligning")):
      infile.next()
      break

  # For now we will suppress printing of the matrix
  S = zeros((len(ids), len(ids)))
  old_row = -1
  for line in infile:
    # Example Sequences line
    #Sequences (1:2) Aligned. Score:  3 
    if(not line.startswith("Sequences")): break
  
    # Don't make any assumptions -- not sure if there are blank or missing 
    # entries if there is no sequence identity or if CLUSTAL ever gets a zero
    # score
    toks = line.split()
    row = int(toks[1].strip("()").split(":")[0]) - 1
    if(not row  == old_row):
      old_row = row
      col = row + 1
  
    # Make this a square matrix
    seq_id = int(toks[4])
    S[row,col] = seq_id
    S[col,row] = seq_id
    col += 1

  return (ids, S)
  
################################################################################

################################################################################


################################################################################

################################################################################
class norm_db_data_entry:

  def __init__(self, line):
    toks = line.rstrip("|\n").split("|")
    self.site_id = toks[0]
    self.long_name = toks[1]
    self.short_name = toks[2]
    self.lig_name = toks[3]
    if(toks[4] == "0"): self.binds_ADE = False
    else: self.binds_ADE = True
    if(toks[5] == "0"): self.binds_PTR = False
    else: self.binds_PTR = True
    if(toks[6] == "0"): self.is_GST_like = False
    else: self.is_GST_like = True

################################################################################

################################################################################
class norm_db_data:

  def __init__(self, fname):

    infile = open(fname, "r")
    self.records = {}
    for line in infile:
      if(line.startswith("site id") or line.startswith("#")): continue

      data = norm_db_data_entry(line)
      self.records[data.site_id] = data
################################################################################

################################################################################
class query_data:

  def __init__(self, query_id):
    (self.best_scores, self.rmsd_of_best_scores, self.best_samp_rmsd) = \
      load_data(query_id)
    (self.best_norm_scores, empty_0, empty_1) = \
      load_data(query_id, norm_db=True)

################################################################################

################################################################################
def load_query_data(res_file, rmsd_file=None, init_mu=0.0, init_sigma=1.0,
                    ICP_mu=0.0, ICP_sigma=1.0):
  init_data = {}
  ICP_data = {}
  for line in res_file:
    if(line.startswith("#")): continue

    toks = line.split("|")
    t_id = toks[0]
    score = float(toks[1])
    rmsd = -1.0
    if(not rmsd_file == None): rmsd = float(rmsd_file.next().split("|")[0])
    if(t_id.startswith("ICP")):
      ICP_data[t_id[4:]] = ((score - ICP_mu)/ICP_sigma, rmsd)
    else: init_data[t_id] = ((score - init_mu)/init_sigma, rmsd)

  return (init_data, ICP_data)
################################################################################

################################################################################
def load_query_norm_data(res_fname, score_tol=-1.5):
  res_file = open(res_fname, "r")
  (init_scores, ICP_scores) = ([], [])
  (ids, ICP_ids) = ([], [])
  for line in res_file:
    if(line.startswith("#")): continue

    toks = line.split("|")
    score = float(toks[1])
    if(line.startswith("ICP")):
      ICP_scores.append(score)
      ICP_ids.append(toks[0][4:])
    else:
      init_scores.append(score)
      ids.append(toks[0])
  init_scores = array(init_scores)
  ICP_scores = array(ICP_scores)

  (mu, sigma) = (mean(init_scores), std(init_scores))
  (ICP_mu, ICP_sigma) = (mean(ICP_scores), std(ICP_scores))

  norm_scores = ( init_scores - mu ) / sigma
  norm_ICP_scores = ( ICP_scores - ICP_mu ) / ICP_sigma

  (sig_ids, ICP_sig_ids) = ([], [])
  for s, id in zip(norm_scores, ids):
    if(s <= score_tol): sig_ids.append(id)
  for s, id in zip(norm_ICP_scores, ICP_ids):
    if(s <= score_tol): sig_ids.append(id)

  return((mu, sigma), (ICP_mu, ICP_sigma), sig_ids, ICP_sig_ids)

################################################################################

################################################################################
def load_data(data_dir, norm_db_data_dir="", rigid=False):
  import os, gzip

  (sig_ids, ICP_sig_ids) = ({}, {})

  if(len(norm_db_data_dir)):
    init_norm_db_data = {}
    ICP_norm_db_data = {}
    for res_dir in os.listdir(norm_db_data_dir):
      if(not res_dir.endswith("results")): continue
      if(res_dir.endswith("_all_results")): continue

      q_id = res_dir[:-8]
      res_fname = "%s/%s/%s_blast.out" % (norm_db_data_dir, res_dir, q_id)
      if(rigid): res_fname += ".rigid"
      (init_norm_db_data[q_id], ICP_norm_db_data[q_id], sig_ids[q_id],
       ICP_sig_ids[q_id]) = load_query_norm_data(res_fname)

  init_data = {}
  ICP_data = {}
  for res_dir in os.listdir(data_dir):
    if(not res_dir.endswith("results")): continue
    if(res_dir.endswith("_all_results")): continue

    q_id = res_dir[:-8]
    res_fname = "%s/%s/%s_blast.out" % (data_dir, res_dir, q_id)
    rmsd_fname = "%s/%s/%s_rmsd.out" % (data_dir, res_dir, q_id)
    if(rigid):
      res_fname += ".rigid"
      rmsd_fname += ".rigid"

    if(os.path.exists("%s.gz" % (res_fname))):
      res_file = gzip.open("%s.gz" % (res_fname), "r")
      rmsd_file = gzip.open("%s.gz" % (rmsd_fname), "r")
    else:
      res_file = open(res_fname, "r")
      rmsd_file = open(rmsd_fname, "r")

    if(len(norm_db_data_dir)):
      (init_data[q_id], ICP_data[q_id]) = \
        load_query_data(res_file, rmsd_file, init_mu=init_norm_db_data[q_id][0],
                        init_sigma=init_norm_db_data[q_id][1],
                        ICP_mu=ICP_norm_db_data[q_id][0],
                        ICP_sigma=ICP_norm_db_data[q_id][1])
    else:
      (init_data[q_id], ICP_data[q_id]) = load_query_data(res_file, rmsd_file)
    
  return (init_data, ICP_data, sig_ids, ICP_sig_ids)
################################################################################

################################################################################
def load_terms(data_dir, rigid=False):
  import os, gzip

  init_terms = {}
  ICP_terms = {}
  for res_dir in os.listdir(data_dir):
    if(not res_dir.endswith("results")): continue
    if(res_dir.endswith("_all_results")): continue

    q_id = res_dir[:-8]
    res_fname = "%s/%s/%s_blast.out" % (data_dir, res_dir, q_id)
    if(rigid): res_fname += ".rigid"

    if(os.path.exists("%s.gz" % (res_fname))):
      res_file = gzip.open("%s.gz" % (res_fname), "r")
    else: res_file = open(res_fname, "r")

    row_ICP_terms = {}
    row_init_terms = {}
    for line in res_file:
      if(line.startswith("#")): continue    

      toks = line.rstrip("|\n").split("|")
      dset_id = toks[0]
      terms = [ float(s) for s in toks[6].split() ]
      for t in terms:
        if(isnan(t)):
          print res_fname, q_id, toks[0], terms
      if(dset_id.startswith("ICP")): row_ICP_terms[dset_id[4:]] = terms
      else: row_init_terms[dset_id] = terms

    ICP_terms[q_id] = row_ICP_terms
    init_terms[q_id] = row_init_terms
    
  return (init_terms, ICP_terms)
################################################################################

################################################################################
def new_load_terms(data_dir):
  import os, gzip
 
  init_terms = {}
  preIK_terms = {}
  final_terms = {}

  for res_dir in os.listdir(data_dir):
    if(not res_dir.endswith("results")): continue
    if(res_dir.endswith("_all_results")): continue

    q_id = res_dir[:-8]
    res_fname = "%s/%s/%s_blast.out" % (data_dir, res_dir, q_id)

    if(os.path.exists("%s.gz" % (res_fname))):
      res_file = gzip.open("%s.gz" % (res_fname), "r")
    else: res_file = open(res_fname, "r")
    #print "fname:", res_fname

    row_init_terms = {}
    row_preIK_terms = {}
    row_final_terms = {}
    for line in res_file:
      if(line.startswith("#")): continue

      toks = line.rstrip("|\n").split("|")
      dset_id = toks[0]
      #print q_id, "line (", line, ")"
      terms = [ float(s) for s in toks[6].split() ]
      if(dset_id.startswith("init_")): row_init_terms[dset_id[5:]] = terms
      elif(dset_id.startswith("preIK_")): row_init_terms[dset_id[6:]] = terms
      else: row_final_terms[dset_id] = terms
      for t in terms:
        if(isnan(t)):
          print res_fname, q_id, toks[0], terms

    init_terms[q_id] = row_init_terms
    preIK_terms[q_id] = row_preIK_terms
    final_terms[q_id] = row_final_terms
    
  return (init_terms, preIK_terms, final_terms)
################################################################################

################################################################################
def new_load_data(data_dir, norm_db_data_dir="", rigid=False):
  import os, gzip

  print "NOT FINISHED"

  (sig_ids, preIK_sig_ids, init_sig_ids) = ({}, {})

  if(len(norm_db_data_dir)):
    init_norm_db_data = {}
    preIK_norm_db_data = {}
    norm_db_data = {}
    for res_dir in os.listdir(norm_db_data_dir):
      if(not res_dir.endswith("results")): continue
      if(res_dir.endswith("_all_results")): continue

      q_id = res_dir[:-8]
      res_fname = "%s/%s/%s_blast.out" % (norm_db_data_dir, res_dir, q_id)
      if(rigid): res_fname += ".rigid"
#      (init_norm_db_data[q_id], ICP_norm_db_data[q_id], sig_ids[q_id],
#       ICP_sig_ids[q_id]) = load_query_norm_data(res_fname)

#  init_data = {}
#  ICP_data = {}
#  for res_dir in os.listdir(data_dir):
#    if(not res_dir.endswith("results")): continue
#    if(res_dir.endswith("_all_results")): continue
#
#    q_id = res_dir[:-8]
#    res_fname = "%s/%s/%s_blast.out" % (data_dir, res_dir, q_id)
#    rmsd_fname = "%s/%s/%s_rmsd.out" % (data_dir, res_dir, q_id)
#    if(rigid):
#      res_fname += ".rigid"
#      rmsd_fname += ".rigid"
#
#    if(os.path.exists("%s.gz" % (res_fname))):
#      res_file = gzip.open("%s.gz" % (res_fname), "r")
#      rmsd_file = gzip.open("%s.gz" % (rmsd_fname), "r")
#    else:
#      res_file = open(res_fname, "r")
#      rmsd_file = open(rmsd_fname, "r")

    if(len(norm_db_data_dir)):
      (init_data[q_id], ICP_data[q_id]) = \
        load_query_data(res_file, rmsd_file, init_mu=init_norm_db_data[q_id][0],
                        init_sigma=init_norm_db_data[q_id][1],
                        ICP_mu=ICP_norm_db_data[q_id][0],
                        ICP_sigma=ICP_norm_db_data[q_id][1])
    else:
      (init_data[q_id], ICP_data[q_id]) = load_query_data(res_file, rmsd_file)
    
  return (init_data, ICP_data, sig_ids, ICP_sig_ids)
