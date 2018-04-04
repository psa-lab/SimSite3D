import gzip
import os
from sys import stderr

# should also bring up a window with the sorted hits listing the score, and id,
# and state number -- this is a makeshift way to display that information
class hits:

  def __init__(self, fname, num_hits_per_pair=1):
    """
    num_hits_per_pair: number of hits to load per query,database pair of
      site maps
    """
    self.hits = [] 

    try:
      if(fname.endswith(".gz")):
        hits_file = gzip.open(fname, "r")
      else:
        hits_file = file(fname, "r")
    except IOError, (errno, strerror):
      print "Unable to open the file", fname
      print "error(%s): %s" % (errno, strerror)
      # this raise is here to allow the exception to proceed to the Pmw widgets
      raise
      return

    self._read_header(hits_file) 
    self._read_hits(hits_file, num_hits_per_pair)

  def _read_header(self, hits_file):

    for line in hits_file:
      if(not line.startswith("#")):
        # return to beginning of file -- buffering is an issue -- can not use
        # tell & seek when using iterators over file
        
        # gzip seek takes only 1 parameter and has no __doc__ string so not
        # sure how to use it.
        if("rewind" in dir(hits_file)): hits_file.rewind() 
        elif("seek" in dir(hits_file)): hits_file.seek(0, 0)
        else:
          print >> stderr, "Unable to rewind the hits file"
        return

      if(line.startswith("# Local start time:")):
        self.start_time = line.split(":")[1].strip()
      elif(line.startswith("# Searchable sitemaps directory:")):
        self.dbase_dir = line.split(":")[1].strip()
      elif(line.startswith("# Corresponding ligands directory:")):
        self.ligs_dir = line.split(":")[1].strip()
      elif(line.startswith("# Query (Model) Sitemap Name:")):
        (self.query_dir, self.query) = \
          os.path.split(line.split(":")[1].strip())
        if(self.query.endswith("_s.csv")): self.query = self.query[:-6]
      elif(line.startswith("# Working directory (via getcwd()):")):
        self.search_dir = line.split(":")[1].strip()

  def _read_hits(self, hits_file, num_hits_per_pair=1):
    (found_init, found_preIK) = (False, False)
   
    ids = {}
    for line in hits_file:
      if(line.startswith("#")):
        continue

      # Assume each record is has a terminating "|"
      toks = line.split("|")[:-1]
      id = toks[0]
      if(id.startswith("init_")): 
        if(not found_init):
          print "NOTE: ignoring initial alignments (lines that start with " + \
                "init_)"
        found_init = True
        continue
      elif(id.startswith("preIK_")):
        if(not found_preIK):
          print "NOTE: assuming IK aligns don't have valid transformations"
          print "\tand the corresponding preIK line follows"
        self.hits[-1] = (toks[0][6:], self.hits[-1][1],
                         [float(s) for s in toks[2].split()],
                         [float(s) for s in toks[3].split()], 
                         self.hits[-1][4], self.hits[-1][5])
        found_preIK = True
        continue

      if(id in ids):
        if(len(ids[id]) >= num_hits_per_pair): continue
        ids[id].append(True)
      else:
        ids[id] = [True]

      self.hits.append( (toks[0], float(toks[1]), 
                         [float(s) for s in toks[2].split()], 
                         [float(s) for s in toks[3].split()], toks[4], toks[5]))

    def tuple_cmp(A,B):
      if(A[1] < B[1]): return -1
      elif(A[1] > B[1]): return 1
      return 0

    self.hits.sort(tuple_cmp)
