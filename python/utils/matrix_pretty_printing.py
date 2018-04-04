from sys import stderr

def load_column_order_file(path, dataset=""):
  fname = path + "/"
  if(len(dataset)):
    fname += dataset + "_column_order.txt"
  else:
    fname += "column_order.txt"

  try:
    bind_res_file = file(fname, "r")
  except IOError, (errno, strerror):
    print "Unable to open the file", fname
    print "error(%s): %s" % (errno, strerror)
    return ([], [], [])

  c_order = []
  short_prot_names = []
  prot_fam = []
  for line in bind_res_file:
    if(line.startswith("#")):
      continue

    toks = line.split("|")[:-1]
    c_order.append( toks[0].strip() )
    short_prot_names.append( toks[2].strip() )
    if(len(toks) > 3):
      prot_fam.append(toks[3].strip())
  bind_res_file.close()
  return (c_order, short_prot_names, prot_fam)
#####################################


class csv:
#####################################

  #####################################
  def __init__(self, rmsd_tol=2.0, score_cutoff=-1.5):
    self.reference_rmsds = {}
    self.sampled_rmsds = {}
    self.scores = {} 
    self.rmsd_of_scores = {}
    self.SSM_Q_score = {}
    self.SSM_seq_id = {}
    self.rmsd_tol = rmsd_tol
    self.score_cutoff = score_cutoff

  #####################################
  def get_pairs(self, SF_lab):
    my_pairs = []
    my_rmsd = self.rmsd_of_scores[SF_lab]
    my_score = self.scores[SF_lab]
    for q in my_score:
      for t in my_score[q]:
        my_pairs.append( (my_rmsd[q][t], my_score[q][t]) )

    return my_pairs

  #####################################
  def load_reference_rmsd(self, fname):
    (lab, tbl) = self.__read_table__(",", file(fname, "r"))
    self.reference_rmsds = tbl
    
#    infile = file(fname, "r")
#    line = infile.readline()
#    col_labs = [ s.strip() for s in line.split(",")[1:-1] ]
#    for i in range(len(col_labs)):
#      line = infile.readline()
#      toks = line.split(",")
#      self.reference_rmsds[toks[0].strip()] = {}
#      # need to be messy since some of the fields could be blank
#      for r,t_id in zip(toks[1:], col_labs):
#        r = r.strip()
#        self.reference_rmsds[toks[0].strip()][t_id] = float(r)
#	      #if(not r == ""):
#          #self.reference_rmsds[toks[0]][t_id] = float(r)
#	      #else:
#          #self.reference_rmsds[toks[0]][t_id] = 1000.0

  #####################################
  def load_best_sampled_rmsd(self, fname):
    (lab, tbl) = self.__read_table__(",", file(fname, "r"))
    self.sampled_rmsds = tbl

#    infile = file(fname, "r")
#    line = infile.readline()
#    col_labs = [ s.strip() for s in line.split(",")[1:-1] ]
#    for i in range(len(col_labs)):
#      line = infile.readline()
#      toks = line.split(",")
#      self.sampled_rmsds[toks[0]] = {}
#      # need to be messy since some of the fields could be blank
#      for r,t_id in zip(toks[1:], col_labs):
#        r = r.strip()
#        self.sampled_rmsds[toks[0]][t_id] = float(r)
#	  #    if (not r == ""):
#    #      self.sampled_rmsds[toks[0]][t_id] = float(r)
#	  #    else:
#    #      self.sampled_rmsds[toks[0]][t_id] = 1000.0

  #####################################
  def load_best_scores(self, fname):
    infile = file(fname, "r")
  
    line = infile.readline()
    SF_lab = ""
    self.error_rate = {}
    while(len(line)):
      if(line.startswith("SF") or line.startswith("two tiered")):
        (SF_lab, score_tbl) = self.__read_table__(line, infile)
        self.scores[SF_lab] = score_tbl
      elif(line.startswith("Average squared error rate,")):
        self.error_rate[SF_lab] = float(line.split(",")[1])
      line = infile.readline()

  #####################################
  def load_rmsd_of_best_scores(self, fname):
    infile = file(fname, "r")

    line = infile.readline()
    SF_lab = ""
    self.error_rate = {}
    while(len(line)):
      if(line.startswith("SF") or line.startswith("two tiered")):
        (SF_lab, rmsd_tbl) = self.__read_table__(line, infile)
        self.rmsd_of_scores[SF_lab] = rmsd_tbl
      line = infile.readline()

  #####################################
  def load_SSM_Q_score(self, fname):
    infile = file(fname, "r")

    line = infile.readline()
    while(len(line)):
      if(line.startswith("database binding residues")):
        (lab, tbl) = self.__read_table__(line, infile)
        self.SSM_Q_score = tbl
      line = infile.readline()

  #####################################
  def load_SSM_seq_id(self, fname):
    infile = file(fname, "r")

    line = infile.readline()
    while(len(line)):
      if(line.startswith("database binding residues")):
        (lab, tbl) = self.__read_table__(line, infile)
        self.SSM_seq_id = tbl
      line = infile.readline()

  #####################################
  def write_rmsd_tbl(self, ofile, col_order=None, table_name="best_sampled",
                     trunc_col_heads=True, remove_underscores=True):
    rmsd_tbl = {}
    if(table_name == "best_sampled"):
      rmsd_tbl = self.sampled_rmsds
    elif(table_name == "top_scored"):
      rmsd_tbl = self.rmsd_of_scores
    elif(table_name == "reference"):
      rmsd_tbl = self.reference_rmsds
    else:
      print >> stderr, "Unsupported table type (%s), unable to write anything" % (table_name)
      return

    if(len(rmsd_tbl) < 1):
      print >> stderr, "Empty %s rmsd table, unable to write anything" % (table_name)
      return

    N = len(rmsd_tbl)
    print >> ofile, "Queries,Database Sitemaps,"

    # Print the column headings
    if(col_order == None):
      first_key = rmsd_tbl.keys()[0]
      print >> ofile, " , " + ", ".join(rmsd_tbl[first_key].keys()) + ","
    elif(trunc_col_heads):
      print >> ofile, " , " + ", ".join([ lab[0:4] for lab in col_order ]) + ","
    else:
      print >> ofile, " , " + ", ".join(col_order) + ","
      
    # print the data
    if(col_order == None):
      for q_id, row in rmsd_tbl.iteritems():
        print >> ofile, q_id.replace("_", " "),
        for t_id, rmsd in row.iteritems():
          print >> ofile, ", %6.2f" % (rmsd),
        print >> ofile, ","
    else:
      rmsd_sum = 0.0
      rmsd_sum_ss2 = 0.0
      for q_id in col_order:
        if(remove_underscores):
          print >> ofile, q_id.replace("_", " "),
        else:
          print >> ofile, q_id,

        row = rmsd_tbl[q_id]
        for t_id in col_order:
          print >> ofile, ", %6.2f" % (row[t_id]),

	  # extremely crude measure to toss out outliers
          if(row[t_id] < 100.0):
            rmsd_sum += row[t_id]
            rmsd_sum_ss2 += row[t_id]*row[t_id]
        print >> ofile, ","

    N = len(col_order)*len(col_order)
    mu = rmsd_sum / N
    from math import sqrt
    sigma = sqrt(rmsd_sum_ss2 / N - mu*mu)
    print >> ofile, "%s rmsd stats (mean; std) (%6.2f; %6.2f)," % (table_name, mu, sigma)

  #####################################
  def diff_rmsds(self, old_mat, old_SF_lab, new_SF_lab):
    old_score_rmsds = old_mat.rmsd_of_scores[old_SF_lab]
    new_score_rmsds = self.rmsd_of_scores[new_SF_lab]
    old_samp_rmsds = old_mat.sampled_rmsds
    new_samp_rmsds = self.sampled_rmsds

    diff_mat = csv(rmsd_tol=0.0)
    for q in self.scores[new_SF_lab]:
      diff_mat.rmsd_of_scores[q] = {}
      diff_mat.sampled_rmsds[q] = {}
      for t in self.scores[new_SF_lab][q]:
        diff_mat.rmsd_of_scores[q][t] = \
          old_score_rmsds[q][t] - new_score_rmsds[q][t]
        diff_mat.sampled_rmsds[q][t] = \
          old_samp_rmsds[q][t] - new_samp_rmsds[q][t]

    return diff_mat

  #####################################

  #####################################
  def __read_table__(self, SF_line, infile):
    SF_lab = SF_line.split(",")[0]
    col_labs = [ s.strip() for s in (infile.readline()).split(",")[1:-1] ] 
    tbl = {}
    for junk in col_labs:
      toks = (infile.readline()).rstrip(",\n").split(",")
      tbl[toks[0]] = {}
      for s,t_id in zip(toks[1:], col_labs):
        s = s.strip()
        #print toks[0], t_id,
        if(s == ""):
          tbl[toks[0]][t_id] = 1e+37
        #  print 1e+37
        else:
          tbl[toks[0]][t_id] = float(s)
        #  print float(s)
    return (SF_lab, tbl)


#####################################
class tex:

  #####################################
  def __init__(self, csv_obj, rmsd_tol=2.0, score_cutoff=-1.5):
    self.sampled_rmsds = csv_obj.sampled_rmsds
    self.scores = csv_obj.scores
    self.rmsd_of_scores = csv_obj.rmsd_of_scores
    self.rmsd_tol = rmsd_tol
    self.score_cutoff = score_cutoff

  #####################################
#  def get_pairs(self, SF_lab):
#    my_pairs = []
#    my_rmsd = self.rmsd_of_scores[SF_lab]
#    my_score = self.scores[SF_lab]
#    for q in my_score:
#      for t in my_score[q]:
#        my_pairs.append( (my_rmsd[q][t], my_score[q][t]) )
#    return my_pairs
#
#  #####################################
#  def load_best_sampled_rmsd(self, fname):
#    infile = file(fname, "r")
#    line = infile.readline()
#    col_labs = [ s.strip() for s in line.split(",")[1:-1] ]
#    for i in range(len(col_labs)):
#      line = infile.readline()
#      toks = line.split(",")
#      self.sampled_rmsds[toks[0]] = {}
#      for r,t_id in zip(toks[1:], col_labs):
#        # need to be messy since some of the fields could be blank
#        r = r.strip()
#	if(not r == ""):
#          self.sampled_rmsds[toks[0]][t_id] = float(r)
#	else:
#          self.sampled_rmsds[toks[0]][t_id] = 1000.0
#  # the next line is now computed by the printing routine
#  #  self.mean_samp_rmsd = float((infile.readline()).split(",")[1])
#
#  #####################################
#  def load_best_scores(self, fname):
#    infile = file(fname, "r")
#  
#    line = infile.readline()
#    SF_lab = ""
#    self.error_rate = {}
#    while(len(line)):
#      if(line.startswith("SF")):
#        (SF_lab, score_tbl) = self.__read_table__(line, infile)
#        self.scores[SF_lab] = score_tbl
#      elif(line.startswith("Average squared error rate,")):
#        self.error_rate[SF_lab] = float(line.split(",")[1])
#      line = infile.readline()
#
#  #####################################
#  def load_rmsd_of_best_scores(self, fname):
#    infile = file(fname, "r")
#
#    line = infile.readline()
#    SF_lab = ""
#    self.error_rate = {}
#    while(len(line)):
#      if(line.startswith("SF")):
#        (SF_lab, rmsd_tbl) = self.__read_table__(line, infile)
#        self.rmsd_of_scores[SF_lab] = rmsd_tbl
#      line = infile.readline()
#
  #####################################
  def write_rmsd_tbl(self, ofile, col_order=None, table_name="best_sampled", 
                     gray_if_above_tol=True):
    rmsd_tbl = {}
    if(table_name == "best_sampled"):
      rmsd_tbl = self.sampled_rmsds
    elif(table_name == "top_scored"):
      rmsd_tbl = self.rmsd_of_scores
    else:
      print >> stderr, "Unsupported table type (%s), unable to write anything" % (table_name)
      return

    if(len(rmsd_tbl) < 1):
      print >> stderr, "Empty %s rmsd table, unable to write anything" % (table_name)
      return

    N = len(rmsd_tbl)
    print >> ofile, "\\begin{sidewaystable}[p]\\tiny"
    print >> ofile, "\\caption{ empty %s rmsd table caption}" % (table_name)
    print >> ofile, "\\begin{center}"
    print >> ofile, "\\resizebox{\\textwidth }{!}{%"
    print >> ofile, "\\begin{tabular}{l*{%d}{r}}" % (N)
    print >> ofile, "\\toprule%"
    print >> ofile, "Queries& \\multicolumn{%d}{c}{Database Sitemaps}\\\\" % (N)
    print >> ofile, "\\midrule"

    # Print the column headings
    if(col_order == None):
      first_key = rmsd_tbl.keys()[0]
      print >> ofile, " & " + "& ".join(rmsd_tbl[first_key].keys())
      print >> ofile, "\\\\"
    else:
      print >> ofile, " & " + "& ".join([ lab[0:4] for lab in col_order ])
      print >> ofile, "\\\\"
    print >> ofile, "\\otoprule"
      
    # print the data
    if(col_order == None):
      for q_id, row in rmsd_tbl.iteritems():
        print >> ofile, q_id.replace("_", "\_"),
        for t_id, rmsd in row.iteritems():
          print >> ofile, "& %6.2f" % (rmsd),
          if(rmsd > self.rmsd_tol and gray_if_above_tol or
             rmsd < self.rmsd_tol and not gray_if_above_tol):
            print >> ofile, "\\cellcolor[gray]{0.8}",
        print >> ofile, "\\\\"
    else:
      rmsd_sum = 0.0
      rmsd_sum_ss2 = 0.0
      for q_id in col_order:
        print >> ofile, q_id.replace("_", "\_"),
        row = rmsd_tbl[q_id]
        for t_id in col_order:
          print >> ofile, "& %6.2f" % (row[t_id]),

	  # extremely crude measure to toss out outliers
          if(row[t_id] < 100.0):
            rmsd_sum += row[t_id]
            rmsd_sum_ss2 += row[t_id]*row[t_id]

          if(row[t_id] > self.rmsd_tol and gray_if_above_tol or
             row[t_id] < self.rmsd_tol and not gray_if_above_tol):
            print >> ofile, "\\cellcolor[gray]{0.8}",
        print >> ofile, "\\\\"

    N = len(col_order)*len(col_order)
    mu = rmsd_sum / N
    from math import sqrt
    sigma = sqrt(rmsd_sum_ss2 / N - mu*mu)
  
    print >> ofile, "\\bottomrule"
    print >> ofile, "\\end{tabular}}"  # end of resize box 
    print >> ofile, "\\label{no_label::%s_rmsd_tbl}" % (table_name)
    print >> ofile, "\\end{center}"
    print >> ofile, "\\end{sidewaystable}"
    print "%s rmsd stats (mean, std) (%6.2f, %6.2f)" % (table_name, mu, sigma)
  

  #####################################
  def write_scores(self, ofile, col_order=None, SF_labs=None):
    if(len(self.scores) < 1):
      print >> stderr, "Empty tables, unable to write anything"
      return

    if(SF_labs == None):
      SF_labs = self.scores.keys()

    for SF_lab in SF_labs:
      N = len(self.scores[SF_lab])
      print >> ofile, "\\begin{sidewaystable}[p]\\tiny"
      print >> ofile, "\\caption{" + SF_lab + " empty score table caption}"
      print >> ofile, "\\begin{center}"
      print >> ofile, "\\resizebox{\\textwidth }{!}{%"
      print >> ofile, "\\begin{tabular}{l*{%d}{r}}" % (N)
      print >> ofile, "\\toprule%"
      print >> ofile, "Queries& \\multicolumn{%d}{c}{Database Sitemaps}\\\\" % (N)
      print >> ofile, "\\midrule"

      # Print the column headings
      if(col_order == None):
        first_key = self.scores[SF_lab].keys()[0]
        print >> ofile, " & " + "& ".join(self.scores[SF_lab][first_key].keys())
        print >> ofile, "\\\\"
      else:
        print >> ofile, " & " + "& ".join([ lab[0:4] for lab in col_order ])
        print >> ofile, "\\\\"
      print >> ofile, "\\otoprule"
      
      # print the data
      my_scores = self.scores[SF_lab]
      my_rmsds = self.rmsd_of_scores[SF_lab]
      rmsd_sum = 0.0
      rmsd_sum_ss2 = 0.0
      for q_id in col_order:
        print >> ofile, q_id.replace("_", "\_"),
        score_row = my_scores[q_id]
        rmsd_row = my_rmsds[q_id]
        for t_id in col_order:
          (score, rmsd) = (score_row[t_id], rmsd_row[t_id])

	  if(score <= self.score_cutoff):
            print >> ofile, "& %6.2f" % (rmsd),
	    if(rmsd > self.rmsd_tol):
	      print >> ofile, "\\cellcolor[rgb]{1.0, 0.0, 0.0}",
	  else:
	    print >> ofile, "& \\cellcolor[gray]{0.8}",

#          print >> ofile, "& %6.2f" % (score),
#
#          if(score > self.score_cutoff):
#            if(rmsd > self.rmsd_tol):
#              print >> ofile, "\\cellcolor[rgb]{0.0, 1.0, 0.0}",
#            else:
#              print >> ofile, "\\cellcolor[rgb]{1.0, 0.0, 1.0}",
#          elif(rmsd > self.rmsd_tol):
#            print >> ofile, "\\cellcolor[rgb]{1.0, 0.0, 0.0}",
          
          rmsd_sum += rmsd
          rmsd_sum_ss2 += rmsd*rmsd
        print >> ofile, "\\\\"
    
      mu = rmsd_sum / (N*N)
      from math import sqrt
      sigma = sqrt(rmsd_sum_ss2 / (N*N) - mu*mu)

      print >> ofile, "\\bottomrule"
      print >> ofile, "\\end{tabular}}"  # end of resize box 
      print >> ofile, "\\label{no_label::" + SF_lab + "}"
      print >> ofile, "\\end{center}"
      print >> ofile, "\\end{sidewaystable}"
      print "%s best scored rmsd stats (mean, std) (%6.2f, %6.2f)" % (SF_lab, mu, sigma)

