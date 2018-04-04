import sys
import os
from xml.dom import minidom
from numpy import array

#####
def my_precision(s, percent=False):
  rv = ""
  if(len(s) > 0):
    if(percent):
    #  rv = "%.f%%" % (100.0 *float(s))
      rv = "%f" % (float(s))
    else:
      rv = "%.2f" % (float(s))
  return rv
#####

#####
def my_float(the_str, empty_val=0):
  if(len(the_str)):
    return float(the_str)
  else:
    return empty_val
#####

#####
class Hit:

  RT_tag_to_idx = { "Rxx" : 0, "Rxy" : 1, "Rxz" : 2, "Tx"  : 3,
                    "Ryx" : 4, "Ryy" : 5, "Ryz" : 6, "Ty"  : 7,
                    "Rzx" : 8, "Rzy" : 9, "Rzz" :10, "Tz"  :11 };

  #####
  def __init__(self, match_node):

    self.q_score = "-1"
    self.p_score = "-1"
    self.z_score = "-1"
    self.rmsd = "-1"
    self.n_align = "-1"
    self.n_gaps = "-1"
    self.seq_id = "-1"
    self.target = None
    self.T_4x4 = [1, 0, 0, 0, \
                  0, 1, 0, 0, \
                  0, 0, 1, 0, \
                  0, 0, 0, 1]

    for node in match_node.childNodes:    
      if(not node.nodeType == node.ELEMENT_NODE):
        continue

      if(node.hasChildNodes()):
        if(node.tagName == "Q-score"):
          self.q_score = node.firstChild.data
        elif(node.tagName == "P-score"):
          self.p_score = node.firstChild.data
        elif(node.tagName == "Z-score"):
          self.z_score = node.firstChild.data
        elif(node.tagName == "RMSD"):
          self.rmsd = node.firstChild.data
        elif(node.tagName == "Nalgn"):
          self.n_align = node.firstChild.data
        elif(node.tagName == "Ngaps"):
          self.n_gaps = node.firstChild.data
        elif(node.tagName == "SeqIdentity"):
          self.seq_id = node.firstChild.data
	elif(node.tagName == "Target"):
	  self.get_target_name(node)
	elif(node.tagName == "RTMatrix"):
	  self.get_rot_trans(node)

  #####
  def get_target_name(self, pnode):
    for node in pnode.childNodes:    
      if(not node.nodeType == node.ELEMENT_NODE):
        continue

      if(node.hasChildNodes() and node.tagName == "name"):
        self.target = node.firstChild.data.split(".")[0]
      elif(node.hasChildNodes() and node.tagName == "chainID"):
        self.chainID = node.firstChild.data.split(".")[0]

  #####
  def get_rot_trans(self, pnode):
    for node in pnode.childNodes:    
      if(not node.nodeType == node.ELEMENT_NODE):
        continue

      if(node.hasChildNodes()):
        self.T_4x4[self.RT_tag_to_idx[node.tagName]] = float(node.firstChild.data)
    
    self.T_4x4 = array(self.T_4x4).reshape((4,4))
    self.R = self.T_4x4[0:3,0:3].T
    self.T = self.T_4x4[0:3,3] 

#####

#####
def getQuery(pnode):
  rv = ""
  for node in pnode.childNodes:
    if(not node.nodeType == node.ELEMENT_NODE):
      continue

    if(node.tagName == "source" and node.hasChildNodes()):
      rv = node.firstChild.data.split(".")[0]
  return rv

#####

#####
def parse_ssm_xml_results(xml_fname):
  xmldoc = minidom.parse(xml_fname)
  
  if(xmldoc.firstChild == None):
    print "Unable to parse the file: " + xml_in
    return None
  
  head_node = xmldoc.childNodes[0]
  if(not head_node.tagName == "SSMResults"):
    print "Expected to find SSMResults as the head node tag; found " + \
      head_node.tagName
    return None
  
  query_node = head_node.getElementsByTagName("Query")[0]
  query_str = query_node.getElementsByTagName("source")[0].firstChild.nodeValue
  hits = {}
  hit_nodes = head_node.getElementsByTagName("Match")
  for node in hit_nodes:
    tmp_hit = Hit(node)
    hits[ tmp_hit.target  + "_" + tmp_hit.chainID] = tmp_hit
    #hits[ tmp_hit.target[:-2] ] = tmp_hit
    #hits.append(Hit(node))

  return query_str[:-9], hits
#####

#####
def generate_CSV_file(struct_ids, seq_id_tbl=False, q_score_tbl=False):
  # process each results file
  hits_tbl = {}
  for q in struct_ids:
    (query_id, hits) = parse_ssm_xml_results(q  + "_results.xml")
  
    if(not query_id == q):
      print "The query id from the results file name does not match the" + \
            "query given in the results file"
      sys.exit(1)
  
    hits_tbl[query_id] = []
    for targ in struct_ids:
      if(not targ in hits):
        hits_tbl[query_id].append("")
      else:
        if(seq_id_tbl):
          hits_tbl[query_id].append(hits[targ].seq_id)
        elif(q_score_tbl):
          hits_tbl[query_id].append(hits[targ].q_score)
  
  # Print results in csv format
  percent = False
  if(seq_id_tbl):
    print "EBI SSM percent of sequence identity,"
    percent = True
  elif(q_score_tbl):
    print "EBI SSM Q Score,"
  print "database binding residues,"
  
  # Put into groups as per groups file and sort based on score to that one
  grps_file = open("groups.txt", "r")
  groups = {}
  ordered_ids = [] 
  for line in grps_file:
    if(line.startswith("#")):
      continue
  
    tmp = line.split()
    codes = dict([(s, 0) for s in tmp])
    all_pairs = dict( [(id, s) for id,s in zip(struct_ids, hits_tbl[tmp[0]]) ] )
    my_pairs = [ (c, my_float(all_pairs[c], empty_val=-1)) for c in codes ]
    my_pairs.sort(cmp = lambda x,y: cmp(x[1], y[1]), reverse=True)
    ordered_ids.extend([ p[0] for p in my_pairs])
  
  # Need to rearrange (sort) the score columns to match the reordered column
  # headings (structure ids)
  orig_idz = dict([(tag,idx) for idx,tag in enumerate(struct_ids)])
  idz = [ orig_idz[tag] for tag in ordered_ids ]
  
  #line = "query," + ",".join([s[:-4] for s in ordered_ids]) + ","
  line = "query," + ",".join(ordered_ids) + ","
  print line
  for id in ordered_ids:
    vals = [ my_precision(hits_tbl[id][i], percent) for i in idz ]
    print id + "," + ",".join(vals) + ","
##########################################################################    

##########################################################################    
if(__name__ == "__main__"):

  # Parse cmdline options
  from optparse import OptionParser
  parser = OptionParser()
  parser.add_option("", "--seq_id", dest="seq_id", action="store_true", 
                    help="Get amount of seq id")
  parser.add_option("", "--Q-score", dest="q_score", action="store_true", 
                    help="Get SSM Q score")
  parser.add_option("", "--xforms", dest="query_id", help="Get SSM transforms")
  (options, args) = parser.parse_args()
  
  if len(args): parser.error("Bare arguments are not supported")
  
  # Get all of the struct_ids
  struct_ids = []
  for f in os.listdir("."):
    if(f.endswith("_results.xml")): struct_ids.append(f[:-12])
  struct_ids.sort()
  
  if(options.seq_id):
    generate_CSV_file(struct_ids, seq_id_tbl=True)
  elif(options.q_score):
    generate_CSV_file(struct_ids, q_score_tbl=True)
  elif(options.query_id != None):
    (query_id, hits) = parse_ssm_xml_results(options.query_id  + "_results.xml")
    for targ_id in hits.keys():
      print "%s|%s|" % (targ_id, " ".join([str(x) for x in hits[targ_id].T_4x4]))
  else:
    parser.error("One of the options must be used")
