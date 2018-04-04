# Issues to address
# Multiple locations
# Modified residues
# hetatoms:
#   Metals
#   Ligands
#   Waters

import os
import sys
from SimSite3DPy.utils import pdb

################################################################################
class PDBsum_entry:

  def __init__(self, line):
    """
    Line:  line from PDB sum protnames.lst
    """

    toks = line.split()
    self.id = toks[0]
    self.dep_date = toks[2]
    self.exp_method = toks[3]
    if(self.exp_method == "X-ray"): self.resolution = float(toks[4])
    else: self.resolution = -1.0
    self.title = line[31:]

  def add_fasta(self, chainID, fasta_seq):
    self.chainID = chainID
    self.fasta_seq = fasta_seq
################################################################################

################################################################################
class PDBsum:
  """
  Open PDB sum flat files and store the info in a dictionary with pdb codes as
  the keys
  """

  def __init__(self, protnames_fname, fastas_fname):
    self.records = {}

    self.load_protnames(protnames_fname)
    self.load_fastas(fastas_fname)

#>101m:A
#MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRVKHLKTEAEMKASEDLKKHGVTVLTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGNFGADAQGAMNKALELFRKDIAAKYKELGYQG
#    self.

  def load_protnames(self, fname):
    
    infile = open(fname, "r")
    for line in infile: self.records[line[0:4]] = PDBsum_entry(line)

  def load_fastas(self, fname):
    infile = open(fname, "r")

    id_line = infile.next()
    fasta = infile.next()
    while(fasta):
      id_line = infile.next()
      fasta = infile.next()
      print id_line

    


################################################################################

################################################################################
class bindingMOAD:

  def __init__(self, every_fname, delim = "|"):
    every_file = file(every_fname, "r")

    pdb_id = ""
    ligs = {}
    for line in every_file:
      if(line.startswith("\"")): continue
      toks = line.split(",")

      # Ignore enzyme class for now -- at least 1 "bad" line
      if(len(toks[0]) or len(toks) < 3): continue
     
      if(len(toks[2])):
        pdb_id = toks[2]
#        print "\n%s%s" % (toks[2], delim),

      if(len(toks[3]) and toks[4] == "valid"):
        lig_id = toks[3]
        if(not lig_id in ligs): ligs[lig_id] = [ pdb_id ]
        else: ligs[lig_id].append(pdb_id)
#        print "%s%s" % (toks[3], delim),

    self.ligs = [ (key, len(ligs[key]), ligs[key]) for key in ligs ]

  def __greater(self, A, B):
    if(A[1] > B[1]): return -1
    elif(A[1] == B[1]): return 0
    return 1

  def sort_by_count(self):
    self.ligs.sort(cmp=self.__greater)
################################################################################

################################################################################
def best_for_each_ligand(my_MOAD, my_PDB):
  print "Ligand|# of occurances in BindingMOAD|PDB rep.|Nominal Resolution|Name|Compound|"

  for i in range(len(my_MOAD.ligs)):
    lig = my_MOAD.ligs[i]

    # for each pdb code in lig[2] get the first with the best nominal rmsd
    best_pdb = None
    best_res = 99.0
    for id in lig[2]:
      if(id.lower() in my_PDB.records):
        my_rec = my_PDB.records[id.lower()]
        if(my_rec.resolution > 0.0 and my_rec.resolution < best_res):
          best_res = my_rec.resolution
          best_pdb = my_rec

    if(not best_pdb == None):
      print "%s|%d|%s|%.2f|%s|%s|" % \
        (lig[0], lig[1], best_pdb.id.lower(), best_res, best_pdb.name,
         best_pdb.compound)
################################################################################

################################################################################
def structs_with_lig(my_MOAD, my_PDB, lig_id):

  lig_info = ("", "", [])
  for LL in my_MOAD.ligs:
    if(LL[0] == lig_id.upper()):
      lig_info = LL
      break

  if(len(LL[0]) == 0):
    print "Could not find \"%s\" in BindingMOAD" % (lig_id)
    return

#  self.ligs = [ (key, len(ligs[key]), ligs[key]) for key in ligs ]
  print "PDB ID|CC ID|Method|Compound|Authors|Resolution|Name|"
  for pdb_id in lig_info[2]:
    my_rec = my_PDB.records[pdb_id.lower()]
    print "%s|%s|%s|%s|%s|%.2f|%s|" % \
      (pdb_id, lig_id.upper(), my_rec.exp_method, my_rec.compound, 
       my_rec.authors, my_rec.resolution, my_rec.name)
 

################################################################################

################################################################################
def count_occurances(my_MOAD):
  my_MOAD.sort_by_count()
  for LL in my_MOAD.ligs:
    print LL[0], LL[1]
  
################################################################################

################################################################################
if(__name__ == "__main__"):
  from optparse import OptionParser
  import RCSB_PDB_records


  data_dir = "/psa/results/SimSite3D_datasets/data_files"

  cmd_parser = OptionParser()
  cmd_parser.add_option("", "--best_example_for_each", default=False,
                        action="store_true",
                        help="Get the structure the best resolution for each ligand in BindingMOAD")
  cmd_parser.add_option("-l", "--structs_containing_lig", 
                        default="", metavar="3 char HET code",
                        help="Get the structures that contain the given ligand")
  cmd_parser.add_option("-c", "--count_ligs", default=False,
                        action="store_true",
                        help="Get the the number of structures that contain each ligand")
#  if(len(sys.argv) < 3):
#    print >> sys.stderr, "\n\t * All flags and arguments are required *"
#    print >> sys.stderr, "\t   (Except help of course)\n"
#    cmd_parser.print_help()
#    sys.exit(-1)
  (cmd_options, cmd_args) = cmd_parser.parse_args()


  my_MOAD = bindingMOAD(data_dir + "/BindingMOAD/every.csv")
  my_MOAD.sort_by_count()

  entries_fname = data_dir + "/RCSB_PDB/entries.idx"
  fastas_fname = data_dir + "/RCSB_PDB/pdb_seqres.txt"
  my_PDB = RCSB_PDB_records.table(entries_fname, fastas_fname)
 
  if(cmd_options.best_example_for_each):
    best_for_each_ligand(my_MOAD, my_PDB)
  elif(len(cmd_options.structs_containing_lig)):
    structs_with_lig(my_MOAD, my_PDB, cmd_options.structs_containing_lig)
  elif(cmd_options.count_ligs):
    count_occurances(my_MOAD)
    

