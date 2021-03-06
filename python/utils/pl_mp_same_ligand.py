#!/usr/bin/env python
import sys
import os
import gzip
import re
################################################################################

################################################################################
def sort_lig_info(lig_tuples, tbl_file):
  max_len = 0
  labs = None
  acts = []
  for idx in range(len(lig_info)):
    if(len(lig_info[idx][5]) > max_len): 
      max_len = len(lig_info[idx][5]) 
      labs = lig_info[idx][4]

    tmp = {}
    for (atom, hit) in zip(lig_info[idx][4], lig_info[idx][5]): tmp[atom] = hit 
    acts.append(tmp)
  
  print >> tbl_file, \
"""PDB id|Lig chain|lig res num|Lig iCode|-- Ligand matchprint --
||||%s|"""  % ("|".join(labs))
  for (t, act) in zip(lig_info, acts):
    toks = []
    
    for lab in labs:
      if(lab in act): 
        toks.append(act[lab])
      else: toks.append(" ")

    print >> tbl_file, "%s|%s|%d|%s|%s|" % \
      (t[0], t[1], t[2], t[3], "|".join(toks))

  for t in lig_info: os.unlink(t[6])
################################################################################

################################################################################

if(__name__ == "__main__"):

  import tempfile
  import cPickle
  from SimSite3DPy import *
  from SimSite3DPy.utils import add_RCSB_hydrogens, misc_utils
  import sys
#  sys.path.append("/psa/results/SimSite3D_datasets/scripts")
  import RCSB_PDB_records
  from optparse import OptionParser

  cmd_parser = OptionParser()
  cmd_parser.add_option("", "--lig_code", help="3 char lig code",
                        metavar="PDB_HET_CODE")
  cmd_parser.add_option("", "--pdb_ids", help="csv file holding pdb ids",
                        metavar="FILE.csv")
  cmd_parser.add_option("", "--hits_tbl", help="file name for hits table",
                        metavar="<fname>.csv") 
  if(len(sys.argv) < 7):
    print >> sys.stderr, "\n\t * All flags and arguments are required *"
    print >> sys.stderr, "\t   (Except help of course)\n"
    cmd_parser.print_help()
    sys.exit(-1)
  (cmd_options, cmd_args) = cmd_parser.parse_args()
 
  pdb_ftp_url = "ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb"

  infile = open(cmd_options.pdb_ids)
  pdb_ids = {}
  for line in infile:
    if(line.startswith("#") or line.startswith("PDB ID|")): continue
    toks = line.split("|")
    pdb_ids[toks[0].rstrip("\n").strip()] = True

  tbl_file = open(cmd_options.hits_tbl, "w+")

  # We must determine which chain(s) interacts with the chosen ligand 
  # Idea: use prot_lig_score
  # This requires: 
  # 1) extracting the ligand(s) from the PDB file 
  # 2) adding hydrogens & charges to the ligand(s)
  # 3) Scoring the ligand versus the protein
  # 4) Parsing the output and getting the percentage of interactions
  #    Use re to find the ( ) parts
  chain_re = re.compile("\([A-Za-z0-9\ ]\)") 

#  print >> tbl_file, "PDB id|Lig chain|lig res num|Lig iCode|-- Ligand matchprint --"
  print "PDB id|Lig chain|Lig res num|Lig iCode|Prot chains interacting with ligand|"
  lig_info = []
  for pdb_id in pdb_ids:
    (prot_fname, prot) = misc_utils.load_pdb_file(pdb_id)

    if(prot.fail): 
      print "ERROR:  could not load the pdb file for", pdb_id
      continue

    for hetgrp in prot.hetgroups:
      if(not hetgrp.hetID == cmd_options.lig_code): continue
      (mol2_fd, mol2_fname) = tempfile.mkstemp(suffix=".mol2", 
                                               prefix=cmd_options.lig_code)
      os.close(mol2_fd)
      add_RCSB_hydrogens.run(hetgrp, mol2_fname, use_lig_expo=True)

      # score it
      my_lig = utils.mol2.molecule(mol2_fname)
      lig_obj = utils._mol2File.mol2(mol2_fname)
      prot_obj = utils.PDBStructure(prot_fname)
      score_info = score.prot_lig_score(prot_obj, lig_obj)
      chainIDs = {}

      hvy_atoms = []
      hit = []
      for lig_atom, str in zip(my_lig.atoms, score_info.lig_act_strings()):
        if(lig_atom.name == "H"): continue

#        #print lig_atom.name_str, len(str) 
        hvy_atoms.append(lig_atom.name_str) 
        hit.append(str)
        #os.unlink(mol2_fname)
      lig_info.append((pdb_id, hetgrp.chainID, hetgrp.seqNum, hetgrp.iCode,
                       hvy_atoms, hit, mol2_fname))
  sort_lig_info(lig_info, tbl_file)
