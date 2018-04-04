import os
import sys
import gzip
from ASCbasePy.utils import pdb, add_RCSB_hydrogens
from ASCbasePy.utils.query_BindingMOAD import *


if(__name__ == "__main__"):
  from optparse import OptionParser

  # At the present these HET groups have missing model and ideal coordinates
  # but the CIF file notes that at least one of the two sets of coordinates
  # do NOT have any missing atoms
  # OXY is molecular oxygen and does not have the typical _loop over a bond
  # section since there is only one double bond
  BAD_HETS = ["GCR", "OXY"]

  # The directory holding RCSB and BindingMOAD directories 
  data_dir = "/psa/results/SimSite3D_datasets/data_files"
  RCSB_ftp_dir = "ftp://ftp.wwpdb.org/pub/pdb/data/biounit/coordinates/all"
  
  cmd_parser = OptionParser()
  cmd_parser.add_option("", "--site_directory", 
                        help="Directory to store the binding sites")
  cmd_parser.add_option("", "--use_full_BindingMOAD", default=False,
                        action="store_true",
                        help="Use the full BindingMOAD instead of the non-redundant set")
  cmd_parser.add_option("", "--overwrite_existing", default=False,
                        action="store_true",
                        help="Overwrite all binding sites in the chosen directory")

  (cmd_options, cmd_args) = cmd_parser.parse_args()
  if(cmd_options.site_directory == None or \
     len(cmd_options.site_directory) == 0):
    print >> sys.stderr, "\n\t--site_directory and its argument is required\n";
    print >> sys.stderr, "Cannot continue!\n"
    sys.exit(-1)

  # create the site directory if needed
  if(not os.path.exists(cmd_options.site_directory)):
    os.mkdir(cmd_options.site_directory)
  sites_dir = cmd_options.site_directory

  if(cmd_options.use_full_BindingMOAD):
    my_MOAD = bindingMOAD(data_dir + "/BindingMOAD/every.csv")
  else:
    my_MOAD = bindingMOAD(data_dir + "/BindingMOAD/nr.csv")

  # get the existing binding sites
  existing_sites = {}
  if(not cmd_options.overwrite_existing):
    junk = os.listdir(cmd_options.site_directory)
    for j in junk:
      if(j.endswith("_s.csv")):
        existing_sites[ j.split("/")[-1][:-6] ] = \
          cmd_options.site_directory + "/" + j
  
  # For each valid PDB ID and HET ID, create the binding site
  for i in range(len(my_MOAD.ligs)):
    (HET_code, num_ligs, pdb_ids) = my_MOAD.ligs[i]
    if(HET_code in BAD_HETS): 
      print "NOTE:", HET_code, "has been flagged as \"BAD\"\n"
      continue

    for pdb_id in pdb_ids:
      print pdb_id, HET_code

      if(len(HET_code) > 3): 
        print "Could not process: ", pdb_id, HET_code
        print "length of HET code was > 3\n"

      # typically if the HET_code has fewer than 3 characters it is
      # filled from the front by spaces
      while(len(HET_code) < 3): HET_code = " " + HET_code

      # struct id for this protein-ligand pair
      if(len(HET_code) == 3):
        struct_id = "%s_%s" % (pdb_id, HET_code.replace(" ", "_"))
      else: struct_id = "%s_%s" % (pdb_id, HET_code.replace(" ", "-"))

      # If a sitemap exists in the database, skip it 
      if(struct_id in existing_sites):
        print "The sitemap for " + struct_id + " already exists"
        continue

      print "\nProcessing structure", struct_id

      # Get the gzipped biounit from the RCSB PDB
      biounit = "%s/%s.pdb1.gz" % (RCSB_ftp_dir, pdb_id.lower())
      os.system("wget " + biounit)
      biounit = "%s.pdb1.gz" % (pdb_id.lower())
      if(not os.path.exists(biounit)): 
#        print >> sys.stderr, "Unable to download ", \
        print "Unable to download %s/%s.pdb1.gz" % \
          (RCSB_ftp_dir, pdb_id.lower())
        continue

      # unzip the biounit since SimSite3D isn't designed to open gzipped files
      # force unzipping here since we are getting bogged down with 
      # user input when redirecting to a file
      os.system("gunzip -f " + biounit)
      biounit = biounit[:-3]
    
      # Load the prot and get the first het group that matches HET_code and
      # whose heavy atoms exactly match those in the corresponding cif file
      # (from Ligand Expo)
      prot = pdb.residues(biounit)
      rv = False

      # If we have 1 het code then just get the first occurence in the
      # first model in the biounit file
      if(len(HET_code) == 3):
        for hetgrp in prot.hetgroups:
          if(hetgrp.hetID == HET_code):
            rv = add_RCSB_hydrogens.run(hetgrp, struct_id + "_l.mol2")
            if(rv): break 
      
      # If we have more than 1 entry in the het code, check if there 
      # is a corresponding chain via the SEQRES entries
      else:
        for seq_chainID, seq in prot.seqres.seqres.iteritems():
          if(not len(seq) == len(HET_code)): continue
          mismatch = False
          for A,B in zip(seq, HET_code):
            if(not A == B):
              mismatch = True
              break
          if(not mismatch):
            print "DID FIND found", struct_id
            # this is not entire difficult to handle, but will take too much
            # time to implement at this poin
            pass
            
      if(rv == False): 
        print "\nNOTE: could not load the hetgrp for", struct_id, "\n"
        continue
  

      # run gen_sitemap from SimSite3D
      # Generate the sitemap without normalizing it with respect to the
      # 140 sitemaps since that step is currently too slow to do for the 
      # database.
      the_site_csv = "%s/%s_s.csv" % (sites_dir, struct_id)
      prot_fname = "%s/%s.pdb1" % (".", pdb_id.lower())
      lig_fname = "%s/%s_l.mol2" % (".", struct_id)
      print "*****\nGenerating the sitemap of\n  " + prot_fname
      print "with respect to the ligand volume of\n  " + lig_fname + "\n"
#    the_cmd = dirs.install_dir + "/bin/gen_points --no_normalization " 
      the_cmd = "/home/vanvoor4/code/SimSite3D_surfaces/src/sitemap/gen_points "
      the_cmd += "--no_normalization "
      the_cmd += " -p " + prot_fname
      the_cmd += " -l " + lig_fname 
#      the_cmd += " --probe_radius 1.39 "
#      print "NOTE: using a slightly smaller probe radius to finesse MSMS"
      the_cmd += " --include_metals --allow_small_site_maps "
      the_cmd += " --msms_surf " + the_site_csv
#      print the_cmd
#      continue
      status = os.system(the_cmd)
      if(status or not os.path.isfile(the_site_csv)):
#        print >> sys.stderr, \
        print \
"""  Warning!
Could not create the sitemap for the protein\n
  %s 
and ligand
  %s
""" % (prot_fname, lig_fname)
        if(status): print "The return value of gen_points was %d" % (status)
      os.unlink(biounit)  
      print "**************************************************************\n"
