import os
import sys
from ASCbasePy import * 
import cPickle

print "-" * 80
print sys.argv[0]
dirs = parameters.parameters_base()
dirs.version()
print "-" * 80
print
dirs.get_params()
if(dirs.fail): sys.exit(1)

print "  Proteins Directory  : " + dirs.dbase_prots
print "  Ligands Directory   : " + dirs.dbase_ligs
print "  Temporary Directory : " + dirs.scratch_dir
print "  Database Directory  : " + dirs.dbase_sites
print "\nGenerating sitemaps:"

# Get the ligands
ligs = {}
junk = os.listdir(dirs.dbase_ligs)
for j in junk:
  if(j.endswith("_l.mol2")):
    ligs[ j.split("/")[-1][:-7] ] = dirs.dbase_ligs + "/" + j

# Get the sitemaps in the database dir
existing_sites = {}
junk = os.listdir(dirs.dbase_sites)
for j in junk:
  if(j.endswith("_s.csv")):
    existing_sites[ j.split("/")[-1][:-6] ] = dirs.dbase_sites + "/" + j


# Get the proteins and run through them
prots = os.listdir(dirs.dbase_prots)
for prot in prots:
  if(not prot.endswith("_p.pdb")):
    continue

  # Cannot use pop at present as we need to support older versions of python
  # lig_fname = ligs.pop(prot[:-6], None) 
  # if(lig_fname == None):
  if(not prot[:-6] in ligs): 
    print >> sys.stderr, "The protein file " + prot + \
      " does not have an associated ligand file"
    print >> sys.stderr, "  Skipping . . .\n"
    continue

  lig_fname = ligs[prot[:-6]]
  del ligs[prot[:-6]]

  # If a sitemap exists in the database, skip it 
  if(prot[:-6] in existing_sites):
    print "The sitemap for " + prot[:-6] + " already exists"
    continue

  # Generate the sitemap without normalizing it with respect to the
  # 140 sitemaps since that step is currently too slow to do for the 
  # database.
  the_site_csv = dirs.dbase_sites + "/" + prot[:-6] + "_s.csv"
  the_prot = dirs.dbase_prots + "/" + prot
  print "*****\nGenerating the sitemap of\n  " + the_prot
  print "with respect to the ligand volume of\n  " + lig_fname
  the_cmd = dirs.install_dir + "/bin/gen_points --no_normalization " + \
    "-l " + lig_fname + " -p " + the_prot + " " + the_site_csv
  status = os.system(the_cmd)
  if(status or not os.path.isfile(the_site_csv)):
    print >> sys.stderr, \
"""  Warning!
Could not create the sitemap for the protein\n
  %s 
and ligand
  %s
""" % (the_prot, lig_fname)
    if(status): print "The return value of gen_points was %d" % (status)
  
 
  # Create the hbond volume file
  site_vol = utils.volumes.box(lig_fname=lig_fname) 
  hbond_caps = sitemaps.hbond_volumes(the_prot, site_vol)
  caps_fname = "%s/%s.pkl" % (dirs.dbase_sites, prot[:-6])
  caps_fout = file(caps_fname, "w+b")
  cPickle.dump(hbond_caps, caps_fout, cPickle.HIGHEST_PROTOCOL)
  caps_fout.close()

# check if any ligands are left
if(len(ligs)):
  print >> sys.stderr, "\nThe following ligands in\n  " + dirs.dbase_ligs
  print >> sys.stderr, "did not have associated proteins in\n  " \
    + dirs.dbase_prots
  for l in ligs:
    print l + "_l.mol2"

