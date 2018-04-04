#!/usr/bin/python
#
# $Source: /psa/share/repository/pfizer_proj/python/auto_gen_sitemaps.py,v $
# $Revision: 1.7 $
# $Author: vanvoor4 $
# $Date: 2008-01-07 15:32:36 $
#
# $Log: not supported by cvs2svn $
# Revision 1.6  2007/11/01 16:44:54  vanvoor4
# Changes to command line flags and names of environment variables
#
# Revision 1.5  2007/10/04 15:01:34  vanvoor4
# Replaced the dictionary pop functionality (requires >=python 2.3) with
# existing python ideas
#
# Revision 1.4  2007/10/02 20:14:07  vanvoor4
# Updated print statements to be "more" correct
#
# Revision 1.3  2007/10/02 19:04:32  vanvoor4
# Reverted to the os.system for the gen_points commmand so that
# the status is caught here
#
# Revision 1.2  2007/10/02 18:42:47  vanvoor4
# Modified to use SimSite3D.utils functions
#
# Revision 1.1  2007/09/26 14:34:50  vanvoor4
# Initial checkin
#
#
#


import os
import sys

print "-" * 80
print sys.argv[0], "(SimSite3D Software 3.3)"
print """Copyright (C) 2006-2009, Michigan State University (MSU) Board of Trustees.
All rights reserved.

Written by Jeffrey R. Van Voorst and Leslie A. Kuhn
"""
print "-" * 80
print

proj_dir = os.getenv("ASCBASE_INSTALL_DIR")
if(proj_dir == None):
  print >> sys.stderr, \
    "\nUnable to get the environment variable $ASCBASE_INSTALL_DIR"
  print >> sys.stderr, "Please set $ASCBASE_INSTALL_DIR to a valid path"
  sys.exit(1)

sys.path.append(proj_dir + "/python")
import utils

dirs = utils.directories()
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
  # the_lig = ligs.pop(prot[:-6], None) 
  # if(the_lig == None):
  if(not prot[:-6] in ligs): 
    print >> sys.stderr, "The protein file " + prot + \
      " does not have an associated ligand file"
    print >> sys.stderr, "  Skipping . . .\n"
    continue

  the_lig = ligs[prot[:-6]]
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
  print "with respect to the ligand volume of\n  " + the_lig
  the_cmd = dirs.install_dir + "/bin/gen_points --no_normalization " + \
    "-l " + the_lig + " -p " + the_prot + " " + the_site_csv
  status = os.system(the_cmd)
  if(status or not os.path.isfile(the_site_csv)):
    print >> sys.stderr, \
"""  Warning!
Could not create the sitemap for the protein\n
  %s 
and ligand
  %s
""" % (the_prot, the_lig)
    if(status): print "The return value of gen_points was %d" % (status)

# check if any ligands are left
if(len(ligs)):
  print >> sys.stderr, "\nThe following ligands in\n  " + dirs.dbase_ligs
  print >> sys.stderr, "did not have associated proteins in\n  " \
    + dirs.dbase_prots
  for l in ligs:
    print l + "_l.mol2"

