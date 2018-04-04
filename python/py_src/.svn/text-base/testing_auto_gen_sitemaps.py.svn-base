#
# h4x0r version to support my multiple datasets
#

import os
import sys
from SimSite3DPy import * 
from numpy import *
import cPickle

print "-" * 80
print sys.argv[0]
dirs = parameters.parameters_base()
dirs.version()
print "-" * 80
print
dirs.get_params()
if(dirs.fail): sys.exit(1)

base_data_dir = os.getcwd()

#datasets = [ "adenines", "pterins", "inhibited_gst"]
datasets = [ "adenines", "pterins"]

for dataset in datasets:

  dirs.dbase_prots = "%s/%s/proteins" % (base_data_dir, dataset)
  dirs.dbase_ligs = "%s/%s/ligands" % (base_data_dir, dataset)
  if(dataset == "enolases"):
    dirs.dbase_sites = "%s/%s/dbase_metals" % (base_data_dir, dataset)
  else: dirs.dbase_sites = "%s/%s/dbase" % (base_data_dir, dataset)

  if(not os.path.exists(dirs.dbase_sites)): os.mkdir(dirs.dbase_sites)

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
    prot_fname = dirs.dbase_prots + "/" + prot
    print "*****\nGenerating the sitemap of\n  " + prot_fname
    print "with respect to the ligand volume of\n  " + lig_fname
#    the_cmd = dirs.install_dir + "/bin/gen_points --no_normalization " 
    the_cmd = "/home/vanvoor4/code/SimSite3D_surfaces/src/sitemap/gen_points "
    the_cmd += "--no_normalization "

#2ak1 & 1ncw have a discrepancy where BEZ is bound - take center of pocket
#(estimated by 2ak1 BEZ H2) with 5(A) radius sphere.
    if(prot.startswith("2ak1") or prot.startswith("1ncw")):
      center = [16.8088,  -7.2212,  36.8497]
      radius = 5.0
      the_cmd += " --sphere \"%f %f %f %f\" " % \
        (center[0], center[1], center[2], radius)
      # don't mess with metal_site_vol since we don't have metals in these 2
      site_vol = utils.volumes.ball(center, radius)
      print "generating sphere -- special case for 2 binding sites"
    else: 
      lig = utils.mol2.molecule(lig_fname)
      centers = array([ a.position for a in lig.atoms ])
      radii = array([ 2.5 for a in lig.atoms ])
      site_vol = utils.volumes.UnionOfBalls(centers, radii)
      the_cmd += " --prune_to_lig " 
      the_cmd += "-l " + lig_fname 

    the_cmd += " -p " + prot_fname
    the_cmd += " --include_metals --allow_small_site_maps "
    the_cmd += " --msms_surf " + the_site_csv
    print the_cmd
    status = os.system(the_cmd)
    if(status or not os.path.isfile(the_site_csv)):
      print >> sys.stderr, \
"""  Warning!
Could not create the sitemap for the protein\n
  %s 
and ligand
  %s
""" % (prot_fname, lig_fname)
      if(status): print "The return value of gen_points was %d" % (status)
  
    # Create the hbond volumes
    print "Determining the hbond cap volumes"
    hbond_cap_vols = sitemaps.hbond_volumes(prot_fname, site_vol)

    # Create the hbond surface caps
    print "Determining the hbond cap surfaces"
    the_prot = utils.pdb.residues(prot_fname)
    hbond_cap_surfs = sitemaps.hbond_surf_caps(the_prot, site_vol,
                                               hbond_cap_vols)

    # Create the metal volumes
    lig = utils.mol2.molecule(lig_fname)
    centers = array([ a.position for a in lig.atoms ])
    radii = array([ 3.0 + 2.5 for a in lig.atoms ])
    metal_site_vol = utils.volumes.UnionOfBalls(centers, radii)
    metal_vols = sitemaps.metal_volumes(prot_fname, metal_site_vol)

    # Create the metal surfs
    metal_surfs = sitemaps.metal_surfs(the_prot, metal_site_vol, metal_vols)

    # Write out the volumes & surface caps
    caps_fname = "%s/%s.pkl" % (dirs.dbase_sites, prot[:-6])
    caps_fout = open(caps_fname, "w+b")

    for s in metal_surfs.spheres:
      print "num circles:", len(s.circles)


    cPickle.dump((hbond_cap_vols, hbond_cap_surfs, metal_vols, metal_surfs),
                 caps_fout)
    print "dumping to:", caps_fout
                 #caps_fout, cPickle.HIGHEST_PROTOCOL)
    caps_fout.close()

# testing this stupid thing
    caps_in = open(caps_fname, "rb")
    (a,b,c,d) = cPickle.load(caps_in)
    for s in d.spheres:
      print "from file, num circles:", len(s.circles)



  
  # check if any ligands are left
  if(len(ligs)):
    print >> sys.stderr, "\nThe following ligands in\n  " + dirs.dbase_ligs
    print >> sys.stderr, "did not have associated proteins in\n  " \
      + dirs.dbase_prots
    for l in ligs:
      print l + "_l.mol2"
