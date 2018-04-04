
def generate_sitemaps(prots_dir, ligs_dir, dbase_dir, msms_surf=False,
                      include_metals=False, normalize_sitemaps=False, 
                      scratch_dir=""):
  """
  Generate site maps for the axis aligned ligand bounding box for each
  ligand & protein pair in prots_dir and ligs_dir.  The proteins and ligands
  are expected to follow the XXXXXX_p.pdb & XXXXXX_l.mol2 naming convention
  (respectively).  Here the XXXXXX is typically an informative structure
  identifer.
  """

  import os
  from sys import stderr
  from SimSite3DPy.sitemaps._parameters import parameters
  from SimSite3DPy.sitemaps.Sitemap import sitemap
  #import cPickle

  print "\nGenerating sitemaps:"
  print "  Proteins Directory  : " + prots_dir
  print "  Ligands Directory   : " + ligs_dir
  print "  Database Directory  : " + dbase_dir
  print "  Temporary Directory : " + scratch_dir + "\n"

   # Get the ligands
  ligs = {}
  junk = os.listdir(ligs_dir)
  for j in junk:
    if(j.endswith("_l.mol2")):
      ligs[ j.split("/")[-1][:-7] ] = ligs_dir + "/" + j

  # Get the sitemaps in the database dir
  existing_sites = {}
  junk = os.listdir(dbase_dir)
  for j in junk:
    if(j.endswith("_s.csv")):
      existing_sites[ j.split("/")[-1][:-6] ] = dbase_dir + "/" + j

  # Get the proteins and run through them
  prots = os.listdir(prots_dir)
  for prot in prots:
    if(not prot.endswith("_p.pdb")):
      continue

    # Cannot use pop at present as we need to support older versions of python
    # lig_fname = ligs.pop(prot[:-6], None) 
    # if(lig_fname == None):
    if(not prot[:-6] in ligs):
      print >> stderr, \
"""The protein file %s does not have an associated ligand file
  Skipping . . .
""" % (prot)
      continue

    # Remove the structure id from the list of ligands so that at the end
    # we know which ligands do not have a corresponding protein
    lig_fname = ligs[prot[:-6]]
    del ligs[prot[:-6]]
 
    # If a sitemap exists in the database, skip it 
    if(prot[:-6] in existing_sites):
      print "The sitemap for " + prot[:-6] + " already exists"
      continue

    # Generate the sitemap
    the_site_csv = dbase_dir + "/" + prot[:-6] + "_s.csv"
    the_prot = prots_dir + "/" + prot

    print \
"""*****\nGenerating the sitemap of
  %s
with respect to the ligand volume of
  %s
""" % (the_prot, lig_fname)

    cc_params = parameters()
    cc_params.pts_fname = the_site_csv
    cc_params.prot_fname = the_prot
    cc_params.lig_fname = lig_fname
    cc_params.normalize = normalize_sitemaps
    cc_params.include_metals = include_metals
    cc_params.call_msms = msms_surf
    if(len(scratch_dir)): cc_params.scratch_dir = scratch_dir
    my_site = sitemap(cc_params)
    my_site.write_files(prot[:-6], dbase_dir)

    # Somethign is amiss with the python verison of Sitemap.fail() -- 
    # it always returns true
    # if(my_site.fail() or not os.path.isfile(the_site_csv)):
    if(not os.path.isfile(the_site_csv)):
      print >> stderr, \
"""  Warning!
Could not create the sitemap for the protein\n
  %s 
and ligand
  %s
""" % (the_prot, lig_fname)

#    # Create the hbond volume file
#    site_vol = utils.volumes.box(lig_fname = lig_fname)
#    hbond_caps = sitemaps.hbond_volumes(the_prot, site_vol)
#    caps_fname = "%s/%s.pkl" % (dbase_dir, prot[:-6])
#    caps_fout = file(caps_fname, "w+b")
#    cPickle.dump(hbond_caps, caps_fout, cPickle.HIGHEST_PROTOCOL)
#    caps_fout.close()

  # check if any ligands are left
  if(len(ligs)):
    print >> stderr, """
The following ligands in
  %s
did not have associated proteins in
  %s""" % (ligs_dir, prots_dir)
    for l in ligs:
      print l + "_l.mol2"

################################################################################

################################################################################
