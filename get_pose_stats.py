#!/usr/bin/python

import os
import sys
import numpy
import tempfile

def RMSD(A,B):
  tmp = A - B
 
  return numpy.sqrt( numpy.sum(tmp*tmp) / A.shape[0])

def setup_directories(base_dir, mode):
  if(not os.path.exists(base_dir + "/moved_ligands")):
    os.makedirs(base_dir + "/moved_ligands", 0770)
  if(not os.path.exists(base_dir + "/ligand_fragments")):
    os.makedirs(base_dir + "/ligand_fragments", 0770)

def copy_search_params(params_py):
  from ASCbasePy import search
  params_cc = search.parameters()

#(params.num_rand_aligns > 0):

  # base parameters
  params_cc.dbase_sites = params_py.dbase_sites
  params_cc.dbase_ligs = params_py.dbase_ligs
  params_cc.dbase_prots = params_py.dbase_prots 
  params_cc.diverse_sites = params_py.diverse_sites
  params_cc.diverse_ligs = params_py.diverse_ligs
  params_cc.proj_output = params_py.proj_output
  params_cc.scratch_dir = params_py.scratch_dir

  # search parameters
  params_cc.model_file_name = \
    params_py.query_sitemap_dir + "/" + params_py.query_site_id
  params_cc.dbase_file_name = \
    params_py.dbase_sitemap_dir + "/" + params_py.dbase_site_id

  params_cc.normalize = params_py.normalize
  params_cc.num_scores_to_keep = params_py.num_scores_to_keep
  params_cc.score_cutoff = params_py.score_cutoff
  params_cc.min_num_atoms = params_py.min_num_atoms
  params_cc.dmetol = params_py.dmetol
  params_cc.lsetol = params_py.lsetol
  params_cc.write_ligands = params_py.write_ligands
  params_cc.allow_hphob_triangles = params_py.hydrophobic_query
  return params_cc

###############################################################################
proj_dir = os.getenv("ASCBASE_INSTALL_DIR")
if(proj_dir == None):
  print >> sys.stderr, \
    "\nUnable to get the environment variable $ASCBASE_INSTALL_DIR"
  print >> sys.stderr, "Please set $ASCBASE_INSTALL_DIR to a valid path"
  sys.exit(1)

sys.path.append(proj_dir + "/python")
from ASCbasePy import *

# Handle command line parameters
params = parameters.search_parameters()
rv = params.cmdline_options(sys.argv)
if(rv == False or params.fail):
  sys.exit(1)

# Check for the necessary directories and create them if they do not exist
print "Results directory is:" , params.proj_output
print "Results file is", params.ofname
setup_directories(params.proj_output, 0770)
params_cc = copy_search_params(params)

timer = None
if(params.time_process): 
  timer = utils.system_timers()
  timer.start();

# Set up the alignment method -- possibly based on the model site
model = search.ModelSitemap(params.query_sitemap_dir, params.query_site_id,
                            params_cc, params.normalize)
model.load_hbond_caps()
hbond_caps_score = search.hbond_volume_score(model, params)
if(params.num_rand_aligns == 0 and params.align_to_query):
  align_method = search.MatchTriangles(model, params.dmetol, params.lsetol,
                                       params.hydrophobic_query)
elif(not params.align_to_query):
  align_method = search.IdentityAlignment()
elif(params.num_rand_aligns > 0):
  align_method = search.RandomAlignments()
else:
  print "Unknown alignment method"
  sys.exit(1)

# Take care of the output -- sent to a file
out_file_thingie = utils.stl_ofstream()
rv = out_file_thingie.open(params.ofname)
if(not rv):
  print "Results file could not be opened"
  sys.exit(1)
params_cc.report(out_file_thingie)
score_method = search.WeightedSumsScore(model, params_cc)
score_method.write_score_header(out_file_thingie)

aligns = search.vector_less__rigid_align_t__greater_()
# Pairwise comparison
if(len(params.dbase_site_id)):
  dbase_site = sitemaps.sitemap(params.dbase_sitemap_dir, params.dbase_site_id,
                                params_cc, params.normalize)
  dbase_site.load_hbond_caps()
  if(params.num_rand_aligns == 0):
    align_method.align(dbase_site, aligns)
    score_method.score_alignments(aligns, dbase_site, out_file_thingie)
  else:
    # Add the Identity alignment as well
    id_align_mthd = search.IdentityAlignment()
    id_align_mthd.align(dbase_site, aligns)

    lig = utils.mol2.molecule("%s/2qx0_ph2_l.mol2" % (params.dbase_ligs))
    print lig
    orig_lig_pos = []
    for atom in lig.atoms:
      orig_lig_pos.append(atom.position)
    orig_lig_pos = numpy.array(orig_lig_pos)   

    out_file_thingie.close()
    results_file = file(params.ofname, "a")
    print "%s/%s_s.csv" % (params.dbase_sitemap_dir, params.dbase_site_id)
    site_pts = utils.pdb.hetgroup("%s/%s_s.pdb" % (params.dbase_sitemap_dir, 
                                                   params.dbase_site_id))
    site_pos = []
    for a in site_pts.atoms: site_pos.append(a.position)
    site_pos = numpy.array(site_pos)
    centroid = numpy.mean(site_pos, 0)
    orient_params = \
      align_method.align(dbase_site, aligns, centroid, N=params.num_rand_aligns,
                         q_width=0.1, half_side_len=0.5)
    orient_params.extend(
      align_method.align(dbase_site, aligns, centroid, N=params.num_rand_aligns,
                         q_width=0.05, half_side_len=0.25)
      )
    orient_params.extend(
      align_method.align(dbase_site, aligns, centroid, N=params.num_rand_aligns,
                         q_width=0.025, half_side_len=0.125)
      )
    orient_params.extend(
      align_method.align(dbase_site, aligns, centroid, N=params.num_rand_aligns,
                         q_width=0.01, half_side_len=0.05)
      )

    orig_site_pos = site_pos.copy()

    (tmp_fd, tmp_fname) = tempfile.mkstemp(".out")
    os.close(tmp_fd)  
    tmp_out = utils.stl_ofstream()
    rv = tmp_out.open(tmp_fname)
    if(not rv): print "Temporary results file could not be opened"
    score_method.score_alignments(aligns, dbase_site, tmp_out)
    tmp_out.close()

    print >> results_file, "# *** NOTE: R is to be premultiplied by the coordinates ***"
    print >> results_file, "# dbase_id | raw SimSite3D score | R | T | matchprint | | quat params | effective translation | hb caps sum | hb caps match print | ligand RMSD | site RMSD|"
    hbond_caps_score.score_alignments(aligns, dbase_site)
    for orient_stuff, align in zip(orient_params, aligns):
      (q_params, trans_shift) = orient_stuff
      R = numpy.array([ align.R[i] for i in range(9)]).reshape((3,3))
      T = numpy.array([ align.T[i] for i in range(3)])
      toks = ["%s|%f" % (params.dbase_site_id, align.score) ]
      toks.append(" ".join(["%f" % (align.R[i]) for i in range(9)]))
      toks.append(" ".join(["%f" % (align.T[i]) for i in range(3)]))
      match_print = ""
      for b in align.match_print:
        if(b): match_print += "1"
        else: match_print += "0"
      toks.append(match_print)
      toks.append("")
      toks.append("%f %f %f" % (q_params[0], q_params[1], q_params[2]))
      toks.append("%f %f %f" % (trans_shift[0], trans_shift[1], trans_shift[2]))
      toks.append("%f" % (align.hb_caps_score))
      toks.append(align.hb_caps_match_print)

      # compute ligand and sitemap rmsd
      lig_pos = numpy.dot(orig_lig_pos, R) 
      lig_pos += numpy.tile(T, (orig_lig_pos.shape[0], 1))
      toks.append("%f" % RMSD(lig_pos, orig_lig_pos))

      site_pos = numpy.dot(orig_site_pos, R) 
      site_pos += numpy.tile(T, (orig_site_pos.shape[0], 1))
      toks.append("%f" % RMSD(site_pos, orig_site_pos))

      print >> results_file, "|".join(toks) + "|"
    os.unlink(tmp_fname)
      
#  score_method.score_alignments(aligns, dbase_site, out_file_thingie)
else:
  # do full db scan 
  if(not len(params.db_index_fname)):
    for site in os.listdir(params.dbase_sites):
      if(not site.endswith("_s.csv")):
        continue

      dbase_site = sitemaps.sitemap(params.dbase_sites, site[:-6], params_cc, 
                                    params.normalize)
      dbase_site.load_hbond_caps()
      aligns.clear()
      align_method.align(dbase_site, aligns)
      score_method.score_alignments(aligns, dbase_site, out_file_thingie)

  # look in the db_index_fname to start with dbstart and finish with
  # dbstop
  else:
    db_index_file = file(params.db_index_fname, "r")
    site_ids = []
    for line in db_index_file:
      if(line.startswith("#")): continue
      site_ids.append(line.strip())
  
    for i in range(params.dbstart, params.dbstop + 1):
      site_id = site_ids[i]
      pos = site_id.find(".")
      if(pos > 0): site_id = site_id[:pos-2]
      dbase_site = sitemaps.sitemap(params.dbase_sites, site_id, params_cc,
                                    params.normalize)
      aligns.clear()
      align_method.align(dbase_site, aligns)
      score_method.score_alignments(aligns, dbase_site, out_file_thingie)

if(timer):
  s = timer.get()
  print "Elapsed time:         %.1f s." % (s[0])
  print "Process compute time: %.1f s." % (s[1])
  print "System compute time:  %.1f s." % (s[2] - s[1])

