#!/usr/bin/python

import os
import sys
from numpy import *
import tempfile

def RMSD(A,B):
  tmp = A - B
 
  return sqrt( sum(tmp*tmp) / A.shape[0])

def setup_directories(base_dir, mode):
  if(not os.path.exists(base_dir + "/moved_ligands")):
    os.makedirs(base_dir + "/moved_ligands", 0770)
  if(not os.path.exists(base_dir + "/ligand_fragments")):
    os.makedirs(base_dir + "/ligand_fragments", 0770)

def add_rand_aligns(dbase_site, align_method, aligns, centroid, 
                    num_rand_aligns):
  from SimSite3DPy import search

  # Add the Identity alignment as well
  id_align_mthd = search.IdentityAlignment()
  id_align_mthd.align(dbase_site, aligns)

  orient_params = \
    align_method.align(dbase_site, aligns, centroid, N=num_rand_aligns,
                       q_width=0.005, half_side_len=0.025)
  orient_params.extend(
    align_method.align(dbase_site, aligns, centroid, N=num_rand_aligns,
                       q_width=0.010, half_side_len=0.05)
    )
  orient_params.extend(
    align_method.align(dbase_site, aligns, centroid, N=num_rand_aligns,
                       q_width=0.025, half_side_len=0.125)
    )
#  orient_params.extend(
#    align_method.align(dbase_site, aligns, centroid, N=num_rand_aligns,
#                       q_width=0.05, half_side_len=0.25)
#    )
#  orient_params.extend(
#    align_method.align(dbase_site, aligns, centroid, N=num_rand_aligns,
#                       q_width=0.1, half_side_len=0.5)
#    )
  return aligns


def copy_search_params(params_py):
  from SimSite3DPy import search
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
  params_cc.load_surf_files = True
  params_cc.require_min_npts = params_py.require_min_npts
  return params_cc

###############################################################################
proj_dir = os.getenv("SIMSITE3D_INSTALL_DIR")
if(proj_dir == None):
  print >> sys.stderr, \
    "\nUnable to get the environment variable $SIMSITE3D_INSTALL_DIR"
  print >> sys.stderr, "Please set $SIMSITE3D_INSTALL_DIR to a valid path"
  sys.exit(1)

sys.path.append(proj_dir + "/python")
from SimSite3DPy import *

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
                            params_cc, params.normalize, True, True)
#model.load_hbond_caps()
hbond_caps_score = search.hbond_volume_score(model, params)
if(params.num_rand_aligns == 0 and params.align_to_query):
  align_method = search.MatchTriangles(model, params.dmetol, params.lsetol,
                                       params.hydrophobic_query)
elif(not params.align_to_query):
  align_method = search.IdentityAlignment()
elif(params.num_rand_aligns > 0):
  align_method = search.RandomAlignments()
  site_pts = utils.pdb.hetgroup("%s/%s_s.pdb" % (params.query_sitemap_dir,
                                                 params.query_site_id))
  site_pos = array([ a.position for a in site_pts.atoms ])
  centroid = mean(site_pos, 0)
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
#score_method = search.WeightedSumsScore(model, params_cc)
score_method = search.HbondSurfacesScore(model, params_cc)
score_method.write_score_header(out_file_thingie)

aligns = search.vector_less__rigid_align_t__greater_()
# Pairwise comparison
if(len(params.dbase_site_id)):
  dbase_site = sitemaps.sitemap(params.dbase_sitemap_dir, params.dbase_site_id,
                                params_cc, False, True, False)
  #dbase_site.load_hbond_caps()
  if(params.num_rand_aligns == 0):
    align_method.align(dbase_site, aligns)
    print "time to score"
    score_method.score_alignments(aligns, dbase_site, out_file_thingie)
  else:
    aligns = add_rand_aligns(dbase_site, align_method, aligns, centroid, 
                             params.num_rand_aligns)

#    # Add the Identity alignment as well
#    id_align_mthd = search.IdentityAlignment()
#    id_align_mthd.align(dbase_site, aligns)
#
#    orient_params = \
#      align_method.align(dbase_site, aligns, centroid, N=params.num_rand_aligns,
#                         q_width=0.1, half_side_len=0.5)
#    orient_params.extend(
#      align_method.align(dbase_site, aligns, centroid, N=params.num_rand_aligns,
#                         q_width=0.05, half_side_len=0.25)
#      )
#    orient_params.extend(
#      align_method.align(dbase_site, aligns, centroid, N=params.num_rand_aligns,
#                         q_width=0.025, half_side_len=0.125)
#      )
#    orient_params.extend(
#      align_method.align(dbase_site, aligns, centroid, N=params.num_rand_aligns,
#                         q_width=0.01, half_side_len=0.05)
#      )

    print "time to score"
    score_method.score_alignments(aligns, dbase_site, out_file_thingie)
    out_file_thingie.close()
      
else:
  # do full db scan 
  if(not len(params.db_index_fname)):
    for site in os.listdir(params.dbase_sites):
      if(not site.endswith("_s.csv")):
        continue

      dbase_site = sitemaps.sitemap(params.dbase_sites, site[:-6], params_cc, 
                                    False, True, False)
      aligns.clear()
      aligns = add_rand_aligns(dbase_site, align_method, aligns, centroid, 
                               params.num_rand_aligns)
      score_method.score_alignments(aligns, dbase_site, out_file_thingie)

  # look in the db_index_fname to start with dbstart and finish with
  # dbstop
  else:
    db_index_file = file(params.db_index_fname, "r")
    site_ids = []
    for line in db_index_file:
      if(line.startswith("#")): continue
      site_ids.append(line.strip())
 
    if(params.dbstart == 0 and params.dbstop == 0):
      params.dbstop = len(site_ids) - 1
 
    for i in range(params.dbstart, params.dbstop + 1):
      site_id = site_ids[i]
      pos = site_id.find(".")
      if(pos > 0): site_id = site_id[:pos-2]
      dbase_site = sitemaps.sitemap(params.dbase_sites, site_id, params_cc,
                                    False, True, False) 
      aligns.clear()
      aligns = add_rand_aligns(dbase_site, align_method, aligns, centroid, 
                               params.num_rand_aligns)
      score_method.score_alignments(aligns, dbase_site, out_file_thingie)

if(timer):
  s = timer.get()
  print "Elapsed time:         %.1f s." % (s[0])
  print "Process compute time: %.1f s." % (s[1])
  print "System compute time:  %.1f s." % (s[2] - s[1])

