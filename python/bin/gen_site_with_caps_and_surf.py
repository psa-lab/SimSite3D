import os
import sys
from numpy import *
import cPickle
from SimSite3DPy import *


################################################################################

################################################################################
def get_options():
  from optparse import OptionParser

  if(len(sys.argv) < 5):
    print >> sys.stderr, "\n\t * A protein and ligand file are required *\n"
    return ([], [])

  cmd_parser.add_option("-i", "--input", metavar="<fname>.csv",
                        help="Pipe delimted csv file with the first field holding the pdb code, the second field holding chainID if --chainIDs is specified")

  opts = OptionParser()
  opts.add_option("", "--proj_output", metavar="<DIR>", default="",
                  help="Directory to save site map files")
  opts.add_option("-p", "--protein", metavar="/path/to/protein.pdb",
                  help="protein pdb file")
  opts.add_option("-l", "--ligand", metavar="/path/to/ligand.mol2",
                  help="ligand mol2 file")

  (flags, args) = opts.parse_args()
  return (flags, args)
################################################################################

################################################################################
def convert_caps_to_Cpp_text(install_dir, prot_fname, lig_fname, proj_output, 
                             site_id):
  params_cc = search.parameters()
  
  prot = utils.pdb.residues(prot_fname)
  ideal_pts_fname = "%s/params/new_optimum_hbonds.dat" % (install_dir)
  ideal_pts = sitemaps.hbond_ideal_pts(ideal_pts_fname)

  lig = utils.mol2.molecule(lig_fname)
  centers = array([ a.position for a in lig.atoms ])
  radii = array([ 2.5 for a in lig.atoms ])
  site_vol = utils.volumes.UnionOfBalls(centers, radii)
  hbond_triplets = ideal_pts.get_polar_atom_triplets(prot, site_vol)

  dbase_site = sitemaps.sitemap(proj_output, site_id, params_cc, False)
  dbase_site.load_hbond_caps()
  ofile = open("%s/%s_surf_caps.csv" % (proj_output, site_id), "w+")

  for cap in dbase_site.hbond_surf_caps.caps:
    toks = []
    # point type
    if(cap.hbond_atom.interact_type== "ACCEPTOR"): toks.append("DONOR")
    elif(cap.hbond_atom.interact_type == "DONOR"): toks.append("ACCEPTOR")
    elif(cap.hbond_atom.interact_type == "DONEPTOR"): toks.append("DONEPTOR")
    else:
      print >> stderr, "expected a polar point, got something else (%s)" % \
        cap.hbond_atom.interact_type 
      sys.exit(-1)

    # hbond atom number
    toks.append("%d" % (cap.hbond_atom.serial))

    for atom_triplet in hbond_triplets:
      (C_nbr, hbond_atom, other_nbr) = atom_triplet
      if(hbond_atom == cap.hbond_atom):
        print hbond_atom
        break
    print C_nbr.serial, other_nbr.serial


    # C nbr atom num
    toks.append("%d" %(C_nbr.serial))
    # 2nd nbr atom number
    toks.append("%d" % (other_nbr.serial))
    # cap number
    toks.append("%d" % (cap.cap_number))
    # cap plane
    val_strs = [ "%f" % (f) for f in cap.ideal_dir ]
    val_strs.extend([ "%f" % (f) for f in cap.cap_plane_P0 ])
    toks.append(" ".join(val_strs))

    # print circles
    for C in cap.circles:
      # sphere
      val_strs = [ "%f %f %f %f" % (C.center[0], C.center[1], C.center[2],
                                    C.radius)]
      # plane
      val_strs.append( "%f %f %f %f %f %f" % (C.N[0], C.N[1], C.N[2],
                                              C.center[0], C.center[1],
                                              C.center[2]))
      # arcs
      for A in C.final_arcs:
        strs = [ "%f" % (f) for f in A.pts[0] ]
        strs.extend([ "%f" % (f) for f in A.mid_pt ])
        strs.extend([ "%f" % (f) for f in A.pts[1] ])
        val_strs.append(" ".join(strs))

      toks.append(",".join(val_strs))

    print >> ofile, "|".join(toks) + "|"
    
  # handle metals
  for metal in dbase_site.metal_surfs.spheres:

    # Metal atoms don't have neighbors per se & have no plane or ideal direction
    # (in our model)
    toks = ["METAL", "%d" % (metal.metal_atom.serial),
            "%d" % (metal.close_atoms[0].serial),
            "%d" % (metal.close_atoms[1].serial), "0", "0 0 0 0 0 0"]

    for C in metal.circles:
      # sphere
      val_strs = [ "%f %f %f %f" % (C.center[0], C.center[1], C.center[2],
                                    C.radius)]
      # plane
      val_strs.append("%f %f %f %f %f %f" % (C.N[0], C.N[1], C.N[2],
                                             C.center[0], C.center[1],
                                             C.center[2]))
      # arcs
      for A in C.final_arcs:
        strs = [ "%f" % (f) for f in A.pts[0] ]
        strs.extend([ "%f" % (f) for f in A.mid_pt ])
        strs.extend([ "%f" % (f) for f in A.pts[1] ])
        val_strs.append(" ".join(strs))
      toks.append(",".join(val_strs))

    print >> ofile, "|".join(toks) + "|"

  ofile.close()
################################################################################

################################################################################


dirs = parameters.parameters_base()
print
print "-" * 80
dirs.version()
print "-" * 80
print

(cmd_flags, cmd_args) = get_options()
if(len(cmd_flags) < 1): sys.exit(1)

dirs.get_params()
if(dirs.fail): sys.exit(1)

if(len(cmd_flags.proj_output)): dirs.proj_output = cmd_flags.proj_output
 
cmd = "%s/bin/gen_points " % (dirs.install_dir)
cmd += " -p %s -l %s " % (cmd_flags.protein, cmd_flags.ligand)
cmd += " --msms_surf --no_normalization "
cmd += " --min_res_per_chain 10 "
cmd += " --proj_output %s " % (dirs.proj_output)
status = os.system(cmd)

site_id = cmd_flags.protein.split("/")[-1]
if(site_id.endswith("_p.pdb")): site_id = site_id[:-6]
elif(site_id.endswith(".pdb")): site_id = site_id[:-4]

if(status or not os.path.isfile("%s/%s_s.csv" % (dirs.proj_output, site_id))):
  print >> sys.stderr, \
"""  Warning!
Could not create the sitemap for the protein\n
  %s 
and ligand
  %s
""" % (cmd_flags.protein, cmd_flags.ligand)
  if(status): print "The return value of gen_points was %d" % (status)

lig = utils.mol2.molecule(cmd_flags.ligand)
centers = array([ a.position for a in lig.atoms ])
radii = array([ 2.5 for a in lig.atoms ])
site_vol = utils.volumes.UnionOfBalls(centers, radii)

# Create the hbond volumes
print "Determining the hbond cap volumes"
hbond_cap_vols = sitemaps.hbond_volumes(cmd_flags.protein, site_vol)

# Create the hbond surface caps
print "Determining the hbond cap surfaces"
the_prot = utils.pdb.residues(cmd_flags.protein)
hbond_cap_surfs = sitemaps.hbond_surf_caps(the_prot, site_vol, hbond_cap_vols)

# Create the metal volumes
lig = utils.mol2.molecule(cmd_flags.ligand)
centers = array([ a.position for a in lig.atoms ])
radii = array([ 3.0 + 2.5 for a in lig.atoms ])
metal_site_vol = utils.volumes.UnionOfBalls(centers, radii)
metal_vols = sitemaps.metal_volumes(cmd_flags.protein, metal_site_vol)

# Create the metal surfs
metal_surfs = sitemaps.metal_surfs(the_prot, metal_site_vol, metal_vols)

# Write out the volumes & surface caps
caps_fname = "%s/%s.pkl" % (dirs.proj_output, site_id)
caps_fout = open(caps_fname, "w+b")

cPickle.dump((hbond_cap_vols, hbond_cap_surfs, metal_vols, metal_surfs),
             caps_fout)
print "dumping to:", caps_fout
caps_fout.close()
convert_caps_to_Cpp_text(dirs.install_dir, cmd_flags.protein, cmd_flags.ligand,
                         dirs.proj_output, site_id) 
