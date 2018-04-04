import os
from numpy import *
from SimSite3DPy import sitemaps
from SimSite3DPy import search
from SimSite3DPy import utils
params_cc = search.parameters()

datasets = [ "adenines", "pterins", "inhibited_gst" ]
pocket_dirs = ["ade_pockets", "ptr_pockets", "hsite_pockets"]

for pocket_dir, dataset in zip(pocket_dirs, datasets):
  data_dir = "../" + dataset

  if(len(pocket_dir)): sitemap_dir = data_dir + "/" + pocket_dir
  else: sitemap_dir = data_dir + "/dbase"
  prots_dir = data_dir + "/proteins"
  ligs_dir = data_dir + "/ligands"

  for fname in os.listdir(sitemap_dir):
    if(not fname.endswith(".pkl")): continue
  
    site_id = fname[:-4]

    lig_fname = ligs_dir + "/" + site_id + "_l.mol2"
    prot_fname = prots_dir + "/" + site_id + "_p.pdb"
    prot = utils.pdb.residues(prot_fname)
    ideal_pts_fname = "/home/vanvoor4/code/SimSite3D/trunk/params/"
    ideal_pts_fname += "new_optimum_hbonds.dat"
    ideal_pts = sitemaps.hbond_ideal_pts(ideal_pts_fname)
  
    lig = utils.mol2.molecule(lig_fname)
    centers = array([ a.position for a in lig.atoms ])
    radii = array([ 2.5 for a in lig.atoms ])
    site_vol = utils.volumes.UnionOfBalls(centers, radii)
    hbond_triplets = ideal_pts.get_polar_atom_triplets(prot, site_vol)
  
    dbase_site = sitemaps.sitemap(sitemap_dir, site_id, params_cc, False)
    dbase_site.load_hbond_caps()
    ofile = open("%s/%s_surf_caps.csv" % (sitemap_dir, site_id), "w+")
  
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
          #['__doc__', '__init__', '__module__', '__repr__', 'center',
          #    'chord_mid_pt', 'contains', 'in_dir', 'intersection', 'mid_pt',
          #    'pts', 'radius']
  
  
        toks.append(",".join(val_strs))
  
      print >> ofile, "|".join(toks) + "|"
  
    # handle metals
    for metal in dbase_site.metal_surfs.spheres:
      print metal
     
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
