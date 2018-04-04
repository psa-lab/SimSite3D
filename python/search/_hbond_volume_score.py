from SimSite3DPy import utils, sitemaps
from numpy import *
import copy
from sys import stdout

class hbond_volume_score:

  def __init__(self, query_site, params):
    self.query_site = query_site
    self.params = params
    # Assume file is always <query_id>_rad.pdb 
    self.query_id = query_site.atoms_file_name().split("/")[-1][:-8]

  def score_alignments(self, aligns, orig_db_site, ofile=stdout):

    # define a copy function for both hbond_volumes and hbond_volume

    orig_db_caps = orig_db_site.hbond_caps.caps
    db_caps = copy.deepcopy(orig_db_site.hbond_caps.caps)

#    for query_cap in self.query_site.hbond_caps.caps:
#      print query_cap.hbond_atom

#    for align in aligns:
    for zzz in range(aligns.size()):
      align = aligns[zzz]
      R = array(align.R).reshape((3,3))
      T = array(align.T)
      db_sphere_rot_flags = [ 0 for i in range(len(db_caps)) ]

      # Transform the positions in the dbase_site
      for i in range(len(db_caps)):
        # Need to make deep copies of the relavent dbase objects?
        # NOTE!!!!!!!!!!!!!!!!!!!! The matrices in the align vector are 
        # for post multiplication
        db_caps[i].center = dot(orig_db_caps[i].center, R) + T
        db_caps[i].ideal_dir = dot(orig_db_caps[i].ideal_dir, R)
        db_caps[i].hbond_atom.position = \
          dot(orig_db_caps[i].hbond_atom.position, R) + T
     
      # for cap in query, check all caps in db
      pt_cnts = []
      for query_cap in self.query_site.hbond_caps.caps:
        pt_flags = [ 0 for i in range(query_cap.cap_pts.shape[0]) ]
        for i in range(len(db_caps)):
          db_cap = db_caps[i]
          if((query_cap.hbond_atom.interact_type == "ACCEPTOR" and
              db_cap.hbond_atom.interact_type == "DONOR") or
             (query_cap.hbond_atom.interact_type == "DONOR" and
              db_cap.hbond_atom.interact_type == "ACCEPTOR")): continue
          
          squared_tol = (query_cap.r1 + db_cap.r1)**2
          sq_dist = utils.squared_dist(query_cap.hbond_atom.position, 
                                       db_cap.hbond_atom.position)
          if(sq_dist >= squared_tol): continue

          # Check if query sphere is too far below db sphere
          # This is done by checking if the query_cap's sphere is too far
          # below the plane defined by the db_cap's direction (axis) and
          # center of db_cap's spheres.  Too far is greater than or equal to
          # the query sphere's radius minus the offset of the circle defined by
          # the nappe and the smaller of the db's concentric spheres
          # (for a cone of 60 degrees this is 0.5*db_radius_0)
          my_dir = query_cap.center - db_cap.center
          signed_dist = sum(my_dir * db_cap.ideal_dir)
          if(signed_dist <= -(db_cap.r1 - 0.5*query_cap.r0)): continue

          my_dir = db_cap.center - query_cap.center
          signed_dist = sum(my_dir * query_cap.ideal_dir)
          if(signed_dist <= -(query_cap.r1 - 0.5*db_cap.r0)): continue
          

## Ignore these checks for now as they could prove to be more expensive
## than just checking the points -- 
          # Here we could check if the caps intersected by computing the 
          # intersection of the planes defined by the 2 cap extents.
          # If the planes intersect, check if the line intersects both of the
          # circles defining the planes -- simple distance from line to radius.
          if(sq_dist >= query_cap.r1):
            pass
          
          if(not db_sphere_rot_flags[i]):
            db_sphere_rot_flags[i] = True
            for j in range(len(db_caps[i].exclusion_spheres)):
              db_caps[i].exclusion_spheres[j].center = \
                dot(orig_db_caps[i].exclusion_spheres[j].center, R) + T

          # Discrete (point based) approximation of percentage of volume overlap
          tmp = query_cap.cap_pts.copy()
          tmp -= tile(db_cap.center, (query_cap.cap_pts.shape[0], 1))
          sq_dists = sum(tmp * tmp, 1)

          # These have been adjusted because self matches missed many of the 
          # points due to numerical errors
          max_sq_dist = (db_cap.r1 + 0.0005)**2
          min_sq_dist = (db_cap.r0 - 0.0005)**2

          # We keep track of the index of the points so that each query point
          # can be satisfied at most once
          kept_pts = []
          idz = []
          for j in range(query_cap.cap_pts.shape[0]):
            if(sq_dists[j] < min_sq_dist or max_sq_dist < sq_dists[j]):
              continue

            add_pt = True
            for ex_sphere in db_cap.exclusion_spheres:
              if(ex_sphere.contains(query_cap.cap_pts[j])):
                add_pt = False
                break
               
            if(add_pt):
              kept_pts.append(query_cap.cap_pts[j].tolist())
              idz.append(j)
          if(not kept_pts): continue  
          kept_pts = array(kept_pts)
        

          # Get unit vectors in dir of db atom center to pt in query cloud
          tmp = kept_pts - tile(db_cap.center, (kept_pts.shape[0], 1))
          dists = sqrt(sum(tmp*tmp, 1))
          pts_vecs = tmp / tile(dists.reshape(dists.shape[0], 1), (1, 3))
 
          cos_angles = pts_vecs * tile(db_cap.ideal_dir, (kept_pts.shape[0], 1))
          cos_angles = sum(cos_angles, 1)
          
          for j in range(kept_pts.shape[0]):
            # We allow a minimum DHA of 120 degrees in SLIDE.
            # The max angle between the two vectors can be estimated as 60
            # degrees.  cos(60 degrees) = 0.5.  Thus any vector with a
            # cos_angle >= 0.5 is kept

            # not the "correct" threshold, but close enough
            if(cos_angles[j] < 0.4995): continue

            # Check if dot prod between the two vectors is close enough
            v2 = utils.unit_vec(query_cap.center, kept_pts[j])
            if(sum(pts_vecs[j] * v2) >= 0.0):
              pt_flags[idz[j]] = 1

        pt_cnts.append(sum(pt_flags))

      ostr = ""
      percent_sum = 0.0
      npts = len(self.query_site.hbond_caps.caps)
      for i, cnt in zip(range(npts), pt_cnts):
        my_percent = 100.0 * cnt / self.query_site.hbond_caps.caps[i].cap_pts.shape[0]
        ostr += "%.2f%% " % (my_percent)
        percent_sum += my_percent

      vals = [ "%f" % (align.R[i]) for i in range(9) ]
      Rstr = " ".join(vals)
      vals = [ "%f" % (align.T[i]) for i in range(3) ]
      Tstr = " ".join(vals)

      db_site_id = orig_db_site.atoms_file_name().split("/")[-1][:-8]
      # Need to do the copying this way as the boost.python index does
      # not support (or maybe cannot easily support) getting an item by
      # index and then assigning a value to one of its members -- so we get
      # the item, adjust it, and then reassign it.
      tmp_align = aligns[zzz]
      tmp_align.hb_caps_score = percent_sum
      tmp_align.hb_caps_match_print = ostr
      aligns[zzz] = tmp_align
      print >> ofile, "%s|%s|%.2f|%s|%s|%s|" % \
        (self.query_id, db_site_id, percent_sum, Rstr, Tstr, ostr)

         
          



# items for testing this module
if __name__ == "__main__":
  from os import path
  import cPickle
  import sys
  dataset = "adenines"
  query = "1b38_atp"
  prot = "1b38_atp_p.pdb"
  lig = "1b38_ade_l.mol2"
  base_data_dir = "/home/vanvoor4/data/new_sampling"
  data_dir = base_data_dir + "/" + dataset

  prot = "%s/proteins/%s" % (data_dir, prot)
  lig = "%s/ligands/%s" % (data_dir, lig)
  csv = "%s/dbase/%s_s.csv" % (data_dir, query)
  rad = "%s/dbase/%s_rad.pdb" % (data_dir, query)
  pkl_fname = "%s/dbase/%s.pkl" % (data_dir, query)
  if(path.exists(pkl_fname)):
    try:
      pkl_in = file(pkl_fname, "rb")
    except IOError, (errno, strerror):
      print "Unable to open the file", fname
      print "error(%s): %s" % (errno, strerror)
      print "in %s" % (self.__module__)

    hbond_caps = cPickle.load(pkl_in)
    pkl_in.close()
  else:
    print  """
Unable to load the hbond caps file 
  (%s) 
because it doesn't exist
""" % (pkl_fname)
    sys.exit(-1)

  #print hbond_caps.caps[0].cap_pts
  from pymol import cmd
  from pymol.cgo import *
  import pymol
  pymol.finish_launching()
  import cgo_items

  cmd.load(rad, "rad_file")

  i = 8
   
  ideal_pt = hbond_caps.caps[i].center + 3.0*hbond_caps.caps[i].ideal_dir
  cgo_items.draw_points(ideal_pt.reshape((1,3)), "center%02d" %(i), color=[1.0,1.0,1.0], radius=0.2) 
  cgo_items.draw_points(hbond_caps.caps[i].cap_pts, "blah", color=[0.8, 0.0, 0.0]) 
  
  cmd.show("sticks", "i. 83")









