from SimSite3DPy import *
import copy as what_is_clobbering_copy
import cPickle
from numpy import *
from pymol import cmd
from pymol.cgo import *
import pymol

def draw_stuff(R, T, q_hbond_vols, db_hbond_vols):
  db_sphere_rot_flags = [ 0 for i in range(len(db_hbond_vols)) ]

  # Transform the positions in the dbase_site
  for i in range(len(db_hbond_vols)):
    # Need to make deep copies of the relavent dbase objects?
    # NOTE!!!!!!!!!!!!!!!!!!!! The matrices in the align vector are 
    # for post multiplication
    db_hbond_vols[i].center = dot(db_hbond_vols[i].center, R) + T
    db_hbond_vols[i].ideal_dir = dot(db_hbond_vols[i].ideal_dir, R)
    db_hbond_vols[i].hbond_atom.position = \
      dot(db_hbond_vols[i].hbond_atom.position, R) + T
    # NOTE we are delaying the rotation of the exclusion spheres till needed

  # for cap in query, check all caps in db
  pt_cnts = []
  for q_idx, query_cap in zip(range(len(q_hbond_vols)), q_hbond_vols):
    pt_flags = [ 0 for i in range(query_cap.cap_pts.shape[0]) ]
    color_flags = [ 0 for i in range(query_cap.cap_pts.shape[0]) ]
    for i in range(len(db_hbond_vols)):
      db_cap = db_hbond_vols[i]
#      if((query_cap.hbond_atom.interact_type == "ACCEPTOR" and
#          db_cap.hbond_atom.interact_type == "DONOR") or
#         (query_cap.hbond_atom.interact_type == "DONOR" and
#          db_cap.hbond_atom.interact_type == "ACCEPTOR")): continue

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
      if(sq_dist >= query_cap.r1): pass

      if(not db_sphere_rot_flags[i]):
        db_sphere_rot_flags[i] = True
        for j in range(len(db_hbond_vols[i].exclusion_spheres)):
          db_hbond_vols[i].exclusion_spheres[j].center = \
            dot(db_hbond_vols[i].exclusion_spheres[j].center, R) + T

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
        # check if point lies outside of the spherical shell
        if(sq_dists[j] < min_sq_dist or max_sq_dist < sq_dists[j]):
          continue

        add_pt = True
        # check if the point is inside one of the exclusion spheres
        for ex_sphere in db_cap.exclusion_spheres:
          if(ex_sphere.contains(query_cap.cap_pts[j])):
            add_pt = False
            break

        if(add_pt):
          kept_pts.append(query_cap.cap_pts[j].tolist())
          idz.append(j)
      kept_pts = array(kept_pts)
      if(kept_pts.shape[0] == 0): continue

      # Need to check if the point lies in the upper nape of the DB
      # cone
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

        # Check if dot prod between the two vectors is close enough -- as with
        # site point scoring, we assume if the dot product is < 0, the
        # points are not complementary
        v2 = utils.unit_vec(query_cap.center, kept_pts[j])
        if(sum(pts_vecs[j] * v2) >= 0.0):
          pt_flags[idz[j]] = 1
          if((query_cap.hbond_atom.interact_type == "ACCEPTOR" and
              db_cap.hbond_atom.interact_type == "DONOR") or
             (query_cap.hbond_atom.interact_type == "DONOR" and
              db_cap.hbond_atom.interact_type == "ACCEPTOR")): 
            if(color_flags[idz[j]] == 0): color_flags[idz[j]] = 1 
          elif((query_cap.hbond_atom.interact_type == "ACCEPTOR" and
                db_cap.hbond_atom.interact_type == "ACCEPTOR") or
               (query_cap.hbond_atom.interact_type == "ACCEPTOR" and
                db_cap.hbond_atom.interact_type == "DONEPTOR") or
               (query_cap.hbond_atom.interact_type == "DONEPTOR" and
                db_cap.hbond_atom.interact_type == "ACCEPTOR")):
            if(color_flags[idz[j]] == 0 or color_flags[idz[j]] == 2): 
              color_flags[idz[j]] = 2 
            elif(color_flags[idz[j]] == 3): color_flags[idz[j]] = 4
          elif((query_cap.hbond_atom.interact_type == "DONOR" and
                db_cap.hbond_atom.interact_type == "DONOR") or
               (query_cap.hbond_atom.interact_type == "DONOR" and
                db_cap.hbond_atom.interact_type == "DONEPTOR") or
               (query_cap.hbond_atom.interact_type == "DONEPTOR" and
                db_cap.hbond_atom.interact_type == "DONOR")):
            if(color_flags[idz[j]] < 2): color_flags[idz[j]] = 3 
            elif(color_flags[idz[j]] == 2): color_flags[idz[j]] = 4
          elif((query_cap.hbond_atom.interact_type == "DONEPTOR" and
                db_cap.hbond_atom.interact_type == "DONEPTOR")):
            color_flags[idz[j]] = 4 
          else: print "missed interaction pairs somehow"

#      # given the current status we could overwrite the color of some points?
#      if(kept_pts.shape[0] > 0):
##        my_cgo = [BEGIN, POINTS]
#        my_cgo = []
#        if((query_cap.hbond_atom.interact_type == "ACCEPTOR" and
#            db_cap.hbond_atom.interact_type == "DONOR") or
#           (query_cap.hbond_atom.interact_type == "DONOR" and
#            db_cap.hbond_atom.interact_type == "ACCEPTOR")): 
#          my_cgo.extend([COLOR, 0.8, 0.2, 0.8])
#        elif((query_cap.hbond_atom.interact_type == "ACCEPTOR" and
#              db_cap.hbond_atom.interact_type == "ACCEPTOR") or
#             (query_cap.hbond_atom.interact_type == "ACCEPTOR" and
#              db_cap.hbond_atom.interact_type == "DONEPTOR") or
#             (query_cap.hbond_atom.interact_type == "DONEPTOR" and
#              db_cap.hbond_atom.interact_type == "ACCEPTOR")):
#          my_cgo.extend([COLOR, 0.2, 0.2, 0.8])
#        elif((query_cap.hbond_atom.interact_type == "DONOR" and
#              db_cap.hbond_atom.interact_type == "DONOR") or
#             (query_cap.hbond_atom.interact_type == "DONOR" and
#              db_cap.hbond_atom.interact_type == "DONEPTOR") or
#             (query_cap.hbond_atom.interact_type == "DONEPTOR" and
#              db_cap.hbond_atom.interact_type == "DONOR")):
#          my_cgo.extend([COLOR, 0.8, 0.2, 0.2])
#        elif((query_cap.hbond_atom.interact_type == "DONEPTOR" and
#              db_cap.hbond_atom.interact_type == "DONEPTOR")):
#          my_cgo.extend([COLOR, 0.9, 0.9, 0.9])
#        else: print "missed interaction pairs somehow"
#
#        for p in kept_pts:
#          my_cgo.append(SPHERE)
#          my_cgo.extend(p)
#          my_cgo.append(0.05)
#        #my_cgo.append(END)
#        cmd.load_cgo(my_cgo, "Q%d_D%d_overlap" % (q_idx, i))

    pt_cnts.append(sum(pt_flags))
    if(sum(pt_flags)):
      my_cgo = []
      for j in range(len(color_flags)):
        if(color_flags[j] == 0): continue
        elif(color_flags[j] == 1):
          my_cgo.extend([COLOR, 0.8, 0.2, 0.8])
        elif(color_flags[j] == 2):
          my_cgo.extend([COLOR, 0.2, 0.2, 0.8])
        elif(color_flags[j] == 3):
          my_cgo.extend([COLOR, 0.8, 0.2, 0.2])
        elif(color_flags[j] == 4):
          my_cgo.extend([COLOR, 0.9, 0.9, 0.9])

        my_cgo.append(SPHERE)
        my_cgo.extend(query_cap.cap_pts[j])
        my_cgo.append(0.05)
      cmd.load_cgo(my_cgo, "Q%d_overlap" % (q_idx))


if(__name__ == "__main__"):
  pymol.finish_launching()

  data_dir = "/psa/results/SimSite3D_datasets/testing/pterins/"
  q_id = "2qx0_ph2"
  #db_id = "2qx0_ph2"
  db_id = "2toh_hbi"
  #db_id = "1cbk_roi"
  #db_id = "1u72_mtx"

  q_fname = "%s/ptr_pockets/%s.pkl" % (data_dir, q_id)
  db_fname = "%s/dbase/%s.pkl" % (data_dir, db_id)
  cmd.load("%s/ptr_pockets/%s_rad.pdb" % (data_dir, q_id), q_id + "_p")
  cmd.load("%s/ptr_pockets/%s_rad.pdb" % (data_dir, db_id), db_id + "_p")
  cmd.load("%s/ligands/%s_l.mol2" % (data_dir, q_id))

  q_pkl_file = open(q_fname, "rb")
  q_hbond_vols = cPickle.load(q_pkl_file)[0].caps
  q_pkl_file.close()

  db_pkl_file = open(db_fname, "rb")
  db_hbond_vols = cPickle.load(db_pkl_file)[0].caps
  db_pkl_file.close()

  #my_vols = blah[0].caps
  db_vol_copy = what_is_clobbering_copy.deepcopy(db_hbond_vols)
  R = array([1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]).reshape((3,3))
  T = array([0.0, 0.0, 0.0])
  draw_stuff(R, T, q_hbond_vols, db_vol_copy)

  
  



