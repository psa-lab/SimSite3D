from ASCbasePy import utils, sitemaps
from numpy import *
from sys import stderr

class hbond_surf_score:

  def __init__(self, query_site, params):
    self.query_site = query_site
    self.params = params
    # Assume file is always <query_id>_rad.pdb 
    self.query_id = query_site.atoms_file_name().split("/")[-1][:-8]

    pos = []
    dir = []
    act = []
    hbond_pts = self.query_site.hbond_points()
    for pt in hbond_pts.fit_pts():
      pos.append([ pt.pos[i] for i in range(3) ])
      dir.append([ pt.dir[i] for i in range(3) ])
      act.append(pt.act_type)

    self.query_pos = array(pos)
    self.query_dir = array(dir)
    self.act_types = act

  def score_alignments(self, aligns, db_site, tol=1.0):
    """
    Inputs:
      aligns:  python interface to std::vector<rigid_align_t> holding the
               alignments to score.
      db_site: dbase hbond surf caps thingie (object, instance, whatever)
      tol: maximum distance for corresponding points to count to score

    Returns:
      correspondances
      some type of scores or features

    Features:
      1) Count of query groups matched -- estimate by counting each group as
         matched at most once.
      2) Sum of dot product weighted matches
      3) Sum of linear error weighted matches (tol - dist)
      4) Sum of squared error weighted matches (tol - dist)**2
      5) Sum of linear error * dot product weighted matches
      6) Sum of squared error * dot product weighted matches
    """
    from ASCbasePy.sitemaps import ACCEPTOR, DONOR, DONEPTOR

    sq_tol = tol*tol

    db_acts = []
    for cap in db_site.hbond_surf_caps.caps:
      print cap.hbond_atom
      if(cap.hbond_atom.interact_type == "ACCEPTOR"): db_acts.append(DONOR)
      elif(cap.hbond_atom.interact_type == "DONOR"): db_acts.append(ACCEPTOR)
      elif(cap.hbond_atom.interact_type == "DONEPTOR"): db_acts.append(DONEPTOR)
      else:
        print >> stderr, "expected a polar point, got something else (%s)" % \
          cap.hbond_atom.interact_type
        return


    # Query sitemap is represented by pts
    # DB sites are represented by caps 
    for zzz in range(aligns.size()):
      align = aligns[zzz]
      R = array(align.R).reshape((3,3))
      T = array(align.T)

      # Move the query points to the dbase reference frame & rotate directions
      # to the dbase ref
      N_q = self.query_pos.shape[0]
      q_pos = dot(self.query_pos - tile(T, (N_q, 1)), R.T)
      q_dir = dot(self.query_dir, R.T)
      q_act = self.act_types

      # For now -- use the site maps as given -- we can increase sampling
      # density at a later date.  Also, treat each point by itself (i.e. don't
      # consider its neighborhood); again, this is reducing the information we
      # have for each alignment, but makes it easier to
      # get something prototyped

# Use brute force to get it going
      
      corr_pt = array([[] for i in range(N_q)])
      terms = zeros((N_q, 6))

      for db_act, db_cap in zip(db_acts, db_site.hbond_surf_caps.caps):
        for i in range(N_q):

# for now we will ignore polar repulsive points
          if((q_act[i] == ACCEPTOR and db_act == DONOR) or \
             (q_act[i] == DONOR and db_act == ACCEPTOR)): continue

          # Compute closest point
          pt = db_cap.closest_point(q_pos[i])
          if(not pt.shape[0]): continue

          # Scoring distance filter
          tmp = pt - q_pos[i]
          sq_dist = sum(tmp*tmp)
          if(sq_dist > sq_tol): continue

          # Compute the direction for pt
          pt_dir = utils.unit_vec(db_cap.center, pt)
          # Dot product filter 
          dot_prod = sum(pt_dir * q_dir[i])
          if(dot_prod <= 0.0): continue

          # 0) Count of query groups matched -- estimate by counting each group
          #    as matched at most once.
          terms[i][0] += 1

          # 1) Sum of dot product weighted matches
          if(terms[i][1] < dot_prod): terms[i][1] = dot_prod

          # 2) Sum of linear error weighted matches (tol - dist) / tol
          # 3) Sum of squared error weighted matches (tol - dist)**2 / tol**2
          linear_w = (tol - sqrt(sq_dist)) / tol
          if(terms[i][2] < linear_w): 
            terms[i][2] = linear_w
            terms[i][3] = linear_w * linear_w

          # 4) Sum of linear error * dot product weighted matches
          # 5) Sum of squared error * dot product weighted matches
          tmp = linear_w * dot_prod
          if(terms[i][4] < tmp):
            terms[i][4] = tmp
            terms[i][5] = tmp * linear_w

      # using the vector_less__rigid_align_t__greater_ construct does not
      # support assignment such as aligns[zzz].blah = 1.0 (since the bracket 
      # operator does not work as one might expect).  So we do it this way.
      align.terms.append( self.count_groups_hit(terms[:,0]) )
      align.terms.extend( terms[:,1:].sum(0).tolist() )
      aligns[zzz] = align

  def count_groups_hit(self, pt_hits):

    hbond_pts = self.query_site.hbond_points()
    groups = {}
    for (pt, num) in zip(hbond_pts.fit_pts(), pt_hits):
      if(num <= 0.0): continue

      ideal_pt = pt.get_ideal_pt()
      atom = pt.get_atom()
      my_key = "%d_%d" % (atom.atom_num, ideal_pt.pt_num)
      groups[my_key] = True
    return len(groups)

    

if(__name__ == "__main__"):
  import sys
  from ASCbasePy import parameters 
  from ASCbasePy.search import ModelSitemap 
  from ASCbasePy.search import vector_less__rigid_align_t__greater_
  from ASCbasePy.search import IdentityAlignment

  params = parameters.search_parameters()
  rv = params.cmdline_options(sys.argv)
  if(rv == False or params.fail): sys.exit(-1)
  params_cc = params.to_C()

  query = ModelSitemap(params.query_sitemap_dir, params.query_site_id,
                       params_cc, params.normalize)

  score_method = hbond_surf_score(query, params)

  aligns = vector_less__rigid_align_t__greater_()
  # Pairwise comparison
  if(len(params.dbase_site_id)):
    dbase_site = sitemaps.sitemap(params.dbase_sitemap_dir, 
                                  params.dbase_site_id, params_cc, False)
    dbase_site.load_hbond_caps()
    id_align_mthd = IdentityAlignment()
    id_align_mthd.align(dbase_site, aligns)
    score_method.score_alignments(aligns, dbase_site)

