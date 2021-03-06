import sys
from numpy import *
from SimSite3DPy import *

# debugging modules
from pymol import cmd
from pymol.cgo import *
import cgo_items

# state number used to draw hbond caps
#s_num = 1
#arc_num = 0

class arc:
  """
  This class has not had checks for numerical stability at this point
  """
  def __init__(self, center, radius, end_pt0, mid_pt, end_pt1):

    self.radius = radius
    self.center = center
    self.pts = array([end_pt0, end_pt1])
    self.mid_pt = mid_pt
    self.in_dir = utils.unit_vec(self.center, mid_pt)
    self.chord_mid_pt = 0.5*(end_pt0 + end_pt1)
 
  def __repr__(self):
    return "%s %f %s %s %s" % (self.center, self.radius, self.pts[0], \
                               self.mid_pt, self.pts[1])

  def contains(self, pt):
     my_dir = pt - self.chord_mid_pt
     if(sum(my_dir * self.in_dir) >= 0): return True
     else: return False  

  def intersection(self, other):
    """
    Assumption:
      Both arcs lie on the same circle

    Compute the intersection of self and the other arc without modifying
    self or the other arc
    """
    self_contains_other_pts = [ self.contains(other.pts[0]), 
                                self.contains(other.pts[1]) ]
    other_contains_self_pts = [ other.contains(self.pts[0]), 
                                other.contains(self.pts[1]) ]

    # Case 1) Other arc is contained in self -- intersection is other_arc
    if(self_contains_other_pts[0] and self_contains_other_pts[1] and \
       not other_contains_self_pts[0] and not other_contains_self_pts[1]): 
      return (other,)

    # Case 2) Other arc contains self -- intersection is self
    elif(other_contains_self_pts[0] and other_contains_self_pts[1] and \
         not self_contains_other_pts[0] and not self_contains_other_pts[1]):
      return (self, )

    # Case 3) Other arc has no overlap with self -- no intersection
    elif(not self_contains_other_pts[0] and not self_contains_other_pts[1] and \
         not other_contains_self_pts[0] and not other_contains_self_pts[1]):
      return (None, )

    # Case 4) Other arc partially overlaps self -- intersection is overlap
    # Case 4a) all points from each arc are contained in the other arc
    elif(self_contains_other_pts[0] and self_contains_other_pts[1] and \
         other_contains_self_pts[0] and other_contains_self_pts[1]):
      
      # idea -- move the chord for self along the in direction until
      # we reach an endpoint of the other arc
      # (y = mx + b)
      sq_dists = array([ utils.squared_dist(self.pts[0], other.pts[0]),
                         utils.squared_dist(self.pts[0], other.pts[1]),
                         utils.squared_dist(self.pts[1], other.pts[0]),
                         utils.squared_dist(self.pts[1], other.pts[1]) ])
      idx = sq_dists.argmin() 
      self_idx = int(floor(idx / 2.0))
      other_idx = idx % 2

      # Because we are taking the two closest points (1 end point from each
      # arc), we know that the in direction is the vector from the center of
      # the circle to midpoint of the corresponding chord
      mid_pt = 0.5*(self.pts[self_idx] + other.pts[other_idx])
      my_dir = utils.unit_vec(self.center, mid_pt)
      mid_pt = self.center + self.radius * my_dir

      if(utils.squared_dist(self.pts[self_idx], other.pts[other_idx]) < 1E-14):
        A0 = None
      else:
        A0 = arc(self.center, self.radius, self.pts[self_idx], mid_pt, 
                 other.pts[other_idx])
      
      # Next, we use the other 2 points do define another chord and compute
      # the direction from the center to the midpoint of the chord.
      self_idx = (self_idx + 1) % 2
      other_idx = (other_idx + 1) % 2
      mid_pt = 0.5*(self.pts[self_idx] + other.pts[other_idx])
      my_dir = utils.unit_vec(self.center, mid_pt)
      mid_pt = self.center + self.radius * my_dir

      if(utils.squared_dist(self.pts[self_idx], other.pts[other_idx]) < 1E-14):
        A1 = None
      else:
        A1 = arc(self.center, self.radius, self.pts[self_idx], mid_pt, 
                 other.pts[other_idx])
     
      # If the second arc contains one of the end points of other_arc, we need 
      # to flip the direction.
      if(A1.contains(A0.pts[0]) or A1.contains(A0.pts[1])):
        my_dir *= -1.0
        A1.mid_pt = A1.center + A1.radius * my_dir
      return (A0,A1)

    # Case 4b) One point from each arc is contained in the other arc
    else:
      end_pts = zeros((2,3))
      if(self_contains_other_pts[0]): end_pts[0] = other.pts[0]
      else: end_pts[0] = other.pts[1]
      if(other_contains_self_pts[0]): end_pts[1] = self.pts[0]
      else : end_pts[1] = self.pts[1]
      
      # We will ignore very small arcs or those that are computed to be
      # a point 
      if(utils.squared_dist(end_pts[0], end_pts[1]) < 1E-14):
        return (None, )
  
      # compute arc mid_pt 
      mid_pt = 0.5*(end_pts[0] + end_pts[1])
      my_dir = utils.unit_vec(self.center, mid_pt)
  
      mid_pt = self.center + self.radius * my_dir
      if(not self.contains(mid_pt) or not other.contains(mid_pt)):
        mid_pt = self.center - self.radius * my_dir
 
      return (arc(self.center, self.radius, end_pts[0], mid_pt, end_pts[1]), )
################################################################################

################################################################################
class iCircle:

  def __init__(self, center, radius, normal, initial_arc=None):
    """
    INPUTS:
      center  Center of the circle
      radius  Radius of the circle
      normal  Normal to the plane containing the circle 
    """
    self.full_circle = True
    self.center = center
    self.radius = radius
    self.squared_rad = radius*radius
    self.N = normal
    # Plane containing the circle is defined by the circle's center and given 
    # normal.  Here self.N _dot_ X + self.d = 0
    self.d = -1.0 * sum(self.N * self.center)
    self.final_arcs = []
    self.removed_arcs = {}
    if(initial_arc == None): return

    self.full_circle = False
    self.final_arcs.append(initial_arc)
    # We may wish to keep track of which neighbors ate parts of our circle
    # then the following would hold not the kept arcs but the eaten arcs
    rm_arc_mid_pt = initial_arc.center - initial_arc.radius * initial_arc.in_dir
    rm_arc = arc(initial_arc.center, initial_arc.radius, initial_arc.pts[0],
                 rm_arc_mid_pt, initial_arc.pts[1])
    self.removed_arcs[rm_arc] = None

  def remove_overlap(self, other, idx=0, jdx=0):
    #print "TEST if self is contained inside the other circle"
#    global s_num
#    s_num += 1
    tmp = []
    for AA in self.final_arcs:
      if(not other.contains(AA.pts[0]) or not other.contains(AA.pts[1])):
        tmp.append(AA)
      else:
        print "ENTRY: entire arc for %d is contained in circle %d" % (idx, jdx)
    self.final_arcs = tmp
    if(not self.full_circle and len(self.final_arcs) == 0): return

    # 1) Check if the two spheres intersect
    D2 = utils.squared_dist(self.center, other.center)
    if(D2 >= (self.radius + other.radius)**2): 
      print "no itersection due to radius constraints %d %d" % (idx, jdx)
      return

    # 1b) Check if self is entirely contained in other
    if(other.contains(self.center)):
      V = utils.unit_vec(other.center, self.center)
      proj_pt = other.center + other.radius * V
      if(utils.squared_dist(self.center, proj_pt) >= self.squared_rad):
        self.full_circle = False
        self.final_arcs = []
        print "ENTRY: entire sphere for %d is contained in sphere %d" % (idx, jdx)
        return
      

    # 2) Determine the line (chord) where the two circles intersect
    # 2a) Plane intersection is the line Y = mX + b
    (m,b) = utils.intersect_planes(self.N, self.d, other.N, other.d)    

#    if(idx == 13 and jdx == 2):
#      my_cgo = [COLOR, 0.8, 0.8, 0.0, SPHERE]
#      my_cgo.extend(self.center)
#      my_cgo.append(self.radius)
#      my_cgo.extend([COLOR, 1.0, 1.0, 1.0, SPHERE])
#      my_cgo.extend(other.center)
#      my_cgo.append(other.radius)
#      #cmd.load_cgo(my_cgo, "circles_%d" % (s_num), state=s_num)
#      #cmd.load_cgo(my_cgo, "circles", state=idx + 1)
#      cmd.load_cgo(my_cgo, "circles", state=2)
#
#      #cgo_items.draw_plane2(self.N, self.center, plane_num=5, state=idx + 1)
#      #cgo_items.draw_plane2(other.N, other.center, plane_num=6, state=idx + 1)
#      cgo_items.draw_plane2(self.N, self.center, plane_num=5, state=2)
#      cgo_items.draw_plane2(other.N, other.center, plane_num=6, state=2)

    # 2b) Determine intersection between line defined by m & b and self
    (p0,p1) = utils.line_intersect_sphere(self.center, self.radius, m, b)
#    if(idx == 13 and jdx == 2):
      #cgo_items.draw_line(m, b, label="line_%d_%d" % (idx, jdx), state=idx + 1)
#      cgo_items.draw_line(m, b, label="line_%d_%d" % (idx, jdx), state=2)

    if(p0.shape[0] == 0 or p1.shape[0] == 0): 
      print "no intersection between the intersection circles %d & %d" % (idx,jdx)
      return


    # Get the 2 arcs for self -- both the remaining or kept arc and the
    # removed arc or arc that falls inside the other sphere
    my_dir = utils.unit_vec(self.center, (p0 + p1)/2.0)
    mid_pts = array([ self.center + self.radius*my_dir,
                      self.center - self.radius*my_dir ])
    sq_dists = [ utils.squared_dist(other.center, mid_pts[i]) for i in range(2)]

    if(sq_dists[0] < other.squared_rad and sq_dists[1] < other.squared_rad):
      # should never get here
      print "Entire circle is contained in the other circle -- ERROR"
    elif(sq_dists[0] < other.squared_rad):
      keep_mid_pt = mid_pts[1]
      rm_mid_pt = mid_pts[0]
    elif(sq_dists[1] < other.squared_rad):
      keep_mid_pt = mid_pts[0]
      rm_mid_pt = mid_pts[1]
    else:
      # should never get here
      print "Circles do not intersect -- ERROR"

    if(utils.squared_dist(p0, p1) < 1E-14):
      # Keep arc is very small, the contribution of self to total difference
      # is negligible
      if(utils.squared_dist(keep_mid_pt, p0) < 1E-14):
        self.full_circle = False
        self.final_arcs = []
        return
      # Remove arc is very small, the contribution of other to total difference
      # is negligible
      elif(utils.squared_dist(rm_mid_pt, p0) < 1E-14):
        return

    # In direction of overlap is the direction from center to midpoint of chord
    keep_arc = arc(self.center, self.radius, p0, keep_mid_pt, p1)
    rm_arc = arc(self.center, self.radius, p0, rm_mid_pt, p1)
    
    # 3) Update the arc fields for self
    self.removed_arcs[rm_arc] = other
    if(self.full_circle):
      self.full_circle = False
      self.final_arcs.append(keep_arc)
    else:
      new_final_arcs = []
      for AA in self.final_arcs:
        if(other.contains(AA.pts[0]) and other.contains(AA.pts[1])):
          print " entire arc for %d is contained in circle %d" % (idx, jdx)
#        print "intersect arcs:" 
#        print "\tfinal arc", AA
        A_tup = keep_arc.intersection(AA)
#        print "\tresult", A_tup
        for i in range(len(A_tup)):
          if(not A_tup[i] == None): 
            print "(%d, %d) kept arc" % (idx, jdx), A_tup[i]
            new_final_arcs.append(A_tup[i])
      self.final_arcs = new_final_arcs
      print "number of final arcs:", len(self.final_arcs)

#    if(idx == 13):
#      #print "in state %d, we have %d arcs" % (s_num, len(self.final_arcs))
#      arc_num = 0
#      for my_arc in self.final_arcs:
#        arc_num += 1
#        print "ARC_%d_%d_%d : %s" % ( idx, jdx, arc_num, my_arc.__repr__())
#        cgo_items.draw_arc(my_arc.center, my_arc.radius, my_arc.pts,
#                           my_arc.mid_pt, 
#                           label="ARC_%d_%d_%d" % (idx, jdx, arc_num),
#                           color=[0.2, 0.8, 0.2], state=2)


  def contains(self, pt):
    """
    Assumption: this method is only to be called for points that lie on the
    surface of the spherical cap

    NOTE: This function is not "arc aware"
    """
    if(utils.squared_dist(self.center, pt) > self.squared_rad): return False
    return True
################################################################################

################################################################################
class hbond_cap:

  def __init__(self, hbond_atom, C_nbr_atom, ideal_pt, prot, site_vol, 
               cap_number, include_metals):
    
    self.cap_number = cap_number
    self.hbond_atom = hbond_atom
    self.center = array(hbond_atom.position)
    self.ideal_dir = utils.unit_vec(self.center, array(ideal_pt))
    self.radius = 3.0

    #debugging
#    print "ideal dir:", self.ideal_dir
#    print "\"point position\":", self.center + self.radius * self.ideal_dir

#    cgo = [COLOR, 0.0, 0.0, 0.8, SPHERE]
#    cgo.extend(self.center)
#    cgo.append(self.radius)
#    cmd.load_cgo(cgo, "LEU 83 N")

    # Set this next variable -- depends on amount of spherical cap we want to
    # use in the model -- 1.5 is only a guess
    print "1.5 here is only a guess"
    h = 1.5
    self.cap_plane_P0 = self.center + h * self.ideal_dir
    self.cap_plane_d = -1.0 * sum(self.ideal_dir * self.cap_plane_P0)
    self.cap_plane_rad = sqrt(self.radius*self.radius - h*h)

    self.circles = self.determine_exclusions(prot, site_vol, include_metals)
    for i in range(len(self.circles)):
      for j in range(len(self.circles)):
        if(i == j): continue
  
        # we can break from the inner loop if a circle becomes irrelavent
        self.circles[i].remove_overlap(self.circles[j], i, j)
        if(not self.circles[i].full_circle and \
           len(self.circles[i].final_arcs) == 0): break

    # remove circles that do not contribute to the final solution
    kept_circles = []
    for C in self.circles:
      if(C.full_circle or len(C.final_arcs) > 0):
        kept_circles.append(C)
    self.circles = kept_circles
# NOTE there there may be references to the non contributing arcs and circles
# in iCircle.removed_arcs

################################################################################

################################################################################
  def determine_exclusions(self, prot, site_vol, include_metals=True, 
                           hbond_sphere_rad=3.0):

    # Determine exclusion spheres by
    # 1) Find all atoms within hbond_sphere_rad + 2.5 (A) of the heavy atom's 
    #    center
    # 2) Removing those atoms that fall outside of the bounding volume
    # 3) Checking if it is an exclusion sphere with respect to the hbonding
    #    model (volume)
    max_sq_dist = (hbond_sphere_rad + 2.5)**2
    min_sq_dist = 2.5*2.5

#    cgo_items.draw_plane2(self.ideal_dir, self.cap_plane_P0, plane_num=2,
#                          state=1, color=[0.8,0.2, 0.2])

    exclusion_spheres = []
    for res in prot:
      for a in res.atoms:
        # omit self
        if(a.serial == self.hbond_atom.serial): continue
        pos = array(a.position)
        # Ignore atoms >= 3.0 + 2.5 (A) of the heavy atom's center
        if(utils.squared_dist(pos, self.center) >= max_sq_dist): continue
        # Ignore atoms that fall too far outside of the site volume
# Turn this off for now        
#        if(not site_vol.intersects_sphere(pos, 2.5)): continue

        my_sphere = self.is_exclusion_sphere(pos)
        if(not my_sphere == None):
          exclusion_spheres.append(my_sphere)
          print a

    # handle metal atoms -- we allow them to get close ...
    if(include_metals):
      for M in prot.metals:
        # Ignore atoms >= 3.0 + min_metal_rad (A) of the heavy atom's center
        max_sq_dist = (M.min_act_rad + hbond_sphere_rad)**2
        pos = array(a.position)
        if(utils.squared_dist(pos, self.center) >= max_sq_dist): continue
        # Ignore atoms that fall too far outside of the site volume
# Turn this off for now        
#        if(not site_vol.intersects_sphere(pos, 2.5)): continue

        my_sphere = self.is_exclusion_sphere(pos, M.min_act_rad)
        if(not my_sphere == None):
          exclusion_spheres.append(my_sphere)

    return exclusion_spheres
################################################################################

################################################################################
  def is_exclusion_sphere(self, center, radius=2.5):
    """
    """
    sq_dist = utils.squared_dist(self.center, center)
    
    tmp = radius + self.radius
    if(sq_dist >= tmp*tmp): return None

    # debugging
#    global s_num
#    s_num += 1
#    cgo = [COLOR, 0.0, 0.8, 0.0, SPHERE]
#    cgo.extend(center)
#    cgo.append(radius)
#    cmd.load_cgo(cgo, "S%d" % (s_num), s_num)

    #########
    # determine the circle (represented by sphere) of intersection between 
    # self and the input sphere
    #########

    # direction from self's center other sphere's center
    I0_normal = utils.unit_vec(self.center, center)
    # distance between two spheres' centers
    r0_rs_dist = sqrt(sq_dist)
    # distance from self's center to center of circle of intersection
    d0 = 0.5 * (sq_dist + self.radius**2 - radius**2) / r0_rs_dist
    # compute center of circle of intersection
    I0_center = self.center + d0 * I0_normal
    I0_radius = sqrt(self.radius**2 - d0**2)

    # Determine if sphere of intersection is above, below or is cut by 
    # plane used to define the spherical cap
    V = I0_center - self.cap_plane_P0
    my_dist = sum(V * self.ideal_dir)
    # No intersection between cap & sphere
    if(-1.0 * my_dist >= I0_radius): return None
    # Cap & sphere intersect and plane defining cap does not cut
    # sphere of intersection
    elif(my_dist >= I0_radius): return iCircle(I0_center, I0_radius, I0_normal)
    # Cap & sphere intersect & plane defining cap is likely to cut 
    # the circle of intersection and does cut sphere of intersection
    else:
      # Get d for the plane equation I0_normal*X + d = 0
      d = -1.0 * sum(I0_normal * I0_center)
      # Plane intersection is the line Y = mX + b
      (m,b) = utils.intersect_planes(self.ideal_dir, self.cap_plane_d, 
                                     I0_normal, d)    
# debugging -- draw planes as well
#      cgo_items.draw_plane2(I0_normal, I0_center, plane_num=1, state=s_num)
#      cgo_items.draw_line(m, b, label="line%d" % (s_num), state=s_num)

      # Determine intersection between line defined by m & b and the sphere
      # of intersection
      (p0,p1) = utils.line_intersect_sphere(I0_center, I0_radius, m, b)
      if(p0.shape[0] == 0 or p1.shape[0] == 0): return None

      initial_arc = None
      # Circle of intersection is cut by cap plane -- find the remaining arc
      if(not p0.shape[0] == 0 and not p1.shape[0] == 0):
        # project midpoint of line intersecting the sphere to the circle
        chord_mid_pt = 0.5*(p0 + p1)    
        my_dir = utils.unit_vec(I0_center, chord_mid_pt)
        arc_mid_pt = I0_center + I0_radius*my_dir
        # Check if the projection is above/below the cap plane
        if(sum( (arc_mid_pt - self.cap_plane_P0) * self.ideal_dir ) <= 0.0):
          arc_mid_pt = I0_center - I0_radius*my_dir
        initial_arc = arc(I0_center, I0_radius, p0, arc_mid_pt, p1)
         
      return iCircle(I0_center, I0_radius, I0_normal, initial_arc)
################################################################################

################################################################################
  def closest_point(self, pt, tol=1.5):
    """
    Project pt onto the closest point on the cap.  If the distance is greater
    than the tolerance return an empty array
    """

    sq_tol = tol * tol
    #####
    # Project pt onto cap
    #####

    # Project pt onto sphere
    my_dir = utils.unit_vec(self.center, pt)
    proj_pt = pt + self.radius * my_dir
    if(utils.squared_dist(proj_pt, pt) > sq_tol): return array([])

    # Determine if sphere of intersection is above, below or is cut by 
    # plane used to define the spherical cap
    V = proj_pt - self.cap_plane_P0
    dist_to_plane = sum(V * self.ideal_dir)
    if(dist_to_plane < -1.0 * tol): return array([])

    # If point is below plane, project pt onto intersection of plane & sphere
    # (i.e. the circle of the cap in the plane)
    if(dist_to_plane < 0.0):
      # project pt to plane
      proj_pt = pt - dist_to_plane * self.ideal_dir
      if(utils.squared_dist(proj_pt, pt) > sq_tol): return array([])

      # project pt in plane to circle (sphere)
      in_plane_dir = utils.unit_vec(self.center, proj_pt)
      proj_pt += self.cap_plane_rad * in_plane_dir
      if(utils.squared_dist(proj_pt, pt) > sq_tol): return array([])
    
    # Here we need to step through the kept circles and determine if
    # the projected point falls inside any of them.
# NOTE: this is not necessarily correct, as in projecting to the closest point
# but it is reasonable enough and will allow us to proceed more rapidly.  In the
# future if we continue with this method we can adjust it, otherwise the extra
# work spent here may provide little additional gain (i.e. if we start allowing
# protein atoms to move we may wish to fall back to representing the cap by
# a triangular mesh rather than analytically)
    arc_pts = []
    sq_dists = []
    for C in self.circles:
      if(not C.contains(proj_pt)): continue

      # do standard projection to circle
      V = proj_pt - C.center
      dist_to_plane = sum(V * C.N)
      pt_in_plane = proj_pt - dist_to_plane * C.N
      U = utils.unit_vec(C.center, pt_in_plane)
      pt_on_C = C.center + C.radius * U

      sq_dist = utils.squared_dist(pt_on_C, pt) 
      if(sq_dist > sq_tol): continue
      elif(C.full_circle): 
        arc_pts.append(pt_on_C)
        sq_dists.append(sq_dist)

      # Test if projected point is inside arc.  If it is not, move it to the
      # closest end point of the arc
      elif(len(C.final_arcs)):
        for A in C.final_arcs:
          if(A.contains(pt_on_C)): 
            arc_pts.append(pt_on_C)
            sq_dists.append(sq_dist)
          else:
            sq_dists = [ utils.squared_dist(pt, A.pts[i]) for i in range(2) ]
            sq_dists = array(sq_dists)
            min_sq_dist = min(sq_dists)
            if(min_sq_dist <= sq_tol): 
              arc_pts.append(A.pts[argmin(sq_dists)])
              sq_dists.append(min_sq_dist)
      else:
        print "Found an empty circle -- not full and no final arcs"

    if(len(arc_pts)): return arc_pts[argmin(sq_dists)]
    else: return proj_pt
################################################################################

################################################################################
class hbond_surf_caps:
  """
  Hydrogen bonding surfaces approximated by points sampled on spherical cap for
  1 site and geometry for the second.
  """

  def __init__(self, prot, site_vol, hb_vols, include_metals=False):
    """
    Just pass in the hbond_volumes class for now
    """

  #def __init__(self, hbond_atom, C_nbr_atom, ideal_pt, prot, site_vol, 
  #             cap_number, include_metals):
    # Basically for each hbond_volume create an hbond cap
    caps = []
    for hb_vol_cap in hb_vols.caps:
      ideal_pt = hb_vol_cap.hbond_atom.position + 3.0 * hb_vol_cap.ideal_dir
      tmp = hbond_cap(hb_vol_cap.hbond_atom, hb_vol_cap.C_nbr_atom, 
                      ideal_pt, prot, site_vol,
                      hb_vol_cap.cap_number, include_metals)
      caps.append(tmp)
    self.caps = caps


#    """
#    for each hbond atom
#    0) Look up the atom in ideal_pts
#    1) get the neighbor atoms from prot
#    2) Get local coordinate system
#    3) compute the offsets in the local coordinate system
#    4) move the ideal_pt to global coordinate system
#    5) Put down points for cap
#    6) Find all neighbor atoms of hbond atom
#    7) Determine if enough points are left for the cap to be of use
#      a) Clip points falling outside of ligand box -- very fast
#      b) Clip points falling inside a neighbor atom
#
#    """
#
#
#    ideal_pts_fname = "/home/vanvoor4/code/SimSite3D/trunk/params/"
#    ideal_pts_fname += "new_optimum_hbonds.dat"
#    ideal_pts = SimSite3DPy.sitemaps.hbond_ideal_pts(ideal_pts_fname)
#    hbond_triplets = ideal_pts.get_polar_atom_triplets(prot, site_vol)

#    my_caps = []
#    for atom_triplet in hbond_triplets:
#      (C_nbr, hbond_atom, other_nbr) = atom_triplet
#      global_pts = \
#        ideal_pts.compute_global_pos(C_nbr, hbond_atom, other_nbr, prot)
#      for i in range(global_pts.shape[0]):
#        if(not site_vol.contains(global_pts[i])): continue
#
#        cap = hbond_volume(hbond_atom, C_nbr, global_pts[i], local_cap_pts, 
#                           prot, site_vol, i, include_metals)
#        if(cap.cap_pts.shape[0] / max_cap_pts >= min_fraction_pts):
#          my_caps.append(cap)
#
#    self.caps = my_caps
#
#
################################################################################

################################################################################

if(__name__ == "__main__"):

  import pymol
  import cPickle
  pymol.finish_launching()

  my_prot = utils.pdb.residues("/home/vanvoor4/data/new_sampling/adenines/proteins/1b38_atp_p.pdb") 
  #cmd.load("/home/vanvoor4/data/new_sampling/adenines/proteins/1b38_atp_p.pdb", "my_prot")

  ideal_pts_fname = "/home/vanvoor4/code/SimSite3D/trunk/params/"
  ideal_pts_fname += "new_optimum_hbonds.dat"
  ideal_pts = sitemaps.hbond_ideal_pts(ideal_pts_fname)
  site_vol = utils.volumes.box(lig_fname = "/home/vanvoor4/data/new_sampling/adenines/ligands/1b38_atp_l.mol2")
  #hbond_triplets = ideal_pts.get_polar_atom_triplets(my_prot, site_vol)

  pkl_fname = "/home/vanvoor4/data/test_SimSite3D_v3.4/adenines/pockets/1b38_atp_ADE.pkl"
  pkl_in = file(pkl_fname, "rb")
  my_hb_vols = cPickle.load(pkl_in)
  pkl_in.close()

  for res in my_prot:
    if(not res.resSeq == 83 and not res.resSeq == 82): continue
    for atom in res:
      if(atom.name == " N  " and res.resSeq == 83): LEU_83_N = atom
      if(atom.name == " O  " and res.resSeq == 82): LEU_82_O = atom
      if(atom.name == " C  " and res.resSeq == 82): LEU_82_C = atom

#  print LEU_83_N
#  global_pts = \
#    ideal_pts.compute_global_pos(LEU_82_C, LEU_83_N, LEU_82_O, my_prot)
#  print "ideal points:", global_pts


  #sphere = [SPHERE]
  #sphere.extend(global_pts[0])
  #sphere.append(0.25)
  #cmd.load_cgo(sphere, "LEU 83 N ideal point")
  #sphere = [SPHERE]
  #sphere.extend([6.60656315, 32.9410327, 6.59851815])
  #sphere.append(0.25)
  #cmd.load_cgo(sphere, "Salsa LEU 83 N ideal point")

#  for pt,i in zip(global_pts, range(len(global_pts))):
#    cap = hbond_cap(LEU_83_N, LEU_82_C, pt, my_prot, site_vol, i, False)


  my_surf_caps = hbond_surf_caps(my_prot, site_vol, my_hb_vols)
  print "number of surface caps:", len(my_surf_caps.caps)

  for cap,i in zip(my_surf_caps.caps, range(2, len(my_surf_caps.caps) + 2) ):
    sphere = [COLOR, 1.0, 1.0, 1.0, SPHERE]
    sphere.extend(cap.center)
    sphere.append(cap.radius)
    cmd.load_cgo(sphere, "cap_%d" % (i), state=i)

    cgo_items.draw_plane2(cap.ideal_dir, cap.cap_plane_P0, plane_num=i,
                          state=i, color=[0.8,0.2, 0.2])

    for C,j in zip(cap.circles, range(1, len(cap.circles) + 1) ):
      my_dir = utils.unit_vec(cap.center, C.center)
      my_dist = utils.dist(cap.center, C.center)
      my_dist += sqrt(2.5*2.5 - C.radius * C.radius)
      other_center = cap.center + my_dist * my_dir
      sphere = [COLOR, 0.0, 0.8, 0.0, SPHERE]
      sphere.extend(other_center)
      sphere.append(2.5)
      cmd.load_cgo(sphere, "nbr_%d_%d" % (i,j), state=i)

      for A,z in zip(C.final_arcs, range(1, len(C.final_arcs) + 1) ):
        cgo_items.draw_arc(A.center, A.radius, A.pts, A.mid_pt, 
                           label="ARC%d_%d_%d" % (i, j, z),
                           color=[0.2, 0.8, 0.2], state=i)


    
# to test we would pickle the caps and then draw the spheres, planes & arcs +
# neighbor atoms as spheres
