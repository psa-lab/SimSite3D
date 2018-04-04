import SimSite3DPy
from SimSite3DPy import utils
from numpy import *

################################################################################
def init_cap_points(spacing=0.5):
  x = array([1,0,0])
  tmp_pts = []
  #tmp_pts.extend(utils.spherical_grid(2.7, 0.3).tolist())
  #tmp_pts.extend(utils.spherical_grid(3.0, 0.3).tolist())
  #tmp_pts.extend(utils.spherical_grid(3.3, 0.3).tolist())
  for radius in linspace(2.5, 3.5, round(1.0 / spacing) + 1):
    tmp_pts.extend(utils.spherical_grid(radius, spacing).tolist())
    
  #tmp_pts.extend(utils.spherical_grid(3.0, 0.5).tolist())
  #tmp_pts.extend(utils.spherical_grid(3.5, 0.5).tolist())
  #tmp_pts.extend(utils.spherical_grid(2.5, 0.5).tolist())
  pts = []
  for i in range(len(tmp_pts)):
    p = array(tmp_pts[i])
    p = p / sqrt(sum(p*p))
    # angle between p and x axis must be <= 60 degrees
    if(sum(p*x) >= 0.5):
      pts.append(tmp_pts[i])
  return array(pts)
################################################################################

################################################################################
class hbond_volume:

  def __init__(self, hbond_atom, C_nbr_atom, other_nbr, ideal_pt, 
               local_cap_pts, residues, site_vol, cap_number, include_metals):
    self.cap_number = cap_number
    self.hbond_atom = hbond_atom
    self.C_nbr_atom = C_nbr_atom
    self.other_nbr = other_nbr
    self.center = array(hbond_atom.position)
    self.ideal_dir = utils.unit_vec(self.center, array(ideal_pt))
    self.r0 = 2.5
    self.r1 = 3.5
    self.cap_pts = self.compute_global_cap_points(hbond_atom, C_nbr_atom, 
                                                  ideal_pt, local_cap_pts)
    self.exclusion_spheres = \
      self.determine_exclusions(residues, site_vol, include_metals, 3.5)
    self.cap_pts = \
      self.cull_cap_points(self.cap_pts, site_vol, self.exclusion_spheres)

################################################################################

################################################################################
  def determine_exclusions(self, residues, site_vol, include_metals=True,
                           hbond_sphere_rad=3.0):

    # Determine exclusion spheres by
    # 1) Find all atoms within hbond_sphere_rad + 2.5 (A) of the heavy atom's 
    #    center
    # 2) Removing those atoms that fall outside of the bounding volume
    # 3) Checking if it is an exclusion sphere with respect to the hbonding
    #    model (volume)
    max_sq_dist = (hbond_sphere_rad + 2.5)**2
    min_sq_dist = 2.5*2.5

    #nearby_atoms = {}
    exclusion_spheres = []
    for res in residues:
      for a in res.atoms:
        # omit self
        if(a.serial == self.hbond_atom.serial): continue
        pos = array(a.position)
        # Ignore atoms >= hbond_sphere_rad + 2.5 (A) of the heavy atom's center
        if(utils.squared_dist(pos, self.center) >= max_sq_dist): continue
        # Ignore atoms that fall too far outside of the site volume
        if(not site_vol.intersects_sphere(pos, 2.5)): continue

        if(self.is_exclusion_sphere(pos)):
          exclusion_spheres.append(utils.volumes.ball(pos, 2.5))

    # handle metal atoms -- we allow them to get close ...
    if(include_metals):
      for M in residues.metals:
        # Ignore atoms >= hbond_sphere_rad + min_metal_rad (A) of the heavy 
        # atom's center
        max_sq_dist = (M.min_act_rad + hbond_sphere_rad)**2
        pos = array(a.position)
        if(utils.squared_dist(pos, self.center) >= max_sq_dist): continue
        # Ignore atoms that fall too far outside of the site volume
        if(not site_vol.intersects_sphere(pos, 2.5)): continue

        if(self.is_exclusion_sphere(pos, M.min_act_rad)):
          exclusion_spheres.append(utils.volumes.ball(pos, M.min_act_rad))

    return exclusion_spheres
################################################################################

################################################################################
  def is_exclusion_sphere(self, center, radius=2.5):
    """
    Determine if the sphere with given center and radius 2.5 (A) impinges on the
    hbonding volume defined by the two concentric spheres (2.5 and 3.5 (A) with 
    the hbond_atom center as the spheres' centers) and nappe with apex at the
    hbond_atom's center and the axis in the direction of the ideal_dir.  Any
    ray origninating at the apex and lying on the nappe's surface makes an
    angle of 60 degrees with the axis.
    """

    sq_dist = utils.squared_dist(self.center, center)

# this is wrong -- we need to reduce the radius from 2.5 to ... and consider
# the metal atom as a possible exclusion sphere -- modeling of metals needs
# to be done differently from hbond spheres since many things will change and
# simplify (at least in a sense).
#    # If the hbonding volume is being used to approximate a metal atom, then
#    # our model uses the entire volume between the 2 concentric spheres
#    if(metal_rad > 0.0): 
#      tmp = metal_rad + self.r1       
#      if(sq_dist >= tmp*tmp): return False
#      return True

    tmp = radius + self.r1
    if(sq_dist >= tmp*tmp): return False
    tmp = max(radius, self.r1) 
    
    # The exclusion sphere impinges upon the hbond volume sphere 
    # (self.center, self.r1) but the center of the exclusion sphere is not 
    # inside the hbond volume sphere
    #
    # In this case we can (within numerical accuracy) exactly compute whether
    # the volume exclude by the exclusion sphere has a nontrival intersection
    # with the hbonding volume.  This is done by determining the circle of
    # intersection between the hbonding volume sphere and the exclusion sphere.
    # If the circle of intersection intersects with the circle defined by the
    # hbonding volume sphere and cone, then the sphere is an exclusion sphere.  
    if(sq_dist >= tmp*tmp):
      # Find the cosine of the angle between the axes of the two cones
      tmp = center - self.center
      radii_dist = sqrt(sq_dist)
      my_dir = tmp / radii_dist
      cos_axes_angle = sum(my_dir * self.ideal_dir)

      # Find the cosine of the sum of the angles the cone surfaces make with
      # their respective axis.
      # For now the hbonding angle is 60 degrees, and we denote the second
      # cone angle as alpha.  Use the cosine sum rule rewrite the consine
      # sum as two terms.  Use the law of cosines to compute cos(alpha)
      cos_alpha = radius*radius - self.r1*self.r1 - sq_dist
      cos_alpha /= -2.0 * self.r1 * radii_dist
      sin_alpha = sqrt(1.0 - cos_alpha*cos_alpha)

      cos_sum_cone_angles = 0.5*cos_alpha - (sqrt(3) / 2.0) * sin_alpha
      if(cos_axes_angle <= cos_sum_cone_angles): return False
     
    # In this case the center of the exclusion sphere is inside the
    # hbond volume sphere (defined by self.center, self.r1).
    #
    # This case is approximated by requiring the exclusion sphere's center
    # to be no more than some signed distance from the plane defined by
    # the normal in the direction of the nappe's axis (pointing away from the
    # apex) and the apex of the cone as a point on the plane.  The distance
    # from the plane depends on the radius of the exclusion sphere and the
    # smaller concentric sphere used to define the hbond volume
    else:
      # plane normal is self.ideal_dir
      # point on plane is self.center

      center = array(center)
      my_dir = center - self.center
      signed_dist = sum(my_dir * self.ideal_dir)
      if(signed_dist <= -(radius - 0.5*self.r0)): return False

    return True
################################################################################

################################################################################
  def compute_global_cap_points(self, hbond_atom, C_nbr_atom, ideal_pt,
                                local_cap_pts):
    """
    Note: this method assumes a symmetric volume representation (i.e.
    it was designed on the model of 2 concentric spheres and a nappe.).

    Purpose:
       Get the rigid body transformation to move the cap from the local
       coordinates to the global coordinates defined by the ideal_pt,
       C_nbr_atom, and hbond_atom.  Do this by mapping the hbond_atom to
       the origin, ideal_pt to the positive x-axis, and the C_nbr_atom to the
       lower half plane.
    """

    # Determine the rigid body transformation
    hb_atom_pos = array(hbond_atom.position)
    C_nbr_pos = array(C_nbr_atom.position)
    ideal_pos = array(ideal_pt)
    (R,T) = utils.get_local_coords(ideal_pos, hb_atom_pos, C_nbr_pos)

    # Shift the cap representation to the global position
    global_pts = dot(local_cap_pts, R)
    global_pts += tile(T, (local_cap_pts.shape[0], 1))

    return global_pts
################################################################################

################################################################################
  def cull_cap_points(self, cap_pts, site_vol, exclusion_spheres):
    npts = cap_pts.shape[0]
    flags = [ 1 for i in range(npts) ]

    # 1) Remove all points falling ouside of the site volume
    for i in range(npts):
      if(not site_vol.contains(cap_pts[i])): flags[i] = 0
    if(sum(flags) <= 0): return array([])

    # 2) Remove all points inside an exclusion sphere
    for i in range(npts):
      if(not flags[i]): continue
      for Xball in exclusion_spheres:
        if(Xball.contains(cap_pts[i])):
          flags[i] = 0
          break
    if(sum(flags) <= 0): return array([])

    tmp_pts = []
    for i in range(npts):
      if(flags[i]): tmp_pts.append(cap_pts[i].tolist())
    return array(tmp_pts)
    
################################################################################

################################################################################
#  def old_compute_global_cap_points(self, local_cap_pts, ideal_pts, hbond_atom, 
#                                    C_nbr_pos, min_pt, max_pt, residues):
#    point_clouds = []
#    neighbors = [] 
#    required_fraction_of_points = 0.1
#  
#    hvy_atom_pos = array(hbond_atom.position)
#    C_nbr_dir = utils.unit_vec(hvy_atom_pos, C_nbr_pos);
#    C_nbr_dist = utils.dist(C_nbr_pos, hvy_atom_pos)
#    for i in range(ideal_pts.shape[0]):
#      # remove any ideal points that do not fall inside the box
#      keep = True
#      for j in range(3):
#        if(ideal_pts[i][j] < min_pt[j] or max_pt[j] < ideal_pts[i][j]):
#          keep = False
#      if(not keep):
#        point_clouds.append(array([]))
#        neighbors.append([])
#        continue
#  
#      pt_pos = ideal_pts[i]
#      pt_dir = utils.unit_vec(hvy_atom_pos, pt_pos)
#  
#      cos_angle = sum(C_nbr_dir*pt_dir)
#      sin_angle = sin(arccos(cos_angle))
#  
#      # Get LSE fit
#      ideal_pt_len = utils.dist(pt_pos, hvy_atom_pos)
#      local_coords = [ [ideal_pt_len, 0, 0], [0, 0, 0],
#                       [cos_angle*C_nbr_dist, sin_angle*C_nbr_dist, 0.0] ]
#      local_coords = array(local_coords)
#      global_coords = array([ pt_pos, hvy_atom_pos, C_nbr_pos ])
#      (R,T) = utils.simple_lse_fit(global_coords, local_coords)
#  
#      # Shift the cap representation to the global position
#      global_pts = dot(local_cap_pts, R)
#      global_pts += tile(T, (local_cap_pts.shape[0], 1))
#  
#      pts_flag = [ 1 for ii in range(global_pts.shape[0]) ]
#      for j in range(global_pts.shape[0]):
#        for k in range(3):
#          if(global_pts[j][k] < min_pt[k] or max_pt[k] < global_pts[j][k]):
#            #print min_pt[k], global_pts[j][k], max_pt[k]
#            pts_flag[j] = 0
#
#      if(sum(pts_flag) == 0):
#        point_clouds.append(array([]))
#        neighbors.append([])
#        continue
#  
#      # Find all atoms within 3.5 + 2.5 (A) of the heavy atom's center
#      max_sq_dist = 36.0
#      min_sq_dist = 2.5*2.5
#  
#      print
#      print hbond_atom.serial, ideal_pts[i]
#      nearby_atoms = {}
#      for res in residues:
#        for a in res.atoms:
#          # omit self
#          if(a.serial == hbond_atom.serial): continue
#          pos = array(a.position)
#          if(utils.squared_dist(pos, hvy_atom_pos) > max_sq_dist): continue
#  
#          # Flag bumped points
#          for j in arange(global_pts.shape[0]):
#            if(pts_flag[j] == 0): continue
#            if(utils.squared_dist(pos, global_pts[j]) < min_sq_dist):
##              nearby_atoms[a.serial] = a
#              pts_flag[j] = 0
#  
#      for a_num in nearby_atoms.keys():
#        print a_num,
#      print
#  
#      tmp_pts = []
#      for j in arange(global_pts.shape[0]):
#        if(pts_flag[j]): tmp_pts.append(global_pts[j].tolist())
#      tmp_pts = array(tmp_pts)
#
#      if(float(tmp_pts.shape[0]) / float(local_cap_pts.shape[0]) >= \
#         required_fraction_of_points):
#        point_clouds.append(array(tmp_pts))
#        neighbors.append([ nearby_atoms[a_num] for a_num in nearby_atoms.keys() ])
#      else:
#        point_clouds.append(array([]))
#        neighbors.append([])
#  
#    return (point_clouds, neighbors)
################################################################################

###############################################################################
class hbond_volumes:
  """
  Hydrogen bonding volumes approximated by pointsets for 1 site and 
  geometry for the second.
  """
  def __init__(self, prot_fname, site_vol, min_fraction_pts=0.10, 
               include_metals=False, cap_point_spacing=0.5):

    """
    Initialize the sitemap either from reading the xml_fname or if it doesn't
    exist, by creating the sitemap from the SimSite3D V3.3 sitemap files

    for each hbond atom
    0) Look up the atom in ideal_pts
    1) get the neighbor atoms from prot
    2) Get local coordinate system
    3) compute the offsets in the local coordinate system
    4) move the ideal_pt to global coordinate system
    5) Put down points for cap
    6) Find all neighbor atoms of hbond atom
    7) Determine if enough points are left for the cap to be of use
      a) Clip points falling outside of ligand box -- very fast
      b) Clip points falling inside a neighbor atom

    """
    prot = utils.pdb.residues(prot_fname)
    ideal_pts_fname = "/home/vanvoor4/code/SimSite3D_surfaces/params/"
    ideal_pts_fname += "new_optimum_hbonds.dat"
    ideal_pts = SimSite3DPy.sitemaps.hbond_ideal_pts(ideal_pts_fname)
    hbond_triplets = ideal_pts.get_polar_atom_triplets(prot, site_vol)

    my_caps = []
    local_cap_pts = init_cap_points(spacing=cap_point_spacing)
    max_cap_pts = float(local_cap_pts.shape[0])
    for atom_triplet in hbond_triplets:
      (C_nbr, hbond_atom, other_nbr) = atom_triplet
      global_pts = \
        ideal_pts.compute_global_pos(C_nbr, hbond_atom, other_nbr, prot)
      for i in range(global_pts.shape[0]):
        if(not site_vol.contains(global_pts[i])): continue

        cap = hbond_volume(hbond_atom, C_nbr, other_nbr, global_pts[i], 
	                   local_cap_pts, prot, site_vol, i, include_metals)
        if(cap.cap_pts.shape[0] / max_cap_pts >= min_fraction_pts):
          my_caps.append(cap)
    if(include_metals):
      print "NEED to implement model for metal atoms"

    if(include_metals):
      print "NEED to implement model for metal atoms"

    self.caps = my_caps
################################################################################

###############################################################################
def lig_bound_vol(lig, tol=2.0):
  """
  returns the bounding box of the ligand -- ([min pt], [max pt])
  """
  min_pt = array(lig.atoms[0].position)
  max_pt = array(lig.atoms[0].position)

  for atom in lig.atoms:
    for i in range(3):
      if(atom.position[i] < min_pt[i]): min_pt[i] = atom.position[i]
      if(max_pt[i] < atom.position[i]): max_pt[i] = atom.position[i]

  return utils.volumes.box(min_pt - tol, max_pt + tol)
###############################################################################
