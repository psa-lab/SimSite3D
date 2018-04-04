import ASCbasePy
from ASCbasePy import utils
from numpy import *

class metal_info:

  def __init__(self, pdb_metal_str, metal_type, min_rad, avg_rad, max_rad):
    self.pdb_metal_str = pdb_metal_str
    self.metal_type = metal_type
    self.min_rad = min_rad
    self.avg_rad = avg_rad
    self.max_rad = max_rad
################################################################################
  
################################################################################
class metals:
 
  def __init__(self):
    self.distances = {}
    self.distances["CA  "] = metal_info("CA  ", 1, 2.0, 2.4, 2.9)
    self.distances["CO  "] = metal_info("CO  ", 2, 1.7, 1.9, 2.6)
    self.distances["CU  "] = metal_info("CU  ", 2, 1.7, 2.1, 2.6)
    self.distances["FE  "] = metal_info("FE  ", 2, 1.7, 2.2, 2.6)
    self.distances[" K  "] = metal_info(" K  ", 1, 2.0, 2.4, 2.9)
    self.distances["MG  "] = metal_info("MG  ", 2, 1.7, 2.1, 2.6)
    self.distances["MN  "] = metal_info("MN  ", 2, 1.7, 2.2, 2.6)
    self.distances["NA  "] = metal_info("NA  ", 1, 2.0, 2.4, 2.9)
    self.distances["NI  "] = metal_info("NI  ", 2, 1.7, 2.2, 2.6)
    self.distances["ZN  "] = metal_info("ZN  ", 2, 1.7, 2.1, 2.6)
################################################################################
  
################################################################################
def init_metal_points(min_rad, avg_rad, max_rad, spacing=0.5):
  tmp_pts = []
  tmp_pts.extend(utils.spherical_grid(max_rad, spacing).tolist())
  tmp_pts.extend(utils.spherical_grid(avg_rad, spacing).tolist())
  tmp_pts.extend(utils.spherical_grid(min_rad, spacing).tolist())
  return array(tmp_pts)
################################################################################

################################################################################
class metal_volume:

  def __init__(self, pdb_metal, prot, site_vol):
#    my_metals = metals()
#    m_dists = my_metals.distances[metal_atom.name]
    #metal_info = []
    #for met in prot.pdb_metals.table:
      #if(met.pdb_metal_name == pdb_metal.name):
        #metal_info = met
        #break;
    pdb_metals = utils.pdb.metal_lookup()
    metal_info = pdb_metals.table[pdb_metal.name]
    local_pts = init_metal_points(metal_info.min_act_rad, metal_info.pt_rad,
                                  metal_info.max_act_rad)
    self.max_npts = float(local_pts.shape[0])

    self.metal_atom = pdb_metal
    self.center = array(pdb_metal.position)
    self.min_rad = metal_info.min_act_rad
    self.max_rad = metal_info.max_act_rad
       
    # translate all the points to be centered at the center of the metal atom
    self.vol_pts = local_pts + tile(self.center, (local_pts.shape[0], 1))
    self.exclusion_spheres = \
      self.determine_exclusions(prot, site_vol)
    self.vol_pts = \
      self.cull_cap_points(self.vol_pts, site_vol, self.exclusion_spheres)

################################################################################

################################################################################
  def determine_exclusions(self, prot, site_vol):

    # Determine exclusion spheres by
    # 1) Find all atoms within 3.5 + 2.5 (A) of the heavy atom's center
    # 2) Removing those atoms that fall outside of the bounding volume
    # 3) Checking if it is an exclusion sphere with respect to the hbonding
    #    model (volume)
    max_sq_dist = (2.5 + self.max_rad)**2

    exclusion_spheres = []
    for res in prot.residues:
      for a in res.atoms:
        # omit self
        if(a.serial == self.metal_atom.serial): continue
        pos = array(a.position)
        # Ignore atoms >= 3.5 + 2.5 (A) of the heavy atom's center
        if(utils.squared_dist(pos, self.center) >= max_sq_dist): continue
        # Ignore atoms that fall too far outside of the site volume
        if(not site_vol.intersects_sphere(pos, 2.5)): continue
   
#        if(self.is_exclusion_sphere(pos, max_metal_rad)): 
        exclusion_spheres.append(utils.volumes.ball(pos, 2.5))
    return exclusion_spheres
    
################################################################################

################################################################################
  def is_exclusion_sphere(self, nbr_center, nbr_rad=2.5):
    """
    Determine if the sphere with given center and radius 2.5 (A) impinges on the
    aetal volume defined by the two concentric spheres (min_metal_rad and 
    max_metal_rad (A) with the metal_atom's center as the spheres' centers) 
    """

    sq_dist = utils.squared_dist(self.center, nbr_center)

    # If the hbonding volume is being used to approximate a metal atom, then
    # our model uses the entire volume between the 2 concentric spheres
    tmp = self.max_rad + nbr_rad
    if(sq_dist >= tmp*tmp): return False
    return True
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
class metal_volumes:
  """
  metal volumes approximated by pointsets for 1 site and 
  geometry for the second.
  """
  def __init__(self, prot_fname, site_vol, min_fraction_pts=0.10):

    """
    For now we assume that site_vol has an additional tolerance of 2.5
    so that we can check the center of the metal (we are assuming metals do
    not have a preferred direction)

    for each metal atom
    0) Look up the atom distances
    1) get the neighbor atoms from prot
    2) Translate sphere to the metal's center
    3) Find all neighbor atoms of hbond atom
    4) Determine if enough points are left for the cap to be of use
      a) Clip points falling outside of ligand box -- very fast
      b) Clip points falling inside a neighbor atom

    """
    prot = utils.pdb.residues(prot_fname)
    my_metals = metals()

    my_spheres = []
    for metal in prot.metals:
#     Need to check to make sure the metal is close enough -- note we want
#     to add the interaction metal radius to check if we keep it. and not
#     use the metal atom center
      if(not site_vol.contains(metal.position)): continue

#      m_dists = my_metals.distances[metal.name]
#      local_pts = init_metal_points(m_dists.min_rad, m_dists.avg_rad, 
#                                    m_dists.max_rad):
#      max_npts = float(local_pts.shape[0])

# translate all the points     
#      global_pts = local_pts + tile(metal.position, (max_npts, 1))

#def __init__(self, metal_atom, vol_pts, residues, site_vol, min_rad, max_rad):
      sphere = metal_volume(metal, prot, site_vol)

      if(sphere.vol_pts.shape[0] / sphere.max_npts >= min_fraction_pts):
        my_spheres.append(sphere)

    self.metals = my_spheres
################################################################################

###############################################################################
