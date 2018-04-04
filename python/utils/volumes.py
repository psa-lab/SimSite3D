from SimSite3DPy.utils._math import *
from SimSite3DPy.utils import mol2
from numpy import *
    
################################################################################
class ball:
  """
  A spherical volume
  """

  def __init__(self, center, radius):
    """
    center : center of sphere
    radius : radius of sphere
    """
    self.center = array(center)
    self.radius = radius
################################################################################
    
################################################################################
  def contains(self, pt, tol=0.0):
    radius = self.radius + tol
    if(squared_dist(self.center, array(pt)) > radius*radius): 
      return False
    return True
################################################################################
    
################################################################################
  def intersects_sphere(self, center, radius):
    """
    Test if another sphere intersects with this sphere
    """
    tmp = self.radius + radius
    if(squared_dist(self.center, array(center)) > tmp*tmp): return False
    return True
################################################################################

################################################################################
class UnionOfBalls:
  """
  Volume occupied by the union of a set of balls (spheres, atoms, whatever)
  """

  def __init__(self, centers, radii):
    self.centers = array(centers)
    self.radii = array(radii)
    self.squared_radii = array(radii) * array(radii)
    #self.bind_radius = 5.0
    #self.RAD_radius = 9.0

  def contains(self, pt, tol=0.0):
    if(tol == 0.0):
      for center, squared_rad in zip(self.centers, self.squared_radii):
        if(squared_dist(center, pt) <= squared_rad): return True
    else:
      for center, radius in zip(self.centers, self.radii):
        tmp = radius + tol
        if(squared_dist(center, pt) <= tmp*tmp): return True

    return False

  def intersects_sphere(self, s_center, s_radius):
    for center, radius in zip(self.centers, self.radii):
      tmp = radius + s_radius
      if(squared_dist(center, s_center) <= tmp*tmp): return True
    return False

################################################################################

################################################################################
class box:
  """
  A rectangular solid (volume)
  """

  def __init__(self, min_pt=array([]), max_pt=array([]), lig_fname=""):
    """
    min_pt : minimum corner of the box
    max_pt : maximum corner of the box
    """
    if(lig_fname.endswith("mol2")): 
      lig = mol2.molecule(lig_fname)
      (min_pt, max_pt) = self.lig_bound_vol(lig)

    self.min_pt = array(min_pt)
    self.max_pt = array(max_pt)
################################################################################

################################################################################
  def contains(self, pt, tol=0.0):
    min_pt = self.min_pt - tol
    max_pt = self.max_pt + tol
    for i in range(3):
      if(pt[i] < min_pt[i] or max_pt[i] < pt[i]): return False
    return True 
################################################################################

################################################################################
  def intersects_sphere(self, center, radius):
    # 0) Check if center of sphere is in the box
    if(self.contains(center)): return True 
    
    # 1) Determine if the center of the sphere is outside of the box plus a 
    #    tolerance of the sphere's radius added to each dimension
    if(not self.contains(center, tol=radius)): return False

    # 2) Determine the face(s) we need to check
    sign = numpy.zeros((3,))
    for i in range(3):
      if(center[i] <= self.min_pt[i]): sign[i] = -1
      elif(self.max_pt[i] <= center[i]): sign[i] = 1

    # 3) Center of the sphere is one of the boxes defined by each face and 
    #    pushing out the face perpendicuarly by the radius of the sphere
    if((sign[0] and not sign[1] and not sign[2]) or
       (not sign[0] and sign[1] and not sign[2]) or
       (not sign[0] and not sign[1] and sign[2])): return True

    # 4) Does the sphere contain the one of the corners of the box?
    if(sign[0] and sign[1] and sign[2]):
      my_dist2 = 0.0
      for i in range(3):
        if(sign[i] > 0): tmp = self.max_pt[i] - center[i]
        else: tmp = self.min_pt[i] - center[i]
        my_dist2 += tmp*tmp
      if(radius*radius >= my_dist2): return True
      else: return False

    # 5) Does the sphere intersect one of the boxes edges?
    if(not sign[0]): return self._line_in_sphere(radius, center, 0, sign)    
    if(not sign[1]): return self._line_in_sphere(radius, center, 1, sign)    
    if(not sign[2]): return self._line_in_sphere(radius, center, 2, sign)    

  def _line_in_sphere(self, radius, center, idx, sign):
    """
    Determine the point on the line segment that is closest to the center of
    the sphere
    """
    pt = zeros((3,))
    RHS = 0.0
    for i in range(3):
      if(sign[i] > 0.0): pt[i] = self.max_pt[i]
      if(sign[i] < 0.0): pt[i] = self.min_pt[i]
      if(sign[i]): RHS -= (pt[i] - center[i])**2

    if(RHS < 0.0): return False

    RHS = sqrt(RHS)
    S1 = center[idx] - RHS
    S2 = center[idx] + RHS
    if((S1 < self.min_pt[idx] and S2 < self.min_pt[idx]) or \
       (S1 > self.max_pt[idx] and S2 > self.max_pt[idx])): return False
    return True
###############################################################################

###############################################################################
  def lig_bound_vol(self, lig, tol=2.0):
    """
    returns the bounding box of the ligand -- ([min pt], [max pt])
    """
    min_pt = array(lig.atoms[0].position)
    max_pt = array(lig.atoms[0].position)
  
    for atom in lig.atoms:
      for i in range(3):
        if(atom.position[i] < min_pt[i]): min_pt[i] = atom.position[i]
        if(max_pt[i] < atom.position[i]): max_pt[i] = atom.position[i]
  
    return (min_pt - tol, max_pt + tol)
###############################################################################




