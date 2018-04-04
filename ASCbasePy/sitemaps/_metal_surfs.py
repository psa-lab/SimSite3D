from numpy import *
from _hbond_surf_caps import arc, iCircle
from _metal_volumes import *


class metal_surf:

  def __init__(self, metal_atom, prot, metal_site_vol):
    pdb_metals = utils.pdb.metal_lookup()
    metal_info = pdb_metals.table[metal_atom.name]

    self.center = array(metal_atom.position)
    self.radius = metal_info.pt_rad
    self.metal_atom = metal_atom

    (close_atoms, self.circles) = \
      self.determine_exclusions(prot, metal_site_vol)

    # metals are typically ligated, we aren't interested in metals floating
    # in the solution
    if(len(close_atoms) < 2):
      self.close_atoms = []
      return

    # Determine the two closest neighbors to the metal to get its orientation
    close_atoms.sort(cmp=self.__atom_dist_cmp)
    self.close_atoms = [close_atoms[0][0], close_atoms[1][0]]

    # determine the overlap between the exclusion circles
    for i in range(len(self.circles)):
      for j in range(len(self.circles)):
        if(i == j): continue
    
        self.circles[i].remove_overlap(self.circles[j], i, j)
        if(not self.circles[i].full_circle and \
           len(self.circles[i].final_arcs) == 0): break

    # remove circles that do not contribute to the final solution
    kept_circles = []
    print "KEPT"
    for C in self.circles:
      if(C.full_circle or len(C.final_arcs) > 0):
        print C.center, C.radius 
        kept_circles.append(C)
    self.circles = kept_circles
    print "num circles", len(self.circles)
################################################################################

################################################################################
  def determine_exclusions(self, prot, site_vol):
    max_sq_dist = (self.radius + 2.5)**2

    close_atoms = []
    exclusion_spheres = []
    for res in prot.residues:
      for a in res.atoms:
        if(a.serial == self.metal_atom.serial): continue
        other_center = array(a.position)

        # Ignore atoms farther than the average interaction radius + 2.5 (A)
        sq_dist = utils.squared_dist(other_center, self.center)
        if(sq_dist >= max_sq_dist): continue
        close_atoms.append((a, sq_dist))

        # compute the center & radius of iCircle, & normal
        print a

        # Direction from self to the exclusion sphere
        I0_normal = utils.unit_vec(self.center, other_center)
        my_dist = sqrt(sq_dist)
        print "dist", my_dist
        # distance from self's center to center of the circle of intersection
        d0 = 0.5 * (sq_dist + self.radius**2 - 2.5**2) / my_dist
        # Compute the center & radius of the circle of intersection
        I0_center = self.center + d0 * I0_normal
        I0_radius = sqrt(self.radius**2 - d0**2)
        exclusion_spheres.append(iCircle(I0_center, I0_radius, I0_normal))
        print "center", I0_center
        print "radius", I0_radius

    for m in prot.metals:
      if(m.serial == self.metal_atom.serial): continue
      other_center = array(m.position)

      # Ignore atoms farther than self's average radius + other metal's
      # min radius
      pdb_metals = utils.pdb.metal_lookup()
      metal_info = pdb_metals.table[m.name]

      sq_dist = utils.squared_dist(other_center, self.center)
      if(sq_dist >= self.radius**2 + metal_info.min_act_rad): continue
      close_atoms.append((m, sq_dist))

      # Direction from self to the exclusion sphere
      I0_normal = utils.unit_vector(self.center, other_center)
      my_dist = sqrt(sq_dist)
      # distance from self's center to center of the circle of intersection
      d0 = 0.5 * (sq_dist + self.radius**2 - 2.5**2) / my_dist
      # Compute the center & radius of the circle of intersection
      I0_center = self.center + d0 * I0_normal
      I0_radius = sqrt(self.radius**2 - d0**2)
      exclusion_spheres.append(iCircle(I0_center, I0_radius, I0_normal))

    return (close_atoms, exclusion_spheres)
################################################################################

################################################################################
  def __atom_dist_cmp(self, A,B):
    if(A[1] < B[1]): return -1
    elif(A[1] == B[1]): return 0
    return 1
################################################################################

################################################################################
class metal_surfs:

  def __init__(self, prot, metal_site_vol, metal_vols):
    self.spheres = []
   
    for metal_vol in metal_vols.metals:
      tmp = metal_surf(metal_vol.metal_atom, prot, metal_site_vol)
      if(len(tmp.close_atoms) > 1): self.spheres.append(tmp)

    for s in self.spheres:
      print "num circles:", len(s.circles)
################################################################################

################################################################################
