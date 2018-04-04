from ASCbasePy import *
import numpy

###############################################################################
class my_ideal_pt:
  def __init__(self, line):
    toks = line.split(",")
    self.atoms = [toks[6], toks[1], toks[7]]
    self.obj = toks[2]
    self.pt_num = int(toks[3])
    self.alpha = float(toks[4])
    self.beta = float(toks[5])
###############################################################################

###############################################################################
# returns ([min pt], [max pt])
def lig_bound_vol( lig, tol=2.0):
  min_pt = lig.atoms[0].position[:]
  max_pt = lig.atoms[0].position[:]

  for atom in lig.atoms:
    for i in range(3):
      if(atom.position[i] < min_pt[i]): min_pt[i] = atom.position[i]
      if(max_pt[i] < atom.position[i]): max_pt[i] = atom.position[i]

  for i in range(3):
    min_pt[i] -= tol 
    max_pt[i] += tol
  
  return (min_pt, max_pt)
###############################################################################

###############################################################################
# get hbond points near the volume
def get_prot_polar_atoms(residues, min_pt, max_pt, tol=3.5):
  p_min_pt = min_pt[:]
  p_max_pt = max_pt[:]
  for i in range(3):
    p_min_pt[i] -= tol
    p_max_pt[i] += tol

  polar_atoms = [] 
  for res in residues:
    for atom in res.atoms:
      if(atom.interact_type == "ACCEPTOR" or atom.interact_type == "DONOR" or \
         atom.interact_type == "DONEPTOR"):

        keep_atom = True
        for i in range(3):
          if(atom.position[i] < p_min_pt[i] or p_max_pt[i] < atom.position[i]):
            keep_atom = False

        if(keep_atom): polar_atoms.append(atom)
  return polar_atoms
###############################################################################

###############################################################################
def read_ideal_pts(fname):
  infile = file(fname, "r")

  pts = {}

  for line in infile:
    if(line.startswith("#")): continue

    toks = line.split(",")
    if(len(toks) == 1): continue

    if(not toks[0] in pts): pts[toks[0]] = {}
    if(not toks[1] in pts[toks[0]]): pts[toks[0]][toks[1]] = []
    pts[toks[0]][toks[1]].append(my_ideal_pt(line))


  # main chain trash
  res_names = [ "ALA", "ARG", "ASP", "ASN", "CYS", "GLU", "GLN", "GLY", "HIS",
                "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TYR", 
                "TRP", "VAL" ]

  for res_name in res_names:
    if(not res_name in pts): pts[res_name] = {}
    pts[res_name]["O"] = pts["MAIN_CHAIN"]["O"]
    if(not res_name == "PRO"):
      pts[res_name]["N"] = pts["MAIN_CHAIN"]["N"]
    
 
  return pts
###############################################################################

###############################################################################
def compute_offset(bond_len, alpha, beta):
  alpha = utils.deg2rad(alpha) 
  beta = utils.deg2rad(beta) 
  cos_alpha = numpy.cos(alpha) 
  cos_beta = numpy.cos(beta) 
  sin_alpha = numpy.sin(alpha) 
  sin_beta = numpy.sin(beta) 
  return [bond_len * cos_alpha, bond_len * cos_beta * sin_alpha, 
          bond_len * sin_beta * sin_alpha]
###############################################################################

###############################################################################
def get_hbond_atom_nbrs(residues, hbond_atom, C_name, other_name):
  # Not sure about the "python" way to do this, but ...
  for i in range(len(residues)):
    if(residues[i]._residue__contains(hbond_atom)):
      if(hbond_atom.name == " O  "):
        #print "carbonyl O"
        for tmp_atom in residues[i]:
          if(tmp_atom.name.strip() == C_name): C_nbr = tmp_atom

        if(i+1 < len(residues)):
          for tmp_atom in residues[i+1]:
            if(tmp_atom.name.strip() == other_name): other_nbr = tmp_atom

           # distance check here
      elif(hbond_atom.name == " N  "):
        #print "amide N"
        #print my_ideal_pts[0].atoms
        if(i > 0):
          for tmp_atom in residues[i-1]:
            if(tmp_atom.name.strip() == C_name): C_nbr = tmp_atom
            elif(tmp_atom.name.strip() == other_name): other_nbr = tmp_atom
           # distance check here
      else:
        #print "other atom"
        for tmp_atom in residues[i]:
          if(tmp_atom.name.strip() == C_name): C_nbr = tmp_atom
          elif(tmp_atom.name.strip() == other_name): other_nbr = tmp_atom
      break
  return (C_nbr, other_nbr)
###############################################################################

###############################################################################
def init_cap_points(spacing):
  x = numpy.array([1,0,0])
  tmp_pts = []
  #tmp_pts.extend(utils.spherical_grid(2.7, 0.3).tolist())
  #tmp_pts.extend(utils.spherical_grid(3.0, 0.3).tolist())
  #tmp_pts.extend(utils.spherical_grid(3.3, 0.3).tolist())
  tmp_pts.extend(utils.spherical_grid(3.0, 0.5).tolist())
  tmp_pts.extend(utils.spherical_grid(3.5, 0.5).tolist())
  tmp_pts.extend(utils.spherical_grid(2.5, 0.5).tolist())
  pts = []
  for i in range(len(tmp_pts)):
    p = numpy.array(tmp_pts[i])
    p = p / numpy.sqrt(sum(p*p))
    if(sum(p*x) >= 0.5):
      pts.append(tmp_pts[i])
  return numpy.array(pts)
###############################################################################

###############################################################################
def compute_global_points(local_cap_pts, ideal_pts, hbond_atom, C_nbr_pos,
                          min_pt, max_pt, residues):
  point_clouds = []
  neighbors = []
  required_fraction_of_points = 0.1

  hvy_atom_pos = numpy.array(hbond_atom.position)
  C_nbr_dir = utils.unit_vec(hvy_atom_pos, C_nbr_pos);
  C_nbr_dist = utils.dist(C_nbr_pos, hvy_atom_pos)
  for i in range(ideal_pts.shape[0]):
    # remove any ideal points that do not fall inside the box
    keep = True
    for j in range(3):
      if(ideal_pts[i][j] < min_pt[j] or max_pt[j] < ideal_pts[i][j]):
        keep = False
    if(not keep):
      point_clouds.append(numpy.array([]))
      neighbors.append([])
      continue

    pt_pos = ideal_pts[i]
    pt_dir = utils.unit_vec(hvy_atom_pos, pt_pos)

    cos_angle = sum(C_nbr_dir*pt_dir)
    sin_angle = numpy.sin(numpy.arccos(cos_angle))

    # Get LSE fit
    ideal_pt_len = utils.dist(pt_pos, hvy_atom_pos)
    local_coords = [ [ideal_pt_len, 0, 0], [0, 0, 0],
                     [cos_angle*C_nbr_dist, sin_angle*C_nbr_dist, 0.0] ]
    local_coords = numpy.array(local_coords)
    global_coords = numpy.array([ pt_pos, hvy_atom_pos, C_nbr_pos ])
    (R,T) = utils.simple_lse_fit(global_coords, local_coords)

    # Shift the cap representation to the global position
    global_pts = numpy.dot(local_cap_pts, R)
    global_pts += numpy.tile(T, (local_cap_pts.shape[0], 1))

    pts_flag = [ 1 for ii in range(global_pts.shape[0]) ]
    for j in range(global_pts.shape[0]):
      for k in range(3):
        if(global_pts[j][k] < min_pt[k] or max_pt[k] < global_pts[j][k]):
          #print min_pt[k], global_pts[j][k], max_pt[k]
          pts_flag[j] = 0 

    if(sum(pts_flag) == 0):
      point_clouds.append(numpy.array([]))
      neighbors.append([])
      continue

    # Find all atoms within 3.5 + 2.5 (A) of the heavy atom's center
    max_sq_dist = 36.0
    min_sq_dist = 2.5*2.5

    print 
    print hbond_atom.serial, ideal_pts[i]
    nearby_atoms = {}
    for res in residues:
      for a in res.atoms:
        # omit self
        if(a.serial == hbond_atom.serial): continue
        pos = numpy.array(a.position)
        if(utils.squared_dist(pos, hvy_atom_pos) > max_sq_dist): continue

        # Flag bumped points
        for j in numpy.arange(global_pts.shape[0]):
          if(pts_flag[j] == 0): continue
          if(utils.squared_dist(pos, global_pts[j]) < min_sq_dist):
            nearby_atoms[a.serial] = a
            pts_flag[j] = 0

    for a_num in nearby_atoms.keys():
      print a_num,
    print

    tmp_pts = []
    for j in numpy.arange(global_pts.shape[0]):
      if(pts_flag[j]): tmp_pts.append(global_pts[j].tolist())
    tmp_pts = numpy.array(tmp_pts)

    if(float(tmp_pts.shape[0]) / float(local_cap_pts.shape[0]) >= \
       required_fraction_of_points): 
      point_clouds.append(numpy.array(tmp_pts))
      neighbors.append([ nearby_atoms[a_num] for a_num in nearby_atoms.keys() ])
    else:
      point_clouds.append(numpy.array([]))
      neighbors.append([])

  return (point_clouds, neighbors)
  
###############################################################################

###############################################################################
def compute_ideal_pts(C_nbr, hbond_atom, other_nbr):

  # 1) Get local coordinate system
  C_pos = numpy.array(C_nbr.position)
  atom_pos = numpy.array(hbond_atom.position)
  other_pos = numpy.array(other_nbr.position)
 
  tmp = C_pos - atom_pos
  C_dist = numpy.sqrt(sum(tmp * tmp))
  tmp = other_pos - atom_pos
  other_dist = numpy.sqrt(sum(tmp*tmp))
  C_dir = utils.unit_vec(atom_pos, C_pos)
  other_dir = utils.unit_vec(atom_pos, other_pos)
  cos_angle = sum(C_dir * other_dir)
  sin_angle = numpy.sin(numpy.arccos(cos_angle))
  local_coords = [ [C_dist, 0.0, 0.0], [0.0, 0.0, 0.0],
                   [cos_angle*other_dist, -1.0*sin_angle*other_dist, 0.0] ]
  local_coords = numpy.array(local_coords)
  global_coords = numpy.array([C_pos, atom_pos, other_pos])
  (R,T) = utils.simple_lse_fit(global_coords, local_coords)
  
  # 2) compute the offset in the local coordinate system
  offsets = []
  for pt in my_ideal_pts:
    offsets.append(compute_offset(3.0, pt.alpha, pt.beta))
  offsets = numpy.array(offsets)
  
  # 3) move the ideal_pt to global coordinate system
  global_pts = numpy.dot(offsets, R)
  global_pts += numpy.tile(atom_pos, (offsets.shape[0], 1))

  return global_pts

###############################################################################

###############################################################################
# returns ([min pt], [max pt])
#def lig_bound_vol(residues, lig, tol=2.0):
#  prot_pos = []
#  lig_pos = []
#
#  for res in residues:
#    for atom in res.atoms:
#      prot_pos.append(atom.position)
#  prot_pos = numpy.array(prot_pos)
#
#  for atom in lig.atoms:
#    lig_pos.append(atom.position)
#  lig_pos = numpy.array(lig_pos)
 

###############################################################################
if __name__ == "__main__":
#  from pymol import cmd
#  from pymol.cgo import *
#  import pymol
#  pymol.finish_launching()
#  import cgo_items

  csv_fname = "/home/vanvoor4/data/new_sampling/pterins/dbase/2qx0_ph2_s.csv"
  prot_fname = "/home/vanvoor4/data/new_sampling/pterins/proteins/2qx0_ph2_p.pdb"
  lig_fname = "/home/vanvoor4/data/new_sampling/pterins/ligands/2qx0_ph2_l.mol2"
  min_fraction_pts = 0.10

  #cmd.load(prot_fname, "2qx0_rad")
  
  prot = utils.pdb.residues(prot_fname)
  lig = utils.mol2.molecule(lig_fname)
  
  (lig_min_pt, lig_max_pt) =  lig_bound_vol(lig, tol=2.0)
  hbond_atoms = get_prot_polar_atoms(prot, lig_min_pt, lig_max_pt)
  ideal_pts = read_ideal_pts("/home/vanvoor4/code/SimSite3D/trunk/params/new_optimum_hbonds.dat")
  
  local_cap_pts = init_cap_points(0.3)
  
  hbond_groups = []
  
  # for each hbond atom
  # 0) Look up the atom in ideal_pts
  # 1) get the neighbor atoms from prot
  # 2) Get local coordinate system
  # 3) compute the offsets in the local coordinate system
  # 4) move the ideal_pt to global coordinate system
  # 5) Put down points for cap
  # 6) Find all neighbor atoms of hbond atom
  # 7) Determine if enough points are left for the cap to be of use
  #   a) Clip points falling outside of ligand box -- very fast
  #   b) Clip points falling inside a neighbor atom
  for hbond_atom in hbond_atoms:
   # if(not hbond_atom.resSeq == 46 or not hbond_atom.name == " N  "): continue

    # 0) Look up the atom in ideal_pts
    my_ideal_pts = ideal_pts[hbond_atom.resName][hbond_atom.name.strip()]  
  
    # 1) get the neighbor atoms from prot
    (C_nbr, other_nbr) = get_hbond_atom_nbrs(prot.residues, hbond_atom,
                                             my_ideal_pts[0].atoms[0],
                                             my_ideal_pts[0].atoms[2])

    global_pts = compute_ideal_pts(C_nbr, hbond_atom, other_nbr)

#    # 2) Get local coordinate system
#    #print C_nbr
#    #print hbond_atom
#    #print other_nbr
#    C_pos = numpy.array(C_nbr.position)
#    atom_pos = numpy.array(hbond_atom.position)
#    other_pos = numpy.array(other_nbr.position)
#   
#    tmp = C_pos - atom_pos
#    C_dist = numpy.sqrt(sum(tmp * tmp))
#    tmp = other_pos - atom_pos
#    other_dist = numpy.sqrt(sum(tmp*tmp))
#    C_dir = utils.unit_vec(atom_pos, C_pos)
#    other_dir = utils.unit_vec(atom_pos, other_pos)
#    #print "NC dir:", C_dir
#    #print "Nother dir:", other_dir
#    #print "C_dir * other_dir", C_dir * other_dir
#    cos_angle = sum(C_dir * other_dir)
#    sin_angle = numpy.sin(numpy.arccos(cos_angle))
#    #print "cos(angle):", sum(C_dir * other_dir)
#    #print "angle:", utils.rad2deg(numpy.arccos(cos_angle))
#    local_coords = [ [C_dist, 0.0, 0.0], [0.0, 0.0, 0.0],
#                     [cos_angle*other_dist, -1.0*sin_angle*other_dist, 0.0] ]
#    local_coords = numpy.array(local_coords)
#    global_coords = numpy.array([C_pos, atom_pos, other_pos])
#    (R,T) = utils.simple_lse_fit(global_coords, local_coords)
#    #print numpy.dot(global_coords - numpy.tile(T, (3,1)), R.T)
#  
#    # 3) compute the offset in the local coordinate system
#    offsets = []
#    for pt in my_ideal_pts:
#      offsets.append(compute_offset(3.0, pt.alpha, pt.beta))
#    offsets = numpy.array(offsets)
#    #print "offsets:", offsets
  
#    # 4) move the ideal_pt to global coordinate system
#    global_pts = numpy.dot(offsets, R)
#    global_pts += numpy.tile(atom_pos, (offsets.shape[0], 1))

#    # remove any global points that do not fall inside the box
#    tmp = []
#    for i in range(global_pts.shape[0]):
#      keep = True
#      for j in range(3):
#        if(global_pts[i][j] < lig_min_pt[j] or \
#           lig_max_pt[j] < global_pts[i][j]):
#          keep = False
#      if(keep): tmp.append(global_pts[i].tolist())
#    global_pts = numpy.array(tmp)
#    if(global_pts.shape[0] == 0): continue
  
    (pt_clouds, nbrs) = compute_global_points(local_cap_pts, global_pts, 
                                              hbond_atom, C_pos, lig_min_pt, 
                                              lig_max_pt, prot.residues)
    hbond_groups.append((hbond_atom, global_pts, pt_clouds, nbrs))

    # check against 2qx0_ph2_s.csv
    #for i in numpy.arange(global_pts.shape[0]):
    #  print global_pts[i]

#    num_ideal_pts = 0
#    cgo_items.draw_points(global_pts, "ideal_pts", color=[0.0, 0.0, 0.8],
#                          radius = 0.2)
#    if(len(hbond_groups)): 
#      cgo_items.draw_plane(C_nbr, hbond_atom, other_nbr)
#      for hb_grp in hbond_groups:
#        print hb_grp[1]
#        print hb_grp[2]
#        cgo_items.draw_points(hb_grp[2], "fit_pts_%d" % (num_ideal_pts))
#        num_ideal_pts += 1
#      break

#  from ASCbasePy.pymol_scripts import SimSite3D
#  my_site = SimSite3D.site_map(csv_fname)
#  my_site.draw("2qx0_ph2_site", sphere_rad=0.2)

  from test_xml import test_xml as test_xml

  my_klas = test_xml(csv_fname)

  cnt = 0
  for hb_grp in hbond_groups:
    (hbond_atom, ideal_pts, clouds, nbrs) = hb_grp

    add = False
    tmp_pts = []
    for i in range(ideal_pts.shape[0]):
      if(clouds[i].shape[0]): 
        tmp_pts.append(ideal_pts[i])
        add = True
      else: tmp_pts.append([])

    if(add): my_klas.add_clouds(hbond_atom, tmp_pts, nbrs)

    #my_klas.add_cloud(hbond_atom, 

  blue = file("blue.xml", "w+")
  my_klas.toxml(blue)
  my_klas.cleanup()

#  color = [1.0, 1.0, 1.0]
#  if(hbond_atom.interact_type == "ACCEPTOR"): color = [0.1, 0.1, 0.8]
#  elif(hbond_atom.interact_type == "DONOR"): color = [0.8, 0.1, 0.1]
#
#  cgo_items.draw_points(cloud, "fit_pts_%d" % cnt, color=color)
#  cnt += 1
#
#
#
#serials = [327, 329, 332, 333, 334, 335, 432, 433, 978, 980, 983, 984, 986]
#for res in prot:
#  if(res.chainID == "A"):
#    for a in res:
#      if(a.serial in serials):
#        cgo = [COLOR, 0.0, 0.8, 0.0]
#        cgo.append(SPHERE)
#        cgo.extend(a.position)
#        cgo.append(2.5)
#        cmd.load_cgo(cgo, "atom_%d" % (a.serial))


