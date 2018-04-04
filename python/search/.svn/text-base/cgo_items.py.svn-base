from ASCbasePy import *
import numpy
from pymol import cmd
from pymol.cgo import *
import pymol

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
def get_hbond_atom_nbrs(residues, hbond_atom, C_name, other_name):
  # Not sure about the "python" way to do this, but ...
  for i in range(len(residues)):
    if(residues[i]._residue__contains(hbond_atom)):
      if(hbond_atom.name == " O  "):
        print "carbonyl O"
        for tmp_atom in residues[i]:
          if(tmp_atom.name.strip() == C_name): C_nbr = tmp_atom

        if(i+1 < len(residues)):
          for tmp_atom in residues[i+1]:
            if(tmp_atom.name.strip() == other_name): other_nbr = tmp_atom

           # distance check here
      elif(hbond_atom.name == " N  "):
        print "amide N"
        print my_ideal_pts[0].atoms
        if(i > 0):
          for tmp_atom in residues[i-1]:
            if(tmp_atom.name.strip() == C_name): C_nbr = tmp_atom
            elif(tmp_atom.name.strip() == other_name): other_nbr = tmp_atom
           # distance check here
      else:
        print "other atom"
        for tmp_atom in residues[i]:
          if(tmp_atom.name.strip() == C_name): C_nbr = tmp_atom
          elif(tmp_atom.name.strip() == other_name): other_nbr = tmp_atom
      break
  return (C_nbr, other_nbr)
###############################################################################

###############################################################################
def draw_plane(C_nbr, hbond_atom, other_nbr):
  # Normal to plane is [0, 0, 1.0]
  # Pt on plane is [0, 0, 0]
  pts = numpy.array([[ 10.0, 10.0,0.0], [-10.0, 10.0, 0.0], [10.0, -10.0, 0.0],
                     [-10.0,-10.0,0.0], [-10.0, 10.0, 0.0], [10.0, -10.0, 0.0]])
  
  # Fit [[0.0, 0.0, 1.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]] to 
  #     [ hbond_atom pos + cros_prod (plane normal), hbond_atom pos, 
  #       hbond_atom pos + hbond_atom to C_nbr]
  C_pos = numpy.array(C_nbr.position)
  hb_pos = numpy.array(hbond_atom.position)
  other_pos = numpy.array(other_nbr.position)

  hb_to_C = utils.unit_vec(hb_pos, C_pos)
  hb_to_other = utils.unit_vec(hb_pos, other_pos)
  N = numpy.cross(hb_to_other, hb_to_C)
  print N
  N = utils.normalize(N)
  print N
  
  local_coords = numpy.array([[0.0, 0.0, 1.0], [0.0, 0.0, 0.0], 
                              [1.0, 0.0, 0.0]])
  global_coords = numpy.array([hb_pos + N, hb_pos, hb_pos + hb_to_C])
  (R,T) = utils.simple_lse_fit(global_coords, local_coords)
  
  pts = numpy.dot(pts, R) + numpy.tile(T, (6,1))

  ## create white plane using 2 triangles
  cgo = []
  cgo.extend([BEGIN, TRIANGLES])
  cgo.extend([COLOR, 1.0, 1.0, 1.0])
  for i in range(6):
   cgo.append(VERTEX)
   cgo.extend(pts[i].tolist())
  cgo.append(END)
  cmd.load_cgo(cgo, "plane")

###############################################################################

###############################################################################
def draw_points(points, pymol_label, color=[1.0, 1.0, 1.0], radius=0.05):
  my_color = [COLOR]
  my_color.extend(color)
  cgo = []
  cgo.extend(my_color)
  for i in range(points.shape[0]):
    cgo.append(SPHERE)
    cgo.extend(points[i].tolist())
    cgo.append(radius)
  cmd.load_cgo(cgo, pymol_label)
###############################################################################

###############################################################################
if __name__ == "__main__":
  pymol.finish_launching()
  prot = "/home/vanvoor4/data/new_sampling/pterins/dbase/2qx0_ph2_rad.pdb"
  #hbond_stuff = (" OD1", 56)
  #pt1 = numpy.array([ 16.17263938,  5.65454786, -0.83748533])
  #pt2 = numpy.array([ 15.53924725,  6.33700816,  3.30172502])
  #old_pt1 = numpy.array([16.585,7.137,4.126 ])
  
  #hbond_stuff = (" O  ", 44)
  #pt1 = numpy.array([15.86395683, 0.71532213, 2.98540801])
  #pt2 = numpy.array([15.39491103, 3.31517718, 6.74658029])
  #old_pt1 = numpy.array([15.395, 3.315, 6.747])
  
  hbond_stuff = (" N  ", 46)
  pt1 =  numpy.array([  15.69736347, 5.82241156, 7.18412877])
  pt2 = []
  old_pt1 = numpy.array([15.697, 5.822, 7.184])
  
  cmd.load(prot, "2qx0_ph2_rad")
  prot = utils.pdb.residues(prot)
  
  # get the hbonding atom
  for res in prot:
    if(res.resSeq == hbond_stuff[1]):
      for a in res:
        if(a.name == hbond_stuff[0]):
          hbond_atom = a
      break
  
  # Get the neighbors
  ideal_pts = read_ideal_pts("/home/vanvoor4/code/SimSite3D_surfaces/params/new_optimum_hbonds.dat")
  my_ideal_pts = ideal_pts[hbond_atom.resName][hbond_atom.name.strip()]
  (C_nbr, other_nbr) = get_hbond_atom_nbrs(prot.residues, hbond_atom,           
                                         my_ideal_pts[0].atoms[0],
                                         my_ideal_pts[0].atoms[2])
  print C_nbr
  print hbond_atom
  print other_nbr
  draw_plane(C_nbr, hbond_atom, other_nbr)


  pt_rad = 0.25
  # add in the "new" points
  cgo = []
  cgo.extend([COLOR, 0.0, 0.0, 0.8])
  cgo.append(SPHERE)
  cgo.extend(pt1.tolist())
  cgo.append(pt_rad)
  cgo.extend([COLOR, 0.0, 0.0, 0.8])
  cgo.append(SPHERE)
  if(len(pt2) == 3):
    cgo.extend(pt2.tolist())
    cgo.append(pt_rad)
  cmd.load_cgo(cgo, "jeff pts")
  
  # add in the simsite point
  cgo = []
  cgo.extend([COLOR, 0.5, 0.5, 0.8])
  cgo.append(SPHERE)
  cgo.extend(old_pt1.tolist())
  cgo.append(pt_rad)
  cmd.load_cgo(cgo, "SimSite pts")


#OD1 = numpy.array([17.688,  7.060,  1.337])
#CG = numpy.array([18.758,  7.603,  1.086])
#ND2 = numpy.array([19.477,  8.243,  2.011])
#
#pt1 = numpy.array([ 16.17263938,  5.65454786, -0.83748533])
#pt2 = numpy.array([ 15.53924725,  6.33700816,  3.30172502])
#old_pt1 = numpy.array([16.585,7.137,4.126 ])
#
#print "dist(OD1, pt1)", utils.dist(OD1, pt1)
#print "dist(OD1, pt2)", utils.dist(OD1, pt2)
#print "dist(OD1, old_pt1)", utils.dist(OD1, old_pt1)
#
#print "angle(pt1, OD1, CG)",  utils.angle(pt1, OD1, CG)
#print "angle(pt2, OD1, CG)",  utils.angle(pt2, OD1, CG)
#print "angle(old_pt1, OD1, CG)", utils.angle(old_pt1, OD1, CG)
#
## check
#CG_to_ND2 = utils.unit_vec(CG, ND2)
#OD1_to_pt1 = utils.unit_vec(OD1, pt1)
#print utils.rad2deg(numpy.arccos(sum(CG_to_ND2 * OD1_to_pt1)))
#
#pts = numpy.array([ND2.tolist(), CG.tolist(), OD1.tolist()])
#print "jeff dist to plane:", utils.dist_to_plane(pts, pt1)
#print "jeff dist to plane2:", utils.dist_to_plane(pts, pt2)
#print "SimSite2 dist to plane:", utils.dist_to_plane(pts, old_pt1)
#
#cgo = []
#
## Normal to plane is [0, 0, 1.0]
## Pt on plane is [0, 0, 0]
#tri1_pts = numpy.array([[10,10,0.0], [-10, 10, 0.0], [10, -10, 0.0]])
#tri2_pts = numpy.array([[-10,-10,0.0], [-10, 10, 0.0], [10, -10, 0.0]])
#
## Fit [[0.0, 0.0, 1.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]] to 
##     [ CG pos + cros_prod (plane normal), CG pos, CG pos + CG_to_ND2]
#CG_to_OD1 = utils.unit_vec(CG,OD1)
#N = numpy.cross(CG_to_ND2, CG_to_OD1)
#local_coords = numpy.array([[0.0, 0.0, 1.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
#global_coords = numpy.array([ CG + N , CG, CG + CG_to_ND2])
#(R,T) = utils.simple_lse_fit(global_coords, local_coords)
#
#tri1_pts = numpy.dot(tri1_pts, R) + numpy.tile(T, (3,1))
#tri2_pts = numpy.dot(tri2_pts, R) + numpy.tile(T, (3,1))
#
## create carboxyamide plane
#cgo.extend([BEGIN, TRIANGLES])
#cgo.extend([COLOR, 1.0, 1.0, 1.0])
#for i in range(3):
# cgo.append(VERTEX)
# cgo.extend(tri1_pts[i].tolist())
#for i in range(3):
# cgo.append(VERTEX)
# cgo.extend(tri2_pts[i].tolist())
#cgo.append(END)
#cmd.load_cgo(cgo, "CNO_plane")
#
#pt_rad = 0.25
## add in the "new" points
#cgo = []
#cgo.extend([COLOR, 0.0, 0.0, 0.8])
#cgo.append(SPHERE)
#cgo.extend(pt1.tolist())
#cgo.append(pt_rad)
#cgo.extend([COLOR, 0.0, 0.0, 0.8])
#cgo.append(SPHERE)
#cgo.extend(pt2.tolist())
#cgo.append(pt_rad)
#cmd.load_cgo(cgo, "jeff pts")
#
## add in the simsite point
#cgo = []
#cgo.extend([COLOR, 0.5, 0.5, 0.8])
#cgo.append(SPHERE)
#cgo.extend(old_pt1.tolist())
#cgo.append(pt_rad)
#cmd.load_cgo(cgo, "SimSite pts")
#
#
#
#
####cmd.load("tmp.pdb", "SC")
