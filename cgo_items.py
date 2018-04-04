from SimSite3DPy import *
from numpy import *
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
  pts = array([[ 10.0, 10.0,0.0], [-10.0, 10.0, 0.0], [10.0, -10.0, 0.0],
               [-10.0,-10.0,0.0], [-10.0, 10.0, 0.0], [10.0, -10.0, 0.0]])
  
  # Fit [[0.0, 0.0, 1.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]] to 
  #     [ hbond_atom pos + cros_prod (plane normal), hbond_atom pos, 
  #       hbond_atom pos + hbond_atom to C_nbr]
  C_pos = array(C_nbr.position)
  hb_pos = array(hbond_atom.position)
  other_pos = array(other_nbr.position)

  hb_to_C = utils.unit_vec(hb_pos, C_pos)
  hb_to_other = utils.unit_vec(hb_pos, other_pos)
  N = cross(hb_to_other, hb_to_C)
  N = utils.normalize(N)
  
  local_coords = array([[0.0, 0.0, 1.0], [0.0, 0.0, 0.0], 
                              [1.0, 0.0, 0.0]])
  global_coords = array([hb_pos + N, hb_pos, hb_pos + hb_to_C])
  (R,T) = utils.simple_lse_fit(global_coords, local_coords)
  
  pts = dot(pts, R) + tile(T, (6,1))

  ## create white plane using 2 triangles
  cgo = []
  cgo.extend([BEGIN, TRIANGLES])
  cgo.extend([COLOR, 1.0, 1.0, 1.0])
  for i in range(6):
   cgo.append(VERTEX)
   cgo.extend(pts[i].tolist())
  cgo.append(END)
  cmd.load_cgo(cgo, "plane")

def draw_plane2(N, p0, plane_num=0, state=0, color=[0.7, 0.7, 0.7]):
  # Assume p0 is in center of "plane"

  # Normal to plane is [0, 0, 1.0]
  # Pt on plane is [0, 0, 0]
  pts = 3.0 * array([[ 10.0, 10.0,0.0], [-10.0, 10.0, 0.0], [10.0, -10.0, 0.0],
                     [-10.0,-10.0,0.0], [-10.0, 10.0, 0.0], [10.0, -10.0, 0.0]])
  #pts = array([[ 3.0, 3.0, 0.0], [-3.0, 3.0, 0.0], [3.0, -3.0, 0.0],
  #             [-3.0,-3.0, 0.0], [-3.0, 3.0, 0.0], [3.0, -3.0, 0.0]])

# determine the rotation of [0.0, 0.0, 1.0] to be in line with N
  Z = array([0.0, 0.0, 1.0])
  cos_phi = sum(N * Z)
  R = eye(3)
  if(cos_phi < 1 - 1E-07):
    C = cross(N,Z)
    sin_phi = sqrt(sum(C*C))
    C /= sin_phi
    Q = utils.Quaternion(cos_theta=cos_phi, sin_theta=sin_phi, V=C)
    R = Q.get_ortho_rot_mat()

  pts = dot(pts,R) + tile(p0, (6,1))

  ## create gray plane using 2 triangles
  cgo = []
  cgo.extend([BEGIN, TRIANGLES])
  cgo.append(COLOR) 
  cgo.extend(color)
  for i in range(6):
    cgo.append(VERTEX)
    cgo.extend(pts[i].tolist())
  cgo.append(END)
  cgo.append(CYLINDER)
  cgo.extend(p0)
  cgo.extend(p0 + 10*N)
  cgo.append(0.25)
  cgo.extend([1.0, 0.0, 0.0, 0.0, 0.0, 1.0])
  if(state): cmd.load_cgo(cgo, "plane%d_%d" % (plane_num, state), state)
  else: cmd.load_cgo(cgo, "plane%d_%d" % (plane_num, state))
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
def draw_line(U, P, label="", color=[1.0, 1.0, 1.0], radius=0.05, state=0):
  """
  Draw a cylinder that denotes the line defined by Ut + P for t in [-1000,1000]
  """
  my_color = color
  #my_color.extend(color)
  my_cgo = []
  U = array(U)
  P = array(P)
  pts = array([-1000.0 * U + P, 1000.0 * U + P])
  my_cgo = [CYLINDER]
# First the 2 points
  my_cgo.extend(pts[0].tolist())
  my_cgo.extend(pts[1].tolist())
# Next radius of cylinder
  my_cgo.append(radius)
# Finally the colors of the 2 end points of the cylinder
  my_cgo.extend(color)
  my_cgo.extend(color)
#  print "cgo:", cgo
#  print "state:", state, "   label:", label
  if(state): cmd.load_cgo(my_cgo, label, state)
  elif(len(label)): cmd.load_cgo(my_cgo, label)
###############################################################################

###############################################################################
def draw_arc(center, radius, end_pts, mid_pt, label="", color=[1.0, 1.0, 1.0],
             state=1):
  my_cgo = [COLOR]
  my_cgo.extend(color)
  sphere_rad = 0.125
  for pt in end_pts:
    my_cgo.append(SPHERE)
    my_cgo.extend(pt.tolist())
    my_cgo.append(sphere_rad)
  my_cgo.append(SPHERE)
  my_cgo.extend(mid_pt.tolist())
  my_cgo.append(sphere_rad)

  arc_pts = array([end_pts[0], mid_pt, end_pts[1]])
  for i in range(5):
    new_arc_pts = []
    for j in range(len(arc_pts) - 1):
      new_arc_pts.append(arc_pts[j])
      mid_pt = (arc_pts[j] + arc_pts[j+1])/2.0
      dir = utils.unit_vec(center, mid_pt)
      new_arc_pts.append(center + dir * radius)
    # must add last point
    new_arc_pts.append(arc_pts[-1])
    arc_pts = array(new_arc_pts)

  for i in range(arc_pts.shape[0] - 1):
    my_cgo.append(CYLINDER)
    my_cgo.extend(arc_pts[i])
    my_cgo.extend(arc_pts[i+1])
    my_cgo.append(sphere_rad/2.0)
    my_cgo.extend(color)
    my_cgo.extend(color)


  if(len(label)): cmd.load_cgo(my_cgo, label, state=state)
###############################################################################

###############################################################################
def draw_msms_surf(vert_fname, face_fname, obj_name, mesh=True,
                   color=[1.0, 1.0, 1.0]):

  try:
    vert_file = open(vert_fname, "r")
  except IOError, (errno, strerror):
    print "Unable to open the file", vert_fname
    print "error(%s): %s" % (errno, strerror)
    return

  try:
    face_file = open(face_fname, "r")
  except IOError, (errno, strerror):
    print "Unable to open the file", face_fname
    print "error(%s): %s" % (errno, strerror)
    return

  # Load verts & normals
  lineno = 0
  verts = []
  normals = []
  for line in vert_file:
    lineno += 1
    if(lineno <= 3): continue

    toks = [ float(s) for s in line.rstrip("\n").split() ]
    verts.append(toks[0:3])
    normals.append(toks[3:6])
  vert_file.close()
  verts = array(verts)


  # Load faces
  lineno = 0
  face_idz = []
  for line in face_file:
    lineno += 1
    if(lineno <= 3): continue

    toks = [ int(s) for s in line.rstrip("\n").split() ]
    face_idz.append(toks[0:3])
  face_file.close()
  # vertices are 1 indexed in the face file -- we need 0 indexed
  face_idz = array(face_idz) - 1

  if(mesh):
    edges = {}
    corr_lines = [LINEWIDTH,1.0, BEGIN, LINES, COLOR]
    corr_lines.extend(color)

    for face in face_idz:
      for i in range(3): 
        end_idx = (i+1)%3
        e0 = "%d_%d" % (face[i], face[end_idx])
        e1 = "%d_%d" % (face[end_idx], face[i])
        if(not e0 in edges and not e1 in edges):
          edges[e0] = True
          edges[e1] = True
          corr_lines.append(VERTEX)
          corr_lines.extend(verts[face[i]])
          corr_lines.append(VERTEX)
          corr_lines.extend(verts[face[end_idx]])
    corr_lines.append(END)
    cmd.load_cgo(corr_lines, obj_name)


#draw_msms_surf("/psa/results/SimSite3D_datasets/testing/pterins/dbase/2qx0_ph2_surf.vert", "/psa/results/SimSite3D_datasets/testing/pterins/dbase/2qx0_ph2_surf.face", "2qx0_mesh", color=[0.2, 0.9, 0.2])
#draw_msms_surf("blah.vert", "blah.face", "blah_mesh", color=[0.2, 0.2, 0.9])
#draw_msms_surf("dbase/YP_HPPK_3ps-0.58_ptr_surf.vert", 
               #"dbase/YP_HPPK_3ps-0.58_ptr_surf.face", "2400ps_mesh",
               #color=[0.9, 0.5, 0.9])
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
  pt1 =  array([  15.69736347, 5.82241156, 7.18412877])
  pt2 = []
  old_pt1 = array([15.697, 5.822, 7.184])
  
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
  ideal_pts = read_ideal_pts("/home/vanvoor4/code/SimSite3D/trunk/params/new_optimum_hbonds.dat")
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
