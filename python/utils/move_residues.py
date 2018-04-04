from numpy import *
from ASCbasePy.utils import *


def add_noise_to_site_residues(residues_fname, out_fname, trans_sigma,
                               rot_sigma):
  my_residues = pdb.residues(residues_fname)
  out_file = open(out_fname, "w+")  

  # For each residue apply a translation along the CA-CB bond and
  # then apply a rotation from the current CA-CB direction,
  # finally use a uniform distribution to rotate new CA-CB direction
  # about the original CA-CB direction

  for res in my_residues:
    if(res.resName == "PRO"): 
      print >> out_file, res
      continue

    for a in res:  
      if(a.name == " CA "): CA_pos = array(a.position)
      elif(a.name == " CB "): CB_pos = array(a.position)
      elif(a.name == " C  "): C_pos = array(a.position)
    orig_dir = unit_vec(CA_pos, CB_pos)

    # Translation in angstroms
    if(trans_sigma > 0.0): t = random.normal(scale=trans_sigma)
    else: t = 0.0

    T = orig_dir * t
    for a in res:  
      if(a.name == " CA " or a.name == " C  " or a.name == " N  " or \
         a.name == " O  "): continue
      a.position = array(a.position) + T
      if(a.name == " CB "): CB_pos = array(a.position)

    (R, T) = get_local_coords(CB_pos, CA_pos, C_pos)

    # Z "offset" angle in radians
    if(rot_sigma > 0.0): angle = random.normal(scale=rot_sigma)
    else: angle = 0.0
    Q0 = Quaternion(cos_theta=cos(angle), sin_theta=sin(angle), 
                    V = array([0.0, 1.0, 0.0]))
    Ry = Q0.get_ortho_rot_mat()

    # X "offset" angle in radians 
    if(rot_sigma > 0.0): angle = random.normal(scale=rot_sigma)
    else: angle = 0.0
    Q0 = Quaternion(cos_theta=cos(angle), sin_theta=sin(angle), 
                    V = array([0.0, 0.0, 1.0]))
    Rz = Q0.get_ortho_rot_mat()

    # Rotate the offset vector about the local X axis to uniformly sample it
#    rot_angle = random.uniform(low=0.0, high=2.0*pi)
#    Q = Quaternion(cos_theta=cos(rot_angle), sin_theta=sin(rot_angle), 
#                   V=array([1.0, 0.0, 0.0]))
#    R2 = Q.get_ortho_rot_mat()

    for a in res:  
      if(a.name == " CA " or a.name == " C  " or a.name == " N  " or \
         a.name == " O  "): continue
      # Get local orientation
      pos = dot(array(a.position) - T, R.T)

      # move vector out of xy plane
      pos = dot(pos, Ry)

      # move vector in xy plane
      #print "R0", R0
      pos = dot(pos, Rz)
      #print dot(pos, R1) + T1


      # Transform back to global coordinates
      a.position = dot(pos, R) + T
      


    print >> out_file, res



#    V = array([cos(angle), sin(angle), 0.0])
#    V = dot(V, R) + T
#    print V
    




# numpy.random.normal(loc=0.0, scale=1.0, size=None)
# numpy.random.uniform(low=0.0, high=1.0, size=1)
