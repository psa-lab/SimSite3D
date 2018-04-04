from ASCbasePy import utils
import numpy

###############################################################################
class hbond_ideal_pt:

  def __init__(self, line):
    """
    Assumes line is from a SimSite3D V3.3 ideal hbonds data file
    """
    toks = line.split(",")
    self.atoms = [toks[6], toks[1], toks[7]]
    self.obj = toks[2]
    self.pt_num = int(toks[3])
    self.alpha = float(toks[4])
    self.beta = float(toks[5])
    self.orbital = toks[8]
###############################################################################

###############################################################################
class hbond_ideal_pts:

  def __init__(self, ideal_pts_fname):
    self.local_ideal_pts = self.read_ideal_pts(ideal_pts_fname)

###############################################################################

###############################################################################
  def read_ideal_pts(self, fname):
    try:
      infile = file(fname, "r")
    except IOError, (errno, strerror):
      print "Unable to open the file", fname
      print "error(%s): %s" % (errno, strerror)
      return {}
  
    pts = {}
  
    for line in infile:
      if(line.startswith("#")): continue
  
      toks = line.split(",")
      if(len(toks) == 1): continue
  
      if(not toks[0] in pts): pts[toks[0]] = {}
      if(not toks[1] in pts[toks[0]]): pts[toks[0]][toks[1]] = []
      pts[toks[0]][toks[1]].append(hbond_ideal_pt(line))
  
  
    # Main chain points are listed once in the file
    res_names = [ "ALA", "ARG", "ASP", "ASN", "CYS", "GLU", "GLN", "GLY", "HIS",
                  "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TYR",
                  "TRP", "VAL" ]
  
    for res_name in res_names:
      if(not res_name in pts): pts[res_name] = {}
      pts[res_name]["O"] = pts["MAIN_CHAIN"]["O"]
      # Proline residues' amide N cannot donate 
      if(not res_name == "PRO"):
        pts[res_name]["N"] = pts["MAIN_CHAIN"]["N"]
    return pts
###############################################################################

###############################################################################
  def compute_global_pos(self, C_nbr, hbond_atom, other_nbr, residues):
    """
    Compute the global location of the ideal point(s) for the polar atom 
    
    It is assumed that C_nbr, hbond_atom, other_nbr all have the position 
    element.  In addition it is assumed that C_nbr is either the PSA lab
    designated preacceptor or donor antecedent and the other nbr is the other
    atom from the same residue that is used to define the acceptor or donor
    plane

    Parameters
    ----------
    C_nbr : atom -- must have position field
    hbond_atom : atom -- the atom capable of forming a hydrogen bond
    other_nbr : The third atom used to define the donor or acceptor plane
    residues : Protein residues containing the 3 protein atoms

    Returns
    -------
    The global position of the ideal hydrogen bond point(s) -- well the
    DHA angle is ideal, but the offset is currently 3.0 (A) 

    """

    # 1) Get local coordinate system
    C_pos = numpy.array(C_nbr.position)
    atom_pos = numpy.array(hbond_atom.position)
    other_pos = numpy.array(other_nbr.position)
    (R,T) = utils.get_local_coords(C_pos, atom_pos, other_pos)
    (R2,T2) = utils.get_local_coords2(C_pos, atom_pos, other_pos)
   
#    tmp = C_pos - atom_pos
#    C_dist = numpy.sqrt(sum(tmp * tmp))
#    tmp = other_pos - atom_pos
#    other_dist = numpy.sqrt(sum(tmp*tmp))
#    C_dir = utils.unit_vec(atom_pos, C_pos)
#    other_dir = utils.unit_vec(atom_pos, other_pos)
#    cos_angle = sum(C_dir * other_dir)
#    sin_angle = numpy.sin(numpy.arccos(cos_angle))
#    local_coords = [ [C_dist, 0.0, 0.0], [0.0, 0.0, 0.0],
#                     [cos_angle*other_dist, -1.0*sin_angle*other_dist, 0.0] ]
#    local_coords = numpy.array(local_coords)
#    global_coords = numpy.array([C_pos, atom_pos, other_pos])
#    (R,T) = utils.simple_lse_fit(global_coords, local_coords)
    
    # 2) compute the offset in the local coordinate system
    offsets = []
    for pt in self.local_ideal_pts[hbond_atom.resName][hbond_atom.name.strip()]:
#      print hbond_atom.resName, hbond_atom.name.strip(), hbond_atom.serial
#      print "alpha: ", pt.alpha, "   beta:", pt.beta, "   orbital:", pt.orbital
      offsets.append(self.compute_offset(3.0, pt.alpha, pt.beta, pt.orbital))
    offsets = numpy.array(offsets)
    
    # 3) move the ideal_pt to global coordinate system
    global_pts = numpy.dot(offsets, R)
    global_pts += numpy.tile(atom_pos, (offsets.shape[0], 1))
  
    return global_pts
  
###############################################################################

###############################################################################
  def compute_offset(self, bond_len, alpha, beta, orbital="SP2"):
    """
    Assumption -- the hbond atom is at the origin and we start with the
    ideal point at [bond_len, 0, 0].  To get the correct orientation we
    first rotate the VECTOR [bond_len, 0.0, 0.0] about the Y axis by beta to 
    set the out of plane angle.  Then we rotate the resulting VECTOR about the 
    Z axis by alpha to set the in plane angle.  
    
    Note: we must be careful when using rotation matrices on whether the 
    objects or coordinate systems are being rotated by the given angles.
    Here we must rotate the VECTOR and not the coordinate system since that
    is the way angles are presented for the hbond points.

    The previous method of doing this was:
      offset[0] = bond_len * cos_alpha;
      offset[1] = bond_len * cos_beta * sin_alpha,
      offset[2] = bond_len * sin_beta * sin_alpha;

    After more deliberation, I have determined that the main issue with the
    previous way of doing things was to try to solve two different problems
    with the same solution.  Such a solution is bound to be wrong in one of the
    two cases.  In particular we have:
      For SP2: the points are defined in plane on an arc and then moving the 
      arc out of plane 20 or so degrees in either direction (rotation about
      the Y axis).  To achieve this we must first rotate to achieve the 
      desired out of plane angle and then rotate about the Z-axis to get the
      desired in plane angle.

      For SP3: the problem is different, we first move the points to have the
      correct angle about the Z axis and then rotate about the X-axis to 
      place the points on the cone with axis given by the X-axis, apex at
      the origin, and the angle of the cone defined by the inplane angle in
      the XY-plane.
    """
    alpha = utils.deg2rad(alpha) 
    beta = utils.deg2rad(beta) 
    cos_alpha = numpy.cos(alpha) 
    cos_beta = numpy.cos(beta) 
    sin_alpha = numpy.sin(alpha) 
    sin_beta = numpy.sin(beta) 

    if(orbital == "SP2"):
      return [bond_len * cos_beta * cos_alpha, bond_len * cos_beta * sin_alpha, 
              -1.0 * bond_len * sin_beta]
    elif(orbital == "SP3"):
      return [bond_len * cos_alpha, bond_len*cos_beta*sin_alpha,
              bond_len*sin_alpha*sin_beta]
    else:
      print "Unknown orbital in compute_offset"
      return []

#    return [bond_len * cos_alpha, bond_len * cos_beta * sin_alpha, 
#            bond_len * sin_beta * sin_alpha]
###############################################################################

###############################################################################
  def get_polar_atom_triplets(self, residues, site_vol, tol=3.5):
    """
    Tolerance is 3.5 since that is the farthest away a polar atom can be 
    from the site volume and hope to contribute to the hydrogen bond part
    of the site map under the current model.
    """
    #obj_to_act = {"H":"DONOR", "LP":"ACCEPTOR", "U":"DONEPTOR" }

    atom_acts = []
    atom_triplets = []
    for i in range(len(residues)):
      for atom in residues[i].atoms:
        if(not atom.interact_type == "ACCEPTOR" and
           not atom.interact_type == "DONOR" and
           not atom.interact_type == "DONEPTOR"): continue
        if(not site_vol.contains(atom.position, tol)): continue
        

        # Just skip OXT for now
        if(atom.name.strip() == "OXT"): continue
    
        my_pts = self.local_ideal_pts[atom.resName][atom.name.strip()]
        C_nbr_name = my_pts[0].atoms[0]
        other_nbr_name = my_pts[0].atoms[2]

        # Carbonyl O
        if(atom.name == " O  "):
          for tmp_atom in residues[i]:
            if(tmp_atom.name.strip() == C_nbr_name): C_nbr = tmp_atom
          if(i+1 >= len(residues)): continue
          for tmp_atom in residues[i+1]:
            if(tmp_atom.name.strip() == other_nbr_name): other_nbr = tmp_atom
          if(2.25 < utils.squared_dist(numpy.array(C_nbr.position), 
                                       numpy.array(other_nbr.position))): continue

        # amide N
        elif(atom.name == " N  "):
          if(i <= 0): continue

          for tmp_atom in residues[i-1]:
            if(tmp_atom.name.strip() == C_nbr_name): C_nbr = tmp_atom
            elif(tmp_atom.name.strip() == other_nbr_name): other_nbr = tmp_atom

          if(2.25 < utils.squared_dist(numpy.array(C_nbr.position), 
                                       numpy.array(atom.position))): continue
        else:
          #print "other atom"
          for tmp_atom in residues[i]:
            if(tmp_atom.name.strip() == C_nbr_name): C_nbr = tmp_atom
            elif(tmp_atom.name.strip() == other_nbr_name): other_nbr = tmp_atom
       
        atom_triplets.append((C_nbr, atom, other_nbr))
    #    atom_acts.append(obj_to_act[my_pts[0].obj])
    #return (atom_triplets, atom_acts)
    return atom_triplets
###############################################################################
