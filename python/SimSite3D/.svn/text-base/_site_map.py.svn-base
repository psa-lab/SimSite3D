from pymol.cgo import *
from pymol import cmd
import numpy
import tempfile
from SimSite3D import xform_obj


class site_map_pt:

  def __init__(self):
    pass

class site_map:

  def __init__(self, csv_fname):
    """
    Input:
      csv_fname : Name of the sitemap csv file

    Load a site map
    """
    print "loading site map: %s" % (csv_fname.split("/")[-1][:-6])
    from os import path

    (self.dir, self.name) = path.split(csv_fname)
    if(len(self.name) < 7):
      print "did not receive a site map name"
      return

    self.name = self.name[:-6] + "_site"
    self.__init_vars()
    self.__load_points(csv_fname)
    self.__load_csv_file(csv_fname)
    self.__load_mesh_verts("%s/%s_surf.vert" % (self.dir, self.name[:-5]))
    self.__load_mesh_faces("%s/%s_surf.face" % (self.dir, self.name[:-5]))

  def __load_csv_file(self, csv_fname):
    csv_file = file(csv_fname, "r")
    
    for line in csv_file:
      if(line.startswith("<point_counts>")):
        line = csv_file.next()
        num_hbond_pts = int(line.split()[0])
        break

    pts_dict = dict([ (pt.serial, pt) for pt in self.points ])
    atoms_dict = dict([ (a.serial, a) for a in self.atoms ])

    polar_pts = {}
    for line in csv_file:
      if(line.startswith("<hydrogen_bond_points>")):
        for i in range(num_hbond_pts):
          line = csv_file.next()
          toks = line.split("|")[:-1]
          pt = site_map_pt()
          #pt.num = int(toks[0])
          pt.chainID = toks[2]
          pt.resName = toks[3]
          pt.resNum = toks[4]
          pt.atomName = toks[5]
          pt.atomNum = int(toks[6])
          pt.hbondPtNum = int(toks[7])
          pt.hbondPtPos = ( float(s) for s in toks[8].split() )
          pt.carbonNum = int(toks[9])
          polar_pts[int(toks[0])] = pt
          pts_dict[int(toks[0])].atom_pos = atoms_dict[pt.atomNum].position
    csv_file.close()

    dirs = []
    for pt in self.points:
      if(pt.atom_pos == None):
        continue
      
      a = numpy.array(pt.atom_pos)
      b = numpy.array(pt.position)
      dirs.append(((a - b) / numpy.tile(numpy.sqrt(sum((a - b)**2)), [1, 3]))[0,:])
    self.polar_dirs = numpy.array(dirs)
    i = 0
    for pt in self.points:
      if(pt.atom_pos == None):
        continue

      pt.dir = self.polar_dirs[i,:]
      i += 1

  ##############################################################################
  # Load the site map points and corresponding atoms
  ##############################################################################
  def __load_points(self, csv_fname):
    from SimSite3D import pdb
    self.atoms = [] 

    try:
      atoms_file = file(csv_fname[:-5] + "a.pdb")
    except IOError, (errno, strerror):
      print "Unable to open the file", csv_fname[:-5] + "a.pdb"
      print "error(%s): %s" % (errno, strerror)
      return False

    for line in atoms_file:
      if(not line.startswith("ATOM") and not line.startswith("HETATM")):
        continue
      self.atoms.append(pdb.atom(line))
    atoms_file.close()

    self.atom_positions = numpy.array([ a.position for a in self.atoms ])
    for a, pos in zip(self.atoms, self.atom_positions):
      a.position = pos

    self.points = []
    self.point_positions = []
    pts_file = file(csv_fname[:-3] + "pdb")
    for line in pts_file:
      if(not line.startswith("HETATM")):
        continue
      self.points.append(pdb.atom(line))
      self.points[-1].atom_pos = None
      self.points[-1].cgo_sphere_rad = -1.0
    pts_file.close()

    self.point_positions = numpy.array([ p.position for p in self.points ])
    for p, pos in zip(self.points, self.point_positions):
      p.position = pos

    return True

  def __load_mesh_verts(self, vert_fname):
    try:
      vert_file = open(vert_fname, "r")
    except IOError, (errno, strerror):
      print "Unable to open the file", vert_fname
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
    self.mesh_verts = numpy.array(verts)
    self.mesh_normals = numpy.array(normals)

  def __load_mesh_faces(self, face_fname):
    try:
      face_file = open(face_fname, "r")
    except IOError, (errno, strerror):
      print "Unable to open the file", face_fname
      print "error(%s): %s" % (errno, strerror)
      return
  
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
    self.mesh_face_idz = numpy.array(face_idz) - 1
  
  def __init_vars(self):
    self.orig_points_positions = None

    cmd.set_color("light_red", [1.0, 0.6, 0.6])
    cmd.set_color("light_blue", [0.6, 0.6, 1.0])
    cmd.set_color("light_green", [0.6, 1.0, 0.6])

    self.num_to_color = {
      '  0.00' : [COLOR, 1.0, 0.0, 0.0],
      ' 50.00' : [COLOR, 0.0, 0.0, 1.0],
      ' 25.00' : [COLOR, 1.0, 1.0, 1.0],
      '100.00' : [COLOR, 0.0, 1.0, 0.0],
      '200.00' : [COLOR, 0.0, 1.0, 1.0],
    }

    self.num_to_color_name = {
      '  0.00' : "red",
      ' 50.00' : "blue",
      ' 25.00' : "white",
      '100.00' : "green",
      '200.00' : "yellow",
    }
    
    self.num_to_pastel_color = {
      '  0.00' : [COLOR, 1.0, 0.6, 0.6],
      ' 50.00' : [COLOR, 0.6, 0.6, 1.0],
      ' 25.00' : [COLOR, 1.0, 1.0, 1.0],
      '100.00' : [COLOR, 0.6, 1.0, 0.6],
      '200.00' : [COLOR, 0.6, 1.0, 1.0],
    }

    self.num_to_pastel_color_name = {
      '  0.00' : "light_red",
      ' 50.00' : "light_blue",
      ' 25.00' : "white",
      '100.00' : "light_green",
      '200.00' : "yellow",
    }
  
  
    # Used for cylinder so no "color" prefix
    self.name_to_color = {
      'white' : [1.0, 1.0, 1.0],
      'yellow': [1.0, 1.0, 0.0],
    }

  def draw(self, name="", sphere_rad = 0.175, hbond_stick_rad = 0.04, 
           stick_color = "white", pastel=False, state=0,
           hphob_stick_rad = 0.0 ):
    surf_color = [0.9, 0.2, 0.9]
    if(pastel): surf_color = [0.2, 0.9, 0.9]
    self.draw_msms_surf("%s_surf" % (self.name[:-5]), color=surf_color, lw=2,
                        state=state)

    my_colors = self.num_to_color
    my_col2 = self.num_to_color_name
    if(pastel): 
      my_colors = self.num_to_pastel_color
      my_col2 = self.num_to_pastel_color_name
    my_cgo = []

    sphere_cgo = []
    N = len(self.points)
    M = numpy.zeros([N, N])

    if(not len(name)): name = self.name
    my_name = "test_"  + name 

    # get a tempfile to load the site points as atoms 
    (site_fd, site_fname) = tempfile.mkstemp(suffix=".pdb")
    site_file = os.fdopen(site_fd, "w+")
    
    idx = 0
    for pt in self.points:
      if(stick_color == "white"):
        for i in range(N):
          tmp = numpy.array(pt.position) - numpy.array(self.points[i].position)
          M[idx][i] = sum(tmp * tmp)
        idx += 1

      # fill in these at later date [ resn, chain, segi, name]
      my_sphere_rad = sphere_rad
      if(pt.cgo_sphere_rad > 0.0): my_sphere_rad = pt.cgo_sphere_rad  

      print >> site_file, "%s%5d %s%s%s %s%4d%s   %8.3f%8.3f%8.3f%6.2f%6.2f" % \
      ("HETATM", pt.serial, pt.name, pt.altLoc, pt.resName,
       pt.chainID, pt.resSeq, pt.iCode, pt.position[0],
       pt.position[1], pt.position[2], pt.occupancy, pt.tempFactor)


#      if(state == 0):
#        cmd.pseudoatom(object=my_name, elem="O", vdw=my_sphere_rad,
#                       color=my_col2["%6.2f" % (pt.tempFactor)],
#                       pos=pt.position.tolist())
#      else:
#        cmd.pseudoatom(object=my_name, elem="O", vdw=my_sphere_rad,
#                       color=my_col2["%6.2f" % (pt.tempFactor)],
#                       pos=pt.position.tolist(), state=state)

      if(pt.tempFactor == 100.0): 
        if(hphob_stick_rad > 0.0):
          my_cgo.append(CYLINDER)
          my_cgo.extend(pt.position)
          my_cgo.extend(pt.atom_pos)
          my_cgo.append(hphob_stick_rad)
          my_cgo.extend(self.name_to_color[stick_color])
          my_cgo.extend(self.name_to_color[stick_color])

      elif(hbond_stick_rad > 0.0):
        my_cgo.append(CYLINDER)
        my_cgo.extend(pt.position)
        my_cgo.extend(pt.atom_pos)
        my_cgo.append(hbond_stick_rad)
        my_cgo.extend(self.name_to_color[stick_color])
        my_cgo.extend(self.name_to_color[stick_color])
 
    site_file.close()
    my_cgo.extend([COLOR, 0.0, 0.0, 0.0])

    # load site points
    cmd.load(site_fname, my_name, state)
    for tfactor, color in my_col2.iteritems():
      cmd.color(color, "%s & b=%s" % (my_name, tfactor.strip()))
    cmd.alter(my_name, "vdw=%f" % (my_sphere_rad))

    for pt in self.points:
      my_sphere_rad = sphere_rad
      if(pt.cgo_sphere_rad > 0.0): my_sphere_rad = pt.cgo_sphere_rad  
      cmd.alter("%s & i. %d" % (my_name, pt.resSeq), "vdw=%f" % (my_sphere_rad))

    cmd.hide("everything", my_name)
    cmd.show("spheres", my_name)



    if(not len(name)): name = self.name
    if(state): cmd.load_cgo(my_cgo, name, state)
    else: cmd.load_cgo(my_cgo, name)

  # Note: this function will not actually transform the coordinates of any
  # drawn objects
  def transform(self, R, T):
    if(self.orig_points_positions == None):
      self.orig_point_positions = numpy.array(self.point_positions)
      self.orig_atom_positions = numpy.array(self.atom_positions)
      self.orig_polar_dirs = numpy.array(self.polar_dirs)
      self.orig_mesh_verts = numpy.copy(self.mesh_verts)
      self.orig_mesh_norms = numpy.copy(self.mesh_normals)

    X = numpy.dot(self.point_positions, R.T) + \
      numpy.tile(T, [self.point_positions.shape[0], 1])
    for i in range(X.shape[0]):
      self.point_positions[i,:] = X[i,:]

    X = numpy.dot(self.atom_positions, R.T) + \
      numpy.tile(T, [self.atom_positions.shape[0], 1])
    for i in range(X.shape[0]):
      self.atom_positions[i,:] = X[i,:]

    X = numpy.dot(self.polar_dirs, R.T)
    for i in range(X.shape[0]):
      self.polar_dirs[i,:] = X[i,:]

    self.mesh_verts = numpy.dot(self.mesh_verts, R.T) + \
      numpy.tile(T, (self.mesh_verts.shape[0], 1))
    self.mesh_normals = numpy.dot(self.mesh_normals, R.T)

  # NOTE: this will not actually revert any drawn objects
  def revert_positions(self):
    for p, orig_p in zip(self.atom_positions, self.orig_atom_positions):
      p[:] = orig_p

    for p, orig_p in zip(self.point_positions, self.orig_point_positions):
      p[:] = orig_p

    for d, orig_d in zip(self.polar_dirs, self.orig_polar_dirs):
      d[:] = orig_d

  def load_points(self, state=0):
    if(state):
      cmd.load("%s/%ss.pdb" % (self.dir, self.name[:-4]), 
               self.name[:-4] + "pts", state)
    else:
      cmd.load("%s/%ss.pdb" % (self.dir, self.name[:-4]), 
               self.name[:-4] + "pts")
    cmd.disable(self.name[:-4] + "pts")

###############################################################################

###############################################################################
  def draw_msms_surf(self, obj_name, mesh=True, color=[1.0, 1.0, 1.0], lw=1,
                     state=0):
    if(mesh):
      edges = {}
      corr_lines = [LINEWIDTH] 
      corr_lines.append(lw)
      corr_lines.extend([ BEGIN, LINES, COLOR])
      corr_lines.extend(color)
#      print corr_lines
    
      for face in self.mesh_face_idz:
        for i in range(3): 
          end_idx = (i+1)%3 
          e0 = "%d_%d" % (face[i], face[end_idx])
          e1 = "%d_%d" % (face[end_idx], face[i])
          if(not e0 in edges and not e1 in edges):
            edges[e0] = True
            edges[e1] = True
            corr_lines.append(VERTEX)
            corr_lines.extend(self.mesh_verts[face[i]])
            corr_lines.append(VERTEX)
            corr_lines.extend(self.mesh_verts[face[end_idx]])
#            print VERTEX, verts[face[i]], VERTEX, verts[face[end_idx]]
#      print END
      corr_lines.append(END)
      if(state): cmd.load_cgo(corr_lines, obj_name, state=state)
      else:  cmd.load_cgo(corr_lines, obj_name)
###############################################################################

###############################################################################

##################
#
##################
class assess_alignments:

  def __init__(self, hits, query_dir, query_id, query_lig, dbase_dir, 
               dbase_ligs, query_C_color="green", dbase_C_color="grey",
               frag_C_color = "purple", hbond_stick_rad=0.015,  
               hphob_stick_rad=0.0, score_tol = -1.5, best_only=True):
    self.q_dir = query_dir
    if(query_dir == ""): query_dir = "."
    self.q_id = query_id

    # Load the query first
    # We have an issue with query having same name as a database entry
    self.q_site = self.load_items(query_dir, query_id, hbond_stick_rad, 
                                  query_lig, query_C_color, query = True,
                                  hphob_stick_rad = hphob_stick_rad)
    self.dbase_sites_states = \
      self.load_dbase_items(hits, dbase_dir, dbase_ligs, dbase_C_color, 
                            frag_C_color = frag_C_color, score_tol = score_tol,
                            best_only=best_only)
    (self.polar_pairs, self.hphob_pairs) = \
      self.compare_sites(self.q_site, self.dbase_sites_states)

    self.draw_sites()
    # this is a bothersome and time consuming process as written
    #self.show_match_distances()
    cmd.set("stick_radius", "0.09", "*_prot")

  # Load all items for a site but do not draw the sitemap
  def load_items(self, site_dir, site_id, hbond_stick_rad=0.015, lig="",
                 C_color="green", frag_C_color = "purple",
                 hphob_stick_rad=0.0, s_num=0, R="", T="", frag_bit_str="", 
                 query = False):
    site_pref = site_dir + "/" + site_id
    site = site_map(site_pref + "_s.csv")
    site.load_points(s_num)
    if(query): prot_name = "query_prot"
    else: prot_name = site_id + "_prot"
    cmd.load(site_pref + "_rad.pdb", prot_name, s_num)
    cmd.color(C_color, "%s & e. C" % (prot_name))

    if(not lig == ""):
      if(query): lig_name = "query_lig"
      else: lig_name = site_id + "_lig"
      cmd.load(lig, lig_name, s_num)
      cmd.color(C_color, "%s & e. C" % (lig_name))
      cmd.show("sticks", lig_name)

    if(len(R) and len(T)):
      R = numpy.array(R)
      T = numpy.array(T)
      R = R.reshape([3,3])
      T = T.reshape([1,3])
      site.transform(R,T)
      xform_obj(prot_name, R, T, state=s_num)
      xform_obj(site_id + "_pts", R, T, state=s_num)

      if(not lig == ""):
        xform_obj(lig_name, R, T, state=s_num)

# this is not well thought out and colors all atoms not just carbons
#        if(len(frag_bit_str)):
#          my_sel = ""
#          # Note: atoms in mol2 files are 1 indexed
#          for i in range(len(frag_bit_str)):
#            if(not frag_bit_str[i] == "0" and len(my_sel)): 
#              my_sel += "| id %d " % (i + 1)
#            elif(not frag_bit_str[i] == "0"): my_sel = "id %d " % (i + 1)
#          my_sel = site_id + "_lig & (" + my_sel + ")"
#          cmd.color(frag_C_color, my_sel)
#          cmd.create(site_id + "_frag", my_sel, s_num, s_num)
  
    return site

  def load_dbase_items(self, hits, dbase_dir, dbase_ligs, hbond_stick_rad=0.015,
                       dbase_C_color="gray", frag_C_color = "purple",
                       hphob_stick_rad=0.0, score_tol = -1.5, starting_state=2,
                       best_only=True):
    sites_info = []
    state = starting_state
    prev_ids = {} 
    for hit in hits:
      if(hit[1] > score_tol):
        break

      # dumb little hack
      if(hit[0].endswith(".mol2")): 
        dbase_id = hit[0][:-13]
        lig = dbase_ligs + "/" + dbase_id + "_l.mol2"
      else: 
        dbase_id = hit[0]
        lig = dbase_ligs + "/" + dbase_id + "_l.mol2"
        frag_C_color = "gray"
        #lig = ""

      # typically we want to only view the best hit for each pair
      if(best_only):
       if(dbase_id in prev_ids): continue
       else: prev_ids[dbase_id] = True
      
      print "Loading the site: %s" % (dbase_id)
      dbase_site = \
        self.load_items(dbase_dir, dbase_id, hbond_stick_rad, lig, 
                        dbase_C_color, frag_C_color, hphob_stick_rad, state, 
                        hit[2], hit[3], hit[5])
      sites_info.append((dbase_site, state))
      state += 1
    return sites_info

  def compare_sites(self, q_site, dbase_sites_states):
    polar_pairs = {}
    hphob_pairs = {}

    for d in dbase_sites_states:
      (dbase_site, state) = d
      print "Compamring query site to: %s" % (dbase_site.name)
      polar_pairs[state] = []
      hphob_pairs[state] = []

      # Calculate all to all pointwise distances
      q_pos = q_site.point_positions
      d_pos = dbase_site.point_positions
      II = q_pos.shape[0]
      JJ = d_pos.shape[0]

      D2 = numpy.empty((II,JJ))
      for i in range(II):
        D2[i,:] = ((d_pos - numpy.tile(q_pos[i,:], (JJ, 1)))**2 ).sum(axis=1)
      min_D2 = [ (x,i) for x,i in zip(numpy.min(D2, axis=1), 
                                      numpy.argmin(D2, axis=1)) ]
  
      # Distinguish between polar & nonpolar and weight the
      # the nonpolar by dot product
      for i, db_tuple in zip(range(II), min_D2):
        if(db_tuple[0] > 2.25): continue
    
        j = db_tuple[1]
        if(q_site.points[i].tempFactor == 100.0):
          if(dbase_site.points[j].tempFactor == 100.0):
            hphob_pairs[state].append((i,j))
        elif(dbase_site.points[j].tempFactor == 100.0):
          continue
        elif(q_site.points[i].tempFactor == dbase_site.points[j].tempFactor
             or q_site.points[i].tempFactor == 25.0 
             or dbase_site.points[j].tempFactor == 25.0):
          dot_prod = sum(q_site.points[i].dir * dbase_site.points[j].dir)
          if(dot_prod > 0.0):
            polar_pairs[state].append((i, j, dot_prod))
            #q_str = "%s & id %d" % (q_pts_name, q_site.points[i].serial)
            #d_str = "%s & id %d" % (d_pts_name, dbase_site.points[j].serial)
            #polar_pairs[state].append((q_str, d_str, dot_prod))

    return (polar_pairs, hphob_pairs)

  def draw_sites(self, no_match_sphere_rad=0.05, max_polar_sphere_rad=0.25,
                 hbond_stick_rad=0.02, hphob_stick_rad=0.0, 
                 hphob_sphere_rad=0.125):
    self.no_match_sphere_rad = no_match_sphere_rad
    self.max_polar_sphere_rad = max_polar_sphere_rad
    self.hbond_stick_rad = hbond_stick_rad
    self.hphob_stick_rad = hphob_stick_rad
    self.hphob_sphere_rad = hphob_sphere_rad

    query_site = self.q_site
    for d in self.dbase_sites_states:
      (dbase_site, state) = d
      print "drawing the site:  %s" % (dbase_site.name)

      # Initialize the spheres' radii to the unmatched radius
      for ddd in dbase_site.points: ddd.cgo_sphere_rad = no_match_sphere_rad
      for qqq in query_site.points: qqq.cgo_sphere_rad = no_match_sphere_rad

      #polar_pairs = self.polar_pairs[state]
      #hphob_pairs = self.hphob_pairs[state]
      for pp in self.hphob_pairs[state]:
        query_site.points[pp[0]].cgo_sphere_rad = hphob_sphere_rad
        dbase_site.points[pp[1]].cgo_sphere_rad = hphob_sphere_rad

      for pp in self.polar_pairs[state]:
        w = no_match_sphere_rad + pp[2] * max_polar_sphere_rad
        query_site.points[pp[0]].cgo_sphere_rad = w
        dbase_site.points[pp[1]].cgo_sphere_rad = w

      query_site.draw(hbond_stick_rad = hbond_stick_rad, state=state,
                      hphob_stick_rad = hphob_stick_rad)
      dbase_site.draw(hbond_stick_rad = hbond_stick_rad, state=state,
                      hphob_stick_rad = hphob_stick_rad, stick_color="yellow",
                      pastel=True)

  def show_match_distances(self):
    q_pts_name = self.q_site.name[:-4] + "pts"
    q_site = self.q_site

    # Set some items to be more amenable to comparing sitemaps
    cmd.set("dash_gap", 0.025)
    cmd.set("dash_length", 0.1)

    num_polar = 0
    num_hphob = 0
    for d in self.dbase_sites_states:
      (dbase_site, state) = d
      d_pts_name = dbase_site.name[:-4] + "pts"

      for pp in self.polar_pairs[state]:
        num_polar += 1
        (i, j, w) = pp
        cmd.distance("polar_dist_%06d" % (num_polar),
                     "%s & id %d" % (q_pts_name, q_site.points[i].serial),
                     "%s & id %d" % (d_pts_name, dbase_site.points[j].serial))

      for pp in self.hphob_pairs[state]:
        num_hphob += 1
        (i,j) = pp
        cmd.distance("hphob_dist_%06d" % (num_hphob),
                     "%s & id %d" % (q_pts_name, q_site.points[i].serial),
                     "%s & id %d" % (d_pts_name, dbase_site.points[j].serial))
#      cmd.delete(d_pts_name)
#    cmd.delete(q_pts_name)

  
def compare(query, dbase, R, T, state=1):
  # Load the points from file into the state,
  # measure the distances

  q_pts_name = query.name[:-4] + "pts"
  dbase_pts_name = dbase.name[:-4] + "pts"
  cmd.load("%s/%ss.pdb" % (dbase.dir, dbase.name[:-4]), dbase_pts_name,
           state=state)
  xform_obj(dbase_pts_name, R, T, state)

  dbase.transform(R, T)

  q_pos = query.point_positions
  d_pos = dbase.point_positions
  II = q_pos.shape[0]
  JJ = d_pos.shape[0]

  D2 = numpy.empty((II,JJ))
  for i in range(II):
    D2[i,:] = ((d_pos - numpy.tile(q_pos[i,:], (JJ, 1)))**2 ).sum(axis=1)
  min_D2 = [ (x,i) for x,i in zip(numpy.min(D2, axis=1), 
                                  numpy.argmin(D2, axis=1)) ]

  # now need to distinguish between polar & nonpolar and weight the
  # the nonpolar by dot product
  hphob_dists = []
  polar_dists = []
  for d in dbase.points:
    d.cgo_sphere_rad = 0.05

  for i, db_tuple in zip(range(II), min_D2):
    query.points[i].cgo_sphere_rad = 0.05

    if(db_tuple[0] > 2.25): continue

    j = db_tuple[1]
    if(query.points[i].tempFactor == 100.0):
      if(dbase.points[j].tempFactor == 100.0):
        dbase.points[j].cgo_sphere_rad = 0.125
        query.points[i].cgo_sphere_rad = 0.125
        hphob_dists.append("hphob_dist_%d" % len(hphob_dists))
        cmd.distance(hphob_dists[-1], 
                     "%s & id %d" % (q_pts_name, query.points[i].serial),
                     "%s & id %d" % (dbase_pts_name, dbase.points[j].serial))

    elif(query.points[i].tempFactor == dbase.points[j].tempFactor
         or query.points[i].tempFactor == 25.0 
         or dbase.points[j].tempFactor == 25.0):
      dot_prod = sum(query.points[i].dir * dbase.points[j].dir)
      if(dot_prod > 0.0):
        w = 0.05 + dot_prod * 0.20
        # dbase points can be matched by more than 1 query
        if(w > dbase.points[j].cgo_sphere_rad):
          dbase.points[j].cgo_sphere_rad = w
        query.points[i].cgo_sphere_rad = w
        polar_dists.append("polar_dist_%d" % len(polar_dists))
        cmd.distance(polar_dists[-1], 
                     "%s & id %d" % (q_pts_name, query.points[i].serial),
                     "%s & id %d" % (dbase_pts_name, dbase.points[j].serial))


  cmd.create("hphob_dists", "hphob_dist_*", state, state)
  cmd.delete("hphob_dist_*")
  cmd.create("polar_dists", "polar_dist_*")#, state, state)
  cmd.delete("polar_dist_*")
  query.draw(state=state, stick_rad = 0.02)
  dbase.draw(stick_color="yellow", pastel=True, state=state, 
             stick_rad = 0.02)
  dbase.revert_positions()
  cmd.delete(dbase_pts_name)

####
# Does not necessarily belong in _site_map
####
#def load_all(site_dir, site, hbond_stick_rad=0.015, lig="",
#             C_color="green", hphob_stick_rad=0.0):
#  site_pref = site_dir + "/" + site
#  site = site_map(site_pref + "_s.csv")
#  site.load_points()
#  if(hbond_stick_rad > 0.0): site.draw(stick_rad=hbond_stick_rad)
#  cmd.load(site_pref + "_rad.pdb", site + "_prot")
#  cmd.color(C_color, "%s & e. C" % (site + "_prot"))
#  if(not lig == ""):
#    cmd.load(lig, site + "_lig")
#    cmd.color(C_color, "%s & e. C" % (site + "_lig"))
#    cmd.show("sticks", site + "_lig")
#
#  return site

#def draw_msms_surf(vert_fname, face_fname, obj_name, mesh=True,
#                   color=[1.0, 1.0, 1.0], lw=1):
#
#  try:
#    vert_file = open(vert_fname, "r")
#  except IOError, (errno, strerror):
#    print "Unable to open the file", vert_fname
#    print "error(%s): %s" % (errno, strerror)
#    return
#
#  try:
#    face_file = open(face_fname, "r")
#  except IOError, (errno, strerror):
#    print "Unable to open the file", face_fname
#    print "error(%s): %s" % (errno, strerror)
#    return
#
#  # Load verts & normals
#  lineno = 0
#  verts = []
#  normals = []
#  for line in vert_file:
#    lineno += 1
#    if(lineno <= 3): continue
#
#    toks = [ float(s) for s in line.rstrip("\n").split() ]
#    verts.append(toks[0:3])
#    normals.append(toks[3:6])
#  vert_file.close()
#  verts = numpy.array(verts)
#
#
#  # Load faces
#  lineno = 0
#  face_idz = []
#  for line in face_file:
#    lineno += 1
#    if(lineno <= 3): continue
#
#    toks = [ int(s) for s in line.rstrip("\n").split() ]
#    face_idz.append(toks[0:3])
#  face_file.close()
#  # vertices are 1 indexed in the face file -- we need 0 indexed
#  face_idz = numpy.array(face_idz) - 1
#
#  if(mesh):
#    edges = {}
#    corr_lines = [LINEWIDTH] 
#    corr_lines.append(lw)
#    corr_lines.extend([ BEGIN, LINES, COLOR])
#    corr_lines.extend(color)
#    print corr_lines
#    
#    for face in face_idz:
#      for i in range(3): 
#        end_idx = (i+1)%3 
#        e0 = "%d_%d" % (face[i], face[end_idx])
#        e1 = "%d_%d" % (face[end_idx], face[i])
#        if(not e0 in edges and not e1 in edges):
#          edges[e0] = True
#          edges[e1] = True
#          corr_lines.append(VERTEX)
#          corr_lines.extend(verts[face[i]])
#          corr_lines.append(VERTEX)
#          corr_lines.extend(verts[face[end_idx]])
#          print VERTEX, verts[face[i]], VERTEX, verts[face[end_idx]]
#    print END
#    corr_lines.append(END)
#    cmd.load_cgo(corr_lines, obj_name)
################################################################################
#
#################################################################################
