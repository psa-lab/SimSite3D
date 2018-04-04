import sys
from numpy import *
from ASCbasePy.utils import *

def load_CA_positions(pdb_fname, res_nums):
  prot = pdb.residues(pdb_fname)
 
  positions = [] 
  for res in prot:
    if(not res.resSeq in res_nums or not res.chainID == "A"): continue

    for a in res:
      if(a.name == " CA "):
        positions.append(a.position)
        break 
       
  return array(positions)
  
################################################################################

################################################################################
def compute_items(cent_pos, top_pos, bottom_pos):
  # Do LSE fit in future

  # Create plane in ref
  U = unit_vec(cent_pos[0], cent_pos[1])
  V = unit_vec(cent_pos[0], cent_pos[2])
  N = cross(U,V)
  N /= sqrt(dot(N,N))
  P0 = cent_pos[0]

  # Drop point to plane -- do distance from pt to plane for all top & bottom
  V = top_pos - tile(P0, (top_pos.shape[0],1))
  signed_D = dot(V, N)
  adj = tile(signed_D, (3,1)).T * tile(N, (cent_pos.shape[0], 1))
  top_pt_on_plane = top_pos - adj
  top_dists = abs(signed_D)

  V = bottom_pos - tile(P0, (bottom_pos.shape[0],1))
  signed_D = dot(V, N)
  adj = tile(signed_D, (3,1)).T * tile(N, (cent_pos.shape[0], 1))
  bottom_pt_on_plane = bottom_pos - adj
  bottom_dists = abs(signed_D)

  return (N, P0, top_dists, bottom_dists, top_pt_on_plane, bottom_pt_on_plane)

################################################################################

################################################################################


bilayer_residues = [43, 108, 149, 199, 241, 285, 326, 356, 393, 424, 465, 502]
top_residues = [56, 92, 157, 186, 259, 272, 335, 342, 403, 414, 477, 488]
bottom_residues = [26, 129, 136, 215, 227, 307, 313, 371, 379, 445, 450, 522]

from optparse import OptionParser
cmd_parser = OptionParser()
cmd_parser.add_option("", "--ref_prot", help="Reference PDB file",
                      metavar="<FILE>.pdb")
cmd_parser.add_option("", "--other_prot", help="Structure to compare with reference",
                      metavar="<FILE>.pdb")
if(len(sys.argv) < 5):
  print >> sys.stderr, "\n\t * All flags and arguments are required *"
  print >> sys.stderr, "\t   (Except help of course)\n"
  cmd_parser.print_help()
  sys.exit(-1)
(cmd_options, cmd_args) = cmd_parser.parse_args()

cent_ref_pos = load_CA_positions(cmd_options.ref_prot, bilayer_residues)
top_ref_pos = load_CA_positions(cmd_options.ref_prot, top_residues)
bottom_ref_pos = load_CA_positions(cmd_options.ref_prot, bottom_residues)
cent_other_pos = load_CA_positions(cmd_options.other_prot, bilayer_residues)
top_other_pos = load_CA_positions(cmd_options.other_prot, top_residues)
bottom_other_pos = load_CA_positions(cmd_options.other_prot, bottom_residues)


(ref_N, ref_P0, ref_top_dists, ref_bottom_dists, ref_top_pts_on_plane,
 ref_bottom_pts_on_plane) = \
  compute_items(cent_ref_pos, top_ref_pos, bottom_ref_pos)

(other_N, other_P0, other_top_dists, other_bottom_dists, other_top_pts_on_plane,
 other_bottom_pts_on_plane) = \
  compute_items(cent_other_pos, top_other_pos, bottom_other_pos)

# Get unit vector for reference helices' tops
ref_top_dirs = top_ref_pos - cent_ref_pos
ref_top_dirs /= tile(sqrt(sum(ref_top_dirs * ref_top_dirs, 1)), (3,1)).T

# Get unit vectors for reference helices' bottoms
ref_bottom_dirs = bottom_ref_pos - cent_ref_pos
ref_bottom_dirs /= \
  tile(sqrt(sum(ref_bottom_dirs * ref_bottom_dirs, 1)), (3,1)).T

# Get unit vectors for other helices' tops relative to ref center
other_top_dirs = top_other_pos - cent_ref_pos
other_top_dirs /= tile(sqrt(sum(other_top_dirs * other_top_dirs, 1)), (3,1)).T
  
# Get unit vectors for other helices' bottoms relative to ref center
other_bottom_dirs = bottom_other_pos - cent_ref_pos
other_bottom_dirs /= \
  tile(sqrt(sum(other_bottom_dirs * other_bottom_dirs, 1)), (3,1)).T

# Helical bending (change in)
print ref_top_dirs
print other_top_dirs
top_cos_of_angle_deltas = sum(ref_top_dirs * other_top_dirs, 1)
bottom_cos_of_angle_deltas = sum(ref_bottom_dirs * other_bottom_dirs, 1)

# Helical distance (change in)
top_dist_deltas = ref_top_dists - other_top_dists
bottom_dist_deltas = ref_bottom_dists - other_bottom_dists

print "helix#|%s|" % ( "|".join([ "%d" % (i) for i in range(1,13) ]) )
print "top delta angle|%s|" % ( "|".join([ "%.01f" % (f) for f in arccos(top_cos_of_angle_deltas) * (180.0/pi) ]))
print "top delta distance|%s|" % \
  ( "|".join([ "%.02f" % (f) for f in top_dist_deltas ]))


print "bottom delta angle|%s|" % ( "|".join([ "%.01f" % (f) for f in arccos(bottom_cos_of_angle_deltas) * (180.0/pi) ]))
print "bottom delta distance|%s|" % \
  ( "|".join([ "%.02f" % (f) for f in bottom_dist_deltas ]))


 




