from numpy import *
import cgo_items
from pymol import cmd
from pymol.cgo import *
import pymol
pymol.finish_launching()

def subdivide(delta):
 
  deltas_out = []
  mid_pts = []
  j=-1
  for i in range(3): 
    mid_pts.append( (delta[i] + delta[j]) / 2.0 )
    j += 1

  # add center triangle
  deltas_out.append(mid_pts)

  for i in range(3):
    deltas_out.append([ delta[i], mid_pts[(i+1) % 3], mid_pts[i] ])

  return array(deltas_out)


#pts = array([ [1,0,0], [-1,0,0], [0,1,0], [0,-1,0], [0,0,1], [0,0,-1] ])

#tris =  \
#[
#  [ [0,0,1], [0,1,0], [-1,0,0] ],
#  [ [0,0,1], [1,0,0], [0,1,0]  ],
#  [ [0,0,1], [0,-1,0], [1,0,0] ],
#  [ [0,0,1], [-1,0,0], [0,-1,0]],
#
#  [ [0,0,-1], [-1,0,0], [0,1,0] ],
#  [ [0,0,-1], [0,1,0], [1,0,0] ],
#  [ [0,0,-1], [1,0,0], [0,-1,0] ],
#  [ [0,0,-1], [0,-1,0], [-1,0,0] ]
#]
tris =  \
[
  [ [0,0,1], [-1,0,0], [0,1,0] ],
  [ [0,0,1], [0,1,0], [1,0,0] ],
  [ [0,0,1], [1,0,0], [0,-1,0]],
  [ [0,0,1], [0,-1,0], [-1,0,0] ],

  [ [0,0,-1], [0,1,0], [-1,0,0] ],
  [ [0,0,-1], [1,0,0], [0,1,0]],
  [ [0,0,-1], [0,-1,0], [1,0,0]],
  [ [0,0,-1], [-1,0,0], [0,-1,0]]
]
tris = array(tris)

for i in range(3):
  new_tris = []
  for delta in tris:
    new_tris.extend(subdivide(delta))
  tris = array(new_tris)

#def draw_points(points, pymol_label, color=[1.0, 1.0, 1.0], radius=0.05):
# Push pts out to surface of sphere -- since it is a unit sphere just normalize
# them


cgo = []
cgo.extend([BEGIN, TRIANGLES])
cgo.extend([COLOR, 1.0, 0.2, 0.2]) 
for delta in tris:
  for pt in delta:
    pt /= sqrt(sum(pt*pt))
    cgo.append(VERTEX)
    cgo.extend(3*pt)
  for pt in delta:
    cgo.append(NORMAL)
    cgo.extend(pt)
   
cgo.append(END)
cmd.load_cgo(cgo, "blue")
