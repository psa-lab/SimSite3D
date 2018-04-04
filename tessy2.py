from numpy import *
import cgo_items
from pymol import cmd
from pymol.cgo import *
import pymol
pymol.finish_launching()

def sort_stuff(A,B):
  if(A[1] < B[1]): return -1
  elif(A[1] == B[1]): return 0
  return 1

def draw_triangles(tris, obj_name="", color=[1.0, 0.2, 0.2], state=0, 
                   do_print=False):
  cgo = [COLOR]
  cgo.extend(color)
  cgo.extend([BEGIN, TRIANGLES])
  for delta in tris:
    for pt in delta:
      pt /= sqrt(sum(pt*pt))
      cgo.append(VERTEX)
      cgo.extend(3.0*pt)
      if(do_print and pt[2] <= 0.55):
        print "pt", pt
    for pt in delta:
      cgo.append(NORMAL)
      cgo.extend(pt)
    
  cgo.append(END)
  cmd.load_cgo(cgo, obj_name, state=state)


def subdivide(delta):
 
  deltas_out = []
  mid_pts = []
  j=-1
  for i in range(3): 
    mid_pts.append( (delta[i] + delta[j]) / 2.0 )
    j += 1
  mid_pts = array(mid_pts)
  for i in range(3): 
    mid_pts[i] /= sqrt(sum(mid_pts[i] * mid_pts[i]))

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

angle = 0.0
level = 1
tris = []
N = 5
S = sqrt(0.75)
#S = 1.0
T = sqrt(1.0 - S*S)
for i in range(N):
  tri = [ [1.0, 0.0, 0.0] ]
  tri.append([T, S*cos(angle), S*sin(angle)])
  angle += 2*pi / N
  tri.append([T, S*cos(angle), S*sin(angle)])
  tris.append(tri)


#tris =  \
#[
#  [ [0,0,1], [-1,0,0], [0,1,0] ],
#  [ [0,0,1], [0,1,0], [1,0,0] ],
#  [ [0,0,1], [1,0,0], [0,-1,0]],
#  [ [0,0,1], [0,-1,0], [-1,0,0] ],
#]
tris = array(tris)

for i in range(level):
  new_tris = []
  for delta in tris:
    new_tris.extend(subdivide(delta))
  tris = array(new_tris)

#def draw_points(points, pymol_label, color=[1.0, 1.0, 1.0], radius=0.05):
# Push pts out to surface of sphere -- since it is a unit sphere just normalize
# them

draw_triangles(tris, "full_hemisphere")
#cgo = []
#cgo.extend([BEGIN, TRIANGLES])
#cgo.extend([COLOR, 1.0, 0.2, 0.2]) 
#for delta in tris:
#  for pt in delta:
#    pt /= sqrt(sum(pt*pt))
#    cgo.append(VERTEX)
#    cgo.extend(3*pt)
#  for pt in delta:
#    cgo.append(NORMAL)
#    cgo.extend(pt)
#   
#cgo.append(END)
#cmd.load_cgo(cgo, "blue", state=1)

# Now remove all triangles with all three points below the cap
kept_tris = []
N = array([1.0, 0.0, 0.0])
p0 = array([0.5, 0.0, 0.0])

import cgo_items
cgo_items.draw_plane2(N, 3.0*p0)#, plane_num=0, state=0, color=[0.7, 0.7, 0.7]):

for tri in tris:
  signed_dists = zeros((3,))
  for pt,i in zip(tri, range(3)):
    V = pt - p0
    signed_dists[i] = sum(N * V)

  if(signed_dists[0] >= 0.0 and signed_dists[1] >= 0.0 and 
     signed_dists[2] >= 0.0):
    kept_tris.append(tri)
  elif(signed_dists[0] < 0.0 and signed_dists[1] < 0.0 and 
       signed_dists[2] < 0.0):
    print "dropped a triangle"
    continue
  else:
    for i in range(3):
      if(signed_dists[i] > 0.0): continue

      # Project pt to plane
      tri[i] -= signed_dists[i] * N

      U = tri[i] - p0
      U /= sqrt(sum(U*U))
      # radius here is not 1 -- it is sqrt(0.75) = 1.0**2 - 0.5**2
      tri[i] = p0 + sqrt(0.75)*U
      print "point", tri[i]
      print "len(pt)", sqrt(sum(tri[i]*tri[i]))

#    for i,pt in zip(range(3), tri):
#      if(signed_dists[i] > 0.0): continue
#
#      # Project pt to plane
#      pt -= signed_dists[i] * N
#
#      U = pt - p0
#      U /= sqrt(sum(U*U))
#      # radius here is not 1 -- it is sqrt(0.75) = 1.0**2 - 0.5**2
#      pt = p0 + sqrt(0.75)*U
#      print "point", pt
#      print "len(pt)", sqrt(sum(pt*pt))

    kept_tris.append(tri)


kept_tris = array(kept_tris)

draw_triangles(kept_tris, "untrimmed_cap", color=[0.2, 0.2, 1.0], state=2,
               do_print=True)

# Need to "spit" out the unique pts & triangles/edges
# make obscene use of a map
V = {}
T = []
for tri in kept_tris:

  t = []
  for pt in tri:
    for i in range(3):
      if(-1E-07 <= pt[i] and pt[i] <= 1E-07): pt[i] = 0.0
    pt *= 3.0
    my_key = "%.7f,%.7f,%.7f," % (pt[0], pt[1], pt[2])
    if(not my_key in V): V[my_key] = len(V)
    t.append(V[my_key])
  T.append(t)

V = [ (k, idx) for k,idx in V.iteritems() ]
V.sort(cmp=sort_stuff)

for v in V:
  #print v[0]
  s = v[0][:]
  x = array([ float(z) for z in s.split(",")[:-1] ])
  if(x[0] > 1.7426771):
    print v[0]
  else:
    x[0] = 1.5
    val = sqrt(6.75/(x[1]*x[1] + x[2]*x[2]))
    x[1] *= val
    x[2] *= val
    print "%f,%f,%f," % (x[0], x[1], x[2])
#1) 1.5*1.5 + sx[1]*sx[1] + sx[2]*sx[2] = 9
#2) sx[1]*sx[1] + sx[2]*sx[3] = 6.75
#3) s*2 = 6.75/(x[1]*x[1] + x[2]*x[2])
#4) s = sqrt(6.75/(x[1]*x[1] + x[2]*x[2]))



print
for t in T:
  print "%d,%d,%d," % (t[0], t[1], t[2])



center = [17.6715,5.29787,3.32608]

pts = array([[18.3656,4.08951,2.24554],
[18.6117,5.61862,2.13245],
[17.7773,6.86699,2.52851],
[16.4229,6.48531,3.54732],
[16.2604,4.94367,3.8691],
[17.0148,3.69079,3.2626]])

my_tris = array([ [center, pts[i], pts[(i + 1) % 6]] for i in range(6)])
print my_tris
#draw_triangles(my_tris, "one_ring_example", color=[0.2, 0.8, 0.2], state=2,
               #do_print=True)
my_cgo = []
#my_cgo.extend([BEGIN, TRIANGLES])
#my_cgo.extend([BEGIN, POINTS])
val = 0.0
for delta in my_tris:
  val += 1.0/6.0
  for pt in delta:
    my_cgo.extend([COLOR, 1.0-val, val, 0.2])
    my_cgo.append(SPHERE)
    print pt
    my_cgo.extend(pt)
    my_cgo.append(0.25)
#    my_cgo.append(VERTEX)
#    my_cgo.extend(pt)
    #my_cgo.append(NORMAL)
    #my_cgo.extend([1.0, 0.0, 0.0])
#my_cgo.append(END)
cmd.load_cgo(my_cgo, "one_ring_example")




