import numpy

def spherical_grid(r, d):
  if(d > r): return []

  # Calculate approximate angle theta to partition pi into partitions of
  # equal size and the distances of the points is close to d.
  # Taken from Wolfram page on triangles
  theta = numpy.arccos(1.0 - 0.5 * (d/r)**2.0)
  n = numpy.round(numpy.pi / theta)
  npts = n + 1
  odd = False
  if(numpy.mod(npts,2)): odd = True

  # Grid points from positive y axis to x axis in the XY plane.
  theta = numpy.pi/2;
  list_x = []
  for ii in range(int(numpy.floor(npts/2))):
    list_x.append([r*numpy.cos(theta), r*numpy.sin(theta), 0.0])
    theta -= numpy.pi/n;

  # Fill in shell
  X = numpy.array(list_x)
  last_jj = int(numpy.floor(npts/2))
  for jj in numpy.arange(1, last_jj):
    d = 4.0*jj
    step = 2.0*numpy.pi / d
    theta = 0
    for kk in numpy.arange(1, d):
      theta += step
      tmp = numpy.array([[ numpy.cos(theta), 0.0, -1.0*numpy.sin(theta)],
                         [ 0.0,              1.0, 0.0                  ],
                         [ numpy.sin(theta), 0.0, numpy.cos(theta)     ]])
      list_x.append(numpy.dot(tmp, X[jj,:].T).tolist())

  if(odd):
    theta = numpy.pi/2.0 - numpy.pi/n;
    for ii in numpy.arange(1, npts-1):
      list_x.append([ r * numpy.sin(theta), 0.0, r*numpy.cos(theta) ])
      theta -= numpy.pi / n

  # Flip points about the X axis
  X = numpy.array(list_x)
  flip_R = [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]]

  list_x.extend(numpy.dot(X, flip_R))

  # add in the poles
  if(odd): list_x.extend([ [r, 0.0, 0.0], [-1.0*r, 0.0, 0.0] ])
  X = numpy.array(list_x)

  #for i in range(X.shape[0]):
    #print "%f %f %f" % (X[i,0], X[i,1], X[i,2])
  return X

# Create demo for testing
if __name__ == '__main__':
  pts = []
  #pts.extend(spherical_grid(3.0, 0.25).tolist())
  #pts.extend(spherical_grid(3.25, 0.25).tolist())
  #pts.extend(spherical_grid(3.5, 0.25).tolist())
  #pts.extend(spherical_grid(2.75, 0.25).tolist())
  #pts.extend(spherical_grid(2.5, 0.25).tolist())
  pts.extend(spherical_grid(3.0, 0.5).tolist())
  pts.extend(spherical_grid(3.5, 0.5).tolist())
  pts.extend(spherical_grid(2.5, 0.5).tolist())



  #import pymol
  #from pymol.cgo import *
  #from pymol import cmd
  #pymol.finish_launching()

  my_cgo = []
  x = numpy.array([1,0,0])
  for i in range(len(pts)):

    p = numpy.array(pts[i])
    p = p / numpy.sqrt(sum(p*p))
    if(sum(p*x) >= 0.5):
      my_cgo.extend([COLOR, 0.0, 0.0, 1.0])
      print "hit"
    else:
      my_cgo.extend([COLOR, 1.0, 0.0, 0.0])
      print "missed"
    my_cgo.append(SPHERE)
    my_cgo.extend(pts[i])
    my_cgo.append(0.05)


  cmd.load_cgo(my_cgo, "sphere_pts")
