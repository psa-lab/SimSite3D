import numpy
from numpy import linalg
from ASCbasePy.utils._quaternion import *

# silly little function to make code easier to read
def between3(min_pos, pos, max_pos):
  """
  min_pos : lower bound
  pos : position to test
  max_pos : upper bound
  """
  for i in range(3):
    if(pos[i] < min_pos[i] or max_pos[i] < pos[i]): return False
  return True

def rad2deg(r):
  return 180.0 * (r / numpy.pi)

def deg2rad(d):
  return (d / 180.0) * numpy.pi

# Assumes numpy arrays
def dist(A,B):
  tmp = numpy.array(A) - numpy.array(B)
  return numpy.sqrt(sum(tmp*tmp))

# Assumes numpy arrays
def squared_dist(A,B):
  tmp = numpy.array(A) - numpy.array(B)
  return sum(tmp*tmp)

# Assumes numpy arrays
def unit_vec(tail, head):
  tmp = head - tail
  return tmp / numpy.sqrt(sum(tmp*tmp))

# Assumes numpy arrays
def angle(A,B,C):
  BA = unit_vec(B,A)
  BC = unit_vec(B,C)
  return rad2deg(numpy.arccos(sum(BA * BC)))

def normalize(U):
  return U / numpy.sqrt(sum(U*U))

def dist_to_line(V, p0, pt):
  """
  Compute the distance from the 3D point "pt" to the 3D line defined by
  the unit vector "V" and the point "p0"
  """

  p0_to_pt = numpy.array(pt) - numpy.array(p0)
  proj = sum(p0_to_pt * numpy.array(V))
  return numpy.sqrt(sum(p0_to_pt * p0_to_pt) - proj*proj)

def dist2line(pt, Vi, Vj, line_segment=False):
  L = Vj - Vi
  edge_len = sqrt(dot(L, L))
  normalize(L)

  # Vector from Vi to pt
  ViP0 = pt - Vi

  # Project pt on the line
  proj_pt = dot(ViP0, L)

  if(line_segment == False): pt_on_line = Vi + proj_pt * L
  else:
    if(proj_pt < 0.0): pt_on_line = Vi
    elif(proj_pt > edge_len): pt_on_line = Vj
    else: pt_on_line = Vi + proj_pt * L

  return (dist(pt, pt_on_line), pt_on_line)

def line_intersect_sphere(c, r, m, b):
  """
  Does the sphere with center c and radius r intersect with the line y = sm + b?
  (Here y, m, and b are vectors and s is a scalar)
  We do this by finding the closest point from the sphere's center to the line
   
  Assumptions:
    c,m,b are of type numpy.array((3,))
    r is a scalar
    m is a unit vector
  """
  # Get vector from point on line to the center (bc)
  b_to_c = numpy.array(c) - numpy.array(b)
  # Get the square of the length of the vector bc
  b_to_c_sq_len = sum(b_to_c * b_to_c)
  # Get the projection of bc onto m
  proj_len = sum(m * b_to_c)
  
  # Compute the squared distance from the center to the line
  h2 = b_to_c_sq_len - proj_len*proj_len
  r2 = r*r
#  print "h2", h2, "   r2", r2
  # if the line is outside of the sphere, return
  if(h2 >= r2): return (numpy.array([]), numpy.array([]))

  # We know the radius and the "height" of the triangle, get the length of the
  # side of the triangle that is lying on the line y = mx + b
  d = numpy.sqrt(r2 - h2)
  pt = b + proj_len * m
  return (pt - d*m, pt + d*m)

#  # quadratic equation x**2 + Bx + C = 0
#  B = sum(m * b_to_c)
#  C = sum(b_to_c * b_to_c) - r**2
#  rad_arg = B*B - C
#
#  # At present we are considering "touching" as a non intersection
#  if(rad_arg <= 0): return (numpy.array([]), numpy.array([]))
#
#  rad_val = sqrt(rad_arg)
#  return (-1.0*B + rad_val, -1.0*B - rad_val)

def dist_to_plane(three_pts, pt):
  BA = unit_vec(three_pts[1], three_pts[0]) 
  BC = unit_vec(three_pts[1], three_pts[2]) 
  N = numpy.cross(BA,BC)
  N = normalize(N)
  
  A_p0 = pt - three_pts[0]
  return sum(A_p0 * N)

def intersect_planes(N1, d1, N2, d2):
  """
  Compute the intersection of planes defined as NX + d = 0
  """
  from numpy import linalg

  m = numpy.cross(N1,N2)
  m /= numpy.sqrt(sum(m*m))
  b = numpy.zeros((3,))

  # Find the largest component of the vector parallel to the line of
  # intersection 
  max_m_i = numpy.max(m)
  idx = numpy.argmax(m) 

  # Find a point on the line
  A = numpy.zeros((4,))
  j = 0
  for i in range(3):
    if(i == idx): continue

    A[j] = N1[i]
    A[j+2] = N2[i]
    j += 1

  det_A = A[0]*A[3] - A[2]*A[1]
  A_inv = numpy.array([ A[3] / det_A, -1.0*A[1] / det_A,
                        -1.0*A[2] / det_A, A[0] / det_A])

  # Determine the values for the 2 kept coordinates
  j = 0
  for i in range(3):
    if(i == idx): continue
  
    b[i] = -1.0 * (A_inv[j] *d1 + A_inv[j+1] * d2)
    j += 2

  # Determine the value of b for the dropped coordinate
  # Assumption: b[idx] is zero
#  b[idx] = (sum(N1 * b) + d1) / (-1.0 * N1[idx])
#  print "b[idx]", b[idx]
#  b[idx] = 0.0
#  b[idx] = (sum(N2 * b) + d2) / (-1.0 * N2[idx])
#  print "b[idx]", b[idx]

  return (m,b)


def gaussian_K(x, mu, sigma):
  return numpy.exp(-0.5 * ((x-mu)/sigma)**2)

# Assume numpy arrays
def simple_lse_fit(A, B):
  """
  Find rigid alignment that aligns B to A
  """
  if(A.shape[0] != B.shape[0]):
    print "A and B must have the same number of points"
    return ([], [])

  A_cent = numpy.mean(A, 0)
  B_cent = numpy.mean(B, 0)

  M_sums = numpy.zeros([3,3])
  for a,b in zip(A,B):
    M_sums += numpy.dot(a[:, numpy.newaxis], b[numpy.newaxis, :])
  M_sums /= numpy.tile([A.shape[0]], [3,3])
  M_sums -= numpy.dot(A_cent[:, numpy.newaxis], B_cent[numpy.newaxis, :])
  M_sums = numpy.reshape(M_sums, [9,])
  
  N = numpy.empty([16])
  N[0] = M_sums[0] + M_sums[4] + M_sums[8]
  N[5] = M_sums[0] - M_sums[4] - M_sums[8]
  N[10] = M_sums[4] - M_sums[8] - M_sums[0]
  N[15] = M_sums[8] - M_sums[0] - M_sums[4]
  N[1] = N[4] = M_sums[5] - M_sums[7]
  N[2] = N[8] = M_sums[6] - M_sums[2]
  N[3] = N[12] = M_sums[1] - M_sums[3]
  N[6] = N[9] = M_sums[1] + M_sums[3]
  N[7] = N[13] = M_sums[6] + M_sums[2]
  N[11] = N[14] = M_sums[5] + M_sums[7]

  [V,D] = linalg.eig(N.reshape([4,4]))
  # Assume linalg.eig will aways return unit eigenvectors
  idx = numpy.argmax(V)
  Q = D[:,idx]
  if(Q[0] < 0.0):
    Q = -1.0 * Q

  QQ = 2 * Q*Q
  xy = 2*Q[1]*Q[2]
  wz = 2*Q[0]*Q[3]
  xz = 2*Q[1]*Q[3]
  wy = 2*Q[0]*Q[2]
  yz = 2*Q[2]*Q[3]
  wx = 2*Q[0]*Q[1]

  R = numpy.array([[ 1 - QQ[2] - QQ[3], xy - wz, xz + wy ],
                   [ xy + wz, 1 - QQ[1] - QQ[3], yz - wx ],
                   [ xz - wy, yz + wx, 1 - QQ[1] - QQ[2] ]])
  T = A_cent - numpy.dot(B_cent, R)
  return (R,T)

# Get the transformation to move from a local coordinate system to global 
# coordinates.  The local coordinate system is defined by B at the origin, A on
# the positive x-axis, and C in the lower half plane (y <= 0.0).
def get_local_coords(A, B, C):
  """
  Returns (R,T) where R(x-T) = y will transform the global coord x in to the
  local coord y.  

  This means to transfrom from local coords y to global coords x we do
  y(R^t) + T = x
  """
  BA_dist = dist(B, A)
  BC_dist = dist(B, C)
  BA_vec = unit_vec(B, A)
  BC_vec = unit_vec(B, C)
  cos_angle = sum(BA_vec * BC_vec)
  sin_angle = numpy.sqrt(1.0 - cos_angle*cos_angle)
  local_coords = [ [BA_dist, 0.0, 0.0], [0.0, 0.0, 0.0],
                   [cos_angle*BC_dist, -1.0*sin_angle*BC_dist, 0.0] ]
  local_coords = numpy.array(local_coords)
#  print "local coords:", local_coords
  global_coords = numpy.array([A, B, C])
#  print "global coords:", global_coords
  return simple_lse_fit(global_coords, local_coords)

def align_planes(X_norm, Y_norm):
  """
  Compute the rotation that brings X_norm parallel to Y_norm
  """
  cos_phi = sum(Y_norm * X_norm)
  Q = Quaternion()
  if(cos_phi < 1 - 1E-07):
    C = numpy.cross(Y_norm, X_norm)
    sin_phi = numpy.sqrt(sum(C*C))
    C /= sin_phi
    Q = Quaternion(cos_theta=cos_phi, sin_theta=sin_phi, V=C)
  return Q

def align_to_X_axis(V):
  """
  Compute the rotation the brings the projection of V onto the postive X axis
  """
  V = numpy.array(V)
  V /= numpy.sqrt(sum(V*V))
  # V[0] = cos(theta), V[1] = sin(theta)
  Q = Quaternion()
  if(V[0] < 1 - 1E-07):
    unit_Z = numpy.array([0.0, 0.0, 1.0])
    Q = Quaternion(cos_theta=V[0], sin_theta=V[1], V=unit_Z)
  return Q

# The alignment method has been updated to be consistent (at least in my mind)
# First the polar atom is placed at the origin.
# Second use the cross product of the vectors originating at the polar atom
# and pointing to its neighbors to determine the "correct" normal direction
# Third rotate the translated coordinates to bring the "normal direction" 
# parallel to the positive Z axis
# Finally rotate about the positive Z axis to bring the polar atom to carbon
# neighbor direction to lie on the positive X axis
def get_local_coords2(A, B, C):
  # Place point B at the origin
  BA_dist = dist(B, A)
  BC_dist = dist(B, C)
  BA_vec = unit_vec(B, A)
  BC_vec = unit_vec(B, C)

  # Compute the right handed normal to the plane defined by A, B, C and 
  # vectors BA and BC
  N = numpy.cross(BA_vec, BC_vec)
  N /= numpy.sqrt(sum(N * N))

  # compute the rotation that brings N parallel to the local positive Z axis
  unit_Z = numpy.array([0.0, 0.0, 1.0])
  Q1 = align_planes(N, unit_Z)
#  print "Q1", Q1
  R1 = Q1.get_ortho_rot_mat()
#  print "R1", R1
#  print "BA", BA_vec
#  print "BC", BC_vec
  BA_vec_prime = numpy.dot(BA_vec, R1)
  BC_vec_prime = numpy.dot(BC_vec, R1)
#  print "BA'", BA_vec_prime
#  print "BC'", BC_vec_prime
   
  # Rotate about the local Z axis to bring BA to lie on the positve X axis
  Q2 = align_to_X_axis(BA_vec_prime)

  # This is dumb -- I still am having some trouble reasoning with quaternions
  # rotation matrices and order of operations -- ah so here it is
  # With my messing around I computed Q1, Q2, and possibly Q3 but used the
  # the conjugate of Q1 & Q2 (Q3 = Q3.conj).  
  # We have (Q1 * Q2 * Q3).conj = Q3.conj * Q2.conj * Q1.conj

  R2 = Q2.get_ortho_rot_mat()
  BC_vec_pp = numpy.dot(BC_vec_prime, R2)
  # require the Y component (in the local coordinates) of the BC vector to be
  # negative
  if(BC_vec_pp[1] > 0.0):
    # flip about X axis
    Q3 = Quaternion(Q=numpy.array([0.0, 1.0, 0.0, 0.0]))
    R = (Q1 * Q2 * Q3).conjugate().get_ortho_rot_mat()
  else: R = (Q1 * Q2).conjugate().get_ortho_rot_mat()

  return (R,B)


if(__name__ == "__main__"):
  center = numpy.array([0, 0, 0])
  radius = 1.0
  m = numpy.array([1.0, 0.0, 0.0])
  b = numpy.array([3.0, 0.0, 0.0])
  print line_intersect_sphere(center, radius, m, b)













