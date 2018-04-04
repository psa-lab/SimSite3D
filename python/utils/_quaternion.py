import numpy
from numpy import random

class Quaternion:
  """
  Simple quaternion class that (in the end) should mirror the
  SimSite3D::Quaternion C++ class
  """


  def __init__(self, Q=numpy.array([]), R=numpy.array([]),
               cos_theta=0.0, sin_theta=0.0, V=numpy.array([])):
    """
    Initializer for Quaternion class

    Parameters
    ----------
    self : object -- passed automagically by python
      Quaternion class 
    Q : numpy.array (4x1), optional
      array holding a quaternion
    R : numpy.array (3x3), optional
      array holding an orthonormal rotation matrix
    """
#    R_is_idnt = False
#    V_is_empty = True
#    if(R.shape[0] == R.shape[1] == 3 and 
#       R[0][0] == R[1][1] == R[2][2] == 1.0 and
#       R[0][1] == R[0][2] == R[1][0] == 0.0 and 
#       R[1][2] == R[2][0] == R[2][1] == 0.0):
#      R_is_idnt = True

    self.q = numpy.array([1.0, 0.0, 0.0, 0.0])

    if(V.shape[0] == 3):
      self.set_using_angles(cos_theta, sin_theta, V)

    elif(R.shape[0] == 3 and R.shape[1] == 3): 
      R = R.reshape((9,))

      tmp = numpy.zeros(4)
      tmp[0] = 0.5 * numpy.sqrt(1 + R[0] + R[4] + R[8])
      tmp[1] = 0.5 * numpy.sqrt(1 + R[0] - R[4] - R[8])
      tmp[2] = 0.5 * numpy.sqrt(1 - R[0] + R[4] - R[8])
      tmp[3] = 0.5 * numpy.sqrt(1 - R[0] - R[4] + R[8])

      idx = 0
      max_val = tmp[0]
      for i in range(1,4):
        if(tmp[i] > max_val):
          max_val = tmp[i]
          idx = i

      q = numpy.zeros(4)
      q[idx] = tmp[idx]
      if(idx == 0):
        q[1] = 0.25 * (R[7] - R[5])/tmp[0]
        q[2] = 0.25 * (R[2] - R[6])/tmp[0]
        q[3] = 0.25 * (R[3] - R[1])/tmp[0]
      elif(idx == 1):
        q[0] = 0.25 * (R[7] - R[5])/tmp[1]
        q[2] = 0.25 * (R[3] + R[1])/tmp[1]
        q[3] = 0.25 * (R[2] + R[6])/tmp[1]
      elif(idx == 2):
        q[0] = 0.25 * (R[2] - R[6])/tmp[2]
        q[1] = 0.25 * (R[3] + R[1])/tmp[2]
        q[3] = 0.25 * (R[7] + R[5])/tmp[2]
      elif(idx == 3):
        q[0] = 0.25 * (R[3] - R[1])/tmp[3]
        q[1] = 0.25 * (R[2] + R[6])/tmp[3]
        q[2] = 0.25 * (R[7] + R[5])/tmp[3]
      self.q = q

    elif(Q.shape[0] == 4): self.q = numpy.copy(Q)

    self.is_unit = False

  def set_using_angles(self, cos_theta, sin_theta, V):
    """
    Rodrigues formula for quaternion -- but parameters are for the full angles
    and not half angles.  Use half angle identities to get the desired values.

    Assumes V is a unit vector
    """
    if(cos_theta > 1 - 1E-07): self.q = numpy.array([1.0, 0.0, 0.0, 0.0])
    elif(cos_theta < -1 + 1E-07): 
      # Assume sin(theta) = 0; this implies sin(theta/2) = 1.0
      self.q[0] = 0.0
      self.q[1:4] = V[:]
    else:
      cos_theta_plus_1 = cos_theta + 1.0
      cos_half_theta = numpy.sqrt(0.5 * cos_theta_plus_1)
      sin_half_theta = sin_theta / numpy.sqrt(2 * cos_theta_plus_1)
      self.q[0] = cos_half_theta
      for i in range(3): self.q[i+1] = sin_half_theta * V[i]

  def set_from_parameters(self, theta0, theta1, s):
    sigma0 = numpy.sqrt(1.0 - s)
    sigma1 = numpy.sqrt(s)
    self.q[0] = numpy.cos(theta1) * sigma1
    self.q[1] = numpy.sin(theta0) * sigma0
    self.q[2] = numpy.cos(theta0) * sigma0
    self.q[3] = numpy.sin(theta1) * sigma1

    # really should be unit, but just check
    self.is_unit = False

  def __repr__(self):
    return "[%f, %f, %f, %f]" % (self.q[0], self.q[1], self.q[2], self.q[3])

  def __mul__(self, other):
    """
    Multiplication of self * other -- not self modifying
    """
    q = self.q
    R = numpy.array([ [ q[0], -q[1], -q[2], -q[3] ],
                      [ q[1],  q[0], -q[3],  q[2] ],
                      [ q[2],  q[3],  q[0], -q[1] ],
                      [ q[3], -q[2],  q[1],  q[0] ] ])
    return Quaternion(Q = numpy.dot(R, other.q))

  def get_ortho_rot_mat(self):
    """
    Get an orthonormal rotation matrix that represents the same rotation as 
    the stored quaternion.  

    Parameters
    ----------
    self : object -- passed automagically by python
      Quaternion class 

    Returns
    R : numpy.array (3,3)
      orthonormal rotation matrix that is equivalent to self.q when used
      to premultiply coordinates (i.e. Rx = y and not xR = y')
    """

    R = numpy.zeros(9)
    if(not self.is_unit): self.normalize()

    q = self.q
    q_2 = q * q
    q0qi = q[0] * q[1:4]
    qxqy = q[1]*q[2]
    qxqz = q[1]*q[3]
    qyqz = q[2]*q[3]

    R[0] = q_2[0] + q_2[1] - q_2[2] - q_2[3]
    R[1] = 2.0*(qxqy - q0qi[2])
    R[2] = 2.0*(qxqz + q0qi[1])
    R[3] = 2.0*(qxqy + q0qi[2])
    R[4] = q_2[0] - q_2[1] + q_2[2] - q_2[3]
    R[5] = 2.0*(qyqz - q0qi[0])
    R[6] = 2.0*(qxqz - q0qi[1])
    R[7] = 2.0*(qyqz + q0qi[0])
    R[8] = q_2[0] - q_2[1] - q_2[2] + q_2[3]

    return R.reshape((3,3))

  def normalize(self):
    # Change the first component to be positive -- multiply quaternion by -1.0
    if(self.q[0] < 0.0): self.q *= -1.0

    # Make the quaternion a unit quaternion
    self.q /= numpy.sqrt(sum(self.q*self.q))
    self.is_unit = True

  def conjugate(self):
    """
    Compute the complex conjugate the the quaterion, does not modify self
    but rather returns the conjugate
    """
    q_conj = self.q.copy()
    q_conj[1:4] *= -1.0
    return Quaternion(Q=q_conj)

  def distance(self, other, W=1.0):
    if(not self.is_unit): self.normalize()
    if(not other.is_unit): other.normalize()

    rv = sum(self.q * other.q)
    if(rv < 0.0): rv *= -1.0
    return W * (1.0 - rv)

  def parameters(self):
    theta0 = numpy.arctan(self.q[1] / self.q[2])
    if(self.q[2] < 0.0): theta0 += numpy.pi
    theta1 = numpy.arctan(self.q[3] / self.q[0])
    if(self.q[2] < 0.0): theta1 += numpy.pi
    s = self.q[3]*self.q[3] + self.q[0]*self.q[0]

    return (theta0, theta1, s)

  def randomize(self):
    """
    Generate a random quaternion as per the method of James J. Kuffner, 
    ICRA 2004.  

    This method is computationally expensive since it calls sqrt, cos and sin.
    If this method is used for something other than testing, it may be useful
    to speed it up.
    """

    s = random.uniform()
    sigma0 = numpy.sqrt(1.0 - s)
    sigma1 = numpy.sqrt(s)
    theta0 = 2*numpy.pi * random.uniform()
    theta1 = 2*numpy.pi * random.uniform()
    self.q[0] = numpy.cos(theta1) * sigma1
    self.q[1] = numpy.sin(theta0) * sigma0
    self.q[2] = numpy.cos(theta0) * sigma0
    self.q[3] = numpy.sin(theta1) * sigma1
    if(self.q[0] < 0.0): self.q *= -1.0
    self.is_unit = True


#  # This function has not been verified to work correctly or as expected
#  def randomize(self, theta_beg = 0.0, theta_end = numpy.pi):
#    v = numpy.array([ random.uniform(-1.0, 1.0), random.uniform(-1.0, 1.0),
#                      random.uniform(-1.0, 1.0) ])
#    v /= numpy.sqrt(sum(v * v))
#
#    tmp = random.uniform(-1.0, 1.0)
#    half_theta = 0
#    if(tmp > 0): half_theta = theta_beg + tmp*(theta_end - theta_beg)
#    else: half_theta = -1.0*theta_beg + tmp*(theta_end - theta_beg)
#    half_theta *= 0.5
#
#    cos_half_theta = numpy.cos(half_theta)
#    sin_half_theta = numpy.sin(half_theta)
#    self.q[0] = cos_half_theta
#    self.q[1] = v[0] * sin_half_theta
#    self.q[2] = v[1] * sin_half_theta
#    self.q[3] = v[2] * sin_half_theta
#
#    self.is_unit = False

def gen_orientations(N, C, half_side_len = 0.5, q_width = 0.1,
                     q_dist_lim = 1.0 - numpy.cos(numpy.pi / 48.0), 
                     t_dist_lim = 0.05, overlap=0.20):
  """
  Generate N random, rigid orientations 
  """


  # not sure how do this at this point
#  if(not C.shape[1] == 3): 
#    print "C needs to be a 1x3 matrix (row vector)"
#    return []

  t_dist_lim_squared = t_dist_lim * t_dist_lim
  half_side_len = (1.0 + overlap) * half_side_len
  pi = numpy.pi

  # Need to be careful with overlap as I haven't checked if (s) is cyclic
  lower_s_bound = 1.0 - q_width - 0.5*q_width*overlap
  upper_s_bound = 1.0
  upper_t_bound = 0.5 * q_width * (1.0 + overlap)
  lower_t_bound = -1.0 * upper_t_bound

  orientations = []
  parameters = []
  for n in range(N):
    Q = Quaternion()
   
    # Need to be careful with overlap as I haven't checked if (s) is cyclic
    # Further checking shows that (s) is not directly cyclic, but there could
    # be a more complex mapping
    vv = numpy.array([ random.uniform(lower_t_bound, upper_t_bound),
                       random.uniform(lower_t_bound, upper_t_bound),
                       random.uniform(lower_s_bound, upper_s_bound) ])

    Q.set_from_parameters(2.0*pi * vv[0], 2.0*pi * vv[1], vv[2])
    R = Q.get_ortho_rot_mat()

    # Note: to get the transformation we need to apply the rotation to the
    # centroid.  The translation compoment is the vector needed to move the
    # translated centroid to its original position
    T = -1.0 * numpy.dot(C, R) + C

    # Adjust the translation by some amount
    too_close = True
    while(too_close):
      T_adj = numpy.array([ random.uniform(-1.0, 1.0), random.uniform(-1.0,1.0),
                            random.uniform(-1.0, 1.0) ])
      T += half_side_len * T_adj
      too_close = False
      for O in orientations:
        (prev_Q, prev_T) = O
        q_dist = Q.distance(prev_Q)
        tmp = T - prev_T
        sq_t_dist = sum(tmp*tmp)
#        if(q_dist < q_dist_lim and sq_t_dist < t_dist_lim_squared):
        if(q_dist_lim < q_dist and sq_t_dist < t_dist_lim_squared):
          too_close = True
          break

    orientations.append((Q,T))
    parameters.append((vv, half_side_len * T_adj))
  return (orientations, parameters)



if __name__ == "__main__":

  R = numpy.array([-0.826285, 0.425134, -0.369480, -0.561396, -0.568387, 
                    0.601474, 0.045700, 0.704413, 0.708317])
  R = R.reshape((3,3))
  print R
  q = Quaternion(R=R)
  print q.get_ortho_rot_mat()

  print q.conjugate()
  print q.randomize()
  print q.q

  centroid = numpy.array([10.0, 20.0, -10.0])
  Oz = gen_orientations(10, centroid, q_width=0.05, half_side_len=0.5, 
                        overlap=0.0)
  for o in Oz:
    (q, t) = o
    params = q.parameters()
    x = [params[i] for i in range(3)]
    x[0] /= 2.0*numpy.pi
    x[1] /= 2.0*numpy.pi
    print x


