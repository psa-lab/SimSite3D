from pymol import cmd
from pymol import stored
import pymol
import numpy

def __init__(self):
  cmd.extend("xform_obj", xform_obj)

def xform_obj(sel, R, T, state = 1):
  M = numpy.array(R)
  M = M.reshape([3,3])

  objName = cmd.identify(sel, state)[0][0]
  stored.pos = []
  cmd.iterate_state(state, objName, "stored.pos.append((x,y,z))")

  nrows = len(stored.pos)
  Y = numpy.dot(stored.pos, M.T) + numpy.tile(T, [nrows, 1])
  stored.pos = Y.tolist()
  cmd.alter_state(state, objName, "(x,y,z)=stored.pos.pop(0)")

#cmd.extend("xform_obj", xform_obj)
