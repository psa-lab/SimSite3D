from numpy import *
from SimSite3DPy.sitemaps._Sitemap import *
import cPickle
from os import path

################################################################################

################################################################################
# define the loader
def load_hbond_caps(self):
  atoms_fname = self.atoms_file_name()
  pkl_fname = "%s.pkl" % (atoms_fname[:-8])
  if(path.exists(pkl_fname)):
    try:
      pkl_in = open(pkl_fname, "rb")
    except IOError, (errno, strerror):
      print "Unable to open the file", fname
      print "error(%s): %s" % (errno, strerror)
      print "in %s" % (self.__module__)
      return

    objs = cPickle.load(pkl_in)
    if(len(objs) == 2):    
      (self.hbond_vols, self.hbond_surf_caps) = objs
    elif(len(objs) == 3):
      (self.hbond_vols, self.hbond_surf_caps, self.metal_vols) = objs
    elif(len(objs) == 4):
      (self.hbond_vols, self.hbond_surf_caps, self.metal_vols, 
       self.metal_surfs) = objs
    pkl_in.close()
    return True
  else:
    print  """
Unable to load the hbond caps file 
  (%s) 
because it doesn't exist
""" % (pkl_fname)
    return False
################################################################################

################################################################################

# Add the functions to the sitemap class
sitemap.load_hbond_caps = load_hbond_caps
