class DIRS:
  pass

import sys
import os
import shutil
sys.path.append(os.getcwd() + "/python")
import py_src.utils as utils

dir_mode = 0754

dirs = DIRS()
#dirs.proj_dir = os.getenv("SIMSITE3D_INSTALL_DIR")
dirs.proj_dir = "/soft/linux64/SimSite3D_v4.5"
if(dirs.proj_dir == None):
  print >> sys.stderr, \
    "\nUnable to get the environment variable $SIMSITE3D_INSTALL_DIR"
  print >> sys.stderr, "Please set $SIMSITE3D_INSTALL_DIR to a valid path"
  sys.exit(1)

# This directory can screw things up between machines
if(os.path.isdir("autom4te.cache")):
  utils.system("/bin/rm -rf autom4te.cache")

print "initializing SimSite3D C/C++ install environment... "
# shouldn't always call flush but the commands won't go in order
sys.stdout.flush()
if(not utils.system("/usr/bin/aclocal")):
  sys.exit(1)
if(not utils.system("/usr/bin/libtoolize --copy --force --automake")):
  sys.exit(1)
if(not utils.system("/usr/bin/automake -a -c")):
  sys.exit(1)
if(not utils.system("/usr/bin/autoconf")):
  sys.exit(1)

print "configuring SimSite3D C/C++ install environment... "
sys.stdout.flush()
if(not utils.system("./configure --prefix=\"" + dirs.proj_dir + "\"")):
  sys.exit(1)
print "building SimSite3D C/C++ executables... "
sys.stdout.flush()
if(not utils.system("/usr/bin/make")):
  sys.exit(1)

cwd = os.getcwd()
os.chdir("python")
os.system("tools/boost-jam-3.1.17-1-linuxx86/bjam build")
os.chdir(cwd)
