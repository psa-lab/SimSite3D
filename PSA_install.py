class DIRS:
  pass

import sys
import os
import shutil
sys.path.append(os.getcwd() + "/python")
import py_src.utils as utils

dir_mode = 0755

dirs = DIRS()
dirs.proj_dir = "/soft/linux64/SimSite3D_v4.5"
dirs.diverse_dir = dirs.proj_dir + "/diverse_sites"
if(not os.path.exists(dirs.diverse_dir)):
  os.mkdir(dirs.diverse_dir, dir_mode)


#if(dirs.proj_dir == None):
#  print >> sys.stderr, \
#    "\nUnable to get the environment variable $ASCBASE_INSTALL_DIR"
#  print >> sys.stderr, "Please set $ASCBASE_INSTALL_DIR to a valid path"
#  sys.exit(1)
#
## This directory can screw things up between machines
#if(os.path.isdir("autom4te.cache")):
#  utils.system("/bin/rm -rf autom4te.cache")
#
#print "initializing SimSite3D C/C++ install environment... "
## shouldn't always call flush but the commands won't go in order
#sys.stdout.flush()
#if(not utils.system("/usr/bin/aclocal")):
#  sys.exit(1)
#if(not utils.system("/usr/bin/libtoolize --copy --force --automake")):
#  sys.exit(1)
#if(not utils.system("/usr/bin/automake -a -c")):
#  sys.exit(1)
#if(not utils.system("/usr/bin/autoconf")):
#  sys.exit(1)
#
#print "configuring SimSite3D C/C++ install environment... "
#sys.stdout.flush()
#if(not utils.system("./configure --prefix=\"" + dirs.proj_dir + "\"")):
#  sys.exit(1)
#print "building SimSite3D C/C++ executables... "
#sys.stdout.flush()
#if(not utils.system("/usr/bin/make install")):
#  sys.exit(1)
#
#print "copying SimSite3D python directory...",
#sys.stdout.flush()
#shutil.copy("python/py_src/non_SimSite3DPy_auto_gen_sitemaps.py", 
#            dirs.proj_dir + "/bin/auto_gen_sitemaps.py");
# Do not want to blow away the python directory if we are installing to the 
# current directory -- if this becomes a problem, we should look at how to
# correctly install python programs
#if(os.getcwd() != dirs.proj_dir):
#  if(os.path.isdir(dirs.proj_dir + "/python")):
#    shutil.rmtree(dirs.proj_dir + "/python")
#  shutil.copytree("python", dirs.proj_dir + "/python")

#if(not os.path.isdir(dirs.proj_dir + "/python")):
#  os.mkdir(dirs.proj_dir + "/python", dir_mode)
#shutil.copy("python/py_src/utils.py", dirs.proj_dir + "/python/utils.py");
#print "finished"
#sys.stdout.flush()

print "Installing the SimSite3D C++ programs...",
if(not utils.system("/usr/bin/make install")):
  sys.exit(1)
print "finished"
sys.stdout.flush()

print "setting up SimSite3D parameters directory...",
sys.stdout.flush()
params_dir = dirs.proj_dir + "/SimSite3DSoftParams"
if(not os.path.isdir(params_dir)):
  os.mkdir(params_dir, dir_mode)
shutil.copy("params/optimum_hbonds.dat", params_dir);
shutil.copy("params/sparse_hbonds.dat", params_dir);
shutil.copy("params/minimal_hbonds.dat", params_dir);
# Copy the example conf file -- do not want to overwrite an existing conf file
shutil.copy("params/ext_prot_lig_score_fcns.conf", 
            params_dir + "/ext_prot_lig_score_fcns.conf.example");
shutil.copy("params/ascbase.conf", 
            params_dir + "/ascbase.conf.example");
print "finished"
sys.stdout.flush()

print "copying SimSite3D examples directory...",
sys.stdout.flush()
if(os.getcwd() != dirs.proj_dir):
  if(os.path.isdir(dirs.proj_dir + "/examples")):
    shutil.rmtree(dirs.proj_dir + "/examples")
  shutil.copytree("examples", dirs.proj_dir + "/examples")
print "finished"
sys.stdout.flush()

print "Installing the normalization database files ... ",
cwd = os.getcwd()
shutil.copy("data/diverse_ligands.tbz2", dirs.diverse_dir)
shutil.copy("data/diverse_sitemaps.tbz2", dirs.diverse_dir)
os.chdir(dirs.diverse_dir)
os.system("tar -xvjf diverse_ligands.tbz2")
os.system("tar -xvjf diverse_sitemaps.tbz2")
os.unlink("diverse_ligands.tbz2")
os.unlink("diverse_sitemaps.tbz2")
os.chdir(cwd)
print "finished"
sys.stdout.flush()

#print "\n\nPlease install the diverse sitemaps and ligands:\n"
#print "  The tar files for the diverse sitemaps and ligands may be found in:"
#print "\t" + dirs.proj_dir + "/data\n"
#print "  The file diverse_sitemaps.tgz should be extracted to " + \
  #"$ASCBASE_DIVERSE_SITES"
#print "  The file diverse_ligands.tgz should be extracted to " + \
  #"$ASCBASE_DIVERSE_LIGS\n\n"

print "Please consult the SimSite3D Quick Guide for help in configuring and",
print "using the\nSimSite3D software tools"
print "To configure SimSite3D, two examples files are included:"
print "  An example SimSite3D configuration file is:"
print "    " + params_dir + "/ascbase.conf.example"
print "  An example external prot-lig scoring functions config file is:"
print "    " + params_dir + "/ext_pro_lig_score_fcns.conf.example\n\n"
print "To get started, you may look at the examples in\n  " + \
  dirs.proj_dir + "/examples\n  and peruse the SimSite3D quick guide\n";
