#!/usr/bin/python
#
# $Source: /psa/share/repository/pfizer_proj/install.py,v $
# $Revision: 1.19 $
# $Author: vanvoor4 $
# $Date: 2008-07-30 16:58:41 $
#
# $Log: not supported by cvs2svn $
# Revision 1.18  2008/01/04 21:39:25  vanvoor4
# Misspelled word
#
# Revision 1.17  2008/01/04 18:14:46  vanvoor4
# Changed the comments printed at the end to be more helpful and
# updated the variable names.
#
# Revision 1.16  2008/01/04 16:12:24  vanvoor4
# Added the copying of simsite3d.conf to simsite3d.conf.example
#
# Revision 1.15  2007/11/14 17:21:52  vanvoor4
# Removed some extra end of line characters
#
# Revision 1.14  2007/11/14 17:00:25  vanvoor4
# Added some comments about installing diverse sitemaps
#
# Revision 1.13  2007/11/07 17:28:47  vanvoor4
# Added comment to ask for installation of diverse sitemaps
#
# Revision 1.12  2007/11/01 19:40:37  vanvoor4
# Removed copying the perl code since it is not updated and has
# really been superceded by the C++ code itself (except for
# auto_gen_sitemaps.pl which is replaced by auto_gen_sitemaps.py).
#
# Revision 1.11  2007/11/01 19:31:22  vanvoor4
# Typo -- pro instead of prot
#
# Revision 1.10  2007/11/01 19:16:02  vanvoor4
# Changed $SIMSITE3D_SOFTWARE_DIR to $SIMSITE3D_INSTALL_DIR
#
# Revision 1.9  2007/11/01 16:47:04  vanvoor4
# Removed the part corresponding to the installing of the
# diverse sitemaps
#
# Revision 1.8  2007/10/22 14:08:21  vanvoor4
# Examples dir must not be blown away either...
#
# Revision 1.7  2007/10/22 13:57:05  vanvoor4
# Do not want to blow away the python dir if the install dir is the same as
# the source dir.
#
# Revision 1.6  2007/10/02 20:13:46  vanvoor4
# Added copying of perl scripts to bin
#
# Revision 1.5  2007/10/02 18:50:08  vanvoor4
# Added python directories to be installed
#
# Revision 1.4  2007/09/26 14:53:56  vanvoor4
# Added a few more section print statements
#
# Revision 1.3  2007/09/26 14:41:52  vanvoor4
# Forgot that aclocal, etc have some output on fresh installs.
#
# Revision 1.2  2007/09/26 14:37:29  vanvoor4
# Moved python dir from ./src to ./
#
# Revision 1.1  2007/09/24 13:52:29  vanvoor4
# Initial checkin.  Used Python instead of Perl because I am now more
# familiar with the string manipulations of Python.
#
#
#

class DIRS:
  pass

import sys
import os
import shutil
sys.path.append(os.getcwd() + "/python")
import py_src.utils as utils

dir_mode = 0754

dirs = DIRS()
dirs.proj_dir = os.getenv("SIMSITE3D_INSTALL_DIR")
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
if(not utils.system("/usr/bin/make install")):
  sys.exit(1)

print "copying SimSite3D python directory...",
sys.stdout.flush()
shutil.copy("python/py_src/non_SimSite3DPy_auto_gen_sitemaps.py", 
            dirs.proj_dir + "/bin/auto_gen_sitemaps.py");
# Do not want to blow away the python directory if we are installing to the 
# current directory -- if this becomes a problem, we should look at how to
# correctly install python programs
#if(os.getcwd() != dirs.proj_dir):
#  if(os.path.isdir(dirs.proj_dir + "/python")):
#    shutil.rmtree(dirs.proj_dir + "/python")
#  shutil.copytree("python", dirs.proj_dir + "/python")
if(not os.path.isdir(dirs.proj_dir + "/python")):
  os.mkdir(dirs.proj_dir + "/python", dir_mode)
shutil.copy("python/py_src/utils.py", dirs.proj_dir + "/python/utils.py");
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
shutil.copy("params/simsite3d.conf", 
            params_dir + "/simsite3d.conf.example");
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

print "\n\nPlease install the diverse sitemaps and ligands:\n"
print "  The tar files for the diverse sitemaps and ligands may be found in:"
print "\t" + dirs.proj_dir + "/data\n"
print "  The file diverse_sitemaps.tgz should be extracted to " + \
  "$SIMSITE3D_DIVERSE_SITES"
print "  The file diverse_ligands.tgz should be extracted to " + \
  "$SIMSITE3D_DIVERSE_LIGS\n\n"
print "Please consult the SimSite3D Quick Guide for help in configuring and",
print "using the\nSimSite3D software tools"
print "To configure SimSite3D, two examples files are included:"
print "  An example SimSite3D configuration file is:"
print "    " + params_dir + "/simsite3d.conf.example"
print "  An example external prot-lig scoring functions config file is:"
print "    " + params_dir + "/ext_pro_lig_score_fcns.conf.example\n\n"
print "To get started, you may look at the examples in\n  " + \
  dirs.proj_dir + "/examples\n  and peruse the SimSite3D quick guide\n";
