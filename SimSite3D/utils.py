#!/usr/bin/python
#
# $Source: /psa/share/repository/pfizer_proj/python/SimSite3D/utils.py,v $
# $Revision: 1.4 $
# $Author: vanvoor4 $
# $Date: 2008-11-18 15:40:32 $
#
# $Log: not supported by cvs2svn $
# Revision 1.3  2008/01/22 16:46:56  vanvoor4
# Changed to reflect the more recent method of handling the SimSite3D
# environment.  At this time the course of action is
# 0) Get $ASCBASE_INSTALL_DIR value
# 1) read the /etc/ascbase/ascbase.conf file if it exists
# 2) read the $HOME/.ascbase/ascbase.conf file if it exists
# 3) Get the ASCBASE environment variables
# 4) use the parameters provided on the command line
#
# Revision 1.2  2007/11/01 19:11:47  vanvoor4
# Changed the order of processing environment and the names
# of some variables
#
# Revision 1.1  2007/09/26 14:36:24  vanvoor4
# Initial checkin
#
#
#

# Is this better than importing os and sys? 
from os import path
from os import getenv
from os import system
from sys import stderr

#####
# Simple class to hold the directories
#####
class directories:
  fail = False

  def __init__(self, check_install_dirs = True):
    self.proj_dir = self.__getenv_dir("ASCBASE_INSTALL_DIR")
    if(self.proj_dir == None):
      return None

    #self.config_dir = self.proj_dir + "/SimSite3DSoftParams"
    #if(check_install_dirs and not os.path.isdir(self.config_dir)):
      #print >> sys.stderr, "The directory \"" + self.config_dir + \
        #"\" does not exist"
      #self.fail = True
      #return None
    self.bin_dir = self.proj_dir + "/bin"
    if(check_install_dirs and not path.isdir(self.bin_dir)):
      print >> stderr, "The directory \"" + self.bin_dir + \
        "\" does not exist"
      self.fail = True
      return None

  
    vars = {"ASCBASE_DBASE_PROTS":None, "ASCBASE_DBASE_LIGS":None, 
            "ASCBASE_DBASE_SITES":None, "ASCBASE_DIVERSE_SITES":None, 
            "ASCBASE_DIVERSE_LIGS":None, "ASCBASE_PROJ_OUTPUT":None, 
            "ASCBASE_SCRATCH_DIR":None}

    # Load the defaults form the system wide environment file
    #sys_env_file = open(self.config_dir + "/ascbase_software.conf")
    #for l in sys_env_file.readlines():
      #if(l.startswith("%") or l.startswith("#")):
        #continue
    #
      ## Assumption is the conf file is a whitespace delimited file
      #toks = l.split()
      #if(toks[0] in vars):
        #vars[toks[0]] = toks[1]

    # Load parameters from conf files if they exist
    if(path.exists("/etc/ascbase/ascbase.conf")):
      self.__load_param_file("/etc/ascbase/ascbase.conf", vars)
    local_conf = self.__getenv_dir("HOME")
    if(path.exists(local_conf + "/.ascbase/ascbase.conf")):
      self.__load_param_file(local_conf + "/.ascbase/ascbase.conf", vars)

    # Load the values set in the environment
    for var in vars:
      tmp_val = getenv(var)
      if(not tmp_val == None):
        vars[var] = tmp_val

    self.prot_dir = vars["ASCBASE_DBASE_PROTS"]
    self.ligs_dir = vars["ASCBASE_DBASE_LIGS"]
    self.dbase_dir = vars["ASCBASE_DBASE_SITES"]
    self.diverse_dir = vars["ASCBASE_DIVERSE_SITES"]
    self.diverse_ligs = vars["ASCBASE_DIVERSE_LIGS"]
    self.results_dir = vars["ASCBASE_PROJ_OUTPUT"]

    if(not vars["ASCBASE_SCRATCH_DIR"] == None):
      self.scratch_dir = vars["ASCBASE_SCRATCH_DIR"]
    else:
      self.scratch_dir = self.results_dir

  ##########

  ##########
  def __load_param_file(self, fname, vars):
    # Load the values from the parameter file
    sys_env_file = open(fname, "r")

    for l in sys_env_file.readlines():
      l_end = len(l)
      comment_chars = "%#"
      for cc in comment_chars:
        pos = l.find(cc)
        if(not pos == -1 and l_end > pos):
          l_end = pos
    
      # Assumption is the conf file is a whitespace delimited file
      toks = l[:l_end].split()
      if(len(toks) == 2 and toks[0] in vars):
        vars[toks[0]] = toks[1]

  ##########

  ##########
  def __getenv_dir(self, var_name):
    val = getenv(var_name)
    if(val == None):
      self.fail = True
      print >> stderr, \
        "\nUnable to get the environment variable $" + var_name
      print >> stderr, "Please set $" + var_name + " to a valid path"
    else:
      if(self.__verify_dir(var_name, val)):
        return val
      else:
        self.fail = True

    return None  
  ##########

  ##########
  def __verify_dir(self, var_name, val):
    if(path.isdir(val)):
      return True
    else:
      print >> stderr, "\nThe value of $" + var_name + " is \"" + val + \
        "\",\nbut that path is not a directory."
      print >> stderr, "Please set $" + var_name + " to a valid path"
      return False 
  ##########

# Currently use os.system since I have not been able to get pipes to work 
# correctly
def system(the_cmd):
  import os
  status = os.system(the_cmd)
  if(status):
    print >> stderr, "\n  ERROR!\n  The return from\n" + the_cmd + \
      "\n  was %d" % (status)
    return False
  return True




