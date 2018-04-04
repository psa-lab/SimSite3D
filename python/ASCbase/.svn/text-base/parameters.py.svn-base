import os

##########
class parameters_base:

  ##########

  ##########
  def __init__(self):
    self.dbase_sites = ""
    self.dbase_ligs = ""
    self.dbase_prots = ""
    self.diverse_sites = ""
    self.diverse_ligs = ""
    self.proj_output = ""
    self.scratch_dir = ""
    
    self.load_surf_files = False
    self.fail = False
    
    self.install_dir = self.__getenv_dir("ASCBASE_INSTALL_DIR")
  ##########

  # Parameters are not verified since they can be overwritten by subclasses
  ##########
  def get_params(self):
    if(self.fail):
      return False  

    vars = {"ASCBASE_DBASE_PROTS":None, "ASCBASE_DBASE_LIGS":None,
            "ASCBASE_DBASE_SITES":None, "ASCBASE_DIVERSE_SITES":None,
            "ASCBASE_DIVERSE_LIGS":None, "ASCBASE_PROJ_OUTPUT":None,
            "ASCBASE_SCRATCH_DIR":None}

    # Load parameters from conf files if they exist
    if(os.path.exists("/etc/ascbase/ascbase.conf")):
      self.__load_param_file("/etc/ascbase/ascbase.conf", vars)
    local_conf = self.__getenv_dir("HOME")
    if(os.path.exists(local_conf + "/.ascbase/ascbase.conf")):
      self.__load_param_file(local_conf + "/.ascbase/ascbase.conf", vars)

    # Load the values set in the environment
    for var in vars:
      tmp_val = os.getenv(var)
      if(not tmp_val == None):
        vars[var] = tmp_val

    self.dbase_prots = vars["ASCBASE_DBASE_PROTS"]
    self.dbase_ligs = vars["ASCBASE_DBASE_LIGS"]
    self.dbase_sites = vars["ASCBASE_DBASE_SITES"]
    self.diverse_sites = vars["ASCBASE_DIVERSE_SITES"]
    self.diverse_ligs = vars["ASCBASE_DIVERSE_LIGS"]
    self.proj_output = vars["ASCBASE_PROJ_OUTPUT"]
    self.scratch_dir = vars["ASCBASE_SCRATCH_DIR"]

    if(self.proj_output == None):
      self.proj_output = os.getcwd()

    if(self.scratch_dir == None):
      self.scratch_dir = self.proj_output
  ##########

  ##########
  def version(self):
    print '''ASCbase Software Development Version
Copyright (C) 2006-2008, Michigan State University (MSU) Board of Trustees.
All rights reserved.
Written by Jeffrey R. Van Voorst and Leslie A. Kuhn

      '''
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
    val = os.getenv(var_name)
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
    if(os.path.isdir(val)):
      return True
    else:
      print >> stderr, "\nThe value of $" + var_name + " is \"" + val + \
        "\",\nbut that path is not a directory."
      print >> stderr, "Please set $" + var_name + " to a valid path"
      return False
  ##########

##########
class search_parameters(parameters_base):

  ##########
  def __init__(self):
    parameters_base.__init__(self)
    parameters_base.get_params(self)

    self.prot_lig_score = "None"
    self.ofname = ""
    self.normalize = True
    self.num_scores_to_keep = 1
    self.score_cutoff = -1.5
    self.min_num_atoms = 5
    self.dmetol = 0.3
    self.lsetol = 0.3
    self.ligand_rmsd = False
    self.sitemap_rmsd = False
    self.write_ligands = True
    self.align_to_query = True
    self.num_rand_aligns = 0
    self.time_process = False
    self.hydrophobic_query = False
    self.query_sitemap = ""
    self.dbase_sitemap = ""
    self.db_index_fname = ""
    self.dbstart = -1
    self.dbstop = -1
  ##########

  ##########
  def cmdline_options(self, argv):
    import getopt
    long_opts = ["version", "usage", "help", "time", "prot_lig_score=", 
                 "proj_output=", "dbase_sites=", "dbase_ligs=", "scratch_dir=", 
                 "lig_frag_size=", "keep_n_scores=", "score_threshold=", 
                 "db_index_file=" "dbstart=", "dbstop=",
                 "no_normalization", "SCORE_ONLY", "hydrophobic_query"]
    try:
      (opts, args) = getopt.getopt(argv[1:], "vuht:o", long_opts)
    except getopt.GetoptError, err:
      self.version()
      print str(err)
      self.usage()
      self.fail = True
      return False

    if(len(args) == 0 or len(args) > 2):
      self.version()
      print "search_sitemaps.py requires 1 or 2 site maps"
      self.usage()
      self.fail = True
      return False
    else:
      self.query_sitemap = args[0]
      if(len(args) == 2):
        self.dbase_sitemap = args[1] 
      

    for opt, arg in opts:
      if(opt in ("-v", "--version") ):
        self.version()
      elif(opt in ("-u", "--usage") ):
        self.version()
        self.usage()
      elif(opt in ("-h", "--help") ):
        self.version()
        self.help()
      elif(opt == "-o"):
        self.ofname = arg
      elif(opt == "--prot_lig_score"):
        self.prot_lig_score = arg
      elif(opt == "proj_output"):
        self.proj_output = arg
      elif(opt == "--dbase_sites"):
        self.dbase_sites = arg
      elif(opt == "--dbase_ligs"):
        self.dbase_ligs = arg
      elif(opt == "--scratch_dir"):
        self.scratch_dir = arg
      elif(opt == "--lig_frag_size"):
        try:
          self.min_num_atoms = int(arg)
        except ValueError:
          print "--lig_frag_size requires an nonnegative integer argument"
          self.fail = True
          return False
      elif(opt == "--keep_n_scores"):
        try:
          self.num_scores_to_keep = int(arg)
        except ValueError:
          print "--keep_n_scores requires a positive integer argument"
          self.fail = True
          return False
      elif(opt == "--score_threshold"):
        try:
          self.score_cutoff = float(arg)
        except ValueError:
          print "--score_threshold requires a real number argument"
          self.fail = True
          return False
      elif(opt == "--db_index_file"):
        self.db_index_fname = arg
      elif(opt == "--dbstart"):
        try:
          self.dbstart = int(arg)
        except ValueError:
          print "--dbstart requires a positive integer argument"
          self.fail = True
          return False
      elif(opt == "--dbstop"):
        try:
          self.dbstop = int(arg)
        except ValueError:
          print "--dbstop requires a positive integer argument"
          self.fail = True
          return False
      elif(opt == "--no_normalization"):
        self.normalize = False
      elif(opt == "--SCORE_ONLY"):
        self.align_to_query = False
      elif(opt == "--hydrophobic_query"):
        self.hydrophobic_query = True

    return self.verify_options() 
  ##########

  ##########
  def verify_options(self):
    if(self.fail):
      return False

    if(not self.test_directory(self.dbase_sites, "database site maps")):
      return False

    if(self.min_num_atoms < 0):
      print "minimum number of ligand fragment atoms may not be negative"
      return False

    if(self.min_num_atoms > 0 and 
       not self.test_directory(self.dbase_ligs, "database ligands")):
      return False

    if(self.normalize and 
       not (self.test_directory(self.diverse_sites, 
                                "diverse site maps directory") and
            self.test_directory(self.diverse_ligs, 
                                "diverse ligands directory"))):
      return False

    if(not self.test_directory(self.proj_output, "project output directory")):
      return False

    if(not self.test_directory(self.scratch_dir, "scratch directory")):
      return False

    if(self.num_scores_to_keep < 1):
      print "Number of hits to keep per database site must be positive"
      return False

    if(len(self.db_index_fname) and not os.path.isfile(self.db_index_fname)):
      print "%s is not a valid file" % (self.db_index_fname)
      print "Please set the db index filename to a valid index file"
      return False

    if(self.dbstart > self.dbstop):
      print "dbstart may not be greater than dbstop"
      return False

    # A bit of kludge but hopefully no one tries to run with -1 for both values
    if(self.dbstart < -1 or self.dbstop < -1):
      print "dbstart and dbstop must be positive numbers"
      return False

    # Need to support these fast
    #self.prot_lig_score = "None"
    # self.ofname = "" -- nothing to check here?

    # support these for development & testing at MSU
    #self.write_ligands = True
    #self.num_rand_aligns = 0

    # these can probably go away
    #self.ligand_rmsd = False
    #self.sitemap_rmsd = False
    return True
  ##########

  ##########
  def test_directory(self, dir_name, dir_string):
    if(not os.path.isdir(dir_name)):
      print "%s is not a valid directory\nPlease set %s to a valid directory" \
        % (dir_name, dir_string)
      self.fail = True 
      return False

    return True
  ##########

  ##########
  def usage(self):
     print '''
search_sitemaps.py [-vuhto] [--version --usage --help --time 
      --prot_lig_score SF_NAME --dbase_sites DIR --dbase_ligs DIR 
      --scratch_dir DIR --lig_frag_size NUM --keep_n_scores NUM 
      --score_threshold REAL --db_index_file FILE --dbstart NUM --dbstop NUM 
      --no_normalization --SCORE_ONLY --hydrophobic_query ]

Database search: 
      search_sitemaps.py [OPTIONS] <query_sitemap>.csv 
Site comparison: 
      search_sitemaps.py [OPTIONS] <query_sitemap>.csv <sitemap2>.csv
        '''
  ##########

  ##########
  def help(self):
     print '''
Name:
      search_sitemaps.py - ASCbase Software tool to compare binding sites

Database search: 
      search_sitemaps.py [OPTIONS] <query_sitemap>.csv 

Site comparison: 
      search_sitemaps.py [OPTIONS] <query_sitemap>.csv <sitemap2>.csv

Description:
      search_sitemaps.py is a tool in the ASCbase Software package that takes
      a query site map and either searches for similar sites in a site map 
      database (directory) or compares the query to a second site map.  For
      more information on generating site maps see create_sitemaps.py and the
      ASCbase documentation.

Note:
      The hits returned by search_sitemaps depends on both the maximum number
      of hits per database site map as well as the score threshold.  In order
      for a hit to be returned it must always meet the score threshold.  If
      more than one hit is found for a database site map, the top N hits for
      that site map are returned (where N is the maximum number of hits per
      database site map).

Directory and File Assumptions:
      To be reasonably flexible, ASCbase requires that the database site maps,
      database ligands, project output, diverse sites, and diverse ligands
      directories be set explicitly.  The typical or default settings for
      these values are taken first from the configuration files and then from
      environment variables.  However, these values, except for the diverse
      site maps and ligands directories, can be set by command line options.

      ASCbase assumes that each protein-ligand binding site has a unique
      identifier represented here by "'XXXXXX".  The files used by ASCbase are:
            XXXXXX_p.pdb		Protein file for XXXXXX
            XXXXXX_l.mol2		Ligand file for XXXXXX
            XXXXXX_s.csv		Site map file for XXXXXX
            XXXXXX_s.pdb		Site map points file for XXXXXX
            XXXXXX_a.pdb		Site map atoms file for XXXXXX
            XXXXXX_rad.pdb		Site map residues file for XXXXXX

Options:
      Most of the options have default values set either in configuration 
      files and/or by environmental variables.

General Options:
      -v,--version Print version information and exit
      -h,--help    Print this help message
      -u,--usage   Print a short help message
      -t,--time    Track the execution time of the search

Search Options:
      The following options affect the search parameters used by 
      search_sitemaps.  Attention is required when using these options as they
      can drastically affect the outcome of a search.

      --lig_frag_size NUMBER [default: 5]
            Set the minimum number of ligand atoms for a ligand fragment.
            This option has two different effects.  If the ligand fragment
            size is set to zero, binding sites that do not have bound ligands
            will also be searched.  If the ligand fragment size is positive,
            it is assumed that you are very interested in the database ligands 
            bound to sites similar to the query site.  In this case, each hit 
            will be guarenteed to have one or more ligand fragments inside the
            query site, and each fragment will have at least the minimum number
            of heavy atoms (i.e. no empty pockets).

      --keep_n_scores NUMBER [default: 1]
            Set the maximum number of hits (results) per database sitemap.

      --score_threshold REAL_NUMBER [default: -1.5]
            Set the score threshold.  During the design of search_sitemaps
            a normalized score of -1.5 tended to reduce the number of false
            positives and return a high number of true positives.  Choosing a
            more stringent threshold (less than -1.5) will result in fewer
            false positives at a cost of missing more remote structural 
            homologs and structurally dissimilar proteins with binding sites
            similar to the query.  Choosing a less stringent threshold 
            (greater than -1.5) tends to increase the number of remotely similar
            binding sites as well as false positives.

            Note: if the --no_normalization flag is used, you will need to 
            adjust this treshold to a higher value (around -0.25).

      --prot_lig_score SF_NAME [default: None]
            Use an external protein-ligand scoring function to compute the
            fit of the database ligand fragments to the query binding site.
             
            Note: search_sitemaps aligns binding sites.  In other words, 
            search_sitemaps seeks to optimize the alignment of the interactions
            of two proteins and does not optimize the orientation of the ligand
            fragments with respect to the binding site.  For this reason, we
            recommend using a protein-ligand scoring function that is 
            softer than you might use for protein-ligand docking. 

      --no_normalization [default: normalize scores]
            Give the raw score for each hit rather than normalizing
            (similar to a Z-score) the scores. 

      --hydrophobic_query
            Give search_sitemaps a hint that the query pocket is highly
            hydrophobic.  Polar points tend to be specificity determinants.  If
            your query has only 1 or 2 hbonding atoms you can run a search
            using this option to focus more on the hydrophobic component and
            have less emphasis on lining up hydrogen bonding atoms.

      --SCORE_ONLY [default: align and score aligned binding sites]
            Assumes the binding sites have been aligned by another method.  
            search_sitemaps is not to align the sites, but score them as given.

Additional Options:
      These options set output locations and the temporary directory.

      -o FILE  
            Set the name of the results file.  The default value is 
            <query_sitemap>_blast.out or <query_sitemap>_to_<sitemap2>.out
            if a database or site comparision is run (respectively).  This
            file is saved in the project output directory.

      --proj_output DIRECTORY
            Set the name of the search output directory.  This value is also
            set by the $ASCBASE_PROJ_OUTPUT variable            

      --scratch_dir DIRECTORY
            Set the name of the temporary directory that ASCbase can use 
            during the course of a search.  This is set by $ASCBASE_SCRATCH_DIR
            the variable.

Site Map Database Options:
      The following options affect the database to be searched.

      --dbase_sites DIRECTORY
            Set the directory holding the database site maps.  This option
            overrides the value set by $ASCBASE_DBASE_SITES and is typically
            used with the --dbase_ligs option.

      --dbase_ligs DIRECTORY
            Set the directory holding the database ligands.  This option
            overrides the value set by $ASCBASE_DBASE_LIGS and is typically
            used with the --dbase_sites option.

      --db_index_file FILE
            Set the index file holding the protein identifiers for the
            database site maps.  Only the protein identifiers in this file 
            will be considered during the search.  This option is required when
            using --dbstart and --dbstop.

      --dbstart NUMBER
            Set the line number to start at in the database index file.

      --dbstop NUMBER
            Set the line number to finish at in the database index file.
        ''' 
  #####

  #####
  def version(self):
    print "\nsearch_sitemaps.py"
    parameters_base.version(self)
  #####
