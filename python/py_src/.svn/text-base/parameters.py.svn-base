import os
from sys import stderr

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
    self.require_min_npts = True
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

    if(self.proj_output == None): self.proj_output = os.getcwd()
    if(self.scratch_dir == None): self.scratch_dir = self.proj_output
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
  def setup_directories(self, mode=0770):
    if(not os.path.exists(self.proj_output)): os.mkdirs(self.proj_output, mode)
    if(not os.path.exists(self.proj_output + "/moved_ligands")):
      os.makedirs(self.proj_output + "/moved_ligands", mode)
    if(not os.path.exists(self.proj_output + "/ligand_fragments")):
      os.makedirs(self.proj_output + "/ligand_fragments", mode)
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
    self.time_process = True
    self.hydrophobic_query = False
    self.query_site_id = ""
    self.query_sitemap_dir = "."
    self.dbase_site_id = ""
    self.dbase_sitemap_dir = "."
    self.db_index_fname = ""
    self.dbstart = 0 
    self.dbstop = 0 
    self.score_method = "WeightedSumsScore"
    self.cmdline_proj_output = ""
    self.internal_prot_lig_scoring = False
  ##########

  ##########
  def cmdline_options(self, argv):
    import getopt
    long_opts = ["version", "usage", "help", "no_timers", "prot_lig_score=", 
                 "proj_output=", "dbase_prots=", "dbase_sites=", "dbase_ligs=",
                 "scratch_dir=", "lig_frag_size=", "keep_n_scores=", 
                 "score_threshold=", "db_index_file=", "dbstart=", "dbstop=",
                 "no_normalization", "SCORE_ONLY", "hydrophobic_query", 
                 "DO_NOT_WRITE_LIGANDS", "num_rand_aligns=",
                 "internal_prot_lig_scoring", "allow_small_site_maps",
                 "no_fine_tuning", "max_corr_surf_pt_dist="]
    try:
      (opts, args) = getopt.getopt(argv[1:], "Vuh:o", long_opts)
    except getopt.GetoptError, err:
      self.version(argv[0])
      print str(err)
      self.usage(argv[0])
      self.fail = True
      return False

    for opt, arg in opts:
      if(opt in ("-V", "--version") ):
        self.version(argv[0])
        return False
      elif(opt in ("-u", "--usage") ):
        self.version(argv[0])
        self.usage(argv[0])
        return False
      elif(opt in ("-h", "--help") ):
        self.version(argv[0])
        self.help(argv[0])
        return False
      elif(opt == "-o"):
        self.ofname = arg
      elif(opt == "--prot_lig_score"):
        self.prot_lig_score = arg
      elif(opt == "--proj_output"):
        self.cmdline_proj_output = arg
      elif(opt == "--dbase_prots"):
        self.dbase_prots = arg
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
        print "got a db index file", arg
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
      elif(opt == "--no_timers"):
        self.time_process = False
      elif(opt == "--DO_NOT_WRITE_LIGANDS"):
        self.write_ligands = False
      elif(opt == "--internal_prot_lig_scoring"):
        self.internal_prot_lig_scoring = True
      elif(opt == "--allow_small_site_maps"):
        self.require_min_npts = False
      elif(opt == "--no_fine_tuning"):
        self.fine_tune_tier2_alignments = False
        self.fine_tune_best_tier2_alignment = False
      elif(opt == "--num_rand_aligns"):
        try:
          self.num_rand_aligns = int(arg)
        except ValueError:
          print "--num_rand_aligns requires a positive integer argument"
          self.fail = True
          return False
      elif(opt == "--max_corr_surf_pt_dist"):
        try:
          self.max_corr_surf_pt_dist = float(arg)
        except ValueError:
          print "--max_corr_surf_pt_dist requires real number argument"
          self.fail = True
          return False

    if(len(args) == 0 or len(args) > 2):
      self.version(argv[0])
      print "* %s requires 1 or 2 site maps *" % (argv[0])
      self.usage(argv[0])
      self.fail = True
      return False
    else:
      self.query_site_id = args[0]
      (self.query_sitemap_dir, self.query_site_id) = \
        os.path.split(self.query_site_id)
      if(not len(self.query_site_id)):
        print "* %s got a directory as the query sitemap *" % (argv[0])
        return False
      if(not len(self.query_sitemap_dir)): self.query_sitemap_dir = "."
      pos = self.query_site_id.find(".")
      if(pos > 0): self.query_site_id = self.query_site_id[:pos-2]

      if(len(args) == 2):
        self.dbase_site_id = args[1] 
        (self.dbase_sitemap_dir, self.dbase_site_id) = \
          os.path.split(self.dbase_site_id)
        if(not len(self.dbase_site_id)):
          print "* Received 2 sitemaps but database sitemap is a directory *"
          return False
        if(not len(self.dbase_sitemap_dir)): self.dbase_sitemap_dir = "."
        pos = self.dbase_site_id.find(".")
        if(pos > 0): self.dbase_site_id = self.dbase_site_id[:pos-2]

    return self.verify_options() 
  ##########

  ##########
  def verify_options(self):
    if(self.fail):
      return False

    if(not self.test_directory(self.dbase_prots, "database protein files")):
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

    if(len(self.cmdline_proj_output)):
      self.proj_output = self.cmdline_proj_output
    else:
      self.proj_output += "/" + self.query_site_id + "_results"

    if(not os.path.isdir(self.proj_output)):
      try:
        os.makedirs(self.proj_output, 0770)
      except IOError, (errno, strerror):
        print >> stderr, "Unable to create the output directory:",
        print >> stderr, self.proj_output
        print >> stderr, "error(%s): %s" % (errno, strerror)
        return False

    if(not self.test_directory(self.scratch_dir, "scratch directory")):
      return False

    if(self.num_scores_to_keep < 1):
      print "Number of hits to keep per database site must be positive"
      return False

    if(self.dbstart > self.dbstop):
      print "dbstart may not be greater than dbstop"
      return False

    if(self.dbstart > 0 and self.dbstop > 0 and \
       len(self.db_index_fname) and not os.path.isfile(self.db_index_fname)):
      print "%s is not a valid file" % (self.db_index_fname)
      print "Please set the db index filename to a valid index file"
      return False

    # shouldn't get here
    if(self.dbstart < 0 or self.dbstop < 0):
      print "dbstart and dbstop must be positive numbers"
      return False

    # Need to support these fast
    #self.prot_lig_score = "None"

    if(len(self.ofname) == 0):
      self.ofname = self.proj_output + "/" + self.query_site_id
      if(len(self.dbase_site_id)): 
        self.ofname += "_to_" + self.dbase_site_id + ".out"
      else: self.ofname += "_blast.out"

    # support these for development & testing at MSU
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
  def usage(self, prog_name):
     print '''
%s [-Vuho] [--version --usage --help
      --prot_lig_score SF_NAME --dbase_prots DIR --dbase_sites DIR 
      --dbase_ligs DIR --scratch_dir DIR --lig_frag_size NUM --keep_n_scores NUM
      --score_threshold REAL --db_index_file FILE --dbstart NUM --dbstop NUM 
      --no_normalization --SCORE_ONLY --hydrophobic_query --no_timers
      --max_corr_surf_pt_dist REAL ]

Database search: 
      %s [OPTIONS] <query_id>
Site comparison: 
      %s [OPTIONS] <query_id> <site_id>
        ''' % (prog_name, prog_name, prog_name)
  ##########

  ##########
  def help(self, prog_name):
     print '''
Name:
      %s - ASCbase Software tool to compare binding sites

Database search: 
      %s [OPTIONS] <query_id>

Site comparison: 
      %s [OPTIONS] <query_id> <site_id>

Description:
      %s is a tool in the ASCbase Software package that takes
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
      -V,--version Print version information and exit
      -h,--help    Print this help message
      -u,--usage   Print a short help message
      --no_timers  Do not use Linux/Posix system timers to track the execution 
                   time of the search

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
            <query_site_id>_blast.out or <query_site_id>_to_<site2_id>.out
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
        ''' % (prog_name, prog_name, prog_name, prog_name)
  #####

  #####
  def version(self, prog_name):
    print "\n%s" % (prog_name)
    parameters_base.version(self)
  #####

  #####
  def to_C(self):
    from ASCbasePy import search
    params_cc = search.parameters()

    # base parameters
    params_cc.dbase_sites = self.dbase_sites
    params_cc.dbase_ligs = self.dbase_ligs
    params_cc.dbase_prots = self.dbase_prots 
    params_cc.diverse_sites = self.diverse_sites
    params_cc.diverse_ligs = self.diverse_ligs
    params_cc.proj_output = self.proj_output
    params_cc.scratch_dir = self.scratch_dir
  
    # search parameters
    params_cc.model_file_name = \
      self.query_sitemap_dir + "/" + self.query_site_id
    params_cc.dbase_file_name = \
      self.dbase_sitemap_dir + "/" + self.dbase_site_id
  
    params_cc.normalize = self.normalize
    params_cc.num_scores_to_keep = self.num_scores_to_keep
    params_cc.score_cutoff = self.score_cutoff
    params_cc.min_num_atoms = self.min_num_atoms
    params_cc.dmetol = self.dmetol
    params_cc.lsetol = self.lsetol
    params_cc.write_ligands = self.write_ligands
    params_cc.allow_hphob_triangles = self.hydrophobic_query
    params_cc.do_internal_prot_lig_score = self.internal_prot_lig_scoring
    params_cc.load_surf_files = True
    params_cc.require_min_npts = self.require_min_npts
    params_cc.fine_tune_tier2_alignments = self.fine_tune_tier2_alignments
    params_cc.fine_tune_best_tier2_alignment = \
      self.fine_tune_best_tier2_alignment

    return params_cc
  #####
