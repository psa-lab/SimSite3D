import os, sys
import SimSite3DPy

################################################################################

################################################################################
class FragGen:
  ##############################################################################

  ##############################################################################
  def __init__(self, cif_ligs_dir, FragGen_dir = "/soft/linux64/fraggen", 
               java_dir="/soft/linux64/jre1.5.0_06/bin/java",
               max_cuts=3, frag_props="HeavyAtomCount 5 15"):
    """
    self: this class
    cif_ligs_dir: Directory holding the cif ligands to screen/search
    FragGen_dir: Directory holding the FragGen .jar files (program)
    java_dir: The dir to the current JRE (needs to be compatible with FragGen)
    max_cuts: Maximum number of bonds that FragGen can cut when generating
              ligand fragments
    frag_props: Properties used by FragGen to determine valid fragments
    """
    self.java_dir = java_dir
    self.cif_ligs_dir = cif_ligs_dir
    self.FragGen_dir = FragGen_dir
    self.max_cuts = max_cuts
    self.frag_props = frag_props

    # GCR is missing at least one model and one ideal coordinate
    # OXY is molecular oxygen and is missing a _loop section for 
    # comp_chem_bond (technically ok, but the current SimSite3D code
    # was not designed for it ).
    self.BAD_HETS = ["GCR", "OXY"]

    # build the "dataset" file
    self.dataset = \
      self.__frag_ligs__(self.cif_ligs_dir, self.max_cuts, self.frag_props)
  ##############################################################################

  ##############################################################################
  def __frag_ligs__(self, cif_ligs_dir, max_cuts, frag_props):
    lig_ids = []
    for fname in os.listdir(cif_ligs_dir):
      if(not fname.endswith(".cif")): continue
      lig_ids.append(fname[0:-4])
    
    dataset = "FragGen_%d_cuts_%s" % (max_cuts, frag_props.replace(" ", "-"))
    sdf_file = open("%s.sdf" % (dataset), "w+")
    for lig_id in lig_ids:
      if(lig_id in self.BAD_HETS): continue

      # Load the ideal cif and use it to convert the pdb hetgrp to sdf format
      my_cif = \
        SimSite3DPy.utils.pdb_cif.molecule("%s/%s.cif" % (cif_ligs_dir, lig_id))
      if(my_cif.fail): 
        print "Could not load the molecule: %s/%s.cif" % (cif_ligs_dir, lig_id)
        continue
      my_chem_comp = my_cif.chemical_component
      my_sdf = my_chem_comp.to_sdf()
      
      print >> sdf_file, my_sdf
    sdf_file.close()

    #
    # Prepare database file
    #
    Prep_cmd = "%s -jar %s/Prep003.jar %s.sdf %s_clean.sdf " % \
	    (self.java_dir, self.FragGen_dir, dataset, dataset)
    Prep_cmd += "-keepLargest -keep PUBCHEM_COMPOUND_CID -unique "
    Prep_cmd += "-properties IsOrganic true"
    print "\nPreparing the ligands\n"
    status = os.system(Prep_cmd)
    print "\nFinished preparing the ligands"
    
    #
    # Generate all fragments, annotate fragments with compoundID, and annotate 
    # database file with fragment keys
    #
    Frag_cmd = "%s -jar %s/FragGen003.jar " % (self.java_dir, self.FragGen_dir)
    Frag_cmd += "%s_clean.sdf %s_all_fragments.sdf %s_keys.sdf " % \
      (dataset, dataset, dataset)
    Frag_cmd += " -cuts %d " % (max_cuts)
    Frag_cmd += " -allFragments "
    Frag_cmd += " -properties %s" % (frag_props)
    print "\nFragmenting the ligands\n"
    os.system(Frag_cmd)
    print "\nFinished fragmenting the ligands\n"

    return "%s_keys.sdf" % (dataset)
  ##############################################################################

  ##############################################################################
  def exact_frag_matches(self, frag_fname, dataset, res_fname):
    """ 
    Given a ligand or fragment in cif or sdf format, use FragGen to find
    all ligand that contain the fragment 
  
    frag_fname: Name of the ligand/fragment file -- i.e. query ligand/fragment
    dataset: Name of the fragmented ligand dataset to screen/search
    res_fname: Name of the results file (exact fragment hits) -- written in .sdf
               format
    """ 

    sdf_frag_fname = frag_fname
    if(frag_fname.endswith(".cif")):
      #
      # Get the query sdf file using the CIF file
      #
      pos = frag_fname.rfind("/")
      if(pos > -1): sdf_frag_fname = frag_fname[pos+1:-3] + "sdf"
      else: sdf_frag_fname = frag_fname[:-3] + "sdf"
     
      q_sdf_file = open(sdf_frag_fname, "w+")
      my_cif = SimSite3DPy.utils.pdb_cif.molecule(frag_fname)
      if(my_cif.fail): 
        print "Warning!  Unable to convert", frag_fname, "to .sdf format"
        return False
      my_chem_comp = my_cif.chemical_component
      my_sdf = my_chem_comp.to_sdf()
      print >> q_sdf_file, my_sdf
      q_sdf_file.close()

    elif(not frag_fname.endswith(".sdf")):
      print "Warning!  Input files must be in PDF .cif format or .sdf format"
      print "\tReceived", frag_fname
      return False

    Screen_cmd = "%s -jar %s/SimScreen003.jar " % \
      (self.java_dir, self.FragGen_dir)
    Screen_cmd += "%s %s %s " % (sdf_frag_fname, dataset, res_fname)
    Screen_cmd += "-fingerPrint FragmentKeys "
    Screen_cmd += "-metric SuperStructure "
    print "\nScreening the ligands for exact fragment matches\n"
    os.system(Screen_cmd)
    print "\nFinished screening the ligands\n"
    return True
  ##############################################################################

  ##############################################################################
  def Search(self, query_lig, dataset, res_fname, metric="SuperStructure"):
    """
    query_lig: path to the query ligand .cif or .sdf file
    dataset: screening dataset
    res_fname: name of the results (hits) file
    metric: A valid FragGen metric for fingerprint searching 
    """
    sdf_lig_fname = self.convert_lig_if_needed(query_lig)
    if(not len(sdf_lig_fname)): return []
 
    Screen_cmd = "%s -jar %s/SimScreen003.jar " % \
      (self.java_dir, self.FragGen_dir)
    Screen_cmd += "%s %s tmp_unsorted.sdf" % (sdf_lig_fname, dataset)
    Screen_cmd += " -minSimilarity 0.10 "
    Screen_cmd += "-metric %s " % (metric)
    print "\nScreening the ligands -- %s\n" % (metric)
    os.system(Screen_cmd)
    
    print "\nSorting the ligands by %s score\n" % (metric)
    hits = SimSite3DPy.utils.sdf.sdf_molecules("tmp_unsorted.sdf")
    hits.sort_by_score(score_type=metric)
    out_file = open(res_fname, "w+")
    for my_mol in hits.mols: print >> out_file, my_mol
    out_file.close()
    os.unlink("tmp_unsorted.sdf")
    return hits
  ##############################################################################

  ##############################################################################
  def convert_lig_if_needed(self, lig_fname):
    sdf_lig_fname = lig_fname
    if(lig_fname.endswith(".cif")):
      #
      # Get the query sdf file using the CIF file
      #
      pos = lig_fname.rfind("/")
      if(pos > -1): sdf_lig_fname = lig_fname[pos:-3] + "sdf"
      else: sdf_lig_fname = lig_fname[:-3] + "sdf"
     
      q_sdf_file = open(sdf_lig_fname, "w+")
      my_cif = SimSite3DPy.utils.pdb_cif.molecule(lig_fname)
      if(my_cif.fail): 
        print "Warning!  Unable to convert", lig_fname, "to .sdf format"
        return ""
      my_chem_comp = my_cif.chemical_component
      my_sdf = my_chem_comp.to_sdf()
      print >> q_sdf_file, my_sdf
      q_sdf_file.close()

    elif(not lig_fname.endswith(".sdf")):
      print "Warning!  Input files must be in PDB .cif format or .sdf format"
      print "\tReceived", lig_fname
      return ""

    return sdf_lig_fname
  ##############################################################################

  ##############################################################################
    
  

#my_java = "/soft/linux64/jre1.5.0_06/bin/java";
#FragGen_dir = "/soft/linux64/fraggen";
#
#max_cuts = 3
#frag_props = "HeavyAtomCount 5 15"
#def run(query_lig_id, cif_ligs_dir, cif_lig_fname, cif_frag_fname=""):
#  """
#  query_lig_id: id to use when naming the output files
#  cif_ligs_dir: directory of cif ligand to screen / search
#  cif_lig_fname: name (path to) the cif ligand or fragment
#  cif_frag_fname: if provided use this file for the fragment screen
#    (exact fragment matching), otherwise use the cif_lig_fname as the
#    fragment for fragment screening
#  """
#
#  lig_ids = []
#  for fname in os.listdir(cif_ligs_dir):
#    if(not fname.endswith(".cif")): continue
#    lig_ids.append(fname[0:3])
#  
#  dataset = "blah"
#  sdf_file = open("%s.sdf" % (dataset), "w+")
#  for lig_id in lig_ids:
#    # Load the ideal cif and use it to convert the pdb hetgrp to sdf format
#    my_cif = \
#      SimSite3DPy.utils.pdb_cif.molecule("%s/%s.cif" % (cif_ligs_dir, lig_id))
#    if(my_cif.fail): 
#      print "Could not load the molecule: %s/%s.cif" % (cif_ligs_dir, lig_id)
#      continue
#    my_chem_comp = my_cif.chemical_component
#    my_sdf = my_chem_comp.to_sdf()
#    
#    print >> sdf_file, my_sdf
#  sdf_file.close()
#  
#  #
#  # Get the query sdf file using the CIF file
#  #
#  frag_fname = "%s.sdf" % (query_lig_id)
#  q_sdf_file = open(frag_fname, "w+")
#  my_cif = SimSite3DPy.utils.pdb_cif.molecule(cif_lig_fname)
#  if(my_cif.fail): sys.exit(-1)
#  my_chem_comp = my_cif.chemical_component
#  my_sdf = my_chem_comp.to_sdf()
#  print >> q_sdf_file, my_sdf
#  q_sdf_file.close()
#
#  #
#  # Get the fragment sdf file using the CIF file
#  #
#  if(len(cif_frag_fname) > 0):
#    pos = cif_frag_fname.rfind("/")
#    if(pos > -1):
#      frag_fname = cif_frag_fname[cif_frag_fname.rfind("/"):-3] + "sdf"
#    else: frag_fname = cif_frag_fname[:-3] + "sdf"
#
#    q_sdf_file = open(frag_fname, "w+")
#    my_cif = SimSite3DPy.utils.pdb_cif.molecule(cif_frag_fname)
#    if(my_cif.fail): sys.exit(-1)
#    my_chem_comp = my_cif.chemical_component
#    my_sdf = my_chem_comp.to_sdf()
#    print >> q_sdf_file, my_sdf
#    q_sdf_file.close()
#  
#  #
#  # Prepare database file
#  #
#  Prep_cmd = "%s -jar %s/Prep003.jar %s.sdf %s_clean.sdf " % \
#	  (my_java, FragGen_dir, dataset, dataset)
#  Prep_cmd += "-keepLargest -keep PUBCHEM_COMPOUND_CID -unique "
#  Prep_cmd += "-properties IsOrganic true"
#  print "\nPreparing the ligands\n"
#  status = os.system(Prep_cmd)
#  print "\nFinished preparing the ligands"
#  
#  #
#  # Generate all fragments, annotate fragments with compoundID, and annotate 
#  # database file with fragment keys
#  #
#  Frag_cmd = "%s -jar %s/FragGen003.jar " % (my_java, FragGen_dir)
#  Frag_cmd += "%s_clean.sdf %s_all_fragments.sdf %s_keys.sdf " % \
#    (dataset, dataset, dataset)
#  Frag_cmd += " -cuts %d " % (max_cuts)
#  Frag_cmd += " -allFragments "
#  Frag_cmd += " -properties %s" % (frag_props)
#  print "\nFragmenting the ligands\n"
#  os.system(Frag_cmd)
#  print "\nFinished fragmenting the ligands\n"
#
#  
#  # Use below for fragment based screening -- exact hash key matching
#  #
#  # Use fragment keys file to run a similarity search for Adenine
#  #
#  Screen_cmd = "%s -jar %s/SimScreen003.jar " % (my_java, FragGen_dir)
#  Screen_cmd += "%s %s_keys.sdf %s_Frag_matches.sdf " % \
#    (frag_fname, dataset, query_lig_id)
#  Screen_cmd += "-fingerPrint FragmentKeys "
#  Screen_cmd += "-metric SuperStructure "
#  print "\nScreening the ligands for exact fragment matches\n"
#  os.system(Screen_cmd)
#  print "\nFinished screening the ligands\n"
#  
#  Screen_cmd = "%s -jar %s/SimScreen003.jar " % (my_java, FragGen_dir)
#  Screen_cmd += "%s.sdf %s_keys.sdf %s_SubStruct_unsorted.sdf " % \
#    (query_lig_id, dataset, query_lig_id)
#  Screen_cmd += " -minSimilarity 0.10 "
#  Screen_cmd += "-metric SuperStructure "
#  print "\nScreening the ligands -- Substructure\n"
#  os.system(Screen_cmd)
#  
#  print "\nSorting the ligands by Substructure score\n"
#  SubStruct_hits = \
#    SimSite3DPy.utils.sdf.sdf_molecules(query_lig_id + "_SubStruct_unsorted.sdf")
#  SubStruct_hits.sort_by_score(score_type="SuperStructure")
#  out_file = open(query_lig_id + "_SubStruct_sorted.sdf", "w+")
#  for my_mol in SubStruct_hits.mols: print >> out_file, my_mol
#  os.unlink(query_lig_id + "_SubStruct_unsorted.sdf")
#  
#  Screen_cmd = "%s -jar %s/SimScreen003.jar " % (my_java, FragGen_dir)
#  Screen_cmd += "%s.sdf %s_keys.sdf %s_Tanimoto_unsorted.sdf " % \
#    (query_lig_id, dataset, query_lig_id)
#  Screen_cmd += " -minSimilarity 0.10 "
#  Screen_cmd += "-metric Tanimoto "
#  print "\nScreening the ligands -- Tanimoto\n"
#  os.system(Screen_cmd)
#  
#  print "\nSorting the ligands by Tanimoto score\n"
#  Tanimoto_hits = \
#    SimSite3DPy.utils.sdf.sdf_molecules(query_lig_id + "_Tanimoto_unsorted.sdf")
#  Tanimoto_hits.sort_by_score(score_type="Tanimoto")
#  out_file = open(query_lig_id + "_Tanimoto_sorted.sdf", "w+")
#  for my_mol in Tanimoto_hits.mols: print >> out_file, my_mol
#  os.unlink(query_lig_id + "_Tanimoto_unsorted.sdf")
#  
#  return (SimSite3DPy.utils.sdf.sdf_molecules(query_lig_id + "_Frag_matches.sdf"),
#          SubStruct_hits, Tanimoto_hits) 
###  print "\nFinished screening the ligands\n"
