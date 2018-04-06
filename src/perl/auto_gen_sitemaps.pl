#!/usr/bin/perl
#
# $Source: /psa/share/repository/pfizer_proj/src/perl/auto_gen_sitemaps.pl,v $
# $Revision: 1.12 $
# $Author: vanvoor4 $
# $Date: 2007-08-21 14:24:53 $
#
# $Log: not supported by cvs2svn $
# Revision 1.11  2007/06/06 20:05:21  vanvoor4
# Changed handling of C++ progs -- don't want to normalize the
# database sitemaps as that currently takes too long.
#
# Revision 1.10  2007/03/23 17:32:54  vanvoor4
# Added print statement for situation if we run out of ligands and
# still have some proteins left
#
# Revision 1.9  2007/03/23 17:18:50  vanvoor4
# Same change as before, just removed the last and moved the
# increment of idx inside the else
#
# Revision 1.8  2007/03/23 17:17:00  vanvoor4
# If a lig_fname > current prot. fname, then we do not want
# to increase the ligand index
#
# Revision 1.7  2007/03/23 16:48:38  vanvoor4
# Added cluster diameter of 3.5 (A) to generation step.
#
# Revision 1.6  2007/02/07 17:24:15  vanvoor4
# Removed all references to the data/query_sitemaps dir
#
# Revision 1.5  2007/02/07 15:16:10  vanvoor4
# Added the necessary parameters to accomodate the changes made to
# search_sitemaps.  Changed the name of the environment variable.
#
# Revision 1.4  2006/11/22 22:55:23  vanvoor4
# Pfizer Name fix
#
# Revision 1.3  2006/11/21 02:23:04  vanvoor4
# Hacked untill lib support.
#
# Revision 1.2  2006/11/20 19:42:58  vanvoor4
# Added CVS header and changed the location of the perl lib dir.
#
#
# The following file format is expected
# Protein:  XXXXXX_p.pdb
# Ligand:   XXXXXX_l.mol2
# Sitemap:  XXXXXX_s.pdb

use Carp;

$dirs = get_directories();
my $proj_bin = $dirs->{project} . "/bin";

# Get a list of all the pdb files in $...PROTEIN...
opendir PROTEINS_DIR, $dirs->{proteins} 
  or die "Unable to open " . $dirs->{proteins} . ": $!";
@files = grep /.*_p.pdb/, readdir PROTEINS_DIR;
close PROTEINS_DIR;
@proteins = sort @files;

# Get a list of all the mol2 files in $...LIGAND... sorted based on the
# protein structure names
opendir LIGANDS_DIR, $dirs->{ligands}
  or die "Unable to open " . $dirs->{ligands} . ": $!";
@files = grep /.*_l.mol2/, readdir LIGANDS_DIR;
close LIGANDS_DIR;
@ligands = sort @files;

$idx = 0;
foreach $protein (@proteins){
  # Assume we only need to strip "_p.pdb" or "_l.mol2" off the end to get the 
  # name
  my $prot_name = substr($protein, 0, length($protein) - 6);
  my $lig_name = "";

  if($idx <= $#ligands){
    $lig_name = substr($ligands[$idx], 0, length($ligands[$idx]) - 7);
  }

  # $protein does not have any associated ligands
  if($lig_name gt $prot_name){
    print "Protein file $protein\n\tdoes not have a corresponding ligand\n";
    print "\tSkipping . . .\n";
  }else{
    # Skip all ligands that do not have an associated protein
    while($idx <= $#ligands && $lig_name lt $prot_name){
      print "Ligand $ligands[$idx]\n";
      print "\tdoes not have a corresponding protein\n\tSkipping . . .\n";
      $idx++;
      if($idx <= $#ligands){
        $lig_name = substr($ligands[$idx], 0, length($ligands[$idx]) - 7);
      }
    }
    
    if($idx > $#ligands){
      carp "Reached the end of sorted ligand names, but still have some " .
        "sorted proteins names left\n";
      last;
    }

    # Assume only one ligand is associated with one protein file.
    my $csv_file = $dirs->{search} . "/" . $prot_name . "_s.csv";

    # Generate the sitemap without normalizing it with respect to the
    # 140 sitemaps since that step is currently too slow to do for the 
    # database.
    print "*****\n* Generating the sitemap of $protein\n";
    print "* with respect to the ligand $ligands[$idx]\n";
    system "$proj_bin/gen_points --no_normalization --hbond_density sparse " .
      " --hphob_baseline --cluster_diameter 3.5 -l " .  $dirs->{ligands} .
      "/$ligands[$idx] -p " . $dirs->{proteins} . "/$protein $csv_file";

    print "* Finished generating the sitemap\n*****\n";
    $idx++;
  }
}


sub get_directories{

  # Directory where SimSite3D is installed
  $dirs->{project} = $ENV{'SIMSITE3D_SOFTWARE_DIR'} or croak "ERROR: " .
    "environment variable SIMSITE3D_SOFTWARE_DIR is not set\n";
  # Added so that it can be changed if needed.
  $dirs->{data} = $dirs->{project} . "/data";
  # Directory holding the Pfizer protein pdb files
  $dirs->{proteins} = $ENV{'PROTEIN_CRYSTAL_STRUCTURES_DIR'} or croak "ERROR: ".    "environment variable PROTEIN_CRYSTAL_STRUCTURES_DIR is not set\n";
  # Directory holding the Pfizer ligand mol2 files
  $dirs->{ligands} = $ENV{'LIGAND_CRYSTAL_STRUCTURES_DIR'} or croak "ERROR: ".
    "environment variable LIGAND_CRYSTAL_STRUCTURES_DIR is not set\n";
  # Directory holding the searchable SimSite3D sitemap pdb files
  $dirs->{search} = $ENV{'SIMSITE3D_SITEMAPS_DIR'} or croak "ERROR: ".
    "environment variable SIMSITE3D_SITEMAPS_DIR is not set\n";

  return $dirs;
}
