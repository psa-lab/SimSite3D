#!/usr/bin/perl -w
#
# $Source: /psa/share/repository/pfizer_proj/src/perl/search_sitemaps.pl,v $
# $Revision: 1.5 $
# $Author: vanvoor4 $
# $Date: 2007-06-06 20:03:02 $
#
# $Log: not supported by cvs2svn $
# Revision 1.4  2007/03/23 16:48:03  vanvoor4
# Changed printing so that it reflects the fact that one search_sitemaps
# operates in 1 of 2 modes
#
# Revision 1.3  2007/02/07 17:29:07  vanvoor4
# Commented out unused variable.
#
# Revision 1.2  2007/02/07 17:24:15  vanvoor4
# Removed all references to the data/query_sitemaps dir
#
# Revision 1.1  2007/02/07 15:11:11  vanvoor4
# Replaces run_site_dock.pl
#
#

use Carp;

$#ARGV == 0 or $#ARGV == 1 or die "USAGE: $0 <query_s.pdb> [search_s.pdb]\n";
my $template_one = $ARGV[0];

my $dirs = get_directories();
my $proj_bin = $dirs->{project} . "/bin";

# If we cannot find the first template file look in the dbase dir. 
if(! -f $template_one){
  $template_one = $dirs->{search} . "/$template_one";
  -f $template_one or die "Unable to find the sitemap file $ARGV[0]";
}
my $q_sitemap = $template_one;

# Look for the second template (if there is one)
my $template_two = "";
if($#ARGV == 1){
  $template_two = $ARGV[1];
  if(! -f "$template_two"){
    $template_two = $dirs->{search} . "/$template_two";
    -f $template_two or die "Unable to find the sitemap file $ARGV[1]";
  }
}

# Need dir to hold the resulting files
$q_name = substr($template_one, rindex($template_one, "/") + 1);
$q_name = substr($q_name, 0, length($q_name) -6);
use File::Path;
mkpath($q_name,0,0744) unless -d $q_name;
mkpath("$q_name/ligands",0,0744) unless -d "$q_name/ligands";

# Do a blast against the database
if($template_two eq ""){
  print "  Blasting query pocket against the database . . . . . . . . . . . . ";
}else{
  print "  Comparing $template_one with $template_two . . . . . . . . . . . . ";
}
system "$proj_bin/search_sitemaps --outdir $q_name $template_one $template_two";
print "finshed\n";

####################
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
