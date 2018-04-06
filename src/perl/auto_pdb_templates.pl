#!/usr/bin/perl
#
# $Source: /psa/share/repository/pfizer_proj/src/perl/auto_pdb_templates.pl,v $
# $Revision: 1.3 $
# $Author: vanvoor4 $
# $Date: 2006-11-21 02:22:50 $
#
# $Log: not supported by cvs2svn $
# Revision 1.2  2006/11/20 19:45:05  vanvoor4
# Added CVS header and changed perl lib dir.
#
#
#
#

use Carp;

my $dirs = get_directories();
my $proj_bin = $dirs->{project} . "/bin";

$where_am_I = $ENV{'HOSTNAME'} or die "ERROR: environment variable ".
  "HOSTNAME is not set\n";
@stuff = split("\\.", $where_am_I);
$stuff[2] eq "msu" && $stuff[3] eq "edu" or 
  die "Development tool, please use auto_gen_sitemaps instead";

print "Running in Ligand Bounding Box Mode\n";

# Get a list of all the pdb files in $...PROTEIN...
opendir PROTEINS_DIR, $dirs->{proteins} 
  or die "Unable to open " . $dirs->{proteins} . ": $!";
@files = grep /.*.pdb/, readdir PROTEINS_DIR;
close PROTEINS_DIR;
@proteins = sort @files;

# Get a list of all the ligand files in $...LIGANDS_DIR...
opendir DIR, $dirs->{ligands} 
  or die "Unable to open " . $dirs->{ligands} . ": $!";
@files = grep /.*.mol2/, readdir DIR;
close DIR;
@ligands = sort @files;

# For each 4 char pdb code in protein file process all ligands with the same
# prefix -- assumes protein pdb codes are unique.
use File::Temp qw/ :mktemp  /;
$lig_idx = 0;
foreach $protein_file (@proteins){
  $pdb_code = substr $protein_file, 0, 4;
  $lig_prefix = substr $ligands[$lig_idx], 0, 4; 
  while($lig_idx <= $#ligands && $lig_prefix lt $pdb_code){
    $lig_idx++;
    $lig_prefix = substr $ligands[$lig_idx], 0, 4; 
  }
  
  # For each ligand with the same prefix (4 char pdb code) as the current pdb
  # file, generate a template of the protein.
  while($lig_idx <= $#ligands && $lig_prefix eq $pdb_code){
    $ligand_file = $dirs->{ligands} . "/$ligands[$lig_idx]";

    # Set output file name 
    @temp = split("\\.", $ligands[$lig_idx]);
    $out_file_name = $dirs->{search} . "/$temp[0]" . "_sitemap.pdb";

    # Generate the template
    system "gen_points " . $dirs->{proteins} . "/$protein_file $ligand_file " . 
      "$out_file_name\n";
    $lig_idx++;
    $lig_prefix = substr $ligands[$lig_idx], 0, 4; 
  } 
}

sub get_directories{

  # Directory where SimSite3D is installed
  $dirs->{project} = $ENV{'SIMSITE3D_SOFTWARE_DIR'} or croak "ERROR: " .
    "environment variable SIMSITE3D_SOFTWARE_DIR is not set\n";
  # Added so that it can be changed if needed.
  $dirs->{data} = $dirs->{project} . "/data";
  # Directory holding the query SimSite3D sitemap pdb file
  $dirs->{queries} = $dirs->{data} . "/query_sitemaps";
  # Directory holding the Pfizer protein pdb files
  $dirs->{proteins} = $ENV{'PROTEIN_CRYSTAL_STRUCTURES_DIR'} or croak "ERROR: ".    "environment variable PROTEIN_CRYSTAL_STRUCTURES_DIR is not set\n";
  # Directory holding the Pfizer ligand mol2 files
  $dirs->{ligands} = $ENV{'LIGAND_CRYSTAL_STRUCTURES_DIR'} or croak "ERROR: ".
    "environment variable LIGAND_CRYSTAL_STRUCTURES_DIR is not set\n";
  # Directory holding the searchable SimSite3D sitemap pdb files
  $dirs->{search} = $ENV{'SIMSITE3D_SEARCHABLE_SITEMAPS_DIR'} or croak "ERROR: ".
    "environment variable SIMSITE3D_SEARCHABLE_SITEMAPS_DIR is not set\n";

  return $dirs;
}
