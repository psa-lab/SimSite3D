#!/usr/bin/perl
# $Source: /psa/share/repository/pfizer_proj/src/perl/gen_sitemap.pl,v $
# $Revision: 1.8 $
# $Author: vanvoor4 $
# $Date: 2007-06-06 20:04:43 $
#
# $Log: not supported by cvs2svn $
# Revision 1.7  2007/03/23 16:48:46  vanvoor4
# Added cluster diameter of 3.5 (A) to generation step.
#
# Revision 1.6  2007/02/07 17:24:15  vanvoor4
# Removed all references to the data/query_sitemaps dir
#
# Revision 1.5  2007/02/07 15:13:33  vanvoor4
# Environment variable name change and parameters to gen_points were
# added.
#
# Revision 1.4  2006/11/22 22:54:57  vanvoor4
# Pfizer Name fix
#
# Revision 1.3  2006/11/21 02:24:05  vanvoor4
# hacked untill lib support.  In addition, untill we decide on how to
# correctly determine which sitemap is being used, we will need to
# normalize each run.
#
# Revision 1.2  2006/11/20 19:45:17  vanvoor4
# Added CVS header and changed perl lib dir.
#
#
#

use Carp;

$#ARGV == 1 or $#ARGV == 2 or die "USAGE: " .
  "$0 <XXXXXX_p.pdb> <XXXXXX_l.mol2> [XXXXXX_s.pdb]\n";

$pdb_file = $ARGV[0];
$ligand_file = $ARGV[1];

my $dirs = get_directories();
my $proj_bin = $dirs->{project} . "/bin";
                                                                                
# If we cannot find the ligand_file, look for it in the ligands_dir
if(! -f "$ligand_file"){
  $ligand_file = $dirs->{ligands} . "/$ligand_file";
  -f $ligand_file or die "Unable to find the ligand file $ARGV[1];  ";
}
                                                                              
# The resulting sitemap file will be in the sitemap format and have the
# same name as the mol2 file (without the .mol2 suffix of course).
$out_file_name = "";
if($#ARGV == 2){
  $out_file_name = $ARGV[2];
}else{
  @temp = split(/\//, $ARGV[$#ARGV]);
  # strip off the "_l.mol2"
  $out_file_name = substr($temp[$#temp], 0, length($temp[$#temp]) - 7);
  $out_file_name .= "_s.pdb";
}

# If we cannot find the pdb_file, look for it in the proteins_dir
if(! -f "$pdb_file"){
  $pdb_file = $dirs->{proteins} . "/$pdb_file";
  -f $pdb_file or die "Unable to find the protein (pdb) file $ARGV[0];  ";
}
                                                                              
system "$proj_bin/gen_points --hbond_density sparse " .
      " --hphob_baseline --cluster_diameter 3.5 -l  $ligand_file " .
      "-p $pdb_file $out_file_name";


#########################
sub get_directories{

  # Directory where SimSite3D is installed
  $dirs->{project} = $ENV{'ASCBASE_SOFTWARE_DIR'} or croak "ERROR: " .
    "environment variable ASCBASE_SOFTWARE_DIR is not set\n";
  # Added so that it can be changed if needed.
  $dirs->{data} = $dirs->{project} . "/data";
  # Directory holding the Pfizer protein pdb files
  $dirs->{proteins} = $ENV{'PROTEIN_CRYSTAL_STRUCTURES_DIR'} or croak "ERROR: ".    "environment variable PROTEIN_CRYSTAL_STRUCTURES_DIR is not set\n";
  # Directory holding the Pfizer ligand mol2 files
  $dirs->{ligands} = $ENV{'LIGAND_CRYSTAL_STRUCTURES_DIR'} or croak "ERROR: ".
    "environment variable LIGAND_CRYSTAL_STRUCTURES_DIR is not set\n";
  # Directory holding the searchable SimSite3D sitemap pdb files
  $dirs->{search} = $ENV{'ASCBASE_SITEMAPS_DIR'} or croak "ERROR: ".
    "environment variable ASCBASE_SITEMAPS_DIR is not set\n";

  return $dirs;
}
