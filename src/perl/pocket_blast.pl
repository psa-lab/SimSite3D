#!/usr/bin/perl
#
# $Source: /psa/share/repository/pfizer_proj/src/perl/pocket_blast.pl,v $
# $Revision: 1.11 $
# $Author: vanvoor4 $
# $Date: 2007-06-06 20:03:30 $
#
# $Log: not supported by cvs2svn $
# Revision 1.10  2007/02/07 17:29:07  vanvoor4
# Commented out unused variable.
#
# Revision 1.9  2007/02/07 17:24:15  vanvoor4
# Removed all references to the data/query_sitemaps dir
#
# Revision 1.8  2007/02/07 15:15:03  vanvoor4
# The assumption is that pockets are relatively small.  This means that
# the mess of storing pocket sitemaps is likely to overshadow the time gained
# by storing pocket sitemaps.
#
# Revision 1.7  2006/11/22 22:55:02  vanvoor4
# Pfizer Name fix
#
# Revision 1.6  2006/11/21 21:55:25  vanvoor4
# Only want to keep the top score below the cutoff (if any).
#
# Revision 1.5  2006/11/21 21:09:59  vanvoor4
# Sitemap header is no longer printed from this script
#
# Revision 1.4  2006/11/21 05:20:42  vanvoor4
# Changed spelling.
#
# Revision 1.2  2006/11/20 19:45:31  vanvoor4
# Added CVS header and changed perl lib dir.
#
#
#
# Given a protein pdb file and a ball specified by a point and radius 
# pocket_blast.pl will search the sitemaps in the directory 
# ASCBASE_SEARCHABLE_SITEMAPS_DIR for matches to the intersection of the
# protein and the ball.
#
# Several key assumptions are needed
# 1) There is only one protein directory
# 2) There is only one ligand directory
# 3) There is only one searchable ASCbase sitemaps directory
# 4) There is only one query ASCbase sitemaps directory -- this directory
#      could be temporary in the future

use Carp;

$#ARGV == 4 or die "USAGE: $0 <XXXXXX_p.pdb> pt_x pt_y pt_z radius\n";

# Process @ARGV
my $q_protein = $ARGV[0];
my @center;
push(@center, $ARGV[1], $ARGV[2], $ARGV[3]);
my $rad = $ARGV[4];

my $dirs = get_directories();
my $proj_bin = $dirs->{project} . "/bin";

# If we cannot find the pdb_file, look for it in the proteins_dir
if(! -f "$q_protein"){
  $q_protein = $dirs->{proteins} . "/$q_protein";
  -f $q_protein or die "Unable to find the protein (pdb) file $ARGV[0];  ";
}

# Need dir to hold the resulting files
$q_name = substr($q_protein, rindex($q_protein, "/") + 1);
$q_name = substr($q_name, 0, length($q_name) -6);
$q_name .= "-$ARGV[1]_$ARGV[2]_$ARGV[3]-$rad";
use File::Path;
mkpath($q_name,0,0744) unless -d $q_name;
mkpath("$q_name/ligands",0,0744) unless -d "$q_name/ligands";

# Generate the sitemap for the pocket
$q_file_path = "$q_name/$q_name" . "_s.pdb";
system "$proj_bin/gen_points --sphere=\"@center $rad\" --hphob_baseline" . 
  " --hbond_density sparse -p $q_protein $q_file_path";

# Do a blast against the database
print "  Blasting query pocket against the database . . . . . . . . . . . . ";
system "$proj_bin/search_sitemaps --ExtScore default --keep_top_scores 1 "
  . " --outdir $q_name $q_file_path";
print "finshed\n";
  
####################
sub get_directories{

  # Directory where ASCbase is installed
  $dirs->{project} = $ENV{'ASCBASE_SOFTWARE_DIR'} or croak "ERROR: " .
    "environment variable ASCBASE_SOFTWARE_DIR is not set\n";
  # Added so that it can be changed if needed.
  $dirs->{data} = $dirs->{project} . "/data";
  # Directory holding the Pfizer protein pdb files
  $dirs->{proteins} = $ENV{'PROTEIN_CRYSTAL_STRUCTURES_DIR'} or croak "ERROR: ".    "environment variable PROTEIN_CRYSTAL_STRUCTURES_DIR is not set\n";
  # Directory holding the Pfizer ligand mol2 files
  $dirs->{ligands} = $ENV{'LIGAND_CRYSTAL_STRUCTURES_DIR'} or croak "ERROR: ".
    "environment variable LIGAND_CRYSTAL_STRUCTURES_DIR is not set\n";
  # Directory holding the searchable ASCbase sitemap pdb files
  $dirs->{search} = $ENV{'ASCBASE_SITEMAPS_DIR'} or croak "ERROR: ".
    "environment variable ASCBASE_SITEMAPS_DIR is not set\n";

  return $dirs;
}

####################
sub calc_norm_stats{
  my($query, $data_dir, $proj_bin) = @_;

  system "$proj_bin/search_sitemaps --dbase $data_dir/diverse_sitemaps/"
    . " --ligs_dir $data_dir/diverse_ligands/ $query";

  # Assume that the query filename ends with the .pdb extension and may have
  # some specified path
  my $infile = substr($query, 0, -6);
  my $idx = rindex($query, "/") + 1;
  if($idx > 0){ $infile = substr($infile, $idx); }
  $infile .= "_blast.out";
  open IN, $infile or die "Unable to open the file $infile: $!\n";

  #$tot_prof_time = 0;
  my @scores;
  while(<IN>){
    if(substr($_,0,1) eq "%" || substr($_,0,1) eq "#"){ next; }
    my @args = split(/\|/, $_);
    push(@scores, $args[1]);
    #$tot_prof_time += $args[4];
  }
  close IN;

  $mu = 0;
  foreach $score (@scores){ $mu += $score; }
  $mu /= ($#scores + 1);

  $sigma = 0;
  foreach my $score (@scores){
    my $tmp = $score - $mu;
    $sigma +=  $tmp * $tmp;
  }

  $sigma = sqrt(1.0 / $#scores * $sigma);

  unlink $infile;
}
