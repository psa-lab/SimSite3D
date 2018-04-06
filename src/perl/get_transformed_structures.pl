#!/usr/bin/perl -w
#
# $Source: /psa/share/repository/pfizer_proj/src/perl/get_transformed_structures.pl,v $
# $Revision: 1.5 $
# $Author: vanvoor4 $
# $Date: 2007-06-07 13:33:56 $
#
# $Log: not supported by cvs2svn $
# Revision 1.4  2007/06/06 20:04:06  vanvoor4
# Changed remark 400 to remark 40
#
# Revision 1.3  2007/03/13 15:08:32  vanvoor4
# Changed to move only the protein structures for the dbase hits.  Added
# the name of the query to the REMARKS section
#
# Revision 1.2  2007/02/07 17:24:15  vanvoor4
# Removed all references to the data/query_sitemaps dir
#
# Revision 1.1  2007/02/07 15:10:48  vanvoor4
# Initial checkin of a potential useless file.  At checkin time, the ligand
# fragments are already moved to the query coordinates by the C++ part of
# SimSite3D
#
#
#

use Carp;

$#ARGV == 1 or die "USAGE: $0 <query_blast>.out <out_dir>\n";
my $results = $ARGV[0];
my $out_dir = $ARGV[1];
use File::Path;
mkpath($out_dir,0,0755) unless -d $out_dir;

my $dirs = get_directories();
my $proj_bin = $dirs->{project} . "/bin";

open RESULTS, $results 
  or die "Unable to open the results file $results for reading: $!\n";

my $query_name = "";
while(<RESULTS>){
  # Skip comments but look for query name
  if(substr($_,0,1) eq "#" || substr($_, 0, 1) eq "%"){
    if(substr($_, 2, 27) eq "Query (Model) Sitemap Name:"){
      my @tmp = split(/ +/, substr($_, 30)); 
      chomp(@tmp);
      $query_name = $tmp[1];
    }
    next;
  }

  my @fields = split(/\|/, $_); 
  # assume only the sitemap file name is stored not the entire or partial path
  my $sitemap = $fields[0];
  my $protein = substr($sitemap, 0, length($sitemap) - 6) . "_p.pdb";
  my @rot_mat = split(" +", $fields[3]);
  my @trans_vec = split(" +", $fields[4]);
  my @trans_mat;
  for($r = 0; $r < 3; $r++){
    for($c = 0; $c < 3; $c++){
      push(@trans_mat, $rot_mat[3*$r + $c]);
    }
    push(@trans_mat, $trans_vec[$r]);
  }
 
  $out_file = $out_dir . "/" . substr($protein, rindex($protein, "/") + 1);
  translate_pdb($query_name, $protein, $out_file, \@trans_mat);
}
close RESULTS;



############################################
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


############################################
sub translate_pdb{

  my($query_name, $pdb_in, $pdb_out, $trans_ref) = @_;

  open IN, $pdb_in 
    or croak "Unable to open the pdb file $pdb_in for reading: $!\n"; 
  open OUT, "+>$pdb_out"
    or croak "Unable to open the pdb file $pdb_out for writing: $!\n"; 

  ############# PRINT HEADER SO WE KNOW WHERE THIS FILE ORIGINATED
  print OUT "REMARK 10\n";
  print OUT "REMARK 10 COMPOUND\n";
  print OUT "REMARK 10 Moved by SimSite3D from original crystallographic coordinates to those\n";
  print OUT "REMARK 10 of $query_name\n";

  while(<IN>){
    # We only want to transform HETATM and ATOM lines
    my $tag = substr($_, 0, 6);
    if($tag ne "HETATM" && $tag ne "ATOM  "){
      print OUT "$_";
    }else{
      # Need to get the cooridnates from the lines
      my $beg = substr($_, 0, 30);
      print OUT $beg;
      my $end = substr($_, 54); 
      my @pt;
      for($i = 0; $i < 3; $i++){
        push(@pt, substr($_, 31 + 8*$i, 8));
      }
      @pt = transform_point($trans_ref, \@pt);
      printf(OUT "%8.3f%8.3f%8.3f", $pt[0], $pt[1], $pt[2]);
      print OUT $end;
    }
  }
  close IN;
  close OUT;

}


############################################
sub translate_mol2{

  my($mol2_in, $mol2_out, $trans_ref) = @_;

  open IN, $mol2_in 
    or croak "Unable to open the pdb file $mol2_in for reading: $!\n"; 
  open OUT, "+>$mol2_out"
    or croak "Unable to open the pdb file $mol2_out for writing: $!\n"; 

  ############# PRINT HEADER SO WE KNOW WHERE THIS FILE ORIGINATED
  print OUT "# Modified by SimSite3D: \n";
  print OUT "# Moved from crystallographic coordinates to those\n";
  print OUT "# of the query from SimSite3D\n";

  # spin till find MOLECULE section
  while(<IN>){
    print OUT $_;
    if(/@<TRIPOS>MOLECULE/){ last; }
  }

  # Get the number of atoms
  $_ = <IN>;
  print OUT $_;
  $_ = <IN>;
  print OUT $_;
## need to handle whitespace issues correctly here
  my @fields = split(" +", $_);
  shift @fields;
  my $num_atoms = $fields[0];
  print "Number of atoms is $num_atoms\n";

  # spin till find ATOM section
  while(<IN>){
    print OUT $_;
    if(/@<TRIPOS>ATOM/){ last; }
  }

  for($n = 0; $n < $num_atoms; $n++){
    $_ = <IN>;
  ## need to handle whitespace issues correctly here
    my @fields = split(/ +/, $_);
    shift @fields;
    my @pt;
    push(@pt, $fields[2], $fields[3], $fields[4]);
    @pt = transform_point($trans_ref, \@pt);
    print OUT "$fields[0] $fields[1] ";
    printf(OUT "%10.4f %10.4f %10.4f ", $pt[0], $pt[1], $pt[2]);
    for($i = 5; $i <= $#fields; $i++){
      print OUT " $fields[$i]";
    }
    #print OUT "\n";
  }

  #  just output rest as is
  while(<IN>){
    print OUT $_;
  }

  close IN;
  close OUT;
}


############################################
# Uses the typical 4x4 transformation matrix but without the bottom row.
# Also, it is typical that the point passed in does not have a value for the
# fourth coordinate (should be 1 in an 3D transformation).
sub transform_point{
  my ($trans_ref, $pt_ref) = @_;
  
  if(scalar(@$pt_ref) == 3){
    push(@$pt_ref, 1);
  }

  my @coords;
  for($r = 0; $r < 3; $r++){
    push(@coords, 0);   
    for($c = 0; $c < 4; $c++){
      $coords[$r] += $trans_ref->[4*$r + $c] * $pt_ref->[$c];
    }
  } 

  return @coords;
}
