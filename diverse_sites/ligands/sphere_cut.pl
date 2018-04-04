#!/usr/bin/perl

$#ARGV == 0 or die "USAGE: $0 <radius>\n";
my $radius = $ARGV[0];

opendir DIR, "." or die "Unable to open . $!\n";
my @files = grep /.*.mol2/, readdir DIR;
close DIR;
my @ligs = sort @files;

opendir DIR, "../dbase" or die "Unable to open ../dbase $!\n";
my @files = grep /.*.pdb/, readdir DIR;
close DIR;
my @sites = sort @files;

my @atoms; 
my $rad_2 = $radius * $radius;
for($ii = 0; $ii <= $#sites; $ii++){
  my $id = substr($ligs[$ii],0, -7);
  #if($id eq substr($sites[$ii], 0, length($id))){

  $#atoms = -1;
  get_mol2_atoms($ligs[$ii], \@atoms);

  open IN, "../dbase/$sites[$ii]" 
    or die "Unable to open ../dbase/$sites[$ii]: $!\n";
  open OUT, "+>../dbase_cut/$id" . "_s.pdb"
    or die "Unable to open ../dbase_cut/$sites[$ii]: $!\n";
  # For each template point if it is not with in $radius of any ligand atom
  # omit
  while(<IN>){
    /(.{6}).{24}(.{8})(.{8})(.{8})/;
    if($1 eq "HETATM" || $1 eq "ATOM  "){ 
      my @pt;
      push(@pt, $2, $3, $4);

      for($a = 0; $a <= $#atoms; $a++){
        my $dist_2 = 0;
        for($jj = 0; $jj < 3; $jj++){
          my $tmp = $pt[$jj] - $atoms[$a]->[$jj];
          $dist_2 += $tmp * $tmp;
        }

        if($dist_2 <= $rad_2){
          print OUT $_;
          last; 
        }
      }

    }else{ print OUT $_; }
  }
  close OUT;
  close IN;
}



###########################################################
sub get_mol2_atoms(){
  my ($fname, $atm_ref) = @_;
  open MOL2, $fname or die "Unable to open $fname: $!\n";
  my $flag = 0;
  while(<MOL2>){
    if(substr($_, 0, 13) eq "@<TRIPOS>ATOM"){ 
      $flag = 1; 
      next;
    }

    if($flag){
      my @tmp = split(/ +/, $_);
      push(@{$atm_ref->[$#atoms+1]}, $tmp[3], $tmp[4], $tmp[5]);
    }
  }
}
