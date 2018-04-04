#!/usr/bin/perl

opendir DIR, "../pdb_ligands" or die "Unable to open ../pdb_ligands: $!\n";
#@files = readdir DIR;
@files = grep /.*.pdb/, readdir DIR;
close DIR;
print "@files\n";

foreach $file (@files){
  $out = substr($file, 0, length($file) - 4) . "_l.mol2";
  system "molcharge -in ../pdb_ligands/$file -out $out";
}
