package ASCbaseUtils;

####################
sub calc_norm_stats{
  my($query, $data_dir, $proj_bin) = @_;

  system "$proj_bin/site_dock --dbase=\"$data_dir/diverse_sitemaps\" $query"; 

  # Assume that the query filename ends with the .pdb extension and may have 
  # some specified path
  my $infile = substr($query, 0, -4);
  my $idx = rindex($query, "/") + 1;
  if($idx > 0){ $infile = substr($infile, $idx); }
  $infile .= "_blast.out";
  open IN, $infile or die "Unable to open the file $infile: $!\n";

  $tot_prof_time = 0;
  my @scores;
  while(<IN>){
    substr($_,0,1) ne "%" or next;
    my @args = split(/\|/, $_);
    push(@scores, $args[1]);
    $tot_prof_time += $args[4];
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

1;
