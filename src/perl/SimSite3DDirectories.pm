package SimSite3DDirectories;
use Carp;

sub new{
  my $type = shift;

  # Directory where SimSite3D is installed
  $self->{project} = $ENV{'ASCBASE_SOFTWARE_DIR'} or croak "ERROR: " .
    "environment variable ASCBASE_SOFTWARE_DIR is not set\n";
  # Added so that it can be changed if needed.
  $self->{data} = $self->{project} . "/data";
  # Directory holding the query SimSite3D sitemap pdb file
  $self->{queries} = $self->{data} . "/query_sitemaps";
  # Directory holding the Pfizer protein pdb files
  $self->{proteins} = $ENV{'PROTEIN_CRYSTAL_STRUCTURES_DIR'} or croak "ERROR: ".
    "environment variable PROTEIN_CRYSTAL_STRUCTURES_DIR is not set\n";
  # Directory holding the Pfizer ligand mol2 files 
  $self->{ligands} = $ENV{'LIGAND_CRYSTAL_STRUCTURES_DIR'} or croak "ERROR: ".
    "environment variable LIGAND_CRYSTAL_STRUCTURES_DIR is not set\n";
  # Directory holding the searchable SimSite3D sitemap pdb files
  $self->{search} = $ENV{'ASCBASE_SEARCHABLE_SITEMAPS_DIR'} or croak "ERROR: ".
    "environment variable ASCBASE_SEARCHABLE_SITEMAPS_DIR is not set\n";

  return bless $self, $type;
}

1;
