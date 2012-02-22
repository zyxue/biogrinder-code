package Grinder::Database;

use strict;
use warnings;
use Bio::SeqIO;

use base qw(Bio::Root::Root); # using throw() and _rearrange() methods

sub new {
   my ($class, @args) = @_;
   my $self = $class->SUPER::new(@args);
   my ($fasta_file, $unidirectional, $forward_reverse_primers, $abundance_file,    
      $delete_chars, $min_len) = $self->_rearrange([qw(FASTA_FILE UNIDIRECTIONAL
      FORWARD_REVERSE_PRIMERS ABUNDANCE_FILE DELETE_CHARS MIN_LEN)], @args);

   #$self->k( defined $k ? $k : 10 );
   #$self->weights($weights)     if defined $weights;
   #$self->add_seqs($seqs, $ids) if defined $seqs;
   #$self->add_file($file)       if defined $file;

   # Defaults
   $unidirectional = 0 if not defined $unidirectional; # bidirectional
   $min_len = 1 if not defined $min_len;

   $self->_create($fasta_file, $unidirectional, $forward_reverse_primers,
      $abundance_file, $delete_chars, $min_len);

   return $self;
}


#sub k {
#   my ($self, $val) = @_;
#   if ($val) {
#      if ($val < 1) {
#         $self->throw("Error: The minimum kmer length is 1 but got $val\n");
#      }
#      $self->{'k'} = $val;
#   }
#   return $self->{'k'};
#}


sub _create {
  # Read and import sequences
  # Parameters:
  #   * FASTA file containing the sequences or '-' for stdin. REQUIRED
  #   * Sequencing unidirectionally? 0: no, 1: yes forward, -1: yes reverse
  #   * Amplicon PCR primers (optional): Should be provided in a FASTA file and
  #     use the IUPAC convention. If a primer sequence is given, any sequence
  #     that does not contain the primer (or its reverse complement for the
  #     reverse primer) is skipped, while any sequence that match is trimmed so
  #     that it is flush with the primer sequence.
  #   * Abundance file (optional): To avoid registering sequences in the database
  #     unless they are needed
  #   * Delete chars (optional): Characters to delete form the sequences.
  #   * Minimum sequence size: Skip sequences smaller than that
  my ($self, $fasta_file, $unidirectional, $forward_reverse_primers,
    $abundance_file, $delete_chars, $min_len) = @_;

  # Input filehandle
  if (not defined $fasta_file) {
    $self->throw("No reference sequences provided\n");
  }
  my $in;
  if ($fasta_file eq '-') {
    $in = Bio::SeqIO->newFh(
      -fh     => \*STDIN,
      -format => 'fasta',
    );
  } else {
    $in = Bio::SeqIO->newFh(
      -file   => $fasta_file,
      -format => 'fasta',
    );
  }
  # Regular expression to catch amplicons present in the database
  my ($forward_regexp, $reverse_regexp);
  if (defined $forward_reverse_primers) {
    # Read primers from FASTA file
    my $primer_in = Bio::SeqIO->newFh(
      -file   => $forward_reverse_primers,
      -format => 'fasta',
    );
    my $first_primer = <$primer_in>;
    if (not defined $first_primer) {
      $self->throw("The file '$forward_reverse_primers' does not contain any primers\n");
    } else {
      # Force the alphabet since degenerate primers can look like protein sequences
      $first_primer->alphabet('dna');
    }
    my $second_primer = <$primer_in>;
    if (defined $second_primer) {
      $second_primer->alphabet('dna');
    }
    undef $primer_in;
    # Regexp for forward primer and optional reverse complement of reverse primer
    $forward_regexp = iupac_to_regexp($first_primer->seq);
    if (defined $second_primer) {
      $second_primer = $second_primer->revcom;
      $reverse_regexp = iupac_to_regexp($second_primer->seq);
    }
  }
  # Get list of all IDs with a manually-specified abundance
  my %ids_to_keep;
  if ($abundance_file) {
    my ($ids) = community_read_abundances($abundance_file);
    for my $comm_num (0 .. $#$ids) {
      for my $gen_num ( 0 .. scalar @{$$ids[$comm_num]} - 1 ) {
        my $id = $$ids[$comm_num][$gen_num];
        $ids_to_keep{$id} = undef;
      }
    }
  }
  # Process database sequences
  my %seq_db;      # hash of BioPerl sequence objects (all amplicons)
  my %seq_ids;     # hash of reference sequence IDs and IDs of their amplicons
  my %mol_types;    # hash of count of molecule types (dna, rna, protein)
  while ( my $ref_seq = <$in> ) {
    # Skip empty sequences
    next if not $ref_seq->seq;
    # Record molecule type
    $mol_types{$ref_seq->alphabet}++;
    # Skip unwanted sequences
    my $ref_seq_id = $ref_seq->id;
    next if (scalar keys %ids_to_keep > 0) && (not exists $ids_to_keep{$ref_seq_id});
    # If we are sequencing from the reverse strand, reverse complement now
    if ($unidirectional == -1) {
      $ref_seq = $ref_seq->revcom;
    }
    # Extract amplicons if needed
    my $amp_seqs;
    if (defined $forward_regexp) {
      $amp_seqs = $self->database_extract_amplicons($ref_seq, $forward_regexp,
        $reverse_regexp, \%ids_to_keep);
      next if scalar @$amp_seqs == 0;
    } else {
      $amp_seqs = [$ref_seq];
    }

    for my $amp_seq (@$amp_seqs) {
      # Remove forbidden chars
      if ( (defined $delete_chars) && (not $delete_chars eq '') ) {
        my $clean_seq = $amp_seq->seq;
        $clean_seq =~ s/[$delete_chars]//gi;
        $amp_seq->seq($clean_seq);
      }
      # Skip the sequence if it is too small
      next if $amp_seq->length < $min_len;
      # Save amplicon sequence and identify them by their unique object reference
      $seq_db{$amp_seq} = $amp_seq;
      $seq_ids{$ref_seq_id}{$amp_seq} = undef;
    }

  }
  undef $in; # close the filehandle (maybe?!)

  # Error if no usable sequences in the database
  if (scalar keys %seq_ids == 0) {
    $self->throw("No genome sequences could be used. If you specified a file ".
      "of abundances for the genome sequences, make sure that their ID match ".
      "the ID in the FASTA file. If you specified amplicon primers, verify ".
      "that they match some genome sequences.\n");
  }

  # Determine database type: dna, rna, protein
  my $db_alphabet = $self->_get_mol_type(\%mol_types);
  $self->{alphabet} = $db_alphabet;

  # Error if using amplicon on protein database
  if ( ($db_alphabet eq 'protein') && (defined $forward_reverse_primers) ) {
    $self->throw("Cannot use amplicon primers with proteic reference sequences\n");
  }

  # Error if using amplicon on protein database
  if ( ($db_alphabet eq 'protein') && ($unidirectional != 1) ) {
    $self->throw("Got <unidirectional> = $unidirectional but can only use ".
      "<unidirectional> = 1 with proteic reference sequences\n");
  }

  my $database = { 'db' => \%seq_db, 'ids' => \%seq_ids };
  return $database;
}


sub _get_mol_type {
  # Given a count of the different molecule types in the database, determine
  # what molecule type it is.
  my ($self, $mol_types) = @_;
  my $max_count = 0;
  my $max_type  = '';
  while (my ($type, $count) = each %$mol_types) {
    if ($count > $max_count) {
      $max_count = $count;
      $max_type  = $type;
    }
  }
  my $other_count = 0;
  while (my ($type, $count) = each %$mol_types) {
    if (not $type eq $max_type) {
      $other_count += $count;
    }
  }
  if ($max_count < $other_count) {
    $self->throw("Cannot determine to what type of molecules the reference ".
        "sequences belong. Got $max_count sequences of type '$max_type' and ".
        "$other_count others.\n");
  }
  if ( (not $max_type eq 'dna') &&
       (not $max_type eq 'rna') &&
       (not $max_type eq 'protein') ) {
    $self->throw("Reference sequences are in an unknown alphabet '$max_type'.\n");
  }
  return $max_type;
}


1;
