package Grinder::Database;

use strict;
use warnings;
use Bio::DB::Fasta;

use base qw(Bio::Root::Root); # using throw() and _rearrange() methods

sub new {
   my ($class, @args) = @_;
   my $self = $class->SUPER::new(@args);
   my ($fasta_file, $unidirectional, $primers, $abundance_file, $delete_chars,
      $minimum_length) = $self->_rearrange([qw(FASTA_FILE UNIDIRECTIONAL PRIMERS
      ABUNDANCE_FILE DELETE_CHARS MINIMUM_LENGTH)], @args);

   $minimum_length = 1  if not defined $minimum_length;
   $self->_set_minimum_length($minimum_length);

   $delete_chars   = '' if not defined $delete_chars;
   $self->_set_delete_chars($delete_chars);

   # Index file, filter sequences and get IDs
   $self->_init_db($fasta_file, $abundance_file, $delete_chars, $minimum_length);

   $unidirectional = 0  if not defined $unidirectional; # bidirectional
   $self->_set_unidirectional($unidirectional);

   # Read amplicon primers
   $self->_set_primers($primers) if defined $primers;

  # Error if using amplicon on protein database
  if ( ($self->get_alphabet eq 'protein') && ($self->get_unidirectional != 1) ) {
    $self->throw("Got <unidirectional> = $unidirectional but can only use ".
      "<unidirectional> = 1 with proteic reference sequences\n");
  }

   return $self;
}


sub get_primers {
   my ($self) = @_;
   return $self->{'primers'};
}


sub _set_primers {
  my ($self, $forward_reverse_primers) = @_;
  # Read primer file and convert primers into regular expressions to catch
  # amplicons present in the database
  if (defined $forward_reverse_primers) {

    # Read primers from FASTA file
    my $primer_in = Bio::SeqIO->newFh(
      -file   => $forward_reverse_primers,
      -format => 'fasta',
    );

    # Mandatory first primer
    my $primer = <$primer_in>;
    if (not defined $primer) {
      $self->throw("The file '$forward_reverse_primers' contains no primers\n");
    }
    $primer->alphabet('dna'); # Force the alphabet since degenerate primers can look like protein sequences
    $self->_set_forward_regexp( iupac_to_regexp($primer->seq) );
    $primer = undef;

    # Take reverse-complement of optional reverse primers
    $primer = <$primer_in>;
    if (defined $primer) {
      $primer->alphabet('dna');
      $primer = $primer->revcom;
      $self->_set_reverse_regexp( iupac_to_regexp($primer->seq) );
    }

  }
  $self->{'primers'} = $forward_reverse_primers;
  return $self->get_primers;
}


sub get_forward_regexp {
   my ($self) = @_;
   return $self->{'forward_regexp'};
}


sub _set_forward_regexp {
   my ($self, $val) = @_;
   $self->{'forward_regexp'} = $val;
   return $self->get_forward_regexp;
}


sub get_reverse_regexp {
   my ($self) = @_;
   return $self->{'reverse_regexp'};
}


sub _set_reverse_regexp {
   my ($self, $val) = @_;
   $self->{'reverse_regexp'} = $val;
   return $self->get_reverse_regexp;
}


sub _set_alphabet {
   my ($self, $val) = @_;
   $self->{'alphabet'} = $val;
   return $self->get_alphabet;
}


sub get_alphabet {
   my ($self) = @_;
   return $self->{'alphabet'};
}


sub _set_ids {
   my ($self, $val) = @_;
   $self->{'ids'} = $val;
   return $self->get_ids;
}


sub get_ids {
   my ($self) = @_;
   return $self->{'ids'};
}


sub _set_unidirectional {
   my ($self, $val) = @_;
   # Error if using wrong direction on protein database
   if ( ($self->get_alphabet eq 'protein') && ($val != 1) ) {
      $self->throw("Got <unidirectional> = $val but can only use ".
         "<unidirectional> = 1 with proteic reference sequences\n");
   }
   $self->{'unidirectional'} = $val;
   return $self->get_unidirectional;
}


sub get_unidirectional {
   my ($self) = @_;
   return $self->{'unidirectional'};
}


sub _set_minimum_length {
   my ($self, $val) = @_;
   $self->{'minimum_length'} = $val;
   return $self->get_minimum_length;
}


sub get_minimum_length {
   my ($self) = @_;
   return $self->{'minimum_length'};
}


sub _set_delete_chars {
   my ($self, $val) = @_;
   $self->{'delete_chars'} = $val;
   return $self->get_delete_chars;
}


sub get_delete_chars {
   my ($self) = @_;
   return $self->{'delete_chars'};
}


sub _set_database {
   my ($self, $val) = @_;
   $self->{'database'} = $val;
   return $self->get_database;
}


sub get_database {
   my ($self) = @_;
   return $self->{'database'};
}


sub _init_db {
   # Read and import sequences
   # Parameters:
   #   * FASTA file containing the sequences or '-' for stdin. REQUIRED
   #   * Abundance file (optional): To avoid registering unwanted sequences
   #   * Delete chars (optional): Characters to delete from the sequences.
   #   * Minimum sequence size: Skip sequences smaller than that
   my ($self, $fasta_file, $abundance_file, $delete_chars, $min_len) = @_;

   # Get list of all IDs with a manually-specified abundance
   my %ids_to_keep;
   my $nof_ids_to_keep = 0;
   if ($abundance_file) {
      my ($ids) = community_read_abundances($abundance_file);
      for my $comm_num (0 .. $#$ids) {
         for my $gen_num ( 0 .. scalar @{$$ids[$comm_num]} - 1 ) {
            my $id = $$ids[$comm_num][$gen_num];
            $ids_to_keep{$id} = undef;
            $nof_ids_to_keep++;
         }
      }
   }

   # Index input file
   my $db = Bio::DB::Fasta->new($fasta_file, -reindex => 1);
   $self->_set_database($db);

   # List sequences that are ok to use
   my @seq_ids;
   my $nof_seqs;
   my %mol_types;
   my $stream = $db->get_PrimarySeq_stream;
   while (my $seq = $stream->next_seq) {

      # Skip empty sequences
      next if not $seq->seq;

      # Record molecule type
      $mol_types{$seq->alphabet}++;

      # Skip unwanted sequences
      my $seq_id = $seq->id;
      next if ($nof_ids_to_keep > 0) && (not exists $ids_to_keep{$seq_id});

      # Remove specified characters
      $seq = $self->_remove_chars($seq, $delete_chars);

      # Skip the sequence if it is too small
      next if $seq->length < $min_len;

      # Keep this sequence
      push @seq_ids, $seq->id;
      $nof_seqs++;
   }

   # Error if no usable sequences in the database
   if ($nof_seqs == 0) {
      $self->throw("No genome sequences could be used. If you specified a file ".
         "of abundances for the genome sequences, make sure that their ID match".
         " the ID in the FASTA file. If you specified amplicon primers, verify ".
         "that they match some genome sequences.\n");
   }

   # Determine database type: dna, rna, protein
   my $db_alphabet = $self->_set_alphabet( $self->_get_mol_type(\%mol_types) );

   # Record the sequence IDs
   $self->_set_ids( \@seq_ids );

   return $db;
}


sub _remove_chars {
   # Remove forbidden chars
   my ($self, $seq, $chars) = @_;
   if ( defined($chars) && not($chars eq '') ) {

      ####
      #print "Removing chars: $chars\n";
      #print "1/ seq->seq: '".$seq->seq."'\n";
      ####

      my $seq_string = $seq->seq;

      ####
      #print "   a) seq_string: '".$seq_string."'\n";
      ####

      $seq_string =~ s/[$chars]//gi;

      ####
      #print "   b) seq_string: '".$seq_string."'\n";
      ####

      $seq->seq( $seq_string );

      ####
      #print "2/ seq->seq: '".$seq->seq."'\n";
      ####

   }
   return $seq;
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



sub next_seq {
   my ($self) = @_;
   # do something about unidirectional
   # them delete chars
   # then fetch amplicons
   # finally remove seqs < min_len

   my $unidirectional = $self->get_unidirectional;
   my $delete_chars   = $self->get_delete_chars;
}

#    # If we are sequencing from the reverse strand, reverse complement now
#    if ($unidirectional == -1) {
#       $ref_seq = $ref_seq->revcom;
#    }
#
#    # Extract amplicons if needed
#    my $amp_seqs;
#    if (defined $forward_regexp) {
#      $amp_seqs = $self->database_extract_amplicons($ref_seq, $forward_regexp,
#        $reverse_regexp, \%ids_to_keep);
#      next if scalar @$amp_seqs == 0;
#    } else {
#      $amp_seqs = [$ref_seq];
#    }
#
#    for my $amp_seq (@$amp_seqs) {
#      # Remove forbidden chars
#      if ( (defined $delete_chars) && (not $delete_chars eq '') ) {
#        my $clean_seq = $amp_seq->seq;
#        $clean_seq =~ s/[$delete_chars]//gi;
#        $amp_seq->seq($clean_seq);
#      }
#      # Skip the sequence if it is too small
#      next if $amp_seq->length < $min_len;
#      # Save amplicon sequence and identify them by their unique object reference
#      $seq_db{$amp_seq} = $amp_seq;
#      $seq_ids{$ref_seq_id}{$amp_seq} = undef;
#    }
#
#


1;
