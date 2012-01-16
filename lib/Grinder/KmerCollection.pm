# This file is part of the Grinder package, copyright 2009,2010,2011
# Florent Angly <florent.angly@gmail.com>, under the GPLv3 license


package Grinder::KmerCollection;

=head1 NAME

Grinder::KmerCollection - A collection of kmers from sequences

=head1 SYNOPSIS

  my $col = Grinder::KmerCollection->new( -k    => 10,
                                          -file => 'seqs.fa' );

=head1 DESCRIPTION

Manage a collection of kmers found in various sequences. Store information about
what sequence a kmer was found in and its starting position on the sequence.

=head1 AUTHOR

Florent Angly <florent.angly@gmail.com>

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


use strict;
use warnings;
use Bio::SeqIO;

use base qw(Bio::Root::Root);

=head2 new

 Title   : new
 Usage   : my $col = Grinder::KmerCollection->new( -k => 10, -file => 'seqs.fa', -revcom => 1 );
 Function: Build a new kmer collection
 Args    : -k        set the kmer length (default: 10 bp)
           -revcom   count kmers before and after reverse-complementing sequences
                     (default: 0)
           -seqs     count kmers in the provided arrayref of sequences (Bio::Seq
                     objects)
           -ids      if specified, index the sequences provided to -seq using the
                     use the IDs in this arrayref instead of using the sequences
                     $seq->id() method
           -file     count kmers in the provided file of sequences
           -weights  if specified, assign the abundance of each sequence from the
                     values in this arrayref

 Returns : Grinder::KmerCollection object

=cut

sub new {
   my ($class, @args) = @_;
   my $self = $class->SUPER::new(@args);
   my($k, $revcom, $seqs, $ids, $file) =
     $self->_rearrange([qw(K REVCOM SEQS IDS FILE)], @args);

   $self->k( defined $k ? $k : 10 );

   $self->add_seqs($seqs, $ids) if defined $seqs;
   $self->add_file($file)       if defined $file;

   ####$self->add_weights($weights)           if defined $weights;
  
   return $self;
}


=head2 k

 Usage   : $col->k;
 Function: Get the length of the kmers
 Args    : None
 Returns : Positive integer

=cut

sub k {
   my ($self, $val) = @_;
   if ($val) {
      if ($val < 1) {
         die("Error: The minimum kmer length is 1 but got $val\n");
      }
      $self->{'k'} = $val;
   }
   return $self->{'k'};
}


=head2 collection_by_kmer

 Usage   : $col->collection_by_kmer;
 Function: Get the collection of kmers, indexed by kmer
 Args    : None
 Returns : A hashref of hashref of arrayref:
              hash->{kmer}->{ID of sequences with this kmer}->[starts of kmer on sequence]

=cut

sub collection_by_kmer {
   my ($self, $val) = @_;
   if ($val) {
      $self->{'collection_by_kmer'} = $val;
   }
   return $self->{'collection_by_kmer'};
}


=head2 collection_by_seq

 Usage   : $col->collection_by_seq;
 Function: Get the collection of kmers, indexed by sequence ID
 Args    : None
 Returns : A hashref of hashref of arrayref:
              hash->{ID of sequences with this kmer}->{kmer}->[starts of kmer on sequence]

=cut

sub collection_by_seq {
   my ($self, $val) = @_;
   if ($val) {
      $self->{'collection_by_seq'} = $val;
   }
   return $self->{'collection_by_seq'};
}


#==============================================================================#


=head2 add_file

 Usage   : $col->add_file('seqs.fa');
 Function: Process the kmers in the given file of sequences.
 Args    : filename
 Returns : Grinder::KmerCollection object

=cut

sub add_file {
   my ($self, $file) = @_;
   my $in = Bio::SeqIO->new( -file => $file );
   while (my $seq = $in->next_seq) {
      $self->add_seqs([ $seq ]);
   }
   $in->close;
   return $self;
}


=head2 add_seqs

 Usage   : $col->add_seqs([$seq1, $seq2]);
 Function: Process the kmers in the given sequences.
 Args    : * arrayref of Bio::Seq objects
           * arrayref of IDs to use for the indexing of the sequences
 Returns : Grinder::KmerCollection object

=cut

sub add_seqs {
   my ($self, $seqs, $ids) = @_;
   my $col_by_kmer = $self->collection_by_kmer || {};
   my $col_by_seq  = $self->collection_by_seq  || {};
   my $i = 0;
   for my $seq (@$seqs) {
      my $kmer_counts = $self->_count_kmers($seq);
      while ( my ($kmer, $positions) = each %$kmer_counts ) {
         my $seq_id;
         if (defined $ids) {
           $seq_id = $$ids[$i];
         } else {
           $seq_id = $seq->id;
         }
         $col_by_kmer->{$kmer}->{$seq_id} = $positions;
         $col_by_seq->{$seq_id}->{$kmer}  = $positions;
      }
      $i++;
   }
   $self->collection_by_kmer($col_by_kmer);
   $self->collection_by_seq($col_by_seq);
   return $self;
}


=head2 add_weights

 Usage   : $col->add_weights({'seq1' => 3, 'seq10' => 0.45});
 Function: Assign weights to sequences
 Args    : hashref where the keys are sequence IDs and the values are the weight
           of the corresponding (e.g. their relative abundance)
 Returns : Grinder::KmerCollection object

=cut

sub add_weights {
   my ($self, $weights) = @_;
   #####
   # TODO
   die "Error: Not implemented\n";
   #####
   return $self;
}


=head2 filter_rare

 Usage   : $col->filter_rare( 2 );
 Function: Remove kmers occurring less than the number of times specified
 Args    : integer
 Returns : Grinder::KmerCollection object

=cut

sub filter_rare {
   my ($self, $min_num) = @_;
   my $changed = 0;
   my $col_by_kmer = $self->collection_by_kmer;
   my $col_by_seq  = $self->collection_by_seq;
   while ( my ($kmer, $sources) = each %$col_by_kmer ) {
      my $count = _sum_from_sources( $sources );
      if ($count < $min_num) {
        # Remove this kmer
        $changed = 1;
        delete $col_by_kmer->{$kmer};
        while ( my ($seq, $seq_kmers) = each %$col_by_seq ) {
          delete $seq_kmers->{$kmer};
          delete $col_by_seq->{$seq} if keys %{$seq_kmers} == 0;
        }
      }
   }
   if ($changed) {
     $self->collection_by_kmer( $col_by_kmer );
     $self->collection_by_seq( $col_by_seq );
   }
   return $self;
}


=head2 filter_shared

 Usage   : $col->filter_shared( 2 );
 Function: Remove kmers occurring in less than the number of sequences specified
 Args    : integer
 Returns : Grinder::KmerCollection object

=cut

sub filter_shared {
   my ($self, $min_num) = @_;
   my $changed = 0;
   my $col_by_kmer = $self->collection_by_kmer;
   my $col_by_seq  = $self->collection_by_seq;
   while ( my ($kmer, $sources) = each %$col_by_kmer ) {
      my $count = scalar keys %$sources;
      if ($count < $min_num) {
        $changed = 1;
        delete $col_by_kmer->{$kmer};
        while ( my ($seq, $seq_kmers) = each %$col_by_seq ) {
          delete $seq_kmers->{$kmer};
          delete $col_by_seq->{$seq} if keys %{$seq_kmers} == 0;
        }
      }
   }
   if ($changed) {
     $self->collection_by_kmer( $col_by_kmer );
     $self->collection_by_seq( $col_by_seq );
   }
   return $self;
}


##### should be able to specify a single kmer to get the count of?

=head2 counts

 Usage   : $col->counts
 Function: Calculate the total count of each kmer
 Args    : 0 to report counts (default), 1 to report frequencies (normalize to 1)
 Returns : * arrayref of the different kmers
           * arrayref of the corresponding total counts

=cut

sub counts {
   my ($self, $freq) = @_;
   my $kmers;
   my $counts;
   my $total = 0;
   my $col_by_kmer = $self->collection_by_kmer;
   while ( my ($kmer, $sources) = each %$col_by_kmer ) {
      push @$kmers, $kmer;
      my $count = _sum_from_sources( $sources );
      push @$counts, $count;
      $total += $count;
   }
   if ($freq && $total) {
     $counts = _normalize($counts, $total);
   }
   return $kmers, $counts;
}


=head2 sources

 Usage   : $col->sources()
 Function: Return the sources of a kmer and their abundance.
 Args    : * kmer to get the sources of
           * sources to exclude from the results
           * 0 to report counts (default), 1 to report frequencies (normalize to 1)
 Returns : * arrayref of the different sources
           * arrayref of the corresponding total counts
           If the kmer requested does not exist, the array will be empty.

=cut

sub sources {
   my ($self, $kmer, $excl, $freq) = @_;
   my $sources = [];
   my $counts = [];
   my $total = 0;
   my $kmer_sources = $self->collection_by_kmer->{$kmer};

   if (defined $kmer_sources) {
      while ( my ($source, $positions) = each %$kmer_sources ) {
         if ( (defined $excl) && ($source eq $excl) ) {
           next;
         }
         push @$sources, $source;
         my $count = scalar @$positions;
         push @$counts, $count;
         $total += $count;
      }
      $counts = _normalize($counts, $total) if $freq; 
   }

   return $sources, $counts;
}


=head2 kmers

 Usage   : $col->kmers('seq1');
 Function: This is the inverse of sources(). Return the kmers found in a sequence
           (given its ID).
 Args    : * sequence ID to get the kmers of
           * 0 to report counts (default), 1 to report frequencies (normalize to 1)
 Returns : * arrayref of sequence IDs
           * arrayref of the corresponding total counts
           If the sequence ID requested does not exist, the arrays will be empty.
=cut

sub kmers {
   my ($self, $seq_id, $freq) = @_;
   my $kmers = [];
   my $counts = [];
   my $total = 0;
   my $seq_kmers = $self->collection_by_seq->{$seq_id};

   if (defined $seq_kmers) {
      while ( my ($kmer, $positions) = each %$seq_kmers ) {
         push @$kmers, $kmer;
         my $count = scalar @$positions;
         push @$counts, $count;
         $total += $count;
      }
      $counts = _normalize($counts, $total) if $freq;
   }

   return $kmers, $counts;
}


=head2 positions

 Usage   : $col->positions()
 Function: Return the positions of the given kmer on a given sequence. An error
           is reported if the kmer requested does not exist
 Args    : * desired kmer
           * desired sequence with this kmer
 Returns : Arrayref of the different positions. The arrays will be empty if the
           desired combination of kmer and sequence was not found.

=cut

sub positions {
   my ($self, $kmer, $source) = @_;
   my $kmer_positions = [];
   my $kmer_sources = $self->collection_by_kmer->{$kmer};
   if (defined $kmer_sources) {
      $kmer_positions = $kmer_sources->{$source} || [];
   }
   return $kmer_positions;
}



#======== Internals ===========================================================#


sub _count_kmers {
   # Count the kmers of size k in a sequence (Bio::Seq) and return a hashref.
   my ($self, $seq) = @_;
   my $k = $self->k;
   my $seq_str = uc $seq->seq; # case-insensitive
   my $seq_len = length $seq_str;
   my $hash = {};
   for (my $i = 0; $i <= $seq_len - $k ; $i++) {
      my $kmer = substr $seq_str, $i, $k;
      push @{$hash->{$kmer}}, $i + 1;
   }
   return $hash;
}


sub _sum_from_sources {
   # Calculate the number of occurences of a kmer.
   my ($sources) = @_;
   my $count = 0;
   while ( my ($source, $positions) = each %$sources ) {
      $count += scalar @$positions;
   }
   return $count;
}


sub _normalize {
   # Normalize an arrayref to 1.
   my ($arr, $total) = @_;
   if (not $total) { # total undef or 0
      die "Error: Need to provide a valid total\n";
   }
   $arr = [ map {$_ / $total} @$arr ];
   return $arr;
}

1;
