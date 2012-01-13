
#### todo:
# revcom
# weights


package KmerCollection;

=head1 NAME

KmerCollection - A collection of kmers from sequences

=head1 SYNOPSIS

    my $col = KmerCollection->new( -k    => 10,
                                   -file => 'seqs.fa' );

=head1 DESCRIPTION

Manage a collection of kmers found in various sequences. Store information about
what sequence a kmer was found in and its starting position on the sequence.

#=head1 FEEDBACK

#=head2 Mailing Lists

#User feedback is an integral part of the evolution of this and other
#Bioperl modules. Send your comments and suggestions preferably to
#the Bioperl mailing list.  Your participation is much appreciated.

#  bioperl-l@bioperl.org                  - General discussion
#  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

#=head2 Support 

#Please direct usage questions or support issues to the mailing list:

#I<bioperl-l@bioperl.org>

#rather than to the module maintainer directly. Many experienced and 
#reponsive experts will be able look at the problem and quickly 
#address it. Please include a thorough description of the problem 
#with code and data examples if at all possible.

#=head2 Reporting Bugs

#Report bugs to the Bioperl bug tracking system to help us keep track
#of the bugs and their resolution. Bug reports can be submitted via the
#web:

#  https://redmine.open-bio.org/projects/bioperl/

=head1 AUTHOR

Florent Angly <florent.angly@gmail.com>

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


use strict;
use warnings;

use base qw(Bio::Root::Root);

=head2 new

 Title   : new
 Usage   : my $col = KmerCollection->new( -k => 10, -file => 'seqs.fa', -revcom => 1 );
 Function: Build a new kmer collection
 Args    : -k       set the kmer length (default: 10 bp)
           -revcom  count kmers before and after reverse-complementing sequences
                    (default: 0)
           -seqs    count kmers in the provided arrayref of sequences (Bio::Seq
                    objects)
           -file    count kmers in the provided file of sequences
 Returns : KmerCollection object

=cut

sub new {
   my ($class, @args) = @_;
   my $self = $class->SUPER::new(@args);
   my($k, $revcom, $seqs, $file) = $self->_rearrange([qw(K REVCOM SEQS FILE)], @args);

   $self->k( defined $k ? $k : 10 );
   $self->revcom( defined $revcom ? $revcom : 0 );

   $self->add_seqs($seqs) if defined $seqs;
   $self->add_file($file) if defined $file;

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
         $self->throw("Error: The minimum kmer length is 1 but got $val\n");
      }
      $self->{'k'} = $val;
   }
   return $self->{'k'};
}


### TODO

=head2 revcom

 Usage   : $col->revcom;
 Function: Get whether or not the kmers are calculated from the sequence and
           their reverse-complement
 Args    : None
 Returns : 1 (yes), 0 (no)

=cut

sub revcom {
   my ($self, $val) = @_;
   if ($val) {
      $self->{'revcom'} = $val;
   }
   return $self->{'revcom'};
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

 Usage   : $col->add_file( 'seqs.fa' );
 Function: Process the kmers in the given file of sequences.
 Args    : filename
 Returns : KmerCollection object

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

 Usage   : $col->add_seqs( [$seq1, $seq2] );
 Function: Process the kmers in the given sequences.
 Args    : arrayref of Bio::Seq objects
 Returns : KmerCollection object

=cut

sub add_seqs {
   my ($self, $seqs) = @_;
   my $col_by_kmer = $self->collection_by_kmer || {};
   my $col_by_seq  = $self->collection_by_seq  || {};
   for my $seq (@$seqs) {
      my $kmer_counts = $self->_count_kmers($seq);
      while ( my ($kmer, $positions) = each %$kmer_counts ) {
         my $seq_id = $seq->id;
         $col_by_kmer->{$kmer}->{$seq_id} = $positions;
         $col_by_seq->{$seq_id}->{$kmer}  = $positions;
      }
   }
   $self->collection_by_kmer($col_by_kmer);
   $self->collection_by_seq($col_by_seq);
   return $self;
}


=head2 filter_rare

 Usage   : $col->filter_rare( 2 );
 Function: Remove kmers occurring less than the number of times specified
 Args    : integer
 Returns : KmerCollection object

=cut

sub filter_rare {
   my ($self, $min_num) = @_;
   my $col_by_kmer = $self->collection_by_kmer;
   while ( my ($kmer, $sources) = each %$col_by_kmer ) {
      my $count = _sum_from_sources( $sources );
      delete $col_by_kmer->{$kmer} if $count < $min_num;
   }
   $self->collection_by_kmer( $col_by_kmer );
   return $self;
}


=head2 filter_shared

 Usage   : $col->filter_shared( 2 );
 Function: Remove kmers occurring in less than the number of sequences specified
 Args    : integer
 Returns : KmerCollection object

=cut

sub filter_shared {
   my ($self, $min_num) = @_;
   my $col_by_kmer = $self->collection_by_kmer;
   while ( my ($kmer, $sources) = each %$col_by_kmer ) {
      my $count = scalar keys %$sources;
      delete $col_by_kmer->{$kmer} if $count < $min_num;
   }
   $self->collection_by_kmer( $col_by_kmer );
   return $self;
}


##### should be able to specify a single kmer to get the count of

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
   $counts = _normalize($counts, $total) if $freq;
   return $kmers, $counts;
}


=head2 sources

 Usage   : $col->sources()
 Function: Return the sources of a kmer and their abundance. An error is reported
           if the kmer requested does not exist
 Args    : * kmer to get the sources of
           * 0 to report counts (default), 1 to report frequencies (normalize to 1)
 Returns : * arrayref of the different sources
           * arrayref of the corresponding total counts

=cut

sub sources {
   my ($self, $kmer, $freq) = @_;
   my $sources;
   my $counts;
   my $total = 0;
   my $kmer_sources = $self->collection_by_kmer->{$kmer};
  
   if (not defined $kmer_sources) {
      $self->throw("Error: kmer $kmer was not found in the collection.\n");
   }

   while ( my ($source, $positions) = each %$kmer_sources ) {
      push @$sources, $source;
      my $count = scalar @$positions;
      push @$counts, $count;
      $total += $count;
   }

   $counts = _normalize($counts, $total) if $freq;
   return $sources, $counts;
}


=head2 kmers

 Usage   : $col->kmers()
 Function: This is the inverse of sources(). Return the kmers found in a sequence
           (given its ID). An error is reported if the sequence ID requested does
           not exist.
 Args    : * sequence ID to get the kmers of
           * 0 to report counts (default), 1 to report frequencies (normalize to 1)
 Returns : * arrayref of sequence IDs
           * arrayref of the corresponding total counts

=cut

sub kmers {
   my ($self, $seq_id, $freq) = @_;
   my $kmers;
   my $counts;
   my $total = 0;
   my $seq_kmers = $self->collection_by_seq->{$seq_id};
  
   if (not defined $seq_kmers) {
      $self->throw("Error: Sequence $seq_id was not found in the collection.\n");
   }

   while ( my ($kmer, $positions) = each %$seq_kmers ) {
      push @$kmers, $kmer;
      my $count = scalar @$positions;
      push @$counts, $count;
      $total += $count;
   }

   $counts = _normalize($counts, $total) if $freq;
   return $kmers, $counts;
}


=head2 positions

 Usage   : $col->positions()
 Function: Return the positions of the given kmer on a given sequence. An error
           is reported if the kmer requested does not exist
 Args    : * desired kmer
           * desired sequence with this kmer
 Returns : arrayref of the different positions

=cut

sub positions {
   my ($self, $kmer, $source) = @_;
   my $kmer_sources = $self->collection_by_kmer->{$kmer};
   if (not defined $kmer_sources) {
      $self->throw("Error: kmer $kmer was not found in the collection.\n");
   }
   my $kmer_positions = $kmer_sources->{$source};
   if (not defined $kmer_sources) {
      $self->throw("Error: kmer $kmer was not found in sequence $source.\n");
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