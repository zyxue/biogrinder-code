
#### todo:
# case-insensitive
# revcom
# weights


package KmerCollection;

=head1 NAME

KmerCollection - A collection of kmers from sequences

=head1 SYNOPSIS

    my $col = KmerCollection->new( -k    => 10,
                                   -file => 'seqs.fa' );

=head1 DESCRIPTION

Manage a collection of kmers and from what position of what sequence they come
from.

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
 Returns : KmerCollection object
 Args    : -k       set the kmer length (default: 10 bp)
           -revcom  count kmers before and after reverse-complementing sequences
                    (default: 0)
           -seqs    count kmers in the provided arrayref of sequences (Bio::Seq
                    objects)
           -file    count kmers in the provided file of sequences

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
 Returns : Positive integer
 Args    : None

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


=head2 revcom

 Usage   : $col->revcom;
 Function: Get whether or not the kmers are calculated from the sequence and
           their reverse-complement
 Returns : 1 (yes), 0 (no)
 Args    : None

=cut

sub revcom {
   my ($self, $val) = @_;
   if ($val) {
      $self->{'revcom'} = $val;
   }
   return $self->{'revcom'};
}


=head2 collection

 Usage   : $col->collection;
 Function: Get the collection of kmers.
 Returns : A hashref of hashref of ...
 Args    : None

=cut

sub collection {
   my ($self, $val) = @_;
   if ($val) {
      $self->{'collection'} = $val;
   }
   return $self->{'collection'};
}


#==============================================================================#


=head2 add_file

 Usage   : $col->add_file( 'seqs.fa' );
 Function: Process the kmers in the given file of sequences.
 Returns : 1
 Args    : string

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
 Returns : 1
 Args    : arrayref of Bio::Seq objects

=cut

sub add_seqs {
   my ($self, $seqs) = @_;
   my $kmer_col = $self->collection || {};
   for my $seq (@$seqs) {
      my $kmer_counts = $self->_count_kmers($seq);
      while ( my ($kmer, $positions) = each %$kmer_counts ) {
         $kmer_col->{$kmer}->{$seq->id} = $positions;
      }
   }
   $self->collection($kmer_col);
   return $self;
}


=head2 filter_rare

 Usage   : $col->filter_rare( 2 );
 Function: Remove kmers occurring less than the number of times specified
 Returns : 1
 Args    : integer

=cut

sub filter_rare {
   my ($self, $min_num) = @_;
   my $collection = $self->collection;
   while ( my ($kmer, $sources) = each %$collection ) {
      my $count = _sum_from_sources( $sources );
      delete $collection->{$kmer} if $count < $min_num;
   }
   $self->collection( $collection );
   return $self;
}


=head2 filter_shared

 Usage   : $col->filter_shared( 2 );
 Function: Remove kmers occurring in less than the number of sequences specified
 Returns : 1
 Args    : integer

=cut

sub filter_shared {
   my ($self, $min_num) = @_;
   my $collection = $self->collection;
   while ( my ($kmer, $sources) = each %$collection ) {
      my $count = scalar keys %$sources;
      delete $collection->{$kmer} if $count < $min_num;
   }
   $self->collection( $collection );
   return $self;
}


=head2 counts

 Usage   : $col->counts
 Function: Calculate the total count of each kmer
 Returns : * arrayref of the different kmers
           * arrayref of the corresponding total counts
 Args    : 0 to report counts (default), 1 to report frequencies (normalize to 1)

=cut

sub counts {
   my ($self, $freq) = @_;
   my $kmers;
   my $counts;
   my $total = 0;
   my $collection = $self->collection;
   while ( my ($kmer, $sources) = each %$collection ) {
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
 Returns : * arrayref of the different sources
           * arrayref of the corresponding total counts
 Args    : * kmer to get the sources of
           * 0 to report counts (default), 1 to report frequencies (normalize to 1)

=cut

sub sources {
   my ($self, $kmer, $freq) = @_;
   my $sources;
   my $counts;
   my $total = 0;
   my $kmer_sources = $self->collection->{$kmer};
  
   if (not defined $kmer_sources) {
      $self->throw("Error: This kmer ($kmer) was not found.\n");
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


=head2 positions

 Usage   : $col->positions()
 Function: Return the positions of the given kmer on a given sequence. An error
           is reported if the kmer requested does not exist
 Returns : arrayref of the different positions
 Args    : * desired kmer
           * desired sequence with this kmer

=cut

sub positions {
   my ($self, $kmer, $source) = @_;
   my $kmer_sources = $self->collection->{$kmer};
   if (not defined $kmer_sources) {
      $self->throw("Error: This kmer ($kmer) was not found in the collection.\n");
   }
   my $kmer_positions = $kmer_sources->{$source};
   if (not defined $kmer_sources) {
      $self->throw("Error: This kmer ($kmer) was not found in sequence $source.\n");
   }
   return $kmer_positions;
}



sub _count_kmers {
   # Count the kmers of size k in a sequence (Bio::Seq) and return a hash
   my ($self, $seq) = @_;
   my $k = $self->k;
   my $seq_str = $seq->seq;
   my $seq_len = length $seq_str;
   my $hash = {};
   for (my $i = 0; $i <= $seq_len - $k ; $i++) {
      my $kmer = substr $seq_str, $i, $k;
      push @{$hash->{$kmer}}, $i + 1;
   }
   return $hash;
}


sub _sum_from_sources {
   my ($sources) = @_;
   my $count = 0;
   while ( my ($source, $positions) = each %$sources ) {
      $count += scalar @$positions;
   }
   return $count;
}


sub _normalize {
   my ($arr, $total) = @_;
   if (not $total) { # total undef or 0
      die "Error: Need to provide a valid total\n";
   }
   $arr = [ map {$_ / $total} @$arr ];
   return $arr;
}

1;
