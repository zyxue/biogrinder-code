#! perl

use strict;
use warnings;
use Test::More;
use t::TestUtils;
use Grinder;

plan tests => 403;


my ($factory, $read, $nof_reads);

ok $factory = Grinder->new(
   -reference_file => data('shotgun_database.fa'),
   -total_reads    => 100                        ,
   -read_dist      => 48                         ,
   -insert_dist    => 250                        ,
), 'Mate pairs';

ok $factory->next_lib;

$nof_reads = 0;
while ( $read = $factory->next_read ) {
   $nof_reads++;
   ok_mate($read, undef, $nof_reads);
};
is $nof_reads, 100;



sub ok_mate {
   my ($read, $req_strand, $nof_reads) = @_;
   is ref($read), 'Bio::Seq::SimulatedRead';
   my $source = $read->reference->id;
   my $strand = $read->strand;
   if (not defined $req_strand) {
      $req_strand = $strand;
   } else {
      is $strand, $req_strand;
   }
   my $letters;
   if ( $source eq 'seq1' ) {
      $letters = 'a';
   } elsif ( $source eq 'seq2' ) {
      $letters = 'c';
   } elsif ( $source eq 'seq3' ) {
      $letters = 'g';
   } elsif ( $source eq 'seq4' ) {
      $letters = 't';
   } elsif ( $source eq 'seq5' ) {
      $letters = 'atg';
   }
   if ( $req_strand == -1 ) { # Take the reverse complement
      $letters = Bio::PrimarySeq->new( -seq => $letters )->revcom->seq;
   };
   ok $read->seq =~ m/[$letters]+/;
   my $id = int($nof_reads/2+0.5).'/'.($nof_reads%2?1:2);
   is $read->id, $id;
   is $read->length, 48;
}

