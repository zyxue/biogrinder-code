#!perl -T

use strict;
use warnings;
use Test::More tests => 105;

use Grinder;
my ($factory, $mate1, $mate2, $read, $nof_reads);

ok $factory = Grinder->new(
   -genome_file => './t/data/shotgun_database.fa',
   -total_reads => 100                           ,
   -insert_dist => 250                            ), 'Mate pairs';

ok $factory->next_lib;

ok $mate1 = $factory->next_read;
ok $mate2 = $factory->next_read;

is ref($mate1), 'Bio::Seq::SimulatedRead';
is ref($mate2), 'Bio::Seq::SimulatedRead';

$nof_reads = 2;
while ( $read = $factory->next_read ) {
   $nof_reads++;
   ok_read($read);
};
is $nof_reads, 100;



sub ok_read {
   my ($read, $req_strand) = @_;
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
}

