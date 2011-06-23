#!perl -T

use strict;
use warnings;
use Test::More tests => 107;
use Bio::Seq;

use Grinder;
my ($factory, $nof_reads, $read);


# Initialization with short argument

ok $factory = Grinder->new(
   -gf => './t/data/shotgun-database.fa',
   -tr => 10                             ), 'Shotgun & short arguments';

ok $factory->next_read;


# Long argument

ok $factory = Grinder->new(
   -genome_file => './t/data/shotgun-database.fa',
   -read_dist   => 52                            ,
   -total_reads => 100                            ), 'Long arguments';

ok $factory->next_lib;

ok $read = $factory->next_read;

is ref($read), 'Bio::Seq::SimulatedRead';
is $read->length, 52;

$nof_reads = 1;
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

