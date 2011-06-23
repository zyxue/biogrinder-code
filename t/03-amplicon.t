#!perl -T

use strict;
use warnings;
use Test::More tests => 811;

use Grinder;
my ($factory, $read, $nof_reads);

# Forward primer only, forward sequencing

ok $factory = Grinder->new(
   -genome_file     => './t/data/amplicon_database.fa',
   -forward_reverse => './t/data/forward_primer.fa'   ,
   -length_bias     => 0                              ,
   -unidirectional  => 1                              ,
   -total_reads     => 100                             ), 'Forward primer only, forward sequencing';

ok $factory->next_lib;

ok $read = $factory->next_read;

is ref($read), 'Bio::Seq::SimulatedRead';

$nof_reads = 1;
while ( $read = $factory->next_read ) {
   $nof_reads++;
   ok_read($read, 1);
};
is $nof_reads, 100;



# Forward and reverse primers

ok $factory = Grinder->new(
   -genome_file     => './t/data/amplicon_database.fa'      ,
   -forward_reverse => './t/data/forward_reverse_primers.fa',
   -length_bias     => 0                                    ,
   -unidirectional  => 1                                    ,
   -total_reads     => 100                                   ), 'Forward then reverse primers, forward sequencing';

ok $factory->next_lib;

ok $read = $factory->next_read;

is ref($read), 'Bio::Seq::SimulatedRead';

$nof_reads = 1;
while ( $read = $factory->next_read ) {
   $nof_reads++;
   ok_read($read, 1);
};
is $nof_reads, 100;


# Reverse primer only, reverse sequencing

ok $factory = Grinder->new(
   -genome_file     => './t/data/amplicon_database.fa',
   -forward_reverse => './t/data/reverse_primer.fa'   ,
   -length_bias     => 0                              ,
   -unidirectional  => -1                             ,
   -total_reads     => 100                             ), 'Reverse primer only, reverse sequencing';

ok $factory->next_lib;

ok $read = $factory->next_read;

is ref($read), 'Bio::Seq::SimulatedRead';

$nof_reads = 1;
while ( $read = $factory->next_read ) {
   $nof_reads++;
   ok_read($read, -1);
};
is $nof_reads, 100;


# Reverse and forward primers, reverse sequencing

ok $factory = Grinder->new(
   -genome_file     => './t/data/amplicon_database.fa'      ,
   -forward_reverse => './t/data/reverse_forward_primers.fa',
   -length_bias     => 0                                    ,
   -unidirectional  => -1                                   ,
   -total_reads     => 100                                   ), 'Reverse then forward primers, reverse sequencing';

ok $factory->next_lib;

ok $read = $factory->next_read;

$nof_reads = 1;
while ( $read = $factory->next_read ) {
   $nof_reads++;
   ok_read($read, -1);
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

