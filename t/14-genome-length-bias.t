#!perl -T

use strict;
use warnings;
use Test::More tests => 10;
use Bio::Seq;

use Grinder;
my ($factory, $nof_reads, $read, %sources);


# Specified genome abundance for a single library

ok $factory = Grinder->new(
   -genome_file    => './t/data/shotgun_database.fa',
   -abundance_file => './t/data/test_abundances.txt',
   -length_bias    => 1                             ,
   -random_seed    => 1910567890                    ,
   -total_reads    => 1000                           ), 'Genome abundance for a single libraries';

while ( $read = $factory->next_read ) {
   my $source = $read->reference->id;
   if (not exists $sources{$source}) {
     $sources{$source} = 1;
   } else {
     $sources{$source}++;
   }
};

ok exists $sources{'seq1'};
ok exists $sources{'seq2'};
ok not exists $sources{'seq3'};
ok exists $sources{'seq4'};
ok exists $sources{'seq5'};


# These tests are quite sensitive to the seed used
ok ( ($sources{'seq1'} > 67 ) && ($sources{'seq1'} < 127) ); # avg = 97
ok ( ($sources{'seq2'} > 99 ) && ($sources{'seq2'} < 159) ); # avg = 129
ok ( ($sources{'seq4'} > 357) && ($sources{'seq4'} < 417) ); # avg = 387
ok ( ($sources{'seq5'} > 357) && ($sources{'seq5'} < 417) ); # avg = 387

