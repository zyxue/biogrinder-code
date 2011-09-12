#! perl

use strict;
use warnings;
use Test::More;
use t::TestUtils;
use Grinder;

plan tests => 14;


my ($factory, $nof_reads, $read, %sources);


# Specified genome abundance for a single library

ok $factory = Grinder->new(
   -reference_file => data('shotgun_database.fa'),
   -abundance_file => data('abundances.txt')     ,
   -length_bias    => 1                          ,
   -random_seed    => 1910567890                 ,
   -total_reads    => 1000                       ,
), 'Genome abundance for a single libraries';

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
between_ok( $sources{'seq1'}, 67 , 127 ); # avg = 97
between_ok( $sources{'seq2'}, 99 , 159 ); # avg = 129
between_ok( $sources{'seq4'}, 357, 417 ); # avg = 387
between_ok( $sources{'seq5'}, 357, 417 ); # avg = 387

