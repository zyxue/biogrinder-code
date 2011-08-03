#! /usr/bin/env perl

use strict;
use warnings;
use Test::More tests => 8;

use Grinder;
my ($factory, $nof_reads, $read, %sources);


# Specified genome abundance for a single library

ok $factory = Grinder->new(
   -abundance_file  => './t/data/abundances2.txt'              ,
   -genome_file     => './t/data/multiple_amplicon_database.fa',
   -forward_reverse => './t/data/forward_reverse_primers.fa'   ,
   -copy_bias       => 1                                       ,
   -unidirectional  => 1                                       ,
   -read_dist       => 48                                      ,
   -random_seed     => 1910567890                              ,
   -total_reads     => 1000                                    ), 'Genome abundance for a single libraries';

while ( $read = $factory->next_read ) {
   my $source = $read->reference->id;
   # Strip amplicon sources of the 'amplicon' part
   $source =~ s/_amplicon.*$//;
   if (not exists $sources{$source}) {
     $sources{$source} = 1;
   } else {
     $sources{$source}++;
   }
};

use Data::Dumper; print Dumper(\%sources);

ok exists $sources{'seq1'};
ok exists $sources{'seq2'};
ok exists $sources{'seq3'};


# These tests are quite sensitive to the seed used. Ideal average answer should
# be 180.45, 812.03 and 7.5188
ok ( ($sources{'seq1'} > 160) && ($sources{'seq1'} < 200) );
ok ( ($sources{'seq2'} > 792) && ($sources{'seq2'} < 832) );
ok ( ($sources{'seq3'} >   0) && ($sources{'seq3'} < 20 ) );

is $factory->next_lib, undef;
%sources = ();
