#! perl

use strict;
use warnings;
use Test::More;
use t::TestUtils;
use Grinder;



my ($factory, $nof_reads, $read, %sources);


# Specified genome abundance for a single library

ok $factory = Grinder->new(
   -abundance_file  => data('abundances2.txt')              ,
   -reference_file  => data('multiple_amplicon_database.fa'),
   -forward_reverse => data('forward_reverse_primers.fa')   ,
   -copy_bias       => 1                                    ,
   -unidirectional  => 1                                    ,
   -read_dist       => 48                                   ,
   -random_seed     => 1910567890                           ,
   -total_reads     => 1000                                 ,
), 'Genome abundance for a single libraries';

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

ok exists $sources{'seq1'};
ok exists $sources{'seq2'};
ok exists $sources{'seq3'};


# These tests are quite sensitive to the seed used. Ideal average answer should
# be 180.45, 812.03 and 7.5188
between_ok( $sources{'seq1'}, 160, 200 );
between_ok( $sources{'seq2'}, 792, 832 );
between_ok( $sources{'seq3'},   0,  20 );

is $factory->next_lib, undef;
%sources = ();

done_testing();
