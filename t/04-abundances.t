#!perl -T

use strict;
use warnings;
use Test::More tests => 23;

use Grinder;
my ($factory, $nof_reads, $read, %sources);


# Specified genome abundance for a single library

ok $factory = Grinder->new(
   -genome_file    => './t/data/shotgun_database.fa',
   -abundance_file => './t/data/abundances.txt',
   -length_bias    => 0                             ,
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

# This tests are quite sensitive to the seed used. Ideal average answer should
# be 250 here
ok ( ($sources{'seq1'} > 230) && ($sources{'seq1'} < 280) );
ok ( ($sources{'seq2'} > 230) && ($sources{'seq2'} < 280) );
ok ( ($sources{'seq4'} > 230) && ($sources{'seq4'} < 280) );
ok ( ($sources{'seq5'} > 230) && ($sources{'seq5'} < 280) );

is $factory->next_lib, undef;
%sources = ();


# Specified genome abundance for multiple libraries

ok $factory = Grinder->new(
   -genome_file    => './t/data/shotgun_database.fa'         ,
   -abundance_file => './t/data/abundances_multiple.txt',
   -length_bias    => 0                                      ,
   -random_seed    => 1232567890                             ,
   -total_reads    => 1000                                    ), 'Genome abundance for multiple libraries';

ok $factory->next_lib;

$nof_reads = 0;
while ( $read = $factory->next_read ) { $nof_reads++ };
is $nof_reads, 1000;

ok $factory->next_lib;

ok $factory->next_lib;

while ( $read = $factory->next_read ) {
   my $source = $read->reference->id;
   if (not exists $sources{$source}) {
     $sources{$source} = 1;
   } else {
     $sources{$source}++;
   }
};

ok not exists $sources{'seq1'};
ok not exists $sources{'seq2'};
ok exists $sources{'seq3'};
ok not exists $sources{'seq4'};
ok not exists $sources{'seq5'};

is $sources{'seq3'}, 1000;

is $factory->next_lib, undef;
