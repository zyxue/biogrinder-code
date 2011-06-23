#!perl -T

use strict;
use warnings;
use Test::More tests => 6;
use Bio::Seq;

use Grinder;
my ($factory, $nof_reads, $read);


# Tracking read information in the read description 

ok $factory = Grinder->new(
   -genome_file => './t/data/shotgun-database.fa',
   -read_dist   => 52                            ,
   -total_reads => 10                            ,
), 'no quality scores';

ok $read = $factory->next_read;
is join(' ',@{$read->qual}), '';


ok $factory = Grinder->new(
   -genome_file => './t/data/shotgun-database.fa',
   -read_dist   => 52                            ,
   -total_reads => 10                            ,
   -qual_levels => '30 10'                       ,
), 'with quality scores';

ok $read = $factory->next_read;
is scalar @{$read->qual}, 52;
