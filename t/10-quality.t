#! perl

use strict;
use warnings;
use Test::More;
use t::TestUtils;
use Grinder;

plan tests => 6;


my ($factory, $nof_reads, $read);


# Outputing basic quality scores

ok $factory = Grinder->new(
   -genome_file => data('shotgun_database.fa'),
   -read_dist   => 52                         ,
   -total_reads => 10                         ,
), 'No quality scores';

ok $read = $factory->next_read;
is join(' ',@{$read->qual}), '';


ok $factory = Grinder->new(
   -genome_file => data('shotgun_database.fa'),
   -read_dist   => 52                         ,
   -total_reads => 10                         ,
   -qual_levels => '30 10'                    ,
), 'With quality scores';

ok $read = $factory->next_read;
is scalar @{$read->qual}, 52;
