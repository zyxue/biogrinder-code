#!perl 

use strict;
use warnings;
use Test::More tests => 6;

use Grinder;

my ($factory, $seed);


# Seed the pseudo-random number generator

ok $factory = Grinder->new(
   -genome_file    => './t/data/shotgun_database.fa' ,
   -random_seed    => 1233567890                     ,
   -total_reads    => 10                              ), 'Set the seed';
ok $seed = $factory->get_random_seed();
is $seed, 1233567890;

ok $factory = Grinder->new(
   -genome_file    => './t/data/shotgun_database.fa',
   -total_reads    => 10                             ), 'Get a seed automatically';
ok $seed = $factory->get_random_seed();
is( ($seed > 0), 1);


