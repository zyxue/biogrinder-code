#!perl -T

use strict;
use warnings;
use Test::More tests => 9;
use Bio::Seq;

use Grinder;
my ($factory, $nof_reads, $read);


# Tracking read information in the read description 

ok $factory = Grinder->new(
   -genome_file => './t/data/shotgun_database.fa',
   -total_reads => 10                            ,
   -desc_track  => 1                             ,
), 'tracking on';

ok $read = $factory->next_read;
is ( $read->desc =~ m/reference=.*position=.*strand=.*/, 1);


ok $factory = Grinder->new(
   -genome_file => './t/data/shotgun_database.fa',
   -total_reads => 10                             ,
   -desc_track  => 0                             ,
), 'tracking off';

ok $read = $factory->next_read;
is $read->desc, undef;


ok $factory = Grinder->new(
   -genome_file => './t/data/shotgun_database.fa',
   -total_reads => 10                            ,
), 'tracking default';

ok $read = $factory->next_read;
is ( $read->desc =~ m/reference=.*position=.*strand=.*/, 1);

