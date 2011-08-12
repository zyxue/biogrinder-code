#!perl -T

use strict;
use warnings;
use Test::More tests => 12;
use Bio::Seq;

use Grinder;
my ($factory, $nof_reads, $read);


# Tracking read information in the read description 

ok $factory = Grinder->new(
   -reference_file => './t/data/shotgun_database.fa',
   -total_reads    => 10                            ,
   -desc_track     => 1                             ,
), 'shotgun tracking';

ok $read = $factory->next_read;
is ( $read->desc =~ m/reference=.*position=.*strand=.*/, 1);


ok $factory = Grinder->new(
   -reference_file  => './t/data/amplicon_database.fa',
   -forward_reverse => './t/data/forward_primer.fa'   ,
   -length_bias     => 0                              ,
   -unidirectional  => 1                              ,
   -total_reads     => 10                             ,
   -desc_track      => 1                              ,
), 'amplicon tracking';

ok $read = $factory->next_read;
is ( $read->desc =~ m/reference=.*amplicon=.*position=.*strand=.*/, 1);


ok $factory = Grinder->new(
   -reference_file => './t/data/shotgun_database.fa',
   -total_reads    => 10                            ,
   -desc_track     => 0                             ,
), 'no tracking';

ok $read = $factory->next_read;
is $read->desc, undef;


ok $factory = Grinder->new(
   -reference_file => './t/data/shotgun_database.fa',
   -total_reads    => 10                            ,
), 'tracking default';

ok $read = $factory->next_read;
is ( $read->desc =~ m/reference=.*position=.*strand=.*/, 1);

