#!perl -T

use strict;
use warnings;
use Test::More tests => 15;
use Bio::Seq;
use File::Spec::Functions;

use Grinder;
my ($factory, $nof_reads, $read);


# Tracking read information in the read description 

ok $factory = Grinder->new(
   -reference_file => catfile(qw{t data shotgun_database.fa}),
   -total_reads    => 10                                     ,
   -desc_track     => 1                                      ,
), 'Shotgun tracking';

ok $read = $factory->next_read;
is ( $read->desc =~ m/reference=.*position=.*strand=.*/, 1);


ok $factory = Grinder->new(
   -reference_file  => catfile(qw{t data amplicon_database.fa}),
   -forward_reverse => catfile(qw{t data forward_primer.fa})   ,
   -length_bias     => 0                                       ,
   -unidirectional  => 1                                       ,
   -total_reads     => 10                                      ,
   -desc_track      => 1                                       ,
), 'Amplicon tracking';

ok $read = $factory->next_read;
is ( $read->desc =~ m/reference=\S+.*amplicon=\d+-\d+.*position=.*strand=.*/, 1);

ok $factory = Grinder->new(
   -reference_file  => catfile(qw{t data amplicon_database.fa}),
   -forward_reverse => catfile(qw{t data forward_primer.fa})   ,
   -length_bias     => 0                                       ,
   -unidirectional  => 1                                       ,
   -total_reads     => 10                                      ,
   -desc_track      => 1                                       ,
   -chimera_perc    => 100                                     ,
), 'Chimeric amplicon tracking';

ok $read = $factory->next_read;
is ( $read->desc =~ m/reference=\S+,\S+.*amplicon=\d+-\d+,\d+-\d+.*position=.*strand=.*/, 1);


ok $factory = Grinder->new(
   -reference_file => catfile(qw{t data shotgun_database.fa}),
   -total_reads    => 10                                     ,
   -desc_track     => 0                                      ,
), 'No tracking';

ok $read = $factory->next_read;
is $read->desc, undef;


ok $factory = Grinder->new(
   -reference_file => catfile(qw{t data shotgun_database.fa}),
   -total_reads    => 10                                     ,
), 'Tracking default';

ok $read = $factory->next_read;
is ( $read->desc =~ m/reference=.*position=.*strand=.*/, 1);

