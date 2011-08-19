#! perl

use strict;
use warnings;
use Test::More;
use t::TestUtils;
use Grinder;

plan tests => 15;


my ($factory, $nof_reads, $read);


# Tracking read information in the read description 

ok $factory = Grinder->new(
   -reference_file => data('shotgun_database.fa'),
   -total_reads    => 10                         ,
   -desc_track     => 1                          ,
), 'Shotgun tracking';

ok $read = $factory->next_read;
is ( $read->desc =~ m/reference=.*position=.*strand=.*/, 1);


ok $factory = Grinder->new(
   -reference_file  => data('amplicon_database.fa'),
   -forward_reverse => data('forward_primer.fa')   ,
   -length_bias     => 0                           ,
   -unidirectional  => 1                           ,
   -total_reads     => 10                          ,
   -desc_track      => 1                           ,
), 'Amplicon tracking';

ok $read = $factory->next_read;
is ( $read->desc =~ m/reference=\S+.*amplicon=\d+-\d+.*position=.*strand=.*/, 1);

ok $factory = Grinder->new(
   -reference_file  => data('amplicon_database.fa'),
   -forward_reverse => data('forward_primer.fa')   ,
   -length_bias     => 0                           ,
   -unidirectional  => 1                           ,
   -total_reads     => 10                          ,
   -desc_track      => 1                           ,
   -chimera_perc    => 100                         ,
), 'Chimeric amplicon tracking';

ok $read = $factory->next_read;
is ( $read->desc =~ m/reference=\S+,\S+.*amplicon=\d+-\d+,\d+-\d+.*position=.*strand=.*/, 1);


ok $factory = Grinder->new(
   -reference_file => data('shotgun_database.fa'),
   -total_reads    => 10                         ,
   -desc_track     => 0                          ,
), 'No tracking';

ok $read = $factory->next_read;
is $read->desc, undef;


ok $factory = Grinder->new(
   -reference_file => data('shotgun_database.fa'),
   -total_reads    => 10                         ,
), 'Tracking default';

ok $read = $factory->next_read;
is ( $read->desc =~ m/reference=.*position=.*strand=.*/, 1);

