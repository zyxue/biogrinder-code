#! perl

use strict;
use warnings;
use Test::More;
use t::TestUtils;
use Grinder;



my ($factory, $nof_reads, $read);


# Prepend a single multiplex identifier (MID): ACGT

ok $factory = Grinder->new(
   -reference_file => data('shotgun_database.fa'),
   -multiplex_ids  => data('mids.fa')            ,
   -num_libraries  => 2                          ,
   -read_dist      => 52                         ,
   -total_reads    => 9                          ,
), 'Single MID';

while ( $read = $factory->next_read ) {
   is $read->length, 52;
   is substr($read->seq, 0, 4), 'ACGT';
};


# Prepend two multiplex identifier: ACGT and AAAATTTT

ok $factory = Grinder->new(
   -reference_file => data('shotgun_database.fa'),
   -multiplex_ids  => data('mids.fa')            ,
   -num_libraries  => 2                          ,
   -read_dist      => 52                         ,
   -total_reads    => 10                         ,
), 'Two MIDs';

while ( $read = $factory->next_read ) {
   is $read->length, 52;
   like $read->id, qr/^1_/;
   is substr($read->seq, 0, 4), 'ACGT';
};

$factory->next_lib;

while ( $read = $factory->next_read ) {
   like $read->id, qr/^2_/;
   is $read->length, 52;
   is substr($read->seq, 0, 8), 'AAAATTTT';
};

done_testing();
