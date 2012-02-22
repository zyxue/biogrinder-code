#! perl

use strict;
use warnings;
use Test::More;
use t::TestUtils;


use_ok 'Grinder::Database';

my ($db);

ok $db = Grinder::Database->new(
   -fasta_file              => data('shotgun_database.fa'),
);

isa_ok $db, 'Grinder::Database';


#$db = Grinder::Database->new(
#   -fasta_file              => data('shotgun_database.fa'),
#   -unidirectional          => 
#   -forward_reverse_primers =>
#   -abundance_file          =>
#   -delete_chars            =>
#   -min_len                 => 1;
#);


done_testing();
