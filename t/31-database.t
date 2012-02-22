#! perl

use strict;
use warnings;
use Test::More;
use t::TestUtils;


use_ok 'Grinder::Database';

my ($db);

ok $db = Grinder::Database->new( );

isa_ok $db, 'Grinder::Database';



done_testing();
