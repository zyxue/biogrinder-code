#! perl

use strict;
use warnings;
use Test::More;


BEGIN {
   use_ok( 'Grinder' );
}

diag( "Testing Grinder $Grinder::VERSION, Perl $], $^X" );

use_ok('t::TestUtils');

done_testing();
