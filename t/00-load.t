#! perl

use strict;
use warnings;
use Test::More tests => 2;


BEGIN {
   use_ok( 'Grinder' );
}

diag( "Testing Grinder $Grinder::VERSION, Perl $], $^X" );

use_ok('t::TestUtils');
