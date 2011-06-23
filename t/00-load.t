#!perl -T

use Test::More tests => 1;

BEGIN {
	use_ok( 'Grinder' );
}

diag( "Testing Grinder $Grinder::VERSION, Perl $], $^X" );
