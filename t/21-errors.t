#! perl

use strict;
use warnings;
use Test::More;
use t::TestUtils;
use Grinder;



my ($factory, $nof_reads, $read, @epositions, $min, $max, $mean, $stddev, $prof,
    $eprof, $coeff, $nof_indels, $nof_substs);


# No errors by default

ok $factory = Grinder->new(
   -reference_file => data('single_seq_database.fa'),
   -unidirectional => 1                             ,
   -read_dist      => 50                            ,
   -total_reads    => 1000                          ,
), 'No errors';

while ( $read = $factory->next_read ) {
   is $read->seq, 'aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa';
   unlike $read->desc, qr/errors/;
}


# Substitutions

ok $factory = Grinder->new(
   -reference_file => data('single_seq_database.fa'),
   -unidirectional => 1                             ,
   -read_dist      => 50                            ,
   -total_reads    => 1000                          ,
   -mutation_ratio => (100, 0)                      ,
   -mutation_dist  => ('uniform', 10)               ,
), 'Substitutions only';

while ( $read = $factory->next_read ) {
   my ($error_str) = ($read->desc =~ /errors=(\S+)/);
   if ($error_str) {
      like   $error_str, qr/%/;
      unlike $error_str, qr/[-+]/;
   } else {
      ok 1;
      ok 1;
   }
}


# Indels

ok $factory = Grinder->new(
   -reference_file => data('single_seq_database.fa'),
   -unidirectional => 1                             ,
   -read_dist      => 50                            ,
   -total_reads    => 1000                          ,
   -mutation_ratio => (0, 100)                      ,
   -mutation_dist  => ('uniform', 10)               ,
), 'Indels only';

while ( $read = $factory->next_read ) {
   my ($error_str) = ($read->desc =~ /errors=(\S+)/);
   if ($error_str) {
      unlike $error_str, qr/%/;
      like   $error_str, qr/[-+]/;
   } else {
      ok 1;
      ok 1;
   }
}


# Indels and substitutions

ok $factory = Grinder->new(
   -reference_file => data('single_seq_database.fa'),
   -unidirectional => 1                             ,
   -read_dist      => 50                            ,
   -total_reads    => 1000                          ,
   -mutation_ratio => (50, 50)                      ,
   -mutation_dist  => ('uniform', 10)               ,
), 'Indels and substitutions';

while ( $read = $factory->next_read ) {
   my ($error_str) = ($read->desc =~ /errors=(\S+)/);
   if ($error_str) {
      like $error_str, qr/[-+%]/;
      $nof_indels += ($error_str =~ tr/-+//);
      $nof_substs += ($error_str =~ tr/%//);
   } else {
      ok 1;
   }
}
between_ok( $nof_substs / $nof_indels, 0.92, 1.08 ); # should be 1


# Uniform distribution (frequent errors)

ok $factory = Grinder->new(
   -reference_file => data('single_seq_database.fa'),
   -unidirectional => 1                             ,
   -read_dist      => 50                            ,
   -total_reads    => 1000                          ,
   -mutation_ratio => (50, 50)                      ,
   -mutation_dist  => ('uniform', 10)               ,
   -random_seed    => 1233567880                    ,
), 'Uniform (frequent errors)';

while ( $read = $factory->next_read ) {
   my @positions = error_positions($read);
   push @epositions, @positions if scalar @positions > 0;
}

$prof = hist(\@epositions, 1, 50);
($min, $max, $mean, $stddev) = stats($prof);
between_ok( $$prof[0] , 70, 130 ); # exp. number of errors at 1st  pos is 100 (10%)
between_ok( $$prof[24], 70, 130 ); # exp. number of errors at 25th pos is 100 (10%)
between_ok( $$prof[-1], 70, 130 ); # exp. number of errors at last pos is 100 (10%)
between_ok( $mean     , 97, 103 ); # exp. mean number is 100 (10%)

SKIP: {
   skip rfit_msg() if not can_rfit();
   test_uniform_dist(\@epositions, 1, 50);
}

@epositions = ();


# Uniform distribution (rare errors)

ok $factory = Grinder->new(
   -reference_file => data('single_seq_database.fa'),
   -unidirectional => 1                             ,
   -read_dist      => 50                            ,
   -total_reads    => 10000                         ,
   -mutation_ratio => (50, 50)                      ,
   -mutation_dist  => ('uniform', 0.1)              ,
   -random_seed    => 1233567880                    ,
), 'Uniform (rare errors)';

while ( $read = $factory->next_read ) {
   my @positions = error_positions($read);
   push @epositions, @positions if scalar @positions > 0;
}

$prof = hist(\@epositions, 1, 50);
($min, $max, $mean, $stddev) = stats($prof);
between_ok( $$prof[0] ,  4,  16 ); # exp. number of errors at 1st  pos is 100 (10%)
between_ok( $$prof[24],  4,  16 ); # exp. number of errors at 25th pos is 100 (10%)
between_ok( $$prof[-1],  4,  16 ); # exp. number of errors at last pos is 100 (10%)
between_ok( $mean     ,  8,  12 ); # exp. mean number is 100 (10%)

SKIP: {
   skip rfit_msg() if not can_rfit();
   test_uniform_dist(\@epositions, 1, 50);
}

@epositions = ();


# Linear distribution

ok $factory = Grinder->new(
   -reference_file => data('single_seq_database.fa'),
   -unidirectional => 1                             ,
   -read_dist      => 50                            ,
   -total_reads    => 1000                          ,
   -mutation_ratio => (50, 50)                      ,
   -mutation_dist  => ('linear', 5, 15)             ,
   -random_seed    => 1233567880                    ,
), 'Linear';

while ( $read = $factory->next_read ) {
   my @positions = error_positions($read);
   push @epositions, @positions if scalar @positions > 0;
}

$prof = hist(\@epositions, 1, 50);
($min, $max, $mean, $stddev) = stats($prof);
between_ok( $$prof[0] ,  30,   70 ); # exp. number of errors at 1st  pos is 50 (5%)
between_ok( $$prof[24],  65,  135 ); # exp. number of errors at 25th pos is 100 (10%)
between_ok( $$prof[-1], 115,  185 ); # exp. number of errors at last pos is 150 (15%)
between_ok( $mean     ,  97,  103 ); # exp. mean number of errors is 100

#SKIP: {
   #skip rfit_msg() if not can_rfit();
   #### TODO
   #TODO: {
   #   $TODO = "Need to implement a linear density distribution in R";
   #   test_linear_dist(\@epositions, 1, 50, 0.0000000001);
   #}
#}

@epositions = ();


# Fourth degree polynomial distribution

ok $factory = Grinder->new(
   -reference_file => data('single_seq_database.fa'),
   -unidirectional => 1                             ,
   -read_dist      => 100                           ,
   -total_reads    => 1000                          ,
   -mutation_ratio => (50, 50)                      ,
   -mutation_dist  => ('poly4', 1, 45)              ,
   -random_seed    => 1233567880                    ,
), 'Polynomial';

while ( $read = $factory->next_read ) {
   my @positions = error_positions($read);
   push @epositions, @positions if scalar @positions > 0;
}

$prof = hist(\@epositions, 1, 100);
($min, $max, $mean, $stddev) = stats($prof);
between_ok( $$prof[0] ,    1,   27 ); # exp. number of errors at 1st  is 10 (1%)
between_ok( $$prof[49],    7,   67 ); # exp. number of errors at 50th is 37.5 (3.75%)
between_ok( $$prof[-1],  405,  495 ); # exp. number of errors at last is 450 (45.0%)
between_ok( $mean     ,   97,  103 ); # exp. mean number of errors is 100 (10.02%)

#SKIP: {
   #skip rfit_msg() if not can_rfit();
   #### TODO
   #TODO: {
   #   $TODO = "Need to implement a polynomial distribution in R";
   #   test_polynomial_dist(\@epositions, 1, 50, 0.0000000001);
   #}
#}

@epositions = ();

done_testing();

