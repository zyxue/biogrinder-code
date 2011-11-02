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
   -mutation_dist  => (10, 'uniform')               ,
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
   -mutation_dist  => (10, 'uniform')               ,
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
   -mutation_dist  => (10, 'uniform')               ,
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


# Uniform distribution

ok $factory = Grinder->new(
   -reference_file => data('single_seq_database.fa'),
   -unidirectional => 1                             ,
   -read_dist      => 50                            ,
   -total_reads    => 1000                          ,
   -mutation_ratio => (50, 50)                      ,
   -mutation_dist  => (10, 'uniform')               ,
), 'Uniform';

while ( $read = $factory->next_read ) {
   my @positions = error_positions($read);
   push @epositions, @positions if scalar @positions > 0;
}

$prof = hist(\@epositions, 1, 50);
($min, $max, $mean, $stddev) = stats($prof);
cmp_ok $min, '>=', 65;
cmp_ok $max, '<=', 135;
between_ok( $mean, 97, 103 ); # should be 100
cmp_ok $stddev, '<', 12;

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
   -mutation_dist  => (10, 'linear', 15)            ,
), 'Linear';

while ( $read = $factory->next_read ) {
   my @positions = error_positions($read);
   push @epositions, @positions if scalar @positions > 0;
}

$prof = hist(\@epositions, 1, 50);
($min, $max, $mean, $stddev) = stats($prof);
between_ok( $$prof[0] ,  30,  70 ); # exp. number of errors at 1st  pos is 50
between_ok( $$prof[-1], 125, 175 ); # exp. number of errors at last pos is 150
between_ok( $mean     ,  97, 103 ); # exp. mean number of errors at each pos is 100

SKIP: {
   skip rfit_msg() if not can_rfit();

   #### TODO
   #TODO: {
   #   $TODO = "Need to implement a linear density distribution in R";
   #   test_linear_dist(\@epositions, 1, 50, 0.0000000001);
   #}
   ####
}

@epositions = ();


# Fourth degree polynomial distribution

ok $factory = Grinder->new(
   -reference_file => data('single_seq_database.fa'),
   -unidirectional => 1                             ,
   -read_dist      => 100                           ,
   -total_reads    => 1000                     ,
   -mutation_ratio => (50, 50)                      ,
   -mutation_dist  => (1, 'poly4', 4.4e-7)         ,
), 'Polynomial';

while ( $read = $factory->next_read ) {
   my @positions = error_positions($read);
   push @epositions, @positions if scalar @positions > 0;
}

$prof = hist(\@epositions, 1, 100);
($min, $max, $mean, $stddev) = stats($prof);

use Data::Dumper; print "prof: ".Dumper($prof);

between_ok( $$prof[0] ,    1,   27 ); # mean number at 1st  pos should be 10 (1%)
between_ok( $$prof[49],   16,   46 ); # mean number at 50th pos should be 37.4 (3.74%)
between_ok( $$prof[-1],  419,  479 ); # mean number at last pos should be 449 (44.9%)
between_ok( $mean     ,   80,  120 ); # mean number at each pos should be 100 (10.02%)

SKIP: {
   skip rfit_msg() if not can_rfit();

   #### TODO
   #TODO: {
   #   $TODO = "Need to implement a linear density distribution in R";
   #   test_linear_dist(\@epositions, 1, 50, 0.0000000001);
   #}
   ####
}

@epositions = ();

done_testing();

