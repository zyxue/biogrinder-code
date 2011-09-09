#! perl

use strict;
use warnings;
use Test::More;
use t::TestUtils;
use Grinder;

plan tests => 7024;


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

cmp_ok $nof_substs / $nof_indels, '>', 0.92; # should be 1
cmp_ok $nof_substs / $nof_indels, '<', 1.08;


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
cmp_ok $mean, '<=', 103; # 100
cmp_ok $mean, '>=', 97;
cmp_ok $stddev, '<', 12;

SKIP: {
   skip rfit_msg(), 5 if not can_rfit();
   test_uniform_dist(\@epositions, 1, 50, 'errors_uniform.txt');
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
cmp_ok $$prof[0] , '>=', 30;  # mean number of errors at 1st position of reads should be 50
cmp_ok $$prof[0] , '<=', 70;
cmp_ok $$prof[-1], '>=', 125; # mean number of errors at last position of read shoul be 150
cmp_ok $$prof[-1], '<=', 175;
cmp_ok $mean     , '<=', 103; # mean number of errors at each position should be 100
cmp_ok $mean     , '>=', 97;

SKIP: {
   skip rfit_msg(), 6 if not can_rfit();
####   test_linear_dist(\@epositions, 150, 10, 'errors_linear.txt');
}

@epositions = ();

