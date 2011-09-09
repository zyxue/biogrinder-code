#! perl

use strict;
use warnings;
use Test::More;
use t::TestUtils;
use Grinder;

plan tests => 7018;


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

#$eprof = uniform(1, 50, 1, 50, 5000); # 1000 reads * 50 bp * 10% error
#$coeff = corr_coeff($prof, $eprof, $mean);
#cmp_ok $coeff, '>', 0.99;

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

#while ( $read = $factory->next_read ) {
#   my @positions = error_positions($read);
#   push @epositions, @positions if scalar @positions > 0;
#}
#
#use Data::Dumper;
#print Dumper($epositions);
#
#$prof = error_profile($epositions, 50);
#($min, $max, $mean, $stddev) = stats($prof);
##print "min = $min, max = $max, mean = $mean, stddev = $stddev\n";
#cmp_ok $mean, '<=', 103; # should be 100
#cmp_ok $mean, '>=', 97;
#cmp_ok $min,  '>=', 30;  # should be 50
#cmp_ok $min,  '<=', 70;
#cmp_ok $max,  '>=', 125; # should be 150
#cmp_ok $max,  '<=', 175;
#$eprof = linear(50, 10, 15, 1000);
#$coeff = corr_coeff($prof, $eprof, $mean);
##print "coeff= $coeff\n";
#cmp_ok $coeff, '>', 0.80;
##cmp_ok $coeff, '>', 0.99;
#@epositions = ();




sub error_positions {
   my ($read) = @_;
   my ($err_str) = ($read->desc =~ /errors=(\S+)/);
   my @error_positions;
   if (defined $err_str) {
      for my $error (split ',', $err_str) {
         my ($pos, $type, $res) = ($error =~ m/(\d+)([%+-])([a-z]*)/i);
         push @error_positions, $pos;
      }
   }
   return @error_positions;
}


sub linear {
   # Evaluate the linear function in the given integer range
   my ($x_max, $mean, $right, $num) = @_;
   my $height = ($right - $mean) * 2;
   my $width  = $x_max - 1;
   my $slope = $height / $width;   
   my $left = $mean - ($right - $mean);
   my @ys;
   for my $x (1 .. $x_max) {
      my $y = $left + ($x-1) * $slope;
      push @ys, $y;
   }
   for (my $x = 0; $x < $x_max; $x++) {
      $ys[$x] *= $num / 100;
   }
   return \@ys;
}


