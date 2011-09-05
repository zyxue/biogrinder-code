#! perl

use strict;
use warnings;
use Test::More;
use t::TestUtils;
use Grinder;

plan tests => 19;


my ($factory, $nof_reads, $read, @reads, $ra, $era, $coeff, $min, $max, $mean,
    $stddev, $struct, $param1, $param2);


# Uniform community structure

ok $factory = Grinder->new(
   -reference_file  => data('shotgun_database.fa'),
   -read_dist       => 48                         ,
   -length_bias     => 0                          ,
   -abundance_model => ('uniform', 0)             ,
   -total_reads     => 1000                       ,
), 'Uniform community structure';

while ( $read = $factory->next_read ) {
   push @reads, $read->reference->id;
}

$ra = rank_abundance(\@reads, 10);
($min, $max, $mean, $stddev) = stats($ra);
$era = uniform(10, 5, 1000);
$coeff = corr_coeff($ra, $era, $mean);
cmp_ok $coeff, '>', 0.97;

@reads = ();


# Linear community structure

ok $factory = Grinder->new(
   -reference_file  => data('shotgun_database.fa'),
   -read_dist       => 48                         ,
   -length_bias     => 0                          ,
   -abundance_model => ('linear', 0)              ,
   -total_reads     => 1000                       ,
), 'Linear community structure';

while ( $read = $factory->next_read ) {
   push @reads, $read->reference->id;
}

$ra = rank_abundance(\@reads, 10);
($min, $max, $mean, $stddev) = stats($ra);
$era = linear(10, 5, 1000);
$coeff = corr_coeff($ra, $era, $mean);
cmp_ok $coeff, '>', 0.97;

@reads = ();


# Power law community structure

ok $factory = Grinder->new(
   -reference_file  => data('shotgun_database.fa'),
   -read_dist       => 48                         ,
   -length_bias     => 0                          ,
   -abundance_model => ('powerlaw', 0.5)          ,
   -total_reads     => 1000                       ,
), 'Power law community structure';

while ( $read = $factory->next_read ) {
   push @reads, $read->reference->id;
}

$ra = rank_abundance(\@reads, 10);
($min, $max, $mean, $stddev) = stats($ra);
$era = powerlaw(10, 5, 0.5, 1000);
$coeff = corr_coeff($ra, $era, $mean);
cmp_ok $coeff, '>', 0.97;

@reads = ();


# Logarithmic community structure

ok $factory = Grinder->new(
   -reference_file  => data('shotgun_database.fa'),
   -read_dist       => 48                         ,
   -length_bias     => 0                          ,
   -abundance_model => ('logarithmic', 0.5)       ,
   -total_reads     => 1000                       ,
), 'Logarithmic community structure';

while ( $read = $factory->next_read ) {
   push @reads, $read->reference->id;
}

$ra = rank_abundance(\@reads, 10);
($min, $max, $mean, $stddev) = stats($ra);
$era = logarithmic(10, 5, 0.5, 1000);
$coeff = corr_coeff($ra, $era, $mean);
cmp_ok $coeff, '>', 0.97;

@reads = ();


# Exponential community structure

ok $factory = Grinder->new(
   -reference_file  => data('shotgun_database.fa'),
   -read_dist       => 48                         ,
   -length_bias     => 0                          ,
   -abundance_model => ('exponential', 0.5)       ,
   -total_reads     => 1000                       ,
), 'Exponential community structure';

$struct = $factory->next_lib;
while ( $read = $factory->next_read ) {
   push @reads, $read->reference->id;
}

$ra = rank_abundance(\@reads, 10);
($min, $max, $mean, $stddev) = stats($ra);
$era = exponential(10, 5, 0.5, 1000);
$coeff = corr_coeff($ra, $era, $mean);
cmp_ok $coeff, '>', 0.97;
is $struct->{param}, 0.5;
@reads = ();


# Communities with random structure parameter value

ok $factory = Grinder->new(
   -reference_file  => data('shotgun_database.fa'),
   -read_dist       => 48                         ,
   -length_bias     => 0                          ,
   -num_libraries   => 2                          ,
   -shared_perc     => 100                        ,
   -abundance_model => ('exponential')            ,
   -total_reads     => 1000                       ,
), 'Communities with random structure parameter value';

$struct = $factory->next_lib;
while ( $read = $factory->next_read ) {
   push @reads, $read->reference->id;
}

$ra = rank_abundance(\@reads, 10);
($min, $max, $mean, $stddev) = stats($ra);
$param1 = $struct->{param};
cmp_ok $param1, '>', 0;
cmp_ok $param1, '<', 1000;
$era = exponential(10, 5, $param1, 1000);
$coeff = corr_coeff($ra, $era, $mean);
cmp_ok $coeff, '>', 0.97;

@reads = ();

$struct = $factory->next_lib;
while ( $read = $factory->next_read ) {
   push @reads, $read->reference->id;
}

$ra = rank_abundance(\@reads, 10);
($min, $max, $mean, $stddev) = stats($ra);
$param2 = $struct->{param};
cmp_ok $param2, '>', 0;
cmp_ok $param2, '<', 1000;
$era = exponential(10, 5, $param2, 1000);
$coeff = corr_coeff($ra, $era, $mean);
cmp_ok $coeff, '>', 0.97;

isnt $param1, $param2;

@reads = ();



sub uniform {
   # Evaluate the uniform function in the given integer range
   my ($x_max, $max, $num) = @_;
   my @ys;
   my $width = $max;
   for my $x (1 .. $x_max) {
      my $y;
      if ( $x <= $max ) {
         $y = $num / $width;
      } else {
         $y = 0;
      }
      push @ys, $y;
   }
   return \@ys;
}


sub linear {
   # Evaluate the linear function in the given integer range
   my ($x_max, $max, $num) = @_;
   my @ys;
   my $sum = 0;
   for (my $x = $max; $x >= 1; $x--) {
      my $y = $x;
      $sum += $y;
      push @ys, $y;
   }
   for (my $x = 0; $x < $max; $x++) {
      $ys[$x] *= $num / $sum;
   }
   push @ys, (0) x ($x_max - scalar @ys);
   return \@ys;
}


sub powerlaw {
   # Evaluate the power function in the given integer range
   my ($x_max, $max, $param, $num) = @_;
   my @ys;
   my $sum = 0;
   for my $x (1 .. $max) {
      my $y = $x**(-$param);
      $sum += $y;
      push @ys, $y;
   }
   for (my $x = 0; $x < $max; $x++) {
      $ys[$x] *= $num / $sum;
   }
   push @ys, (0) x ($x_max - scalar @ys);
   return \@ys;
}


sub logarithmic {
   # Evaluate the logarithmic function in the given integer range
   my ($x_max, $max, $param, $num) = @_;
   my @ys;
   my $sum = 0;
   for my $x (1 .. $max) {
      my $y = (log($x+1))**(-$param);
      $sum += $y;
      push @ys, $y;
   }
   for (my $x = 0; $x < $max; $x++) {
      $ys[$x] *= $num / $sum;
   }
   push @ys, (0) x ($x_max - scalar @ys);
   return \@ys;
}


sub exponential {
   # Evaluate the exponential function in the given integer range
   my ($x_max, $max, $param, $num) = @_;
   my @ys;
   my $sum = 0;
   for my $x (1 .. $max) {
      my $y = exp(-$x*$param);
      $sum += $y;
      push @ys, $y;
   }
   for (my $x = 0; $x < $max; $x++) {
      $ys[$x] *= $num / $sum;
   }
   push @ys, (0) x ($x_max - scalar @ys);
   return \@ys;
}


sub rank_abundance {
   my ($data, $max) = @_;
   # Put a data series into bins
   my %hash;
   for my $val (@$data) {
      $hash{$val}++;
   }
   my @y_data = sort { $b <=> $a } (values %hash);
   push @y_data, (0) x ($max - scalar @y_data);
   return \@y_data;
}


