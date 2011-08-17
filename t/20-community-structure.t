#!perl -T

use strict;
use warnings;
use Test::More tests => 19;
use Bio::Seq;
use File::Spec::Functions;

use Grinder;
my ($factory, $nof_reads, $read, @reads, $ra, $era, $coeff, $min, $max, $mean,
    $stddev, $struct, $param1, $param2);


# Uniform community structure

ok $factory = Grinder->new(
   -genome_file     => catfile(qw{t data shotgun_database.fa}),
   -read_dist       => 48                                     ,
   -length_bias     => 0                                      ,
   -abundance_model => ('uniform', 0)                         ,
   -total_reads     => 1000                                   ,
), 'Uniform community structure';

while ( $read = $factory->next_read ) {
   push @reads, $read->reference->id;
}

$ra = rank_abundance(\@reads, 10);
($min, $max, $mean, $stddev) = stats($ra);
$era = uniform(10, 5, 1000);
$coeff = corr_coeff($ra, $era, $mean);
ok ($coeff > 0.97);

@reads = ();


# Linear community structure

ok $factory = Grinder->new(
   -genome_file     => catfile(qw{t data shotgun_database.fa}),
   -read_dist       => 48                                     ,
   -length_bias     => 0                                      ,
   -abundance_model => ('linear', 0)                          ,
   -total_reads     => 1000                                   ,
), 'Linear community structure';

while ( $read = $factory->next_read ) {
   push @reads, $read->reference->id;
}

$ra = rank_abundance(\@reads, 10);
($min, $max, $mean, $stddev) = stats($ra);
$era = linear(10, 5, 1000);
$coeff = corr_coeff($ra, $era, $mean);
ok ($coeff > 0.97);

@reads = ();


# Power law community structure

ok $factory = Grinder->new(
   -genome_file     => catfile(qw{t data shotgun_database.fa}),
   -read_dist       => 48                                     ,
   -length_bias     => 0                                      ,
   -abundance_model => ('powerlaw', 0.5)                      ,
   -total_reads     => 1000                                   ,
), 'Power law community structure';

while ( $read = $factory->next_read ) {
   push @reads, $read->reference->id;
}

$ra = rank_abundance(\@reads, 10);
($min, $max, $mean, $stddev) = stats($ra);
$era = powerlaw(10, 5, 0.5, 1000);
$coeff = corr_coeff($ra, $era, $mean);
ok ($coeff > 0.97);

@reads = ();


# Logarithmic community structure

ok $factory = Grinder->new(
   -genome_file     => catfile(qw{t data shotgun_database.fa}),
   -read_dist       => 48                                     ,
   -length_bias     => 0                                      ,
   -abundance_model => ('logarithmic', 0.5)                   ,
   -total_reads     => 1000                                   ,
), 'Logarithmic community structure';

while ( $read = $factory->next_read ) {
   push @reads, $read->reference->id;
}

$ra = rank_abundance(\@reads, 10);
($min, $max, $mean, $stddev) = stats($ra);
$era = logarithmic(10, 5, 0.5, 1000);
$coeff = corr_coeff($ra, $era, $mean);
ok ($coeff > 0.97);

@reads = ();


# Exponential community structure

ok $factory = Grinder->new(
   -genome_file     => catfile(qw{t data shotgun_database.fa}),
   -read_dist       => 48                                     ,
   -length_bias     => 0                                      ,
   -abundance_model => ('exponential', 0.5)                   ,
   -total_reads     => 1000                                   ,
), 'Exponential community structure';

$struct = $factory->next_lib;
while ( $read = $factory->next_read ) {
   push @reads, $read->reference->id;
}

$ra = rank_abundance(\@reads, 10);
($min, $max, $mean, $stddev) = stats($ra);
$era = exponential(10, 5, 0.5, 1000);
$coeff = corr_coeff($ra, $era, $mean);
ok ($coeff > 0.97);
is $struct->{param}, 0.5;
@reads = ();


# Communities with random structure parameter value

ok $factory = Grinder->new(
   -genome_file     => catfile(qw{t data shotgun_database.fa}),
   -read_dist       => 48                                     ,
   -length_bias     => 0                                      ,
   -num_libraries   => 2                                      ,
   -shared_perc     => 100                                    ,
   -abundance_model => ('exponential')                        ,
   -total_reads     => 1000                                   ,
), 'Communities with random structure parameter value';

$struct = $factory->next_lib;
while ( $read = $factory->next_read ) {
   push @reads, $read->reference->id;
}

$ra = rank_abundance(\@reads, 10);
($min, $max, $mean, $stddev) = stats($ra);
$param1 = $struct->{param};
ok $param1 > 0;
ok $param1 < 1000;
$era = exponential(10, 5, $param1, 1000);
$coeff = corr_coeff($ra, $era, $mean);
ok ($coeff > 0.97);

@reads = ();

$struct = $factory->next_lib;
while ( $read = $factory->next_read ) {
   push @reads, $read->reference->id;
}

$ra = rank_abundance(\@reads, 10);
($min, $max, $mean, $stddev) = stats($ra);
$param2 = $struct->{param};
ok $param2 > 0;
ok $param2 < 1000;
$era = exponential(10, 5, $param2, 1000);
$coeff = corr_coeff($ra, $era, $mean);
ok ($coeff > 0.97);

ok $param1 != $param2;

@reads = ();



sub stats {
   # Calculates min, max, mean, stddev
   my ($vals) = @_;
   my ($min, $max, $mean, $sum, $sqsum, $stddev) = (1E99, 0, 0, 0, 0, 0);
   my $num = scalar @$vals;
   for my $val (@$vals) {
      $min = $val if $val < $min;
      $max = $val if $val > $max;
      $sum += $val;
      $sqsum += $val**2
   }
   $mean = $sum / $num;
   $stddev = sqrt( $sqsum / $num - $mean**2 );
   return $min, $max, $mean, $stddev;
}


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


sub corr_coeff {
   # The correlation coefficient R2 is
   #    R2 = 1 - ( SSerr / SStot )
   # where
   #    SSerr = sum( (y - f)**2 )
   # and
   #    SStot = sum( (y - mean)**2 )
   my ($y, $f, $mean) = @_;
   my $SSerr = 0;
   my $SStot = 0;
   for my $i ( 0 .. scalar @$y - 1 ) {
      $SSerr += ($$y[$i] - $$f[$i])**2;
      $SStot += ($$y[$i] - $mean)**2;
   }
   my $R2 = 1 - ($SSerr / $SStot);
   return $R2;
}
