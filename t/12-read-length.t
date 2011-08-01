#!perl -T

use strict;
use warnings;
use Test::More tests => 18;
use Bio::Seq;
use constant PI => 4 * atan2(1, 1);

use Grinder;
my ($factory, $nof_reads, $read, @reads, $min, $max, $mean, $stddev, $hist,
    $ehist, $coeff);


# All sequences the same length

ok $factory = Grinder->new(
   -genome_file => './t/data/shotgun_database.fa' ,
   -read_dist   => 50                             ,
   #-random_seed => 1910567890                     ,
   -total_reads => 1000                            ), 'Same length reads';

while ( $read = $factory->next_read ) {
   push @reads, $read->length;
};
($min, $max, $mean, $stddev) = stats(\@reads);
is $min, 50;
is $max, 50;
is $mean, 50;
is $stddev, 0;
@reads = ();


# Uniform distribution
ok $factory = Grinder->new(
   -genome_file => './t/data/shotgun_database.fa' ,
   -read_dist   => (50, 'uniform', 10)            ,
   #-random_seed => 1910567890                     ,
   -total_reads => 1000                            ), 'Uniform distribution';

while ( $read = $factory->next_read ) {
   push @reads, $read->length;
};
($min, $max, $mean, $stddev) = stats(\@reads);
is $min, 40;
is $max, 60;
is int($mean+0.5), 50;
ok $stddev < 6.3; # should be 5.79
ok $stddev > 5.3;

$hist = hist(\@reads, 1, 100);
$ehist = uniform(1, 100, 40, 60, 1000);
$coeff = corr_coeff($hist, $ehist, $mean);
ok ($coeff > 0.99);

@reads = ();


# Normal distribution
ok $factory = Grinder->new(
   -genome_file => './t/data/shotgun_database.fa' ,
   -read_dist   => (50, 'normal', 10)             ,
   #-random_seed => 191057890                     ,
   -total_reads => 1000                            ), 'Normal distribution';

while ( $read = $factory->next_read ) {
   push @reads, $read->length;
};
($min, $max, $mean, $stddev) = stats(\@reads);
ok $mean > 49; # should be 50.0
ok $mean < 51;
ok $stddev < 11; # should be 10.0
ok $stddev > 9;

$hist = hist(\@reads, 1, 100);
$ehist = normal(1, 100, $mean, $stddev**2, 1000);
$coeff = corr_coeff($hist, $ehist, $mean);
ok ($coeff > 0.99);

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


sub hist {
   my ($data, $min, $max) = @_;
   # Put a data series into bins
   my %hash;
   for my $val (@$data) {
      $hash{$val}++;
   }
   my @x_data = ($min .. $max);
   my @y_data;
   for my $x (@x_data) {
      my $y = $hash{$x} || 0;
      push @y_data, $y;
   }
   return \@y_data;
}


sub normal {
   # Evaluate the normal function in the given integer range
   my ($x_min, $x_max, $mean, $variance, $num) = @_;
   my @ys;
   for my $x ($x_min .. $x_max) {
      my $proba = 1 / sqrt(2 * PI * $variance) * exp( - ($x - $mean)**2 / (2 * $variance));
      my $y = $proba * $num;
      push @ys, $y;
   }
   return \@ys;
}


sub uniform {
   # Evaluate the uniform function in the given integer range
   my ($x_min, $x_max, $min, $max, $num) = @_;
   my @ys;
   my $width = $max - $min + 1;
   for my $x ($x_min .. $x_max) {
      my $y;
      if ( ($x >= $min) and ($x <= $max) ) {
         $y = $num / $width;
      } else {
         $y = 0;
      }
      push @ys, $y;
   }
   return \@ys;
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
