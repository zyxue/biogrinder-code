#! perl

use strict;
use warnings;
use Test::More;
use t::TestUtils;
use Grinder;

plan tests => 30;


my ($factory, $nof_reads, $mate1, $mate2, @inserts, $min, $max, $mean, $stddev,
    $hist, $ehist, $coeff);


# All inserts the same length

ok $factory = Grinder->new(
   -reference_file => data('single_seq_database.fa'),
   -total_reads    => 1000                          ,
   -read_dist      => 50                            ,
   -insert_dist    => 150                           ,
), 'Same size inserts';

while ( $mate1 = $factory->next_read ) {
   $mate2 = $factory->next_read;
   # insert size includes mate1 + spacer + mate2
   my $insert_length = abs($mate2->end - $mate1->start + 1);
   push @inserts, $insert_length;
};

($min, $max, $mean, $stddev) = stats(\@inserts);
is $min, 150;
is $max, 150;
is $mean, 150;
is $stddev, 0;
@inserts = ();


# Uniformly distributed inserts

ok $factory = Grinder->new(
   -reference_file => data('single_seq_database.fa'),
   -total_reads    => 1000                          ,
   -read_dist      => 50                            ,
   -insert_dist    => (150, 'uniform', 30)          ,
), 'Uniform distribution';

while ( $mate1 = $factory->next_read ) {
   $mate2 = $factory->next_read;
   # insert size includes mate1 + spacer + mate2
   my $insert_length = abs($mate2->end - $mate1->start + 1);
   push @inserts, $insert_length;
};
write_data(\@inserts, 'insert_uniform.txt');

($min, $max, $mean, $stddev) = stats(\@inserts);
cmp_ok $min, '>=', 120;
cmp_ok $max, '<=', 180;
cmp_ok $mean, '<', 152;
cmp_ok $mean, '>', 148;
cmp_ok $stddev, '<', 19;
cmp_ok $stddev, '>', 15;

$hist = hist(\@inserts, 50, 250);
$ehist = uniform(50, 250, 120, 180, 1000);
$coeff = corr_coeff($hist, $ehist, $mean);
cmp_ok $coeff, '>', 0.99;

if ( can_rfit() ) {
   test_uniform_dist(\@inserts, 120, 180);
} else {
   SKIP: {
      skip "Cannot use the fitdistrplus R module on this system", 5;
   }
}

@inserts = ();


# Normally distributed inserts

ok $factory = Grinder->new(
   -reference_file => data('single_seq_database.fa'),
   -total_reads    => 1000                          ,
   -read_dist      => 50                            ,
   -insert_dist    => (150, 'normal', 10)           ,
), 'Normal distribution';

while ( $mate1 = $factory->next_read ) {
   $mate2 = $factory->next_read;
   # insert size includes mate1 + spacer + mate2
   my $insert_length = abs($mate2->end - $mate1->start + 1);
   push @inserts, $insert_length;
};
write_data(\@inserts, 'insert_normal.txt');

($min, $max, $mean, $stddev) = stats(\@inserts);
cmp_ok $mean, '<', 151; # should be 150
cmp_ok $mean, '>', 149;
cmp_ok $stddev, '<', 11; # should be 10
cmp_ok $stddev, '>', 9;

$hist = hist(\@inserts, 50, 250);
$ehist = normal(50, 250, $mean, $stddev**2, 1000);
$coeff = corr_coeff($hist, $ehist, $mean);
cmp_ok $coeff, '>', 0.99;

if ( can_rfit() ) {
   test_normal_dist(\@inserts, 150, 10);
} else {
   SKIP: {
      skip "Cannot use the fitdistrplus R module on this system", 6;
   }
}

@inserts = ();




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


