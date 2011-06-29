#!perl -T

use strict;
use warnings;
use Test::More tests => 17;
use Bio::Seq;

use Grinder;
my ($factory, $nof_reads, $read, @reads, $min, $max, $mean, $stddev);


# All sequences the same length

ok $factory = Grinder->new(
   -genome_file => './t/data/shotgun_database.fa' ,
   -read_dist   => 50                             ,
   -random_seed => 1910567890                     ,
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
   -random_seed => 1910567890                     ,
   -total_reads => 1000                            ), 'Uniform distribution';

while ( $read = $factory->next_read ) {
   push @reads, $read->length;
};
($min, $max, $mean, $stddev) = stats(\@reads);
ok $min >= 40;
ok $max <= 60;
ok $mean < 51;
ok $mean > 49;
ok $stddev < 6;
ok $stddev > 5;
@reads = ();

# Normal distribution
ok $factory = Grinder->new(
   -genome_file => './t/data/shotgun_database.fa' ,
   -read_dist   => (50, 'normal', 10)             ,
   -random_seed => 191057890                     ,
   -total_reads => 1000                            ), 'Normal distribution';

while ( $read = $factory->next_read ) {
   push @reads, $read->length;
};
($min, $max, $mean, $stddev) = stats(\@reads);
ok $mean < 51;
ok $mean > 49;
ok $stddev < 11;
ok $stddev > 9;
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
