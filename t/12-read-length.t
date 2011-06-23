#!perl -T

use strict;
use warnings;
use Test::More tests => 13;
#use Test::More tests => 17;
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
   push @reads, $read;
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
   push @reads, $read;
};
($min, $max, $mean, $stddev) = stats(\@reads);
ok $min >= 40;
ok $max <= 60;
ok $mean < 51;
ok $mean > 49;
####ok $stddev < XXX;
####ok $stddev > XXX;
@reads = ();

# Normal distribution
ok $factory = Grinder->new(
   -genome_file => './t/data/shotgun_database.fa' ,
   -read_dist   => (50, 'normal', 10)             ,
   -random_seed => 1910567890                     ,
   -total_reads => 1000                            ), 'Normal distribution';

while ( $read = $factory->next_read ) {
   push @reads, $read;
};
($min, $max, $mean, $stddev) = stats(\@reads);
ok $mean < 51;
ok $mean > 49;
####ok $stddev < XXX;
####ok $stddev > XXX;
@reads = ();


sub stats {
   # Calculates min, max, mean, stddev
   my ($reads) = @_;
   my ($min, $max, $mean, $sum, $sqsum, $stddev) = (1E99, 0, 0, 0, 0, 0);
   my $num = scalar @$reads;
   for my $read (@$reads) {
      my $length = $read->length;
      $min = $length if $length < $min;
      $max = $length if $length > $max;
      $sum += $length;
      $sqsum += $length**2
   }
   $mean = $sum / $num;
  
   #### calculate stddev

   return $min, $max, $mean, $stddev;
}
