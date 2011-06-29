#!perl -T

use strict;
use warnings;
use Test::More tests => 17;
use Bio::Seq;

use Grinder;
my ($factory, $nof_reads, $mate1, $mate2, @inserts, $min, $max, $mean, $stddev);


# All inserts the same length

ok $factory = Grinder->new(
   -genome_file => './t/data/single_seq_database.fa',
   -total_reads => 1000                             ,
   -random_seed => 1910567890                       ,
   -read_dist   => 50                               ,
   -insert_dist => 150                               ), 'Same size inserts';

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
   -genome_file => './t/data/single_seq_database.fa',
   -total_reads => 1000                             ,
   -random_seed => 1910567890                       ,
   -read_dist   => 50                               ,
   -insert_dist => (150, 'uniform', 30)              ), 'Uniform distribution';

while ( $mate1 = $factory->next_read ) {
   $mate2 = $factory->next_read;
   # insert size includes mate1 + spacer + mate2
   my $insert_length = abs($mate2->end - $mate1->start + 1);
   push @inserts, $insert_length;
};

($min, $max, $mean, $stddev) = stats(\@inserts);
ok $min >= 120;
ok $max <= 180;
ok $mean < 151;
ok $mean > 149;
ok $stddev < 18;
ok $stddev > 16;
@inserts = ();


# Normally distributed inserts

ok $factory = Grinder->new(
   -genome_file => './t/data/single_seq_database.fa',
   -total_reads => 1000                             ,
   -random_seed => 1910567890                       ,
   -read_dist   => 50                               ,
   -insert_dist => (150, 'normal', 10)              ), 'Normal distribution';

while ( $mate1 = $factory->next_read ) {
   $mate2 = $factory->next_read;
   # insert size includes mate1 + spacer + mate2
   my $insert_length = abs($mate2->end - $mate1->start + 1);
   push @inserts, $insert_length;
};

($min, $max, $mean, $stddev) = stats(\@inserts);
ok $mean < 151;
ok $mean > 149;
ok $stddev < 11;
ok $stddev > 9;
@inserts = ();



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
