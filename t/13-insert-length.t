#! perl

use strict;
use warnings;
use Test::More;
use t::TestUtils;
use Grinder;

plan tests => 29;


my ($factory, $nof_reads, $mate1, $mate2, @ilengths, $min, $max, $mean, $stddev,
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
   push @ilengths, $insert_length;
};

($min, $max, $mean, $stddev) = stats(\@ilengths);
is $min, 150;
is $max, 150;
is $mean, 150;
is $stddev, 0;
@ilengths = ();


# Uniformly distributed inserts

ok $factory = Grinder->new(
   -reference_file => data('single_seq_database.fa'),
   -total_reads    => 1000                          ,
   -read_dist      => 50                            ,
   -insert_dist    => (150, 'uniform', 15)          ,
), 'Uniform distribution';

while ( $mate1 = $factory->next_read ) {
   $mate2 = $factory->next_read;
   # insert size includes mate1 + spacer + mate2
   my $insert_length = abs($mate2->end - $mate1->start + 1);
   push @ilengths, $insert_length;
};

($min, $max, $mean, $stddev) = stats(\@ilengths);
cmp_ok $min, '>=', 135;
cmp_ok $max, '<=', 165;
cmp_ok $mean, '<', 152;
cmp_ok $mean, '>', 148;
cmp_ok $stddev, '<', 10;
cmp_ok $stddev, '>', 7;

$hist = hist(\@ilengths, 50, 250);
$ehist = uniform(50, 250, 135, 165, 1000);
$coeff = corr_coeff($hist, $ehist, $mean);
cmp_ok $coeff, '>', 0.99;

SKIP: {
   skip rfit_msg(), 5 if not can_rfit();
   test_uniform_dist(\@ilengths, 135, 165, 'inserts_uniform.txt');
}

@ilengths = ();


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
   push @ilengths, $insert_length;
};

($min, $max, $mean, $stddev) = stats(\@ilengths);
cmp_ok $mean, '<', 151; # should be 150
cmp_ok $mean, '>', 149;
cmp_ok $stddev, '<', 11; # should be 10
cmp_ok $stddev, '>', 9;

$hist = hist(\@ilengths, 50, 250);
$ehist = normal(50, 250, $mean, $stddev**2, 1000);
$coeff = corr_coeff($hist, $ehist, $mean);
cmp_ok $coeff, '>', 0.99;

SKIP: {
   skip rfit_msg(), 6 if not can_rfit();
   test_normal_dist(\@ilengths, 150, 10, 'inserts_normal.txt');
}

@ilengths = ();

