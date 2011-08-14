#!perl -T

use strict;
use warnings;
use Test::More tests => 6013;
use Bio::Seq;

use Grinder;
my ($factory, $nof_reads, $read, $errors, $min, $max, $mean, $stddev, $prof,
    $eprof, $coeff);


# No errors by default

ok $factory = Grinder->new(
   -genome_file    => './t/data/single_seq_database.fa',
   -unidirectional => 1                                ,
   -read_dist      => 50                               ,
   -total_reads    => 1000                              ), 'No errors';

while ( $read = $factory->next_read ) {
   is $read->seq, 'aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa';
   ok $read->desc !~ /errors/;
}


# Substitutions

ok $factory = Grinder->new(
   -genome_file    => './t/data/single_seq_database.fa',
   -unidirectional => 1                                ,
   -read_dist      => 50                               ,
   -total_reads    => 1000                             ,
   -mutation_ratio => (100, 0)                         ,
   -mutation_dist  => (10, 'uniform')                     ), 'Substitutions only';

while ( $read = $factory->next_read ) {
   my ($error_str) = ($read->desc =~ /errors=(\S+)/);
   if ($error_str) {
      ok $error_str =~ m/%/g;
      ok $error_str !~ m/[+-]/g;
   } else {
      ok 1;
      ok 1;
   }
}


# Indels

ok $factory = Grinder->new(
   -genome_file    => './t/data/single_seq_database.fa',
   -unidirectional => 1                                ,
   -read_dist      => 50                               ,
   -total_reads    => 1000                             ,
   -mutation_ratio => (0, 100)                         ,
   -mutation_dist  => (10, 'uniform')                     ), 'Indels only';

while ( $read = $factory->next_read ) {
   my ($error_str) = ($read->desc =~ /errors=(\S+)/);
   if ($error_str) {
      ok $error_str !~ m/%/g;
      ok $error_str =~ m/[+-]/g;
   } else {
      ok 1;
      ok 1;
   }
}


# Uniform distribution

ok $factory = Grinder->new(
   -genome_file    => './t/data/single_seq_database.fa',
   -unidirectional => 1                                ,
   -read_dist      => 50                               ,
   -total_reads    => 1000                             ,
   -mutation_ratio => (50, 50)                         ,
   -mutation_dist  => (10, 'uniform')                   ), 'Uniform';

while ( $read = $factory->next_read ) {
   my ($error_str) = ($read->desc =~ /errors=(\S+)/);
   $errors = add_error_position($error_str, $errors);
}
$prof = error_profile($errors, 50);
($min, $max, $mean, $stddev) = stats($prof);
#print "min = $min, max = $max, mean = $mean, stddev = $stddev\n";


### should be closer to 100 ideally, e.g. (98, 102)
ok $min >= 60;
ok $max <= 140;
ok $mean < 105; # 100
ok $mean > 90;
ok $stddev < 15;
#$eprof = uniform(50, 10, 1000); # 1000 reads * 50 bp * 10% error
#$coeff = corr_coeff($prof, $eprof, $mean);
#print "coeff = $coeff\n";
#ok ($coeff > 0.99);
$errors = {};


# Linear distribution

ok $factory = Grinder->new(
   -genome_file    => './t/data/single_seq_database.fa',
   -unidirectional => 1                                ,
   -read_dist      => 50                               ,
   -total_reads    => 1000                             ,
   -mutation_ratio => (50, 50)                         ,
   -mutation_dist  => (10, 'linear', 15)                 ), 'Linear';

while ( $read = $factory->next_read ) {
   my ($error_str) = ($read->desc =~ /errors=(\S+)/);
   $errors = add_error_position($error_str, $errors);
}
$prof = error_profile($errors, 50);
($min, $max, $mean, $stddev) = stats($prof);
#print "min = $min, max = $max, mean = $mean, stddev = $stddev\n";

### mean should be closer to 100 ideally, e.g. (98, 102) and coeff closer to 0.99
ok $mean < 105; # 100
ok $mean > 90;
$eprof = linear(50, 10, 15, 1000);
$coeff = corr_coeff($prof, $eprof, $mean);
#print "coeff= $coeff\n";
ok ($coeff > 0.75);
#ok ($coeff > 0.99);
###

$errors = {};



sub error_profile {
  my ($errors, $max) = @_;
  my @profile;
  for my $position ( 1 .. $max ) {
    push @profile, $$errors{$position} || 0;
  }
  return \@profile;
}


sub add_error_position {
   my ($err_str, $err_h) = @_;
   if (defined $err_str) {
      my @errors = split ',', $err_str;
      for my $error (@errors) {
         my ($pos, $type, $res) = ($error =~ m/(\d+)([%+-])([a-z]*)/i);
         $$err_h{$pos}++
      }
   }
   return $err_h;
}


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
   my ($x_max, $mean, $num) = @_;
   my @ys;
   for my $x (1 .. $x_max) {
      my $y = $num * $mean / 100;
      push @ys, $y;
   }
   return \@ys;
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
      ####
      #print "".($i+1)."  ".$$y[$i]."  ".$$f[$i]."\n";
      ####
      $SSerr += ($$y[$i] - $$f[$i])**2;
      $SStot += ($$y[$i] - $mean)**2;
   }
   ####
   #print "SSerr = $SSerr\n";
   #print "SStot = $SStot\n";
   ####
   my $R2 = 1 - ($SSerr / $SStot);
   return $R2;
}
