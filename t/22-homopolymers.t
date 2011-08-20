#! perl

use strict;
use warnings;
use Test::More;
use t::TestUtils;
use Grinder;

plan tests => 6648;


my ($factory, $nof_reads, $read, $hpols, $min, $max, $mean, $stddev,
    $expected_mean, $expected_stddev, $hist, $ehist, $coeff);
my $delta_perc = 0.05;


# Balzer homopolymer distribution

ok $factory = Grinder->new(
   -genome_file      => data('homopolymer_database.fa'),
   -unidirectional   => 1                              ,
   -read_dist        => 220                            ,
   -total_reads      => 1000                           ,
   -homopolymer_dist => 'balzer'                       ,
), 'Balzer';

while ( $read = $factory->next_read ) {
   my ($error_str) = ($read->desc =~ /errors=(\S+)/);
   $hpols = add_homopolymers($error_str, $read->reference->seq, $hpols);
#   if ($error_str) {
#      ok $error_str !~ m/%/g;
#      ok $error_str =~ m/[+-]/g;
#   } else {
#      ok 1;
#      ok 1;
#   }
}


for my $homo_len ( 10 ) {
#for my $homo_len ( sort {$b <=> $a} (keys %$hpols) ) {
   ### make tests work for < 4 !
#   next if $homo_len < 4; 
   ###

   my $values = $$hpols{$homo_len};
   print "Homopolymer length: $homo_len\n";
   ($min, $max, $mean, $stddev) = stats($values);
   print "   min = $min, max = $max, mean = $mean, stddev = $stddev\n";
   ($expected_mean, $expected_stddev) = balzer($homo_len);
   print "   expected mean = $expected_mean, expected stddev = $expected_stddev\n";

#   ok $mean > (1-$delta_perc)*$expected_mean;
#   ok $mean < (1+$delta_perc)*$expected_mean;
#   ok $stddev > (1-$delta_perc)*$expected_stddev;
#   ok $stddev < (1+$delta_perc)*$expected_stddev; 
#   $hist = hist($$hpols{$homo_len}, 1, 20);
#   $ehist = normal(1, 20, $mean, $stddev**2, 4000); # 4 homopolymers of each size in the 1000 reads
#   $coeff = corr_coeff($hist, $ehist, $mean);
#   print "   coeff = $coeff\n";
#   ok ($coeff > 0.80);
#   #ok ($coeff > 0.99);

    test_normal_dist($values, $expected_mean, $expected_stddev);

}
$hpols = {};


# Richter homopolymer distribution

#ok $factory = Grinder->new(
#   -genome_file      => data('homopolymer_database.fa'),
#   -unidirectional   => 1                              ,
#   -read_dist        => 220                            ,
#   -total_reads      => 1000                           ,
#   -homopolymer_dist => 'richter'                      ,
#), 'Richter';

#while ( $read = $factory->next_read ) {
#   my ($error_str) = ($read->desc =~ /errors=(\S+)/);
#   if ($error_str) {
#      ok $error_str !~ m/%/g;
#      ok $error_str =~ m/[+-]/g;
#   } else {
#      ok 1;
#      ok 1;
#   }
#   $hpols = add_homopolymers($error_str, $read->reference->seq, $hpols);
#}

#for my $homo_len ( sort {$b <=> $a} (keys %$hpols) ) {
#   ### make tests work for < 4 !
#   next if $homo_len < 4; 
#   ###
#   print "Homopolymer length: $homo_len\n";
#   ($min, $max, $mean, $stddev) = stats($$hpols{$homo_len});
#   print "   min = $min, max = $max, mean = $mean, stddev = $stddev\n";
#   ($expected_mean, $expected_stddev) = richter($homo_len);
#   print "   expected mean = $expected_mean, expected stddev = $expected_stddev\n";
##   ok $mean > (1-$delta_perc)*$expected_mean;
##   ok $mean < (1+$delta_perc)*$expected_mean;
##   ok $stddev > (1-$delta_perc)*$expected_stddev;
##   ok $stddev < (1+$delta_perc)*$expected_stddev; 
#   $hist = hist($$hpols{$homo_len}, 1, 20);
#   $ehist = normal(1, 20, $mean, $stddev**2, 4000); # 4 homopolymers of each size in the 1000 reads
#   $coeff = corr_coeff($hist, $ehist, $mean);
#   print "   coeff = $coeff\n";
#   ok ($coeff > 0.80);
#   #ok ($coeff > 0.99);
#}
#$hpols = {};

####
#die "Temp exit\n";
####

# Margulies homopolymer distribution

#ok $factory = Grinder->new(
#   -genome_file      => data('homopolymer_database.fa'),
#   -unidirectional   => 1                              ,
#   -read_dist        => 220                            ,
#   -total_reads      => 1000                           ,
#   -homopolymer_dist => 'margulies'                    ,
#), 'Margulies';

#while ( $read = $factory->next_read ) {
#   my ($error_str) = ($read->desc =~ /errors=(\S+)/);
#
#   print "$error_str\n";
#
#   if ($error_str) {
#      ok $error_str !~ m/%/g;
#      ok $error_str =~ m/[+-]/g;
#   } else {
#      ok 1;
#      ok 1;
#   }
#   $hpols = add_homopolymers($error_str, $read->reference->seq, $hpols);
#}

#use Data::Dumper; print Dumper($hpols);

#for my $homo_len ( sort {$b <=> $a} (keys %$hpols) ) {
#   ### make tests work for < 4 !
#   next if $homo_len < 4; 
#   ###
#   print "Homopolymer length: $homo_len\n";
#   ($min, $max, $mean, $stddev) = stats($$hpols{$homo_len});
#   print "   min = $min, max = $max, mean = $mean, stddev = $stddev\n";
#   ($expected_mean, $expected_stddev) = margulies($homo_len);
#   print "   expected mean = $expected_mean, expected stddev = $expected_stddev\n";
#   ok $mean > (1-$delta_perc)*$expected_mean;
#   ok $mean < (1+$delta_perc)*$expected_mean;
#   ok $stddev > (1-$delta_perc)*$expected_stddev;
#   ok $stddev < (1+$delta_perc)*$expected_stddev; 
#   $hist = hist($$hpols{$homo_len}, 1, 20);
#
#   use Data::Dumper; print Dumper($hist);
#
#   $ehist = normal(1, 20, $mean, $stddev**2, 4000); # 4 homopolymers of each size in the 1000 reads
#   $coeff = corr_coeff($hist, $ehist, $mean);
#   print "   coeff = $coeff\n";
#   ok ($coeff > 0.80);
#   #ok ($coeff > 0.99);
#}
#$hpols = {};



sub add_homopolymers {
   my ($err_str, $ref_seq, $err_h) = @_;

   # Record position and length of homopolymer errors
   my %errors;
   if (defined $err_str) {
      $err_str = combine_dels($err_str);
      my @errors = split ',', $err_str;
      for my $error (@errors) {
         # Record homopolymer error
         my ($pos, $type, $repl) = ($error =~ m/(\d+)([%+-])(.*)/i);
         my $elen; # error length
         if ($type eq '-') {
            $repl .= '-';
            $elen = - length($repl);
         } elsif ($type eq '+') {
            $elen = + length($repl);
         }
         $errors{$pos} = $elen;
         # Test that proper residue was added to homopolymer
         my $hres    = substr $ref_seq, $pos-1, 1; # first residue of homopolymer
         my $new_res = substr $repl, 0, 1;         # residue to add to homopolymer
         if ( $type eq '+' ) { # in case of insertion (not deletion)
            ok $hres eq $new_res;
         }
      }
   }

   # Record all homopolymer lengths
   while ( $ref_seq =~ m/(.)(\1+)/g ) {
      # Found a homopolymer
      my $hlen = length($2) + 1;               # length of the error-free homopolymer
      my $pos  = pos($ref_seq) - $hlen + 1;    # start of the homopolymer (residue no.)
      my $elen = $hlen + ($errors{$pos} || 0); # length of the error-containing homopolymer
      push @{$$err_h{$hlen}}, $elen;
   }

   return $err_h;
}



sub combine_dels {
   # Put deletions in adjacent position into a single error entry.
   # Ex: 45-,46-,47-    becomes 45---
   my ($err_str) = @_;
   my %errors;
   for my $error (split ',', $err_str) {
      my ($pos, $type, $repl) = ($error =~ m/(\d+)([%+-])([a-z]*)/i);
      push @{$errors{$pos}{$type}}, $repl;
      if ( $type eq '-' ) {
         $errors{$pos}{'-'} = [ undef ]; # this prevents redundant deletions
      }
   }

   for my $pos (sort {$b <=> $a} (keys %errors)) {
      if ( exists $errors{$pos}{'-'} && exists $errors{$pos-1} && exists $errors{$pos-1}{'-'} ) {
         push @{$errors{$pos-1}{'-'}}, @{$errors{$pos}{'-'}};
         delete $errors{$pos}{'-'};
         delete $errors{$pos} if scalar keys %{$errors{$pos}} == 0;
      }
   }

   $err_str = '';
   for my $pos (sort {$a <=> $b} (keys %errors)) {
      while ( my ($type, $repls) = each %{$errors{$pos}} ) {
         my $repl;
         if ($type eq '-') {
            $repl = '-' x (scalar @$repls - 1);
         } else {
            $repl = join '', @$repls;
         }
         $err_str .= $pos.$type.$repl.',';
      }
   }
   $err_str =~ s/,$//;
   return $err_str;
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

      ####
      #$y = int($y + 0.5);
      ####

      push @ys, $y;
   }
   return \@ys;
}


sub margulies {
   my ($homo_len) = @_;
   my $mean     = $homo_len;
   my $stddev = $homo_len * 0.15;
   return $mean, $stddev;
}


sub richter {
   my ($homo_len) = @_;
   my $mean     = $homo_len;
   my $stddev = sqrt($homo_len) * 0.15;
   return $mean, $stddev;
}


sub balzer {
   my ($homo_len) = @_;
   my $mean     = $homo_len;
   my $variance = 0.03494 + $homo_len * 0.06856;
   return $mean, $stddev;
}
