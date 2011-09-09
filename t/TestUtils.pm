package t::TestUtils;


use strict;
use warnings;
use Test::More;
use Statistics::R;
use File::Spec::Functions;
use List::Util qw( min max );

use vars qw{@ISA @EXPORT};
BEGIN {
   @ISA     = 'Exporter';
   @EXPORT  = qw{
      PI
      data
      stats
      hist
      corr_coeff
      write_data
      can_rfit
      test_normal_dist
      test_uniform_dist      
   };
}

our $can_rfit;


#------------------------------------------------------------------------------#


# The Pi mathematical constant
use constant PI => 4 * atan2(1, 1);


sub data {
   # Get the complete filename of a test data file
   return catfile('t', 'data', @_);
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


sub hist {
   # Count the number of occurence of each integer:
   # E.g. given the arrayref:
   #    [ 1, 1, 1, 3, 3, 4 ]
   # Return the arrayref:
   #    [ 3, 0, 2, 4 ]
   # The min and the max of the range to consider can be given as an option
   my ($data, $min, $max) = @_;
   if (not defined $data) {
      die "Error: no data provided to hist()\n";
   }
   my %hash;
   for my $val (@$data) {
      $hash{$val}++;
   }
   $min = min(@$data) if not defined $min;
   $max = max(@$data) if not defined $max;
   my @y_data;
   for my $x ($min .. $max) {
      push @y_data, $hash{$x} || 0;
   }
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
      #print "  ".($i+1)."  ".$$y[$i]."  ".$$f[$i]."\n";
      $SSerr += ($$y[$i] - $$f[$i])**2;
      $SStot += ($$y[$i] - $mean)**2;
   }
   my $R2 = 1 - ($SSerr / $SStot);
   return $R2;
}


sub write_data {
   # Write a data series (array reference) to a file with the specified name, or
   # 'data.txt' by default
   my ($data, $filename) = @_;
   $filename = 'data.txt' if not defined $filename;
   open my $out, '>', $filename or die "Error: Could not write file $filename\n$!\n";
   for my $datum (@$data) {
      print $out "$datum\n";
   }
   close $out;
   return $filename;
}


sub can_rfit {
   # Test if a system can run the fitdistrplus R module through the Statistics::R
   # Perl interface. Return 1 if it can, 0 otherwise.
   if (not defined $can_rfit) {
      eval {
         require Statistics::R;
         my $R = Statistics::R->new();
         my $ret = $R->run(q`library(fitdistrplus)`);
         $R->stop();
      };
      if ($@) {
         $can_rfit = 0;
         my $msg = "Warning: The Statistics::R module for Perl, R (R-Project) ".
            "or the fitdistrplus module for R could not be found on this system.".
            " Some tests will be skipped...\n";
         warn $msg;
      } else {
         $can_rfit = 1;
      }
   }
   return $can_rfit;
}


#sub test_linear_dist {
#   # Test that the datapoints provided follow a linear distribution
#   my ($values, $want_min, $want_max, $want_slope) = @_;
#   my ($min, $max, $slope, $chisqtest) = fit_linear($values);
#
#   # Assume a 5% standard deviation 
#   my $min_sd   = 0.05 * ($max - $min);
#   my $max_sd   = $min_sd;
#   my $slope_sd = 0.05 * $slope;
#
#   # Test
#   cmp_ok $want_min  , '>', $min   - $min_sd  ;
#   cmp_ok $want_min  , '<', $min   + $min_sd  ;
#   cmp_ok $want_max  , '>', $max   - $max_sd  ;
#   cmp_ok $want_max  , '<', $max   + $max_sd  ;
#   cmp_ok $want_slope, '>', $slope - $slope_sd;
#   cmp_ok $want_slope, '<', $slope + $slope_sd;
#   is $chisqtest, 'not rejected';
#
#   # a beta distribution with shape1=1 & shape2=2 is linearly decreasing
#   # a beta distribution with shape1=2 & shape2=1 is linearly increasing
#
#   return 1;
#}


sub test_uniform_dist {
   # Test that the datapoints provided follow a uniform distribution with the 
   # specified minimum and maximum. Note that you probably need over 30-100
   # values for the statistical test to work!
   my ($values, $want_min, $want_max, $filename) = @_;
   my ($min, $max, $chisqtest) = fit_uniform($values, $want_min, $want_max);
   # Assume a standard deviation that is 5% of the min - max interval
   my $min_sd = 0.05 * ($max - $min);
   my $max_sd = $min_sd;
   # Test now
   #print "min = $min +- $min_sd, max = $max +- $max_sd\n"; 
   cmp_ok $want_min, '>', $min - $min_sd, 'fitdist() uniform';
   cmp_ok $want_min, '<', $min + $min_sd;
   cmp_ok $want_max, '>', $max - $max_sd;
   cmp_ok $want_max, '<', $max + $max_sd;
   is( $chisqtest, 'not rejected', 'Chi square') or write_data($values, $filename);
   return 1;
}


sub test_normal_dist {
   # Test that the datapoints provided follow a normal distribution with the 
   # specified mean and standard deviation. Note that you probably need over
   # 30-100 values for the statistical test to work!
   my ($values, $want_mean, $want_sd, $filename) = @_;
   my ($mean, $mean_sd, $sd, $sd_sd, $cvmtest, $adtest) = fit_normal($values, $want_mean, $want_sd);
   #print "mean = $mean +- $mean_sd, sd = $sd +- $sd_sd\n";
   # A interval of mean+-1.96*sd corresponds to the 2.5 and 97.5 percentiles of 
   # the normal distribution, i.e. the 95% confidence interval
   cmp_ok $want_mean, '>', $mean - 1.96 * $mean_sd, 'fitdist() normal';
   cmp_ok $want_mean, '<', $mean + 1.96 * $mean_sd;
   cmp_ok $want_sd  , '>',   $sd - 1.96 * $sd_sd  ;
   cmp_ok $want_sd  , '<',   $sd + 1.96 * $sd_sd  ;
   # Cramer-von Mises test
   is( $cvmtest, 'not rejected', 'Cramer-von Mises' ) or write_data($values, $filename);
   # Anderson-Darling test (emphasizes the tails of a distribution)
   is( $adtest , 'not rejected', 'Anderson-Darling' ) or write_data($values, $filename);
   return 1;
}


#------------------------------------------------------------------------------#


#sub fit_linear {
#   my ($values, $ascending) = @_;
#
#   # Find min and max of series
#   my ($min, $max) = (1E100, 0);
#   for my $value (@$values) {
#      $min = $value if $value < $min;
#      $max = $value if $value > $max;
#   }
#
#   # Rescale values between in 0 and 1 instead of min and max
#   my $rescaled_values;
#   if ( ($min == 0) and ($max == 1) ) {
#      $rescaled_values = $values;
#   } else {
#      for my $value (@$values) {
#         push @$rescaled_values, ($value - $min) / $max;
#      }
#   }
#
#   # Now we can run fit_beta()
#   my ($shape1, $shape1_sd, $shape2, $shape2_sd, $chisqtest) = fit_beta($values);
#
#   ###
#   print "xmin = $min, xmax = $max, shape1 = $shape1 +- $shape1_sd, shape2 = $shape2 +- $shape2_sd\n";
#   ###
#
#   ### what is the slope??
#   return $min, $max;
#}


sub fit_beta {
   # Try to fit a beta distribution to a series of data points using a maximum
   # goodness of fit method. Return the shape1 parameter, its standard error,
   # the shape2 parameter, its standard error, and the results of Chi square
   # statistics. Optional optimization starting values for shape1 and shape2 can
   # be specified.
   my ($values, $start_shape1, $start_shape2) = @_;
   my $start_params = '';
   if ( (defined $start_shape1) or (defined $start_shape2) ) {
      $start_params .= ', start=list(';
      if (defined $start_shape1) {
         $start_params .= "shape1=$start_shape1, ";
      }
      if (defined $start_shape2) {
         $start_params .= "shape2=$start_shape2, ";
      }
      $start_params =~ s/, $//;
      $start_params .= ')';
   }
   my $out = '';
   my $R = Statistics::R->new();
   $R->run(q`library(fitdistrplus)`);
   $R->set('x', $values);
   $R->run(qq`f <- fitdist(x, distr="beta", method="mle"$start_params)`);
   $R->run(q`g <- gofstat(f)`);
   $R->run(q`shape1 <- f$estimate[1]`);
   my $shape1 = $R->get('shape1');
   $R->run(q`shape2 <- f$estimate[2]`);
   my $shape2 = $R->get('shape2');
   $R->run(q`shape1_sd <- f$sd[1]`);
   my $shape1_sd = $R->get('shape1_sd');
   $R->run(q`shape2_sd <- f$sd[2]`);
   my $shape2_sd = $R->get('shape2_sd');
   $R->run(q`chisqpvalue <- g$chisqpvalue`);
   my $chisqtest = test_result( $R->get('chisqpvalue') );
   $R->stop();
   return $shape1, $shape1_sd, $shape2, $shape2_sd, $chisqtest;
}


sub fit_uniform {
   # Try to fit a uniform distribution to a series of integers using a maximum
   # goodness of fit method. Return the minimum, maximum and the results of
   # the Chi square statistics. Optional optimization starting values for the
   # mean and standard deviation can be specified. If used, these starting values
   # also defined the breaks to use for the Chi square test.
   my ($values, $start_min, $start_max) = @_;
   my $start_p = '';
   my $breaks_p = '';
   if ( (defined $start_min) and (defined $start_max) ) {
      $start_p  .= ", start=list(min=$start_min, max=$start_max)";
      $breaks_p .= ', chisqbreaks=seq('.($start_min-0.5).', '.($start_max+0.5).')';
   }
   my $out = '';
   my $R = Statistics::R->new();
   $R->set('x', $values);
   $R->run(qq`library(fitdistrplus)`);
   $R->run(qq`f <- fitdist(x, distr="unif", method="mge", gof="CvM"$start_p)`);
   $R->run(qq`g <- gofstat(f$breaks_p)`);
   my $min         = $R->get('f$estimate[1]');
   my $max         = $R->get('f$estimate[2]');
   my $chisqpvalue = $R->get('g$chisqpvalue');
   my $chisqtest   = test_result($chisqpvalue);
   $R->stop();
   return $min, $max, $chisqtest;
}


sub fit_normal {
   # Try to fit a normal distribution to a series of data points using a maximum
   # likelihood method. Return the mean and its standard deviation, the standard
   # deviation and its standard deviation, and the results of Cramer-von Mises
   # and Anderson-Darling statistics. Optional optimization starting values for
   # the mean and standard deviation can be specified.
   my ($values, $start_mean, $start_sd) = @_;
   my $start_params = '';
   if ( (defined $start_mean) or (defined $start_sd) ) {
      $start_params .= ', start=list(';
      if (defined $start_mean) {
         $start_params .= "mean=$start_mean, ";
      }
      if (defined $start_sd) {
         $start_params .= "sd=$start_sd, ";
      }
      $start_params =~ s/, $//;
      $start_params .= ')';
   }
   my $out = '';
   my $R = Statistics::R->new();
   $R->run(q`library(fitdistrplus)`);
   $R->set('x', $values);
   $R->run(qq`f <- fitdist(x, distr="norm", method="mle"$start_params)`);
   $R->run(q`g <- gofstat(f)`);
   $R->run(q`mean <- f$estimate[1]`);
   my $mean = $R->get('mean');
   $R->run(q`sd <- f$estimate[2]`);
   my $sd = $R->get('sd');
   $R->run(q`mean_sd <- f$sd[1]`);
   my $mean_sd = $R->get('mean_sd');
   $R->run(q`sd_sd <- f$sd[2]`);
   my $sd_sd = $R->get('sd_sd');
   $R->run(q`adtest <- g$adtest`);
   my $adtest = $R->get('adtest');
   $R->run(q`cvmtest <- g$cvmtest`);
   my $cvmtest = $R->get('cvmtest');
   $R->stop();
   return $mean, $mean_sd, $sd, $sd_sd, $cvmtest, $adtest;
}


sub test_result {
   # Reject a statistical test if the p value is less than 0.05
   my ($p_value) = @_;
   my $test_result;
   if ( lc $p_value eq 'nan' ) {
      $p_value = 1; # probably a very large p value
   }
   my $thresh = 0.05;
   if ($p_value <= $thresh) {
      $test_result = 'rejected';
   } elsif ($p_value > $thresh) {
      $test_result = 'not rejected';
   } else {
      die "Error: '$p_value' is not a supported p value\n";
   }
   return $test_result;
}


1;
