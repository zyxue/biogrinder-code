package t::TestUtils;


use strict;
use warnings;
use Test::More;
use Statistics::R;
use File::Spec::Functions;

use vars qw{@ISA @EXPORT};
BEGIN {
   @ISA     = 'Exporter';
   @EXPORT  = qw{
      PI
      data
      can_rfit
      test_normal_dist
      test_uniform_dist

      
   };
}


#------------------------------------------------------------------------------#

# The Pi mathematical constant
use constant PI => 4 * atan2(1, 1);


sub data {
   # Get the complete filename of a test data file
   return catfile('t', 'data', @_);
}


sub can_rfit {
   # Test if a system can run the fitdistrplus R module through the Statistics::R
   # Perl interface. Return 1 if it can, 0 otherwise.
   my $can_run = 1;
   eval {
      require Statistics::R;
      my $R = Statistics::R->new();
      $R->startR;
      $R->send(q`library(fitdistrplus)`);
      my $ret = $R->read;
      die "$ret\n" if $ret;
      $R->stopR();
   };
   if ($@) {
      $can_run = 0;
      my $msg = "Skipping these tests because one of the following is not installed:\n".
                "   Statistics::R module for Perl, R (R-Project) or fitdistrplus module for R\n";
      warn $msg;
   }
   return $can_run;
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
#   ok $want_min > $min - $min_sd;
#   ok $want_min < $min + $min_sd;
#   ok $want_max > $max - $max_sd;
#   ok $want_max < $max + $max_sd;
#   ok $want_slope > $slope - $slope_sd;
#   ok $want_slope < $slope + $slope_sd;
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
   my ($values, $want_min, $want_max) = @_;
   my ($min, $max, $chisqtest) = fit_uniform($values, $want_min, $want_max);
   # Assume a standard deviation that is 5% of the min - max interval
   my $min_sd = 0.05 * ($max - $min);
   my $max_sd = $min_sd;
   # Test now
   #print "min = $min +- $min_sd, max = $max +- $max_sd\n"; 
   ok $want_min > $min - $min_sd;
   ok $want_min < $min + $min_sd;
   ok $want_max > $max - $max_sd;
   ok $want_max < $max + $max_sd;
   is $chisqtest, 'not rejected';
   return 1;
}


sub test_normal_dist {
   # Test that the datapoints provided follow a normal distribution with the 
   # specified mean and standard deviation. Note that you probably need over
   # 30-100 values for the statistical test to work!
   my ($values, $want_mean, $want_sd) = @_;
   my ($mean, $mean_sd, $sd, $sd_sd, $cvmtest, $adtest) = fit_normal($values, $want_mean, $want_sd);
   #print "mean = $mean +- $mean_sd, sd = $sd +- $sd_sd\n"; 
   ok $want_mean > $mean - $mean_sd;
   ok $want_mean < $mean + $mean_sd;
   ok $want_sd > $sd - $sd_sd;
   ok $want_sd < $sd + $sd_sd;
   is $cvmtest, 'not rejected';
   is $adtest , 'not rejected';
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
   my $x = join ', ', @$values;
   my $out = '';
   my $R = Statistics::R->new();
   $R->startR;
   $R->run(q`library(fitdistrplus)`);
   $R->run(qq`x <- c($x)`);
   $R->run(qq`f <- fitdist(x, distr="beta", method="mle"$start_params)`);
   $R->run(q`g <- gofstat(f)`);
   $R->run(q`shape1 <- f$estimate[1]`);
   my $shape1 = $R->value('shape1');
   $R->run(q`shape2 <- f$estimate[2]`);
   my $shape2 = $R->value('shape2');
   $R->run(q`shape1_sd <- f$sd[1]`);
   my $shape1_sd = $R->value('shape1_sd');
   $R->run(q`shape2_sd <- f$sd[2]`);
   my $shape2_sd = $R->value('shape2_sd');
   $R->run(q`chisqpvalue <- g$chisqpvalue`);
   my $chisqpvalue = $R->value('chisqpvalue');
   my $chisqtest = $chisqpvalue < 0.05 ? 'rejected' : 'not rejected';
   return $shape1, $shape1_sd, $shape2, $shape2_sd, $chisqtest;
}


sub fit_uniform {
   # Try to fit a uniform distribution to a series of data points using a maximum
   # goodness of fit method. Return the minimum, maximum and the results of
   # the Chi square statistics. Optional optimization starting values for the
   # mean and standard deviation can be specified.
   my ($values, $start_min, $start_max) = @_;
   my $start_params = '';
   if ( (defined $start_min) or (defined $start_max) ) {
      $start_params .= ', start=list(';
      if (defined $start_min) {
         $start_params .= "min=$start_min, ";
      }
      if (defined $start_max) {
         $start_params .= "max=$start_max, ";
      }
      $start_params =~ s/, $//;
      $start_params .= ')';
   }
   my $x = join ', ', @$values;
   my $out = '';
   my $R = Statistics::R->new();
   $R->startR;
   $R->run(q`library(fitdistrplus)`);
   $R->run(qq`x <- c($x)`);
   $R->run(qq`f <- fitdist(x, distr="unif", method="mge", gof="CvM"$start_params)`);
   $R->run(q`g <- gofstat(f)`);
   $R->run(q`min <- f$estimate[1]`);
   my $min = $R->value('min');
   $R->run(q`max <- f$estimate[2]`);
   my $max = $R->value('max');
   $R->run(q`chisqpvalue <- g$chisqpvalue`);
   my $chisqpvalue = $R->value('chisqpvalue');
   my $chisqtest = $chisqpvalue < 0.05 ? 'rejected' : 'not rejected';
   $R->stopR();
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
   my $x = join ', ', @$values;
   my $out = '';
   my $R = Statistics::R->new();
   $R->startR;
   $R->run(q`library(fitdistrplus)`);
   $R->run(qq`x <- c($x)`);
   $R->run(qq`f <- fitdist(x, distr="norm", method="mle"$start_params)`);
   $R->run(q`g <- gofstat(f)`);
   $R->run(q`mean <- f$estimate[1]`);
   my $mean = $R->value('mean');
   $R->run(q`sd <- f$estimate[2]`);
   my $sd = $R->value('sd');
   $R->run(q`mean_sd <- f$sd[1]`);
   my $mean_sd = $R->value('mean_sd');
   $R->run(q`sd_sd <- f$sd[2]`);
   my $sd_sd = $R->value('sd_sd');
   $R->run(q`adtest <- g$adtest`);
   my $adtest = $R->value('adtest');
   $R->run(q`cvmtest <- g$cvmtest`);
   my $cvmtest = $R->value('cvmtest');
   return $mean, $mean_sd, $sd, $sd_sd, $cvmtest, $adtest;
}


#------------------------------------------------------------------------------#

package Statistics::R;


sub run {
   my ($self, $cmd) = @_;
   # Execute command
   $self->send($cmd);
   # Read command output
   my $output = $self->read;
   # Check for errors
   if ($output =~ m/^<.*>$/) {
      die "Error: There was a problem running the R command\n".
          "   $cmd\n".
          "Because\n".
          "   $output\n";
   }
   return $output;
}


sub value {
   # Get the value of a variable through a Statistics::R object
   my ($self, $varname) = @_;
   $self->send(qq`print($varname)`);
   my $string = $self->read;
   my $value;
   if ($string eq 'NULL') {
      $value = undef;
   } elsif ($string =~ m/^\s*\[\d+\]/) {
      # String look like: ' [1]  6.4 13.3  4.1  1.3 14.1 10.6  9.9  9.6 15.3
      # [16]  5.2 10.9 14.4'
      my @lines = split /\n/, $string;
      for (my $i = 0; $i < scalar @lines; $i++) {
         $lines[$i] =~ s/^\s*\[\d+\] //;
      }
      $value = join ' ', @lines;
      # may need to split array and remove quotes
      $value =~ s/^"(.*)"$/$1/;
   } else {
      my @lines = split /\n/, $string;
      if (scalar @lines == 2) {
         # String looks like: '    mean 
         # 10.41111 '
         # Extract value from second line
         $value = $lines[1];
         $value =~ s/^\s*(\S+)\s*$/$1/;
      } else {
         die "Error: Don't know how to handle this R output\n$string\n";
      }
   }
   return $value;
}


1;
