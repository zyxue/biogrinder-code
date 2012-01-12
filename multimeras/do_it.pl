#! /usr/bin/env perl

use base 'Grinder';

use strict;
use warnings;
use Data::Dumper;
$Data::Dumper::Indent = 1;
#$Data::Dumper::Maxdepth = 5;
use Grinder;



use lib '.';
use KmerCollection;

main();


sub main {

   my $file = 'multimeras.fa';
   my $k = 8;

   my $self = Grinder->new( -reference_file => $file );
   my @seqs = values %{$self->{database}->{db}};
   my $kmer_col = KmerCollection->new(-k=>$k,-seqs=>\@seqs)->filter_shared(2);
   #print Dumper($kmer_col->collection);
   my ($kmers, $freqs) = $kmer_col->counts(1);
   #_dump_arrays($kmers, $freqs);
   my $cdf = $self->proba_cumul($freqs);

   for (1 .. 10) {

      # Size of multimera
      my $max_m = 3;

      my $multimera = rand_multimera($self, $max_m, $k, $kmers, $cdf, $kmer_col );

   }

   return 1;
}


sub rand_multimera {
   
   my ($self, $max_m, $k, $kmers, $cdf, $kmer_col) = @_;

   # Pick a random kmer
   my $kmer = $$kmers[Grinder::rand_weighted($cdf)];
   print "Got random kmer $kmer\n";

   # Pick two sequences with that kmer
   my ($seq1, $pos1, $seq2, $pos2) = rand_bimera( $self, $k, $kmer, $kmer_col );

   for (my $m = 3; $m <= $max_m; $m++) {
        # Given a starting sequence, find another one (with a matching kmer) to
        # add to the multimera
        
   }
  
   ### sort positions in order
   ### join sequences
   ###return $seq;
}


sub rand_kmer_source {
   # 
   my ($self, $kmer, $kmer_col) = @_;
   my ($sources, $freqs) = $kmer_col->sources($kmer, 1);
   my $cdf = $self->proba_cumul($freqs);
   my $source = $$sources[Grinder::rand_weighted($cdf)];
   
}


sub rand_kmer_start {
   #
   my ($self, $kmer, $source, $kmer_col) = @_;
   my $kmer_starts = $kmer_col->positions($kmer, $source);
   return $kmer_starts->[ int rand scalar @$kmer_starts ];
}


sub rand_bimera {
   my ($self, $k, $kmer, $kmer_col) = @_;

   # Pick two sequences with that kmer
   my $seq1 = rand_kmer_source( $self, $kmer, $kmer_col );
   my $pos1 = rand_kmer_start( $self, $kmer, $seq1, $kmer_col );

   my $seq2 = rand_kmer_source( $self, $kmer, $kmer_col );
   my $pos2 = rand_kmer_start( $self, $kmer, $seq2, $kmer_col ) + $k;

   return $seq1, $pos1, $seq2, $pos2;
}


sub _dump_arrays {
   my ($arr1, $arr2) = @_;
   for (my $i = 0; $i < scalar @$arr1; $i++) {
      my $val1 = $$arr1[$i];
      my $val2 = $$arr2[$i];
      print "$val1  $val2\n";
   }
   return 1;
}






























