#! perl

use strict;
use warnings;
use Test::More;
use Test::Warn;
use t::TestUtils;
use Grinder;



my ($factory, $read, $nof_reads);

ok $factory = Grinder->new(
   -reference_file => data('shotgun_database_extended.fa'),
   -total_reads    => 101                                 ,
   -read_dist      => 48                                  ,
   -insert_dist    => 250                                 ,
), 'Mate pairs';

{
  my @warnings;
  local $SIG{__WARN__} = sub {
     push @warnings, @_;
  };
  $factory->next_lib;
  like $warnings[0], qr{.*added a read.*}i;
}

$nof_reads = 0;
while ( $read = $factory->next_read ) {
   $nof_reads++;
   ok_mate($read, undef, $nof_reads);
};
is $nof_reads, 102;



# Coverage fold

ok $factory = Grinder->new(
   -reference_file => data('shotgun_database_extended.fa'),
   -read_dist      => 48                                  ,
   -coverage_fold  => 6.04                                ,
   -insert_dist    => 250                                 ,
), 'Coverage fold';

ok $factory->next_lib;

$nof_reads = 0;
while ( $read = $factory->next_read ) {
   $nof_reads++;
};
is $nof_reads, 112;


done_testing();



sub ok_mate {
   my ($read, $req_strand, $nof_reads) = @_;
   isa_ok $read, 'Bio::Seq::SimulatedRead';
   my $source = $read->reference->id;
   my $strand = $read->strand;
   if (not defined $req_strand) {
      $req_strand = $strand;
   } else {
      is $strand, $req_strand;
   }
   my $letters;
   if ( $source eq 'seq1' ) {
      $letters = 'a';
   } elsif ( $source eq 'seq2' ) {
      $letters = 'c';
   } elsif ( $source eq 'seq3' ) {
      $letters = 'g';
   } elsif ( $source eq 'seq4' ) {
      $letters = 't';
   } elsif ( $source eq 'seq5' ) {
      $letters = 'atg';
   } elsif ( $source eq 'seq7' ) {
      $letters = 'a';
   }
   if ( $req_strand == -1 ) { # Take the reverse complement
      $letters = Bio::PrimarySeq->new( -seq => $letters )->revcom->seq;
   };
   like $read->seq, qr/[$letters]+/;
   my $id = round($nof_reads/2).'/'.($nof_reads%2?1:2);
   is $read->id, $id;
   if ($source eq 'seq7') {
      is $read->length, 1;
   } else {
      is $read->length, 48;
   }
}

