#! perl

use strict;
use warnings;
use Test::More;
use t::TestUtils;
use Grinder;

plan tests => 1204;


my ($factory, $read, $nof_reads);

# Template with several matching amplicons and forward and reverse primers

ok $factory = Grinder->new(
   -reference_file  => data('multiple_amplicon_database.fa'),
   -forward_reverse => data('forward_reverse_primers.fa')   ,
   -length_bias     => 0                                    ,
   -unidirectional  => 1                                    ,
   -read_dist       => 100                                  ,
   -total_reads     => 100                                  ,
), 'Forward and reverse primers';

$nof_reads = 0;
while ( $read = $factory->next_read ) {
   $nof_reads++;
   ok_read_forward_reverse($read, 1, $nof_reads);
};
is $nof_reads, 100;


# Template with several matching amplicons and forward primer

ok $factory = Grinder->new(
   -reference_file  => data('multiple_amplicon_database.fa'),
   -forward_reverse => data('forward_primer.fa')            ,
   -length_bias     => 0                                    ,
   -unidirectional  => 1                                    ,
   -read_dist       => 100                                  ,
   -total_reads     => 100                                  ,
), 'Forward primer only';

$nof_reads = 0;
while ( $read = $factory->next_read ) {
   $nof_reads++;
   ok_read_forward_only($read, 1, $nof_reads);
};
is $nof_reads, 100;


sub ok_read_forward_reverse {
   my ($read, $req_strand, $nof_reads) = @_;
   isa_ok $read, 'Bio::Seq::SimulatedRead';
   like $read->reference->id, qr/^seq\d+$/;
   my $strand = $read->strand;
   if (not defined $req_strand) {
      $req_strand = $strand;
   } else {
      is $strand, $req_strand;
   }
   my $readseq = $read->seq;
   ok (($readseq eq 'AAACTUAAAGGAATTGACGGaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaGTACACACCGCCCGT')
     or ($readseq eq 'AAACTTAAAGGAATTGRCGGttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttGTACACACCGCCCGT'));
   is $read->id, $nof_reads;
   is $read->length, 95;
}

sub ok_read_forward_only {
   my ($read, $req_strand, $nof_reads) = @_;
   isa_ok $read, 'Bio::Seq::SimulatedRead';
   like $read->reference->id, qr/^seq\d+$/;
   my $strand = $read->strand;
   if (not defined $req_strand) {
      $req_strand = $strand;
   } else {
      is $strand, $req_strand;
   }
   my $readseq = $read->seq;
   ok ( ($readseq eq 'AAACTUAAAGGAATTGACGGaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaGTACACACCGCCCGTccccc' )
     or ($readseq eq 'AAACTTAAAGGAATTGRCGGttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttGTACACACCGCCCGTggggg')
     or ($readseq eq 'AAACTTAAAGGAATTGRCGGttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttGTACACACCGCCCGT'     ));
   is $read->id, $nof_reads;
   my $readlength = $read->length;
   ok ( ($readlength == 95) or ($readlength == 100) );
}
