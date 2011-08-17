#!perl -T

use strict;
use warnings;
use Test::More tests => 1204;
use File::Spec::Functions;

use Grinder;
my ($factory, $read, $nof_reads);

# Template with several matching amplicons and forward and reverse primers

ok $factory = Grinder->new(
   -genome_file     => catfile(qw{t data multiple_amplicon_database.fa}),
   -forward_reverse => catfile(qw{t data forward_reverse_primers.fa})   ,
   -length_bias     => 0                                                ,
   -unidirectional  => 1                                                ,
   -read_dist       => 100                                              ,
   -total_reads     => 100                                              ,
), 'Forward and reverse primers';

$nof_reads = 0;
while ( $read = $factory->next_read ) {
   $nof_reads++;
   ok_read_forward_reverse($read, 1, $nof_reads);
};
is $nof_reads, 100;


# Template with several matching amplicons and forward primer

ok $factory = Grinder->new(
   -genome_file     => catfile(qw{t data multiple_amplicon_database.fa}),
   -forward_reverse => catfile(qw{t data forward_primer.fa})            ,
   -length_bias     => 0                                                ,
   -unidirectional  => 1                                                ,
   -read_dist       => 100                                              ,
   -total_reads     => 100                                              ,
), 'Forward primer only';

$nof_reads = 0;
while ( $read = $factory->next_read ) {
   $nof_reads++;
   ok_read_forward_only($read, 1, $nof_reads);
};
is $nof_reads, 100;


sub ok_read_forward_reverse {
   my ($read, $req_strand, $nof_reads) = @_;
   is ref($read), 'Bio::Seq::SimulatedRead';
   my $source = $read->reference->id;
   ok ($source =~ m/^seq\d+$/);
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
   is ref($read), 'Bio::Seq::SimulatedRead';
   my $source = $read->reference->id;
   ok ($source =~ m/^seq\d+$/);
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
