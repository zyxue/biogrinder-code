#! perl

use strict;
use warnings;
use Test::More;
use t::TestUtils;
use Grinder;



my ($factory, $read, $nof_reads, $nof_chimeras, $nof_regulars);


# No Chimeras

ok $factory = Grinder->new(
   -reference_file  => data('amplicon_database.fa')      ,
   -forward_reverse => data('forward_reverse_primers.fa'),
   -length_bias     => 0                                 ,
   -unidirectional  => 1                                 ,
   -chimera_perc    => 0                                 ,
   -total_reads     => 100                               ,
), 'No chimeras';

while ( $read = $factory->next_read ) {
   is nof_references($read->desc), 1;
   # Remove forward and reverse primer
   my $seq = $read->seq;
   $seq = remove_primers($seq, 'AAACT.AAA.GAATTG.CGG', 'G.ACACACCGCCCGT');
   # Now the amplicon is simply long homopolymeric sequences
   like $seq, qr/^(a+|c+|g+|t+)+$/;
}


# 50% chimeras

ok $factory = Grinder->new(
   -reference_file  => data('amplicon_database.fa')      ,
   -forward_reverse => data('forward_reverse_primers.fa'),
   -length_bias     => 0                                 ,
   -unidirectional  => 1                                 ,
   -chimera_perc    => 50                                ,
   -total_reads     => 100                               ,
), '50% chimeras';

while ( $read = $factory->next_read ) {
   # Remove forward and reverse primer
   $nof_chimeras += nof_references($read->desc);
   $nof_regulars += nof_references($read->desc);
}
between_ok( $nof_chimeras / $nof_regulars, 0.9, 1.1 );


# 100% chimeras

ok $factory = Grinder->new(
   -reference_file  => data('amplicon_database.fa')      ,
   -forward_reverse => data('forward_reverse_primers.fa'),
   -length_bias     => 0                                 ,
   -unidirectional  => 1                                 ,
   -chimera_perc    => 100                               ,
   -total_reads     => 100                               ,
), '100% chimeras';

while ( $read = $factory->next_read ) {
   # Remove forward and reverse primer
   is nof_references($read->desc), 2;
}

done_testing();





sub remove_primers {
   my ($seq, $forward_re, $reverse_re) = @_;
   $seq =~ s/$forward_re//i;
   $seq =~ s/$reverse_re//i;
   return $seq;
}


sub nof_references {
   my ($desc) = @_;
   $desc =~ m/reference=(\S+)/;
   my $refs = $1;
   my $nof_refs = scalar split ',', $refs;
   return $nof_refs;
}


sub matches_ref {
   my ($read) = @_;
   my $read_seq = $read->seq;
   my $ref_seq  = $read->reference->seq;
   my $matches  = 0;
   if ($ref_seq =~ m/$read_seq/) {
      $matches = 1;
      print "$read_seq\nmatches\n$ref_seq\n\n";
   }
   return $matches;
}
