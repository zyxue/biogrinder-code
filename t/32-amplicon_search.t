#! perl

use strict;
use warnings;
use Test::More;
use t::TestUtils;
use Bio::PrimarySeq;

use_ok 'Grinder::AmpliconSearch';

my ($search);

ok $search = Grinder::AmpliconSearch->new();
isa_ok $search, 'Grinder::AmpliconSearch';

my $seq = Bio::PrimarySeq->new(
   -seq => 'AAACTTAAAGGAATTGACGGaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaGTACACACCGCCCGT'
);

my $forward = Bio::PrimarySeq->new(
   -seq => 'AAACTTAAAGGAATTGACGG'
);

my $reverse = Bio::PrimarySeq->new(
   -seq => 'GTACACACCGCCCGT'
);




ok $search = Grinder::AmpliconSearch->new(
   -template       => $seq,
   -forward_primer => $forward,
);
is $search->get_template->seq, 'AAACTTAAAGGAATTGACGGaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaGTACACACCGCCCGT';
is $search->get_forward_primer->seq, 'AAACTTAAAGGAATTGACGG';
is $search->get_reverse_primer, undef;




ok $search = Grinder::AmpliconSearch->new(
   -template       => $seq,
   -forward_primer => $forward,
   -reverse_primer => $reverse,
);
is $search->get_template->seq, 'AAACTTAAAGGAATTGACGGaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaGTACACACCGCCCGT';
is $search->get_forward_primer->seq, 'AAACTTAAAGGAATTGACGG';
is $search->get_reverse_primer->seq, 'GTACACACCGCCCGT';




ok $search = Grinder::AmpliconSearch->new(
   -template    => $seq,
   -primer_file => data('forward_reverse_primers.fa'),
);
is $search->get_template->seq, 'AAACTTAAAGGAATTGACGGaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaGTACACACCGCCCGT';
is $search->get_forward_primer->seq, 'AAACTYAAAKGAATTGRCGG';
is $search->get_reverse_primer->seq, 'ACGGGCGGTGTGTRC';

### use Data::Dumper; print Dumper($search);

done_testing();
