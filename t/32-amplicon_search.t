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

ok $search = Grinder::AmpliconSearch->new(
   -template    => $seq,
   -primer_file => data('amplicon_database.fa'),
);



done_testing();
