#! perl

use strict;
use warnings;
use Test::More;
use t::TestUtils;
use Bio::PrimarySeq;


use_ok 'Grinder::KmerCollection';
my ($col, $seq1, $seq2, $by_kmer, $by_seq, $file, $sources, $counts, $kmers, $pos);


# 
$seq1 = Bio::PrimarySeq->new(
  -id => 'seq1',
  -seq => 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
);

$seq2 = Bio::PrimarySeq->new(
  -id => 'seq4',
  -seq => 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGCCCCCCCC'
);


ok $col = Grinder::KmerCollection->new( -k => 8 );

isa_ok $col, 'Grinder::KmerCollection';
is $col->k, 8;

ok $col->add_seqs([$seq1]);
ok $col->add_seqs([$seq2]);

ok $by_kmer = $col->collection_by_kmer;
ok exists $by_kmer->{'AAAAAAAA'}->{'seq1'};
ok exists $by_kmer->{'AAAAAAAA'}->{'seq4'};
ok exists $by_kmer->{'CCCCCCCC'}->{'seq4'};
ok exists $by_kmer->{'CCCCGGGG'}->{'seq4'};
ok exists $by_kmer->{'ACCCCCCC'}->{'seq4'};

ok $by_kmer = $col->collection_by_seq;
ok exists $by_kmer->{'seq1'}->{'AAAAAAAA'};
ok exists $by_kmer->{'seq4'}->{'AAAAAAAA'};
ok exists $by_kmer->{'seq4'}->{'CCCCCCCC'};
ok exists $by_kmer->{'seq4'}->{'CCCCGGGG'};
ok exists $by_kmer->{'seq4'}->{'ACCCCCCC'};

ok $col = $col->filter_rare(2);
isa_ok $col, 'Grinder::KmerCollection';

ok $by_kmer = $col->collection_by_kmer;
ok exists $by_kmer->{'AAAAAAAA'}->{'seq1'};
ok exists $by_kmer->{'AAAAAAAA'}->{'seq4'};
ok exists $by_kmer->{'CCCCCCCC'}->{'seq4'};
ok not exists $by_kmer->{'CCCCGGGG'};
ok not exists $by_kmer->{'ACCCCCCC'};

ok $by_kmer = $col->collection_by_seq;
ok exists $by_kmer->{'seq1'}->{'AAAAAAAA'};
ok exists $by_kmer->{'seq4'}->{'AAAAAAAA'};
ok exists $by_kmer->{'seq4'}->{'CCCCCCCC'};
ok not exists $by_kmer->{'seq4'}->{'CCCCGGGG'};
ok not exists $by_kmer->{'seq4'}->{'ACCCCCCC'};


ok $col = Grinder::KmerCollection->new( -k => 8, -seqs => [$seq1, $seq2] );

ok $col = $col->filter_shared(2);
isa_ok $col, 'Grinder::KmerCollection';

ok $by_kmer = $col->collection_by_kmer;
ok exists $by_kmer->{'AAAAAAAA'}->{'seq1'};
ok exists $by_kmer->{'AAAAAAAA'}->{'seq4'};
ok not exists $by_kmer->{'CCCCCCCC'};
ok not exists $by_kmer->{'CCCCGGGG'};
ok not exists $by_kmer->{'ACCCCCCC'};

ok $by_kmer = $col->collection_by_seq;
ok exists $by_kmer->{'seq1'}->{'AAAAAAAA'};
ok exists $by_kmer->{'seq4'}->{'AAAAAAAA'};
ok not exists $by_kmer->{'seq4'}->{'CCCCCCCC'};
ok not exists $by_kmer->{'seq4'}->{'CCCCGGGG'};
ok not exists $by_kmer->{'seq4'}->{'ACCCCCCC'};

($kmers, $counts) = $col->counts();
is_deeply $kmers , ['AAAAAAAA'];
is_deeply $counts, [       74 ];

($sources, $counts) = $col->sources('AAAAAAAA');
is_deeply $sources, ['seq4', 'seq1'];
is_deeply $counts , [    1 ,    73 ];

($sources, $counts) = $col->sources('ZZZZZZZZ');
is_deeply $sources, [];
is_deeply $counts , [];

($kmers, $counts) = $col->kmers('seq1');
is_deeply $kmers , ['AAAAAAAA'];
is_deeply $counts, [       73 ];

($kmers, $counts) = $col->kmers('seq4');
is_deeply $kmers , ['AAAAAAAA'];
is_deeply $counts, [        1 ];

($kmers, $counts) = $col->kmers('asdf');
is_deeply $kmers , [];
is_deeply $counts, [];

$pos = $col->positions('AAAAAAAA', 'seq1');
is_deeply $pos, [1..73];

$pos = $col->positions('AAAAAAAA', 'seq4');
is_deeply $pos, [34];

$pos = $col->positions('CCCCGGGG', 'seq4');
is_deeply $pos, [];

$pos = $col->positions('AAAAAAAA', 'seq3');
is_deeply $pos, [];

$file = data('kmers.fa');
ok $col = Grinder::KmerCollection->new( -k => 8, -file => $file );

ok $col = Grinder::KmerCollection->new( -k => 8, -file => $file )->filter_rare(2);
isa_ok $col, 'Grinder::KmerCollection';


done_testing;
