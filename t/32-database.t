#! perl

use strict;
use warnings;
use t::TestUtils;
use Test::More;
use_ok 'Grinder::Database';

my ($db, $seq);


# Test minium sequence length and forbidden characters

ok $db = Grinder::Database->new(
   -fasta_file => data('shotgun_database.fa'),
);
isa_ok $db, 'Grinder::Database';
is $db->get_minimum_length, 1;
is $db->get_delete_chars, '';
is_deeply [sort @{$db->get_ids}], ['seq1', 'seq2', 'seq3', 'seq4', 'seq5'];


ok $db = Grinder::Database->new(
   -fasta_file     => data('shotgun_database.fa'),
   -minimum_length => 200,
);
is $db->get_minimum_length, 200;
is_deeply $db->get_ids, ['seq1', 'seq2'];


ok $db = Grinder::Database->new(
   -fasta_file   => data('shotgun_database.fa'),
   -delete_chars => 'ac',
);
is $db->get_delete_chars, 'ac';
is_deeply [sort @{$db->get_ids}], ['seq3', 'seq4', 'seq5'];


# Test retrieving sequences and subsequences

is $db->get_seq('zzz'), undef;

ok $seq = $db->get_seq('seq5');
is $seq->id, 'seq5';
is $seq->seq, 'aaaaaaaaaattttttttttttttttttttttttttttttttttttttttttttttttttttttttttttgggggggggg';

ok $seq = $db->get_seq('seq5:2..11');
is $seq->id, 'seq5';
is $seq->seq, 'aaaaaaaaat';

ok $seq = $db->get_seq('seq5/-1');
is $seq->id, 'seq5';
is $seq->seq, 'ccccccccccaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaatttttttttt';

ok $seq = $db->get_seq('seq5:2..11/-1');
is $seq->id, 'seq5';
is $seq->seq, 'attttttttt';


#ok $db = Grinder::Database->new(
#   -fasta_file   => data('shotgun_database.fa'),
#   -unidirectional => -1,
#);
#is $db->get_unidirectional, -1;


#use Data::Dump qw(dump);
#print Data::Dump::dump($db);


### test alphabet like in t/25-molecule-type.t


#$db = Grinder::Database->new(
#   -fasta_file              => data('shotgun_database.fa'),
#   -unidirectional          => 
#   -forward_reverse_primers =>
#   -abundance_file          =>
#   -delete_chars            =>
#   -min_len                 => 1;
#);


# next seq for shotgun

done_testing();
