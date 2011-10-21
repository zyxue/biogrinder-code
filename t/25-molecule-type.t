#! perl

use strict;
use warnings;
use Test::More;
use t::TestUtils;
use Grinder;

plan tests => 305;


my ($factory, $read, $nof_reads, $nof_chimeras, $nof_regulars);


# DNA database

ok $factory = Grinder->new(
   -reference_file  => data('database_dna.fa')   ,
   -total_reads     => 100                       ,
), 'DNA';

is $factory->{alphabet}, 'dna';


# RNA

ok $factory = Grinder->new(
   -reference_file  => data('database_rna.fa')   ,
   -total_reads     => 100                       ,
), 'RNA';

is $factory->{alphabet}, 'rna';


# Protein

ok $factory = Grinder->new(
   -reference_file  => data('database_protein.fa')   ,
   -total_reads     => 100                       ,
), 'Protein';

is $factory->{alphabet}, 'protein';


# Protein

ok $factory = Grinder->new(
   -reference_file  => data('database_mixed.fa')   ,
   -total_reads     => 100                       ,
), 'Protein';

is $factory->{alphabet}, 'protein';
