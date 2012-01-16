#! perl

use strict;
use warnings;
use Test::More;
use t::TestUtils;
use Grinder;



my ($factory, $read, $nof_reads, $nof_chimeras, $nof_regulars);
my %chim_sizes;


# Bimeras

ok $factory = Grinder->new(
   -reference_file  => data('kmers.fa'),
   -length_bias     => 0               ,
   -unidirectional  => 1               ,
   -chimera_perc    => 100             ,
   -chimera_dist    => (1)             ,
   -chimera_kmer    => 8               ,
   -total_reads     => 100             ,
), 'Bimeras';

while ( $read = $factory->next_read ) {
   is nof_references($read), 2;
}


# Trimeras

ok $factory = Grinder->new(
   -reference_file  => data('kmers.fa'),
   -length_bias     => 0               ,
   -unidirectional  => 1               ,
   -chimera_perc    => 100             ,
   -chimera_dist    => (0, 1)          ,
   -chimera_kmer    => 8               ,
   -total_reads     => 100             ,
), 'Trimeras';

while ( $read = $factory->next_read ) {
#   is nof_references($read->desc), 3;
}

ok $factory = Grinder->new(
   -reference_file  => data('kmers.fa'),
   -length_bias     => 0               ,
   -unidirectional  => 1               ,
   -chimera_perc    => 100             ,
   -chimera_dist    => (1)             ,
   -chimera_kmer    => 8               ,
   -total_reads     => 100             ,
), 'Bimeras';

while ( $read = $factory->next_read ) {
   is nof_references($read), 2;
}


# Quadrameras

ok $factory = Grinder->new(
   -reference_file  => data('kmers.fa'),
   -length_bias     => 0               ,
   -unidirectional  => 1               ,
   -chimera_perc    => 100             ,
   -chimera_dist    => (0, 0, 1)       ,
   -chimera_kmer    => 8               ,
   -total_reads     => 100             ,
), 'Quadrameras';

while ( $read = $factory->next_read ) {
#   is nof_references($read->desc), 4;
}

done_testing();
