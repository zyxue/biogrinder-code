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
   is nof_references($read), 3;
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
   is nof_references($read), 4;
}


# 100% chimeras (bimeras, trimeras, quadrameras)

ok $factory = Grinder->new(
   -reference_file  => data('kmers.fa'),
   -length_bias     => 0                                 ,
   -unidirectional  => 1                                 ,
   -chimera_perc    => 100                               ,
   -chimera_dist    => (1, 1, 1)                         ,
   -chimera_kmer    => 8                                 ,
   -total_reads     => 1000                              ,
), '100% chimeras (bimeras, trimeras, quadrameras)';

while ( $read = $factory->next_read ) {
   # Remove forward and reverse primer
   my $nof_refs = nof_references($read);
   $chim_sizes{$nof_refs}++;
   between_ok( $nof_refs, 2, 4 );
}
between_ok( $chim_sizes{2}, 333.3 * 0.9, 333.3 * 1.1 );
between_ok( $chim_sizes{3}, 333.3 * 0.9, 333.3 * 1.1 );
between_ok( $chim_sizes{4}, 333.3 * 0.9, 333.3 * 1.1 );


done_testing();
