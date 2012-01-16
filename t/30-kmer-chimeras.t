#! perl

use strict;
use warnings;
use Test::More;
use t::TestUtils;
use Grinder;



my ($factory, $read, $nof_reads);
my %chim_sizes;
my %refs;


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

%refs = ();
while ( $read = $factory->next_read ) {
   my @refs = get_references($read);
   is scalar @refs, 2;
   for my $ref (@refs) {
     $refs{$ref}++;
   }
}
ok exists $refs{'seq1'};
ok exists $refs{'seq2'};
ok exists $refs{'seq3'};
ok exists $refs{'seq4'};
ok not exists $refs{'seq5'};


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

%refs = ();
while ( $read = $factory->next_read ) {
   my @refs = get_references($read);
   is scalar @refs, 3;
   for my $ref (@refs) {
     $refs{$ref}++;
   }
}
ok exists $refs{'seq1'};
ok exists $refs{'seq2'};
ok exists $refs{'seq3'};
ok exists $refs{'seq4'};
ok not exists $refs{'seq5'};


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

%refs = ();
while ( $read = $factory->next_read ) {
   my @refs = get_references($read);
   is scalar @refs, 4;
   for my $ref (@refs) {
     $refs{$ref}++;
   }
}
ok exists $refs{'seq1'};
ok exists $refs{'seq2'};
ok exists $refs{'seq3'};
ok exists $refs{'seq4'};
ok not exists $refs{'seq5'};


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

%refs = ();
while ( $read = $factory->next_read ) {
   # Remove forward and reverse primer
   my @refs = get_references($read);
   my $nof_refs = scalar @refs;
   $chim_sizes{$nof_refs}++;
   between_ok( $nof_refs, 2, 4 );
   for my $ref (@refs) {
     $refs{$ref}++;
   }
}
between_ok( $chim_sizes{2}, 333.3 * 0.9, 333.3 * 1.1 );
between_ok( $chim_sizes{3}, 333.3 * 0.9, 333.3 * 1.1 );
between_ok( $chim_sizes{4}, 333.3 * 0.9, 333.3 * 1.1 );
ok exists $refs{'seq1'};
ok exists $refs{'seq2'};
ok exists $refs{'seq3'};
ok exists $refs{'seq4'};
ok not exists $refs{'seq5'};





# using sequence 2, 3, 4, 

done_testing();
