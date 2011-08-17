#!perl 

use strict;
use warnings;
use Test::More tests => 27;
use File::Spec::Functions;

use Grinder;

my ($factory, $read);


# Exclude forbidden characters

ok $factory = Grinder->new(
   -genome_file    => catfile(qw{t data dirty_database.fa}),
   -read_dist      => 80                                   ,
   -random_seed    => 1233567890                           ,
   -total_reads    => 10                                   ,
), 'With dubious chars';

ok $read = $factory->next_read;
ok $read->seq =~ m/[N-]/i;


ok $factory = Grinder->new(
   -genome_file    => catfile(qw{t data dirty_database.fa}),
   -exclude_chars  => 'n-'                                 , # case independent
   -read_dist      => 30                                   ,
   -random_seed    => 1233567890                           ,
   -total_reads    => 10                                   ,
), 'Exclude chars';

while ( $read = $factory->next_read ) {
  ok $read->seq !~ m/[N-]/i;
}


ok $factory = Grinder->new(
   -genome_file    => catfile(qw{t data dirty_database.fa}),
   -exclude_chars  => 'N-'                                 ,
   -read_dist      => 71                                   ,
   -random_seed    => 1233567890                           ,
   -total_reads    => 10                                   ,
), 'Cannot generate read';

eval { $read = $factory->next_read };
ok $@ =~ m/error/i;


# Remove forbidden characters

ok $factory = Grinder->new(
   -genome_file    => catfile(qw{t data dirty_database.fa}),
   -delete_chars   => 'N-'                                 ,
   -read_dist      => 70                                   ,
   -random_seed    => 1233567890                           ,
   -total_reads    => 10                                   ,
), 'Delete chars';

while ( $read = $factory->next_read ) {
  ok $read->seq !~ m/[N-]/i;
}

