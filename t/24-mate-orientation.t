#! perl

use strict;
use warnings;
use Test::More;
use t::TestUtils;
use Grinder;



my ($factory, $read, $nof_reads);


# FR-oriented mates

ok $factory = Grinder->new(
   -reference_file   => data('oriented_database.fa'),
   -total_reads      => 100                         ,
   -read_dist        => 80                          ,
   -insert_dist      => 240                         ,
   -unidirectional   => +1                          ,
   -mate_orientation => 'FR'                        ,
), 'FR-oriented mates';

$nof_reads = 0;
while ( $read = $factory->next_read ) {
   $nof_reads++;
   if ( is_first_mate($nof_reads) ) {
      ok_mate_1_F($read);
   } else {
      ok_mate_2_R($read);
   }
};
is $nof_reads, 100;

ok $factory = Grinder->new(
   -reference_file   => data('oriented_database.fa'),
   -total_reads      => 100                         ,
   -read_dist        => 80                          ,
   -insert_dist      => 240                         ,
   -unidirectional   => -1                          ,
   -mate_orientation => 'FR'                        ,
);

$nof_reads = 0;
while ( $read = $factory->next_read ) {
   $nof_reads++;
   if ($nof_reads % 2 == 1) {
      ok_mate_2_R($read);
   } else {
      ok_mate_1_F($read);
   }
};
is $nof_reads, 100;


# FF-oriented mates

ok $factory = Grinder->new(
   -reference_file   => data('oriented_database.fa'),
   -total_reads      => 100                         ,
   -read_dist        => 80                          ,
   -insert_dist      => 240                         ,
   -unidirectional   => +1                          ,
   -mate_orientation => 'FF'                        ,
), 'FF-oriented mates';

$nof_reads = 0;
while ( $read = $factory->next_read ) {
   $nof_reads++;
   if ( is_first_mate($nof_reads) ) {
      ok_mate_1_F($read);
   } else {
      ok_mate_2_F($read);
   }
};
is $nof_reads, 100;

ok $factory = Grinder->new(
   -reference_file   => data('oriented_database.fa'),
   -total_reads      => 100                         ,
   -read_dist        => 80                          ,
   -insert_dist      => 240                         ,
   -unidirectional   => -1                          ,
   -mate_orientation => 'FF'                        ,
);

$nof_reads = 0;
while ( $read = $factory->next_read ) {
   $nof_reads++;
   if ( is_first_mate($nof_reads) ) {
      ok_mate_2_R($read);
   } else {
      ok_mate_1_R($read);
   }
};
is $nof_reads, 100;


# RF-oriented mates

ok $factory = Grinder->new(
   -reference_file   => data('oriented_database.fa'),
   -total_reads      => 100                         ,
   -read_dist        => 80                          ,
   -insert_dist      => 240                         ,
   -unidirectional   => +1                          ,
   -mate_orientation => 'RF'                        ,
), 'RF-oriented mates';

$nof_reads = 0;
while ( $read = $factory->next_read ) {
   $nof_reads++;
   if ( is_first_mate($nof_reads) ) {
      ok_mate_1_R($read);
   } else {
      ok_mate_2_F($read);
   }
};
is $nof_reads, 100;

ok $factory = Grinder->new(
   -reference_file   => data('oriented_database.fa'),
   -total_reads      => 100                         ,
   -read_dist        => 80                          ,
   -insert_dist      => 240                         ,
   -unidirectional   => -1                          ,
   -mate_orientation => 'RF'                        ,
);

$nof_reads = 0;
while ( $read = $factory->next_read ) {
   $nof_reads++;
   if ( is_first_mate($nof_reads) ) {
      ok_mate_2_F($read);
   } else {
      ok_mate_1_R($read)
   }
};
is $nof_reads, 100;


# RR-oriented mates

ok $factory = Grinder->new(
   -reference_file   => data('oriented_database.fa'),
   -total_reads      => 100                         ,
   -read_dist        => 80                          ,
   -insert_dist      => 240                         ,
   -unidirectional   => +1                          ,
   -mate_orientation => 'RR'                        ,
), 'RR-oriented mates';

$nof_reads = 0;
while ( $read = $factory->next_read ) {
   $nof_reads++;
   if ( is_first_mate($nof_reads) ) {
      ok_mate_1_R($read);
   } else {
      ok_mate_2_R($read);
   }
};
is $nof_reads, 100;

ok $factory = Grinder->new(
   -reference_file   => data('oriented_database.fa'),
   -total_reads      => 100                         ,
   -read_dist        => 80                          ,
   -insert_dist      => 240                         ,
   -unidirectional   => -1                          ,
   -mate_orientation => 'RR'                        ,
);

$nof_reads = 0;
while ( $read = $factory->next_read ) {
   $nof_reads++;
   if ( is_first_mate($nof_reads) ) {
      ok_mate_2_F($read);
   } else {
      ok_mate_1_F($read);
   }
};
is $nof_reads, 100;


done_testing();



sub ok_mate_1_F {
   my ($read) = @_;
   is $read->strand, 1;
   is $read->seq, 'CCCaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa';
}


sub ok_mate_2_F {
   my ($read) = @_;
   is $read->strand, 1;
   is $read->seq, 'aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaTTT';
}


sub ok_mate_1_R {
   my ($read) = @_;
   is $read->strand, -1;
   is $read->seq, 'tttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttGGG';
}

sub ok_mate_2_R {
   my ($read) = @_;
   is $read->strand, -1;
   is $read->seq, 'AAAttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt';
}


sub ok_mate {
   my ($read, $req_strand, $nof_reads) = @_;
   isa_ok $read, 'Bio::Seq::SimulatedRead';

   my $source = $read->reference->id;
   my $strand = $read->strand;
   if (not defined $req_strand) {
      $req_strand = $strand;
   } else {
      is $strand, $req_strand;
   }
   my $letters;
   if ( $source eq 'seq1' ) {
      $letters = 'a';
   } elsif ( $source eq 'seq2' ) {
      $letters = 'c';
   } elsif ( $source eq 'seq3' ) {
      $letters = 'g';
   } elsif ( $source eq 'seq4' ) {
      $letters = 't';
   } elsif ( $source eq 'seq5' ) {
      $letters = 'atg';
   }
   if ( $req_strand == -1 ) { # Take the reverse complement
      $letters = Bio::PrimarySeq->new( -seq => $letters )->revcom->seq;
   };
   like $read->seq, qr/[$letters]+/;
   my $id = round($nof_reads/2).'/'.($nof_reads%2?1:2);
   is $read->id, $id;
   is $read->length, 48;
}


sub is_first_mate {
   # Given the number of a read, is this one the first mate in a pair?
   # This assume that your count is 1-based and that first and second mate
   # alternate, starting with a first mate.
   my ($read_num) = @_;
   return ($read_num % 2) ? 1 : 0;
}
