package Grinder;

use 5.006;
use strict;
use File::Spec;
use Bio::SeqIO;
use Bio::Seq::SimulatedRead;
use Getopt::Euclid qw( :minimal_keys :defer );
use Math::Random::MT::Perl qw(srand rand);
our $VERSION = '0.3.6';

#---------- GRINDER POD DOC ---------------------------------------------------#

=head1 NAME

Grinder - a simulator of random shotgun and amplicon sequence libraries

=head1 DESCRIPTION

Grinder is a program to create artificial random shotgun and amplicon sequence
libraries based on reference sequences in a FASTA file. Features include:

=over

=item *

shotgun library or amplicon library

=item *

arbitrary read length distribution and number of reads

=item *

simulation of PCR and sequencing errors (chimeras, point mutations, homopolymers)

=item *

support for creating paired-end (mate pair) datasets

=item *

specific rank-abundance settings or manually given abundance for each genome

=item *

creation of datasets with a given richness (alpha diversity)

=item *

independent datasets can share a variable number of genomes (beta diversity)

=item *

modeling of the bias created by varying genome lengths or gene copy number

=item *

profile mechanism to store preferred options

=item *

API to automate the creation of a large number of simulated dataset

=back

Grinder can thus produce metagenomic, amplicon or shotgun sequence datasets
which can be used to test the accuracy of bioinformatic tools or help
decide between alternative sequencing methods in an experiment.

=head1 CITATION

If you use Grinder in your research, please cite:

   Angly FE, Willner D, Prieto-Davó A, Edwards RA, Schmieder R, et al. (2009) The
   GAAS Metagenomic Tool and Its Estimations of Viral and Microbial Average Genome
   Size in Four Major Biomes. PLoS Comput Biol 5(12): e1000593.
   L<http://dx.doi.org/10.1371/journal.pcbi.1000593>

=head1 VERSION

0.3.6

=head1 AUTHOR

Florent Angly <florent.angly@gmail.com>

=head1 INSTALLATION

You need to install these dependencies first:

=over

=item *

Perl (L<http://www.perl.com/download.csp>)

=back

The following Perl modules are dependencies that will be installed automatically
for you:

=over


=item *

Bio::SeqIO (from Bioperl)

=item *

Bio::Seq::SimulatedRead (from Bioperl but included here because it is so recent)

=item *

Getopt::Euclid

=item *

Math::Random::MT::Perl

=back

To install Grinder globally on your system, run the following commands in a
terminal or command prompt:

On Linux, Unix, MacOS:

   perl Makefile.PL
   make
   make test

And finally, with administrator privileges:

   make install

On Windows, run the same commands but with nmake instead of make.

If you do not have administrator rights and want to install the module locally,
try something along these lines:

   perl Makefile.PL INSTALL_BASE=/home/fangly/bin/perl


=head1 RUNNING GRINDER

After installation, you can run Grinder using a command-line interface (CLI) or
an application programming interface (API). To get the usage of the CLI, type:

  Grinder --help

If you are interested in running Grinder from within other Perl programs, see
the documentation of the Grinder API:

  perldoc Grinder

The 'utils' folder included in the Grinder package contains utilities:

=over

=item average genome size:

This calculates the average genome size (in bp) of a simulated random library
produces by Grinder.

=item change_paired_read_orientation:

This reverses the orientation of each second mate-pair read (ID endind in /2)
in a FASTA file.

=back

=head1 REFERENCE SEQUENCE DATABASE

A variety of FASTA databases can be used as input for Grinder. For example, the
GreenGenes database (L<http://greengenes.lbl.gov/Download/Sequence_Data/Fasta_data_files/Isolated_named_strains_16S_aligned.fasta>)
contains over 180,000 16S rRNA clone sequences from various species which would
be appropriate to produce a 16S amplicon dataset. A set of over 41,000 OTU
representative sequences and their affiliation in seven different taxonomic
sytems can also be used for the same purpose (L<http://greengenes.lbl.gov/Download/OTUs/gg_otus_6oct2010/rep_set/gg_97_otus_6oct2010.fasta>
and L<http://greengenes.lbl.gov/Download/OTUs/gg_otus_6oct2010/taxonomies/>).
While 16S rRNA is a popular gene, datasets containing any type of gene could be used
in the same fashion to generate simulated amplicon datasets, provided appropriate
primers are used.

The >2,400 curated microbial genome sequences in the NCBI RefSeq collection
(L<ftp://ftp.ncbi.nih.gov/refseq/release/microbial/>) would also be suitable for
producing 16S rRNA simulated datasets (using the adequate primers). However, the
lower diversity of this database compared to the previous two makes it more
appropriate for producing artificial microbial metagenomes. Individual genomes
from this database are also very suitable for the simulation of single or
double-barreled shotgun libraries. Similarly, the RefSeq database contains
over 3,100 curated viral sequences (<ftp://ftp.ncbi.nih.gov/refseq/release/viral/>)
which can be used to produce artificial viral metagenomes.

Quite a few eukaryotic organisms have been sequenced and their genome and the
genes it contains can be the basis for simulating genomic and transcriptomic
(RNA-seq) datasets. For example, the human genome is available at
L<ftp://ftp.ncbi.nih.gov/refseq/H_sapiens/RefSeqGene/> and its transcripts can be
downloaded from L<ftp://ftp.ncbi.nih.gov/refseq/H_sapiens/mRNA_Prot/human.rna.fna.gz>

=head1 CLI EXAMPLES

=over

=item *

A shotgun library with a coverage of 0.1X

   Grinder -reference_file genomes.fna -coverage_fold 0.1

=item *

Same thing but save the result files in a specific folder and with a specific name

   Grinder -reference_file genomes.fna -coverage_fold 0.1 -base_name my_name -output_dir my_dir

=item *

A shotgun library with 1000 reads

   Grinder -reference_file genomes.fna -total_reads 1000

=item *

A shotgun library where species are distributed according to a power law

   Grinder -reference_file genomes.fna -abundance_model powerlaw 0.1

=item *

A shotgun library with 123 species

   Grinder -reference_file genomes.fna -diversity 123

=item *

Two shotgun libraries that have 50% of the species in common

   Grinder -reference_file genomes.fna -num_libraries 2 -shared_perc 50

=item *

A shotgun library where species relative abundances are manually specified

   Grinder -reference_file genomes.fna -abundance_file my_abundances.txt

=item *

A shotgun library with Sanger reads

   Grinder -reference_file genomes.fna -read_dist 800 -mutation_dist 1.5 linear 2 -mutation_ratio 4

=item *

A shotgun library with first-generation 454 reads

   Grinder -reference_file genomes.fna -read_dist 100 normal 10 -homopolymer_dist balzer

=item *

A paired-end shotgun library (insert size normally distributed around 2.5 kbp
with 0.2 kbp standard deviation)

   Grinder -reference_file genomes.fna -insert_dist 2500 normal 200

=item *

A 16S amplicon library

   Grinder -reference_file 16Sgenes.fna -forward_reverse 16Sprimers.fna -length_bias 0 -unidirectional 1

=item *

The same amplicon library with 20% of chimeric reads

   Grinder -reference_file 16Sgenes.fna -forward_reverse 16Sprimers.fna -length_bias 0 -unidirectional 1 -chimera_perc 20

=item *

Three 16S amplicon libraries with specified MIDs

   Grinder -reference_file 16Sgenes.fna -forward_reverse 16Sprimers.fna -length_bias 0 -unidirectional 1 -num_libraries 3 -multiplex_ids MIDs.fna

=item *

Reading reference sequences from the standard input, which allows you to decompress FASTA files on the fly

   zcat microbial_db.fna.gz | Grinder -total_reads 100

=back

=head1 OPTIONS

=head2 Basic parameters

=over

=item -rf <reference_file> | -reference_file <reference_file> | -gf <reference_file> | -genome_file <reference_file>

FASTA file that contains the input reference sequences (full genomes, 16S rRNA
genes, transcripts, ...) or '-' to read them from the standard input. See the
README file for examples of databases you can use and their location.
Default: reference_file.default

=for Euclid:
   reference_file.type: readable
   reference_file.default: '-'

=item -cf <coverage_fold> | -coverage_fold <coverage_fold>

Desired fold coverage of the input reference sequences (the output FASTA length divided
by the input FASTA length).
Default: coverage_fold.default x

=for Euclid:
   coverage_fold.type: +number
   coverage_fold.default: 0.1
   coverage_fold.excludes: total_reads

=item -tr <total_reads> | -total_reads <total_reads>

Number of shotgun or amplicon reads to generate for each library. Do not specify
this if you specify the coverage.

=for Euclid:
   total_reads.type: +integer

=back

=head2 Advanced shotgun and amplicon parameters

=over

=item -rd <read_dist>... | -read_dist <read_dist>...

Desired shotgun or amplicon read length distribution specified as:
   average length, distribution ('uniform' or 'normal') and standard deviation
Only the first element is required.
Examples:

  All sequences exactly 250 bp long: 250
  Uniform distribution around 100+-10 bp: 100 uniform 10
  Read normally distributed with an average of 800 and a standard deviation of 100
    bp: 800 normal 100

Genomes smaller than the specified length are not used. Default: read_dist.default

=for Euclid:
   read_dist.type: string
   read_dist.default: [100]

=item -id <insert_dist>... | -insert_dist <insert_dist>...

Create shotgun paired end reads (mate pairs) spanning the given insert length
(the reads are interior to the insert):
   0 : off,
   or: insert size distribution in bp, in the same format as the read length
       distribution (a typical value is 2,500 bp)
Two distinct reads are generated whether or not the mate pair overlaps.
Default: insert_dist.default

=for Euclid:
   insert_dist.type: string
   insert_dist.default: [0]

=item -ec <exclude_chars> | -exclude_chars <exclude_chars>

Do not create reads containing any of the specified characters (case 
insensitive), e.g. 'N-' to prevent reads with gaps (-) or ambiguities (N).
Default: 'exclude_chars.default'

=for Euclid:
   exclude_chars.type: string
   exclude_chars.default: ''

=item -dc <delete_chars> | -delete_chars <delete_chars>

Remove the specified characters from the reference sequences (case
insensitive), e.g. 'N-' to renove gaps (-) and ambiguities (N).
Default: delete_chars.default

=for Euclid:
   delete_chars.type: string
   delete_chars.default: ''

=item -fr <forward_reverse> | -forward_reverse <forward_reverse>

Use amplicon sequencing using the given forward and reverse PCR primer sequences
(in a FASTA file). It is recommended to use the <length_bias> and <unidirectional>
options to generate amplicon reads. To sequence from the forward strand
(<unidirectional> = 1), put the forward primer first and reverse primer second.
To sequence from the reverse strand, invert the primers in the FASTA file and
use <unidirectional> = -1. The second primer sequence in the FASTA file is always
optional. The sequences should use the IUPAC convention for degenerate residues.
Example: AAACTYAAAKGAATTGRCGG and ACGGGCGGTGTGTRC for the 926F and 1392R primers
respectively (primers that target the v6 to v9 region of the 16S rRNA gene).
Genome sequences that do not match the specified primers are excluded.

=for Euclid:
   forward_reverse.type: readable

=item -un <unidirectional> | -unidirectional <unidirectional>

Instead of producing reads bidirectionally, i.e. from the reference strand and
its reverse complement, proceed unidirectionally, i.e. from one strand only
(forward or reverse). Values: 0 (off), 1 (forward), -1 (reverse)
Default: unidirectional.default

=for Euclid:
   unidirectional.type: integer, unidirectional >= -1 && unidirectional <= 1
   unidirectional.type.error: <unidirectional> must be 0, 1 or -1 (not unidirectional)
   unidirectional.default: 0

=item -lb <length_bias> | -length_bias <length_bias>

In shotgun libraries, sample species proportionally to their genome length:
at the same relative abundance, larger genomes contribute more reads than smaller
genomes. 0 = no, 1 = yes.
Default: length_bias.default

=for Euclid:
   length_bias.type: integer, length_bias == 0 || length_bias == 1
   length_bias.type.error: <length_bias> must be 0 or 1 (not length_bias)
   length_bias.default: 1

=item -cb <copy_bias> | -copy_bias <copy_bias>

In amplicon libraries, sample species proportionally to the number of copies of
the target gene: at equal relative abundance, genomes that have multiple copies
of the target gene contribute more amplicon reads than genomes that have a
single copy. Note: you should use full genomes in <reference_file> to make use
of this option. 0 = no, 1 = yes.
Default: copy_bias.default

=for Euclid:
   copy_bias.type: integer, copy_bias == 0 || copy_bias == 1
   copy_bias.type.error: <copy_bias> must be 0 or 1 (not copy_bias)
   copy_bias.default: 1

=back

=head2 Aberrations and sequencing errors

=over

=item -md <mutation_dist>... | -mutation_dist <mutation_dist>...

Introduce sequencing errors in the reads, under the form of mutations
(substitutions, insertions and deletions) using a specified frequency
distribution:
   average probability (%),
   model (uniform, linear),
   value at 3' end (not applicable for uniform model).
For example, for Sanger-type errors, use:
   1.5 linear 2
Use the <mutation_ratio> option to alter how many of these mutations are
substitutions or indels.
Default: mutation_dist.default

=for Euclid:
   mutation_dist.type: string
   mutation_dist.default: [0, 'uniform', 0]

=item -mr <mutation_ratio> | -mutation_ratio <mutation_ratio>

Indicate the ratio of the number of substitutions to the number of indels
(insertions and deletions). For example, use 4 (4 substitutions for 1 indel)
for Sanger reads.
Default: mutation_ratio.default

=for Euclid:
   mutation_ratio.type: num, mutation_ratio >= 0
   mutation_ratio.default: 0

=item -hd <homopolymer_dist> | -homopolymer_dist <homopolymer_dist>

Introduce sequencing errors in the reads under the form of homopolymeric
stretches (e.g. AAA, CCCCC) using a specified model where the homopolymer length
follows a normal distribution N(mean, standard deviation) that is function of
the homopolymer length n:

  Margulies: N(n, 0.15 * n)             ,  Margulies et al. 2005.
  Richter  : N(n, 0.15 * sqrt(n))       ,  Richter et al. 2008.
  Balzer   : N(n, 0.03494 + n * 0.06856),  Balzer et al. 2010.

Default: homopolymer_dist.default

=for Euclid:
   homopolymer_dist.type: string
   homopolymer_dist.default: 0

=item -cp <chimera_perc> | -chimera_perc <chimera_perc>

Specify the percent of reads in amplicon libraries that should be chimeric
sequences. The 'reference' field in the description of chimeric reads will
contain the ID of all the reference sequences forming the chimeric template. A
typical value is 10%. Default: chimera_perc.default %

=for Euclid:
   chimera_perc.type: number, chimera_perc >= 0 && chimera_perc <= 100
   chimera_perc.type.error: <chimera_perc> must be a number between 0 and 100 (not chimera_perc)
   chimera_perc.default: 0

=back

=head2 Community structure and diversity

=over

=item -af <abundance_file> | -abundance_file <abundance_file>

Specify the relative abundance of the genomes manually in an input file. Each
line of the file should contain a sequence name and its relative abundance (%),
e.g. 'seqABC 82.1' or 'seqABC 82.1 10.2' if you are specifying 2 different
communities.

=for Euclid:
   abundance_file.type: readable

=item -am <abundance_model>... | -abundance_model <abundance_model>...

Relative abundance model for the input genomes: uniform, linear, powerlaw,
logarithmic or exponential. The uniform and linear models do not require a
parameter, but the other models take a parameter in the range [0, infinity). If
this parameter is not specified, then it is randomly picked.
Examples:

  uniform distribution: uniform
  powerlaw distribution with parameter 0.1: powerlaw 0.1
  exponential distribution with automatically chosen parameter: exponential

Default: abundance_model.default

=for Euclid:
   abundance_model.type: string
   abundance_model.default: ['uniform', 1]

=item -nl <num_libraries> | -num_libraries <num_libraries>

Number of independent libraries to create. Specify how diverse and similar they
should be with <diversity>, <shared_perc> and <permuted_perc>. Assign them
different MID tags with <multiplex_mids>.
Default: num_libraries.default

=for Euclid:
   num_libraries.type: +integer
   num_libraries.default: 1

=item -mi <multiplex_ids> | -multiplex_ids <multiplex_ids>

Specify an optional FASTA file that contains multiplex sequence identifiers
(a.k.a MIDs or barcodes) to add to the sequences (one per library). The MIDs
are included in the length specified with the -read_dist option.

=for Euclid:
   multiplex_ids.type: readable

=item -di <diversity>... | -diversity <diversity>...

Richness, or number of genomes to include in the shotgun libraries. Use 0 for
the maximum diversity possible (based on the number of reference sequences
available). Provide one value to make all libraries have the same diversity, or
one diversity value per library otherwise.
Default: diversity.default

=for Euclid:
   diversity.type: 0+integer
   diversity.default: [ 0 ]

=item -sp <shared_perc> | -shared_perc <shared_perc>

For multiple libraries, percent of genomes they should have in common (relative
to the diversity of the least diverse library).
Default: shared_perc.default %

=for Euclid:
   shared_perc.type: number, shared_perc >= 0 && shared_perc <= 100
   shared_perc.type.error: <shared_perc> must be a number between 0 and 100 (not shared_perc)
   shared_perc.default: 0

=item -pp <permuted_perc> | -permuted_perc <permuted_perc>

For multiple libraries, percent of the most-abundant genomes to permute in
rank-abundance.
Default: permuted_perc.default %

=for Euclid:
   permuted_perc.type: number, permuted_perc >= 0 && permuted_perc <= 100
   permuted_perc.type.error: <permuted_perc> must be a number between 0 and 100 (not permuted_perc)
   permuted_perc.default: 0

=back

=head2 Miscellaneous

=over

=item -rs <random_seed> | -random_seed <random_seed>

Seed number to use for the pseudo-random number generator.

=for Euclid:
   random_seed.type: +integer

=item -dt <desc_track> | -desc_track <desc_track>

Track read information (reference sequence, position, errors, ...) by writing
it in the read description.
Default: desc_track.default

=for Euclid:
   desc_track.type: number, desc_track == 0 || desc_track == 1
   desc_track.type.error: <desc_track> must be 0 or 1 (not desc_track)
   desc_track.default: 1

=item -ql <qual_levels>... | -qual_levels <qual_levels>...

Generate very basic quality scores for the simulated reads. Good residues are
given a specified good score (e.g. 30) and residues that are the result of an
insertion or substitution are given a specified bad score (e.g. 10). Specify
first the good score and then the bad score on the command-line, e.g.: 30 10
Default: qual_levels.default

=for Euclid:
   qual_levels.type: 0+integer
   qual_levels.default: [ ]

=item -bn <base_name> | -base_name <base_name>

Prefix of the output files.
Default: base_name.default

=for Euclid:
   base_name.type: string
   base_name.default: 'grinder'

=item -od <output_dir> | -output_dir <output_dir>

Directory where the results should be written. This folder will be created if
needed.
Default: output_dir.default

=for Euclid:
   output_dir.type: writable
   output_dir.default: '.'

=item -pf <profile_file> | -profile_file <profile_file>

A file that contains Grinder arguments. This is useful if you use many options
or often use the same options. Lines with comments (#) are ignored. Consider the
profile file, 'simple_profile.txt':

  # A simple Grinder profile
  -read_dist 105 normal 12
  -total_reads 1000

Running: Grinder -reference_file viral_genomes.fa -profile_file simple_profile.txt

Translates into: Grinder -reference_file viral_genomes.fa -read_dist 105 normal 12 -total_reads 1000

Note that the arguments specified in the profile should not be specified again on the command line.

=back

=head1 API EXAMPLES

  use Grinder;

  # Set up a new factory (see the OPTIONS section for a complete list of parameters)
  my $factory = Grinder->new( -reference_file => 'genomes.fna' );

  # Process all shotgun libraries requested
  while ( my $struct = $factory->next_lib ) {

    # The ID and abundance of the 3rd most abundant genome in this community
    my $id = $struct->{ids}->[2];
    my $ab = $struct->{abs}->[2];

    # Create shotgun reads
    while ( my $read = $factory->next_read) {

      # The read is a Bioperl sequence object with these properties:
      my $read_id     = $read->id;     # read ID given by Grinder
      my $read_seq    = $read->seq;    # nucleotide sequence
      my $read_mid    = $read->mid;    # MID or tag attached to the read
      my $read_errors = $read->errors; # errors that the read contains
 
      # Where was the read taken from? The reference sequence refers to the
      # database sequence for shotgun libraries, amplicon obtained from the
      # database sequence, or could even be a chimeric sequence
      my $ref_id     = $read->reference->id; # ID of the reference sequence
      my $ref_start  = $read->start;         # start of the read on the reference
      my $ref_end    = $read->end;           # end of the read on the reference
      my $ref_strand = $read->strand;        # strand of the reference
      
    }
  }

  # Similarly, for shotgun mate pairs
  my $factory = Grinder->new( -reference_file => 'genomes.fna',
                              -insert_dist    => 250            );
  while ( $factory->next_lib ) {
    while ( my $read = $factory->next_read ) {
      # The first read is the first mate of the mate pair
      # The second read is the second mate of the mate pair
      # The third read is the first mate of the next mate pair
      # ...
    }
  }

  # To generate an amplicon library
  my $factory = Grinder->new( -reference_file  => 'genomes.fna',
                              -forward_reverse => '16Sgenes.fna',
                              -length_bias     => 0,
                              -unidirectional  => 1              );
  while ( $factory->next_lib ) {
    while ( my $read = $factory->next_read) {
      # ...
    }
  }

=head1 API METHODS

The rest of the documentation details the available Grinder API methods.

=head2 new

Title   : new
Function: Create a new Grinder factory initialized with the passed arguments.
          Available parameters described in the OPTIONS section.
Usage   : my $factory = Grinder->new( -reference_file => 'genomes.fna' );
Returns : a new Grinder object

=head2 next_lib

Title   : next_lib
Function: Go to the next shotgun library to process.
Usage   : my $struct = $factory->next_lib;
Returns : Community structure to be used for this library, where $struct->{ids}
          is an array reference containing the IDs of the genome making up the
          community (sorted by decreasing relative abundance) and $struct->{abs}
          is an array reference of the genome abundances (in the same order as
          the IDs).

=head2 next_read

Title   : next_read
Function: Create a amplicon or shotgun read  for the current library.
Usage   : my $read  = $factory->next_read; # for single read
          my $mate1 = $factory->next_read; # for mate pairs
          my $mate2 = $factory->next_read; 
Returns : A sequence represented as a Bio::Seq::SimulatedRead object

=head2 get_random_seed

Title   : get_random_seed
Function: Return the number used to seed the pseudo-random number generator
Usage   : my $seed = $factory->get_random_seed;
Returns : seed number


=head1 COPYRIGHT

Copyright 2009,2010,2011 Florent ANGLY <florent.angly@gmail.com>

Grinder is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License (GPL) as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
Grinder is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with Grinder.  If not, see <http://www.gnu.org/licenses/>.

=head1 BUGS

All complex software has bugs lurking in it, and this program is no exception.
If you find a bug, please report it on the SourceForge Tracker for Grinder:
L<http://sourceforge.net/tracker/?group_id=244196&atid=1124737>

Bug reports, suggestions and patches are welcome. Grinder's code is developed
on Sourceforge (L<https://sourceforge.net/scm/?type=git&group_id=244196>) and is
under Git revision control. To get started with a patch, do:

   git clone git://biogrinder.git.sourceforge.net/gitroot/biogrinder/biogrinder

=cut

#---------- GRINDER FUNCTIONAL API --------------------------------------------#

sub Grinder {
  my (@args) = @_;

  # Create Grinder object
  my $factory = Grinder->new(@args);

  # Print diversity and percent shared and permuted
  diversity_report( $factory->{num_libraries}, $factory->{shared_perc},
    $factory->{permuted_perc}, $factory->{overall_diversity} );

  # Create the output directory if needed
  if ( not -d $factory->{output_dir} ) {
    mkdir $factory->{output_dir} or die "Error: Could not create output folder ".
      $factory->{output_dir}."\n$!\n";
  }

  # Generate sequences
  while ( my $c_struct = $factory->next_lib ) {
    my $cur_lib = $factory->{cur_lib};

    # Output filenames
    my $lib_str = '';
    if ($factory->{num_libraries} > 1) {
      $lib_str = '-'.sprintf('%0'.length($factory->{num_libraries}).'d', $cur_lib);
    }
    my $out_fasta_file = File::Spec->catfile($factory->{output_dir}, 
        $factory->{base_name}.$lib_str."-reads.fa");
    my $out_qual_file;
    if (scalar @{$factory->{qual_levels}} > 0) {
      $out_qual_file = File::Spec->catfile($factory->{output_dir}, 
        $factory->{base_name}.$lib_str."-reads.qual");
    }
    my $out_ranks_file = File::Spec->catfile($factory->{output_dir},
      $factory->{base_name}.$lib_str."-ranks.txt");

    # Write community structure file
    $factory->write_community_structure($c_struct, $out_ranks_file);

    # Prepare output FASTA file
    my $out_fasta = Bio::SeqIO->new( -format => 'fasta',
                                     -flush  => 0,
                                     -file   => ">$out_fasta_file" );

    my $out_qual;
    if ( defined $out_qual_file ) {
       $out_qual = Bio::SeqIO->new( -format => 'qual',
                                    -flush  => 0,
                                    -file   => ">$out_qual_file" );
    }

    # Library report
    my $diversity = $factory->{diversity}[$cur_lib-1];
    library_report( $cur_lib, $out_ranks_file, $out_fasta_file, $out_qual_file,
      $factory->{cur_coverage_fold}, $factory->{cur_total_reads}, $diversity);

    # Generate shotgun or amplicon reads and write them to a file
    while ( my $read = $factory->next_read ) {
      $out_fasta->write_seq($read);
      $out_qual->write_seq($read) if defined $out_qual
    }
    $out_fasta->close;
    $out_qual->close if defined $out_qual;
  }

  return 1;
}


sub diversity_report {
  my ($num_libraries, $perc_shared, $perc_permuted, $overall_diversity) = @_;
  my $format = '%.1f';
  print "Overall diversity = $overall_diversity genomes\n";
  if ($num_libraries > 1) {
    my $nof_shared  = $perc_shared / 100 * $overall_diversity;
    $perc_shared = sprintf($format, $perc_shared);
    print "Percent shared   = $perc_shared % ($nof_shared genomes)\n";
    my $nof_permuted  = $perc_permuted / 100 * $overall_diversity;
    $perc_permuted = sprintf($format, $perc_permuted);
    print "Percent permuted = $perc_permuted % ($nof_permuted top genomes)\n";
  }
  return 1;
}


sub write_community_structure {
  my ($self, $c_struct, $filename) = @_;
  open(OUT, ">$filename") || die("Error: Could not write in file $filename: $!\n");
  print OUT "# rank\tseqID\trel. abundance\n";
  my $diversity = scalar @{$c_struct->{ids}};
  for my $rank ( 1 .. $diversity ) {
    my $oid   = $c_struct->{'ids'}->[$rank-1];
    my $seqid = $self->database_get_parent_id($oid);
    my $relab = $c_struct->{'abs'}->[$rank-1];
    print OUT "$rank\t$seqid\t$relab\n";
  }
  close OUT;
  return 1;
}


sub library_report {
  my ($cur_lib, $ranks_file, $fasta_file, $qual_file, $coverage, $nof_seqs, $diversity) = @_;
  my $format = '%.3f';
  $coverage = sprintf($format, $coverage);
  print "Shotgun library $cur_lib:\n";
  print "   Community structure  = $ranks_file\n";
  print "   FASTA file           = $fasta_file\n";
  print "   QUAL file            = $qual_file\n" if defined $qual_file;
  print "   Library coverage     = $coverage x\n";
  print "   Number of reads      = $nof_seqs\n";
  print "   Diversity (richness) = $diversity\n";
  return 1;
}


#---------- GRINDER OO API ----------------------------------------------------#


sub new {
  my ($class, @args) = @_;
  my $self = {};
  bless $self, ref($class) || $class;
  $self->argparse(\@args);
  $self->initialize();
  return $self;
}


sub next_lib {
  my ($self) = @_;
  $self->{cur_lib}++;
  $self->{cur_read} = 0;
  $self->{cur_total_reads} = 0;
  $self->{cur_coverage_fold} = 0;
  $self->{next_mate} = undef;
  $self->{positions} = undef;
  my $c_struct = $self->{c_structs}[$self->{cur_lib}-1];
  if ( defined $c_struct ) {
    # Create probabilities of picking genomes from community structure
    $self->{positions} = $self->proba_create($c_struct, $self->{length_bias},
      $self->{copy_bias});
    # Calculate needed number of sequences based on desired coverage
    ($self->{cur_total_reads}, $self->{cur_coverage_fold}) = $self->lib_coverage($c_struct);
  }
  return $c_struct;
}


sub next_read {
  my ($self) = @_;
  $self->next_lib if not $self->{cur_lib};
  $self->{cur_read}++;
  my $read;
  if ( $self->{cur_read} <= $self->{cur_total_reads} ) {
    # Generate the next read
    if ($self->{mate_length}) {
      # Generate a mate pair read
      if ( not $self->{next_mate} ) {
        # Generate a new pair of reads
        ($read, my $read2) = $self->next_mate_pair( );
        # Save second read of the pair for later
        $self->{next_mate} = $read2;        
      } else {
        # Use saved read
        $read = $self->{next_mate};
        $self->{next_mate} = undef;
      }
    } else {
      # Generate a single shotgun or amplicon read
      $read = $self->next_single_read( );
    }
  }
  return $read;
}


sub get_random_seed {
  my ($self) = @_;
  return $self->{random_seed};
}


#---------- GRINDER INTERNALS -------------------------------------------------#


sub argparse {
  # Process arguments
  my ($self, $args) = @_;
  # Read profile file
  $args = process_profile_file($args);
  # Parse and validate arguments with Getopt::Euclid
  Getopt::Euclid->process_args($args);
  # Get parsed arguments from %ARGV and put them in $self
  for my $arg (keys %ARGV) {
    # Skip short argument names (they are also represented with long names)
    next if length($arg) <= 2;
    # Process long argument names. Copy their value into $self
    my $ref = ref($ARGV{$arg});
    if (not $ref) {
      $self->{$arg} = $ARGV{$arg};
    } elsif ($ref eq 'ARRAY') {
      @{$self->{$arg}} = @{$ARGV{$arg}};
    } else {
      die "Error: unsupported operation on argument '$arg' which is a reference".
         "of type $ref\n";
    }
  }
  return 1;
}


sub process_profile_file {
  # Find profile file in arguments and read the profiles. The profile file
  # only contains Grinder arguments, and lines starting with a '#' are comments.
  my ($args) = @_;
  my $file;
  for (my $i = 0; $i < scalar @$args; $i++) {
    my $arg = $$args[$i];
    if ($arg =~ m/^-profile_file/ || $arg =~ m/-pf/) {
      $file = $$args[$i+1];
      if ( (not defined $file) || ($file =~ m/^-/) ) {
        die "Error: no value was given to --profile_file\n";
      }
    }
  }
  if (defined $file) {
    open my $in, '<', $file or die "Error: Could not read file '$file'\n$!\n";
    my $profile = '';
    while (my $line = <$in>) {
      chomp $line;
      next if $line =~ m/^\s*$/;
      next if $line =~ m/^\s*#/;
      $profile .= "$line ";
    }
    close $in;
    push @$args, split /\s+/, $profile;
  }
  return $args;
}


sub initialize {
  my ($self) = @_;
  # Returns:

  # Parameter processing: read length distribution
  if ( (not ref $self->{read_dist}) or (ref $self->{read_dist} eq 'SCALAR') ){
    $self->{read_dist}   = [$self->{read_dist}];
  }
  $self->{read_length} = $self->{read_dist}[0] || 100;
  $self->{read_model}  = $self->{read_dist}[1] || 'uniform';
  $self->{read_delta}  = $self->{read_dist}[2] || 0;
  delete $self->{read_dist};

  # Parameter processing: mate insert length distribution
  if ( (not ref $self->{insert_dist}) or (ref $self->{insert_dist} eq 'SCALAR') ){
    $self->{insert_dist} = [$self->{insert_dist}];
  }
  $self->{mate_length} = $self->{insert_dist}[0] || 0;
  $self->{mate_model}  = $self->{insert_dist}[1] || 'uniform';
  $self->{mate_delta}  = $self->{insert_dist}[2] || 0;
  delete $self->{insert_dist};

  # Parameter processing: genome abundance distribution
  if ( (not ref $self->{abundance_model}) or (ref $self->{abundance_model} eq 'SCALAR') ){
    $self->{abundance_model} = [$self->{abundance_model}];
  }
  $self->{distrib} = $self->{abundance_model}[0] || 'uniform';
  $self->{param}   = $self->{abundance_model}[1];
  delete $self->{abundance_model};

  # Parameter processing: point sequencing error distribution
  if ( (not ref $self->{mutation_dist}) or (ref $self->{mutation_dist}  eq 'SCALAR') ) {
    $self->{mutation_dist} = [$self->{mutation_dist}];
  }
  $self->{mutation_freq}  = $self->{mutation_dist}[0] || 0;
  $self->{mutation_model} = $self->{mutation_dist}[1] || 'uniform';
  $self->{mutation_end}   = $self->{mutation_dist}[2] || 0;
  delete $self->{mutation_dist};

  # Parameter processing: homopolymer model
  $self->{homopolymer_dist} = lc $self->{homopolymer_dist};

  # Seed the random number generator (manually or automatically)
  $self->{random_seed} = srand( $self->{random_seed} );

  # Sequence length check
  my $max_read_length = $self->{read_length} + $self->{read_delta}; # approximation
  if ($self->{mate_length}) {
    my $min_mate_length = $self->{mate_length} - $self->{mate_delta};
    if ($max_read_length > $min_mate_length) {
      die("Error: The mate insert length cannot be smaller than read length. ".
        "Try increasing the mate insert length or decreasing the read length\n");
    }
  }
  
  # Read MIDs
  $self->{multiplex_ids} = $self->read_multiplex_id_file($self->{multiplex_ids}, 
    $self->{num_libraries}) if defined $self->{multiplex_ids};

  # Import genome sequences, skipping genomes too short
  $self->{database} = $self->database_create( $self->{reference_file},
    $self->{unidirectional}, $self->{forward_reverse}, $self->{abundance_file},
    $self->{delete_chars} );

  # Genome relative abundance in the different independent libraries to create
  $self->{c_structs} = $self->community_structures( $self->{database}->{ids},
    $self->{abundance_file}, $self->{distrib}, $self->{param},
    $self->{num_libraries}, $self->{shared_perc}, $self->{permuted_perc},
    $self->{diversity}, $self->{forward_reverse} );

  # Markers to keep track of computation progress
  $self->{cur_lib}  = 0;
  $self->{cur_read} = 0;

  return $self;
}


sub read_multiplex_id_file {
  my ($self, $file, $nof_indep) = @_;
  my @mids;
  # Read FASTA file containing the MIDs
  my $in = Bio::SeqIO->newFh(
    -file   => $file,
    -format => 'fasta',
  );
  while (my $mid = <$in>) {
    push @mids, $mid->seq;
  }
  undef $in;
  # Sanity check
  my $nof_mids = scalar @mids;
  if ($nof_mids < $nof_indep) {
    die "Error: $nof_indep communities were requested but the MID file ".
      "had only $nof_mids sequences.\n"; 
  } elsif ($nof_mids > $nof_indep) {
    warn "Warning: $nof_indep communities were requested but the MID file ".
      "contained $nof_mids sequences. Ignoring extraneous MIDs...\n";
  }
  return \@mids;
}


sub community_structures {
  # Create communities with a specified structure, alpha and beta-diversity
  my ($self, $seq_ids, $abundance_file, $distrib, $param, $nof_indep,
    $perc_shared, $perc_permuted, $diversities, $forward_reverse) = @_;

  # Calculate community structures
  my $c_structs;
  if ($abundance_file) {
    # Sanity check
    if ( (scalar @$diversities > 1) || $$diversities[0] ) {
      warn "Warning: Diversity cannot be specified when an abundance file is specified. Ignoring it...\n";
    }
    if ( $perc_shared || $perc_permuted ) {
      warn "Warning: Percent shared and percent permuted cannot be specified when an abundance file is specified. Ignoring them...\n";
    }
    if ( $nof_indep != 1 ) { # 1 is the default value
      warn "Warning: The number of libraries cannot be specified when an abundance file is specified. Ignoring it...\n";
    }
    # One or several communities with specified rank-abundances
    $c_structs = community_given_abundances($abundance_file, $seq_ids);
    # Calculate number of libraries
    $nof_indep = scalar @$c_structs;
    $self->{num_libraries} = $nof_indep;
    # Calculate diversities based on given community abundances
    ($self->{diversity}, $self->{overall_diversity}, $self->{shared_perc},
      $self->{permuted_perc}) = community_calculate_diversities($c_structs);
  } else {
    # One or several communities with rank-abundance to be calculated
    # Sanity check
    if ($nof_indep == 1) { # 1 is the default value
      $nof_indep = scalar @$diversities;
    }
    if ($nof_indep != scalar @$diversities) {
      if (scalar @$diversities == 1) {
        # Use same diversity for all libraries
        my $diversity = $$diversities[0];
        for my $i (1 .. $nof_indep-1) {
          push @$diversities, $diversity;
        }
      } else {
        die "Error: The number of diversities provided does not match the requested number of libraries.\n";
      }
    }
    $self->{num_libraries} = $nof_indep;
    # Select shared species
    my $c_ids;
    my $overall_diversity = 0;
    ($c_ids, $overall_diversity, $diversities, $perc_shared) = community_shared(
      $seq_ids, $nof_indep, $perc_shared, $diversities );

    # Shuffle the abundance-ranks of the most abundant genomes
    ($c_ids, $perc_permuted) = community_permuted($c_ids, $perc_permuted);
    # Update values in $self object
    $self->{overall_diversity} = $overall_diversity;
    $self->{diversity} = $diversities;
    $self->{shared_perc} = $perc_shared;
    $self->{permuted_perc} = $perc_permuted;
    # Put results in a community structure "object"
    for my $c (1 .. $nof_indep) {
      # Assign a random parameter if needed
      my $comm_param = defined $param ? $param : randig(1,0.05);
      # Calculate relative abundance of the community members
      my $diversity = $self->{diversity}[$c-1];
      my $c_abs = community_calculate_species_abundance($distrib, $comm_param,
         $diversity);
      my $c_ids = $$c_ids[$c-1];
      my $c_struct;
      $c_struct->{'ids'}   = $c_ids;
      $c_struct->{'abs'}   = $c_abs;
      $c_struct->{'param'} = $comm_param;
      $c_struct->{'model'} = $distrib;
      push @$c_structs, $c_struct;
    }
  }

  # Convert sequence IDs to object IDs
  for my $c_struct (@$c_structs) {
    my ($c_abs, $c_ids) = community_calculate_amplicon_abundance(
      $c_struct->{'abs'}, $c_struct->{'ids'}, $seq_ids );
  }

  return $c_structs;
}

sub community_calculate_diversities {
  my ($c_structs) = @_;
  my ($diversities, $overall_diversity, $perc_shared, $perc_permuted) = (0, 0, 0, 0);

  # Calculate diversity (richness) based on given community abundances
  my $nof_libs = scalar @$c_structs;
  my %all_ids;
  my @richnesses;
  for my $c_struct (@$c_structs) {
    my $richness = 0;
    for my $i (0 .. scalar @{$$c_struct{ids}} - 1) {
      my $id = $$c_struct{ids}[$i];
      my $ab = $$c_struct{abs}[$i];
      next if not $ab;
      $richness++;
      
      if (defined $all_ids{$id}) {
        $all_ids{$id}++;
      } else {
        $all_ids{$id} = 1;
      }
    }
    push @richnesses, $richness;
  }
  $overall_diversity = scalar keys %all_ids;


  # Calculate percent shared
  my $nof_non_shared = 0;
  for my $id (keys %all_ids) {
    $nof_non_shared++ if $all_ids{$id} < $nof_libs;
  }
  $perc_shared = ($overall_diversity - $nof_non_shared) * 100 / $overall_diversity;

  #####
  # TODO
  # Calculate percent permuted
  # ...
  ##### 

  return \@richnesses, $overall_diversity, $perc_shared, $perc_permuted;
}

sub community_given_abundances {
  # Read a file of genome abundances. The file should be space or tab-delimited. 
  # The first column should be the IDs of genomes, and the subsequent columns is
  # for their relative abundance in different communities. An optional list of
  # valid IDs can be provided. Then the abundances are normalized so that their
  # sum is 1.
  my ($file, $seq_ids) = @_;

  # Read abundances
  my ($ids, $abs) = community_read_abundances($file);
  # Remove genomes with unknown IDs and calculate cumulative abundance
  my $totals;
  for my $comm_num (0 .. $#$ids) {
    my $i = 0;
    while ( $i < scalar @{$$ids[$comm_num]} ) {
      my $id = $$ids[$comm_num][$i];
      my $ab = $$abs[$comm_num][$i];
      if ( (scalar keys %$seq_ids == 0) || (exists $$seq_ids{$id}) ) {
        $$totals[$comm_num] += $ab;
        $i++;
      } else {
        die "Error: Requested reference sequence '$id' in file '$file' does not".
          " exist in the input database.\n";
        splice @{$$ids[$comm_num]}, $i, 1;
        splice @{$$abs[$comm_num]}, $i, 1;
      }
    }
  }
  # Process the communities
  my @c_structs;
  for my $comm_num (0 .. scalar @$ids - 1) {
    my $comm_ids   = $$ids[$comm_num];
    my $comm_abs   = $$abs[$comm_num];
    my $comm_total = $$totals[$comm_num];
    if ($comm_total == 0) {
      warn "Warning: The abundance of all the genomes for community ".($comm_num+1)." was zero. Skipping this community...\n";
      next;
    }
    # Normalize the abundances
    $comm_abs = [ map { $_ / $comm_total } (@$comm_abs) ];
    # Sort relative abundances by decreasing 
    ($comm_abs, $comm_ids) = two_array_sort($comm_abs, $comm_ids);
    $comm_abs = [reverse(@$comm_abs)];
    $comm_ids = [reverse(@$comm_ids)];
    # Save community structure
    my $c_struct = { 'ids' => $comm_ids, 'abs' => $comm_abs };
    push @c_structs, $c_struct;
  }
  return \@c_structs;
}


sub community_read_abundances {
  my ($file) = @_;
  # Read abundances of genomes from a file
  my $ids; # genome IDs
  my $abs; # genome relative abundance
  open my $io, '<', $file or die "Error: Could not read file '$file'\n$!\n";
  while ( my $line = <$io> ) {
    # Ignore comment or empty lines
    if ( $line =~ m/^\s*$/ || $line =~ m/^#/ ) {
      next;
    }
    # Read abundance info from line
    my ($id, @rel_abs) = ($line =~ m/(\S+)/g);
    if (defined $id) {
      for my $comm_num (0 .. $#rel_abs) {
        my $rel_ab = $rel_abs[$comm_num];
        push @{$$ids[$comm_num]}, $id;
        push @{$$abs[$comm_num]}, $rel_ab;
      }
    } else {
      warn "Warning: Line $. of file '$file' has an unknown format. Skipping it...\n";
    }
  }
  close $io;
  return $ids, $abs;
}


sub community_permuted {
  # Change the abundance rank of species in all but the first community.
  # The number of species changed in abundance is determined by the percent
  # permuted, i.e. a given percentage of the most abundant species in this community.
  my ($c_ids, $perc_permuted) = @_;
  my $nof_indep = scalar @$c_ids;

  # Leave the first community alone, but permute the ones after
  for my $c ( 2 .. $nof_indep ) {
    my $ids = $$c_ids[$c-1];
    my $diversity = scalar @$ids;
    # Number of top genomes to permute
    # Percent permuted is relative to diversity in this community
    my $nof_permuted = $perc_permuted / 100 * $diversity;
    $nof_permuted = int($nof_permuted + 0.5); # round number

    # Method published in Angly et al 2006 PLOS Biology    
    # Take the $nof_permuted first ranks (most abundant genomes) and shuffle
    # (permute) their ranks amongst the $nof_permuted first ranks.
    # Caveat: cannot permute only 1 genome
    my $idxs;
    if ($nof_permuted > 0) {
      # Add shuffled top genomes
      my $permuted_idxs = randomize( [0 .. $nof_permuted-1] );
      push @$idxs, @$permuted_idxs;
    }
    if ($diversity - $nof_permuted > 0) {
      # Add other genomes in same order
      my $non_permuted_idxs = [$nof_permuted .. $diversity-1];
      push @$idxs, @$non_permuted_idxs;
    }
    @$ids = @$ids [ @$idxs ];

  }

  return $c_ids, $perc_permuted;
}


sub community_shared {
  # Randomly split a library of sequences into a given number of groups that
  # share a specified percent of their genomes.
  # The % shared is the number of species shared / the total diversity in all communities
  # Input: arrayref of sequence ids
  #        number of communities to produce
  #        percentage of genomes shared between the communities
  #        diversity (optional, will use all genomes if not specified) 
  # Return: arrayref of IDs that are shared
  #         arrayref of arrayref with the unique IDs for each community
  my ($seq_ids, $nof_indep, $perc_shared, $diversities) = @_;

  # If diversity is not specified (is '0'), use the maximum value possible
  my $nof_refs = scalar keys %$seq_ids;
  my $min_diversity = 1E99;
  for my $i (0 .. scalar @$diversities - 1) {
    if ($$diversities[$i] == 0) {
      $$diversities[$i] = $nof_refs / ( $perc_shared/100 + $nof_indep*(1-$perc_shared/100) );
      $$diversities[$i] = int( $$diversities[$i] );
      if ( ($i > 0) && ($$diversities[$i-1] != $$diversities[$i]) ) {
        die "Error: Define either all the diversities or none.\n";
      }
    }
    if ($$diversities[$i] < $min_diversity) {
      $min_diversity = $$diversities[$i];
    }
  }

  if ($min_diversity == 0) {
    die "Error: Cannot make $nof_indep libraries sharing $perc_shared % species".
      " from $nof_refs references\n";
  }

  # Calculate the number of sequences to share, noting that the percent shared
  # is relative to the diversity of the least abundant library
  my $nof_shared = int($min_diversity * $perc_shared / 100);
  $perc_shared = $nof_shared * 100 / $min_diversity;

  # Unique sequences
  my @nof_uniques;
  my $sum_not_uniques = 0;
  for my $diversity (@$diversities) {
    my $nof_unique = $diversity - $nof_shared;
    $sum_not_uniques += $nof_unique;
    push @nof_uniques, $nof_unique;
  }

  # Overall diversity
  my $overall_diversity = $nof_shared + $sum_not_uniques;
  if ($nof_refs <  $overall_diversity) {
    die "Error: The number of reference sequences available ($nof_refs) is not".
      " large enough to support the requested diversity ($overall_diversity ".
      "genomes overall with $perc_shared % genomes shared between $nof_indep ".
      "libraries)\n";
  }

  # Add shared sequences
  my @ids = keys %$seq_ids;
  my @shared_ids;
  for (0 .. $nof_shared - 1) {
    # Pick a random sequence
    my $rand_offset = int(rand($nof_refs));
    my $rand_id = splice @ids, $rand_offset, 1;
    $nof_refs = scalar(@ids);
    # Add this sequence in all independent libraries
    push @shared_ids, $rand_id;
  }

  # Add sequences not shared
  my @unique_ids;
  for my $lib_num (0 .. $nof_indep-1) {
    my $nof_unique = $nof_uniques[$lib_num];
    for (0 .. $nof_unique - 1) {
      # Pick a random sequence
      my $rand_offset = int(rand($nof_refs));
      my $rand_id = splice @ids, $rand_offset, 1;
      $nof_refs = scalar(@ids);
      # Add this sequence in this independent library only
      push @{$unique_ids[$lib_num]}, $rand_id;
    }
  }

  # Randomly pick the rank of the shared IDs
  my $shared_ranks = randomize( [1 .. $min_diversity] );
  @$shared_ranks = splice @$shared_ranks, 0, $nof_shared;

  # Construct community ranks
  my @c_ranks;
  for my $lib_num (0 .. $nof_indep-1) {
    my $diversity = $$diversities[$lib_num];
    my @ranks = (undef) x $diversity;
    # Add shared IDs
    for my $i (0 .. $nof_shared-1) {
      my $id   = $shared_ids[$i];
      my $rank = $$shared_ranks[$i];
      $ranks[$rank-1] = $id;
    }
    # Add unique IDs
    my $ids = $unique_ids[$lib_num];
    for my $rank (1 .. $diversity) {
      next if defined $ranks[$rank-1];
      $ranks[$rank-1] = pop @$ids;   
    }
    push @c_ranks, \@ranks;
  }

  return \@c_ranks, $overall_diversity, $diversities, $perc_shared;
}


sub community_calculate_species_abundance {
  # Calculate relative abundance based on a distribution and its parameters.
  # Input is a model, its 2 parameters, and the number of values to generate
  # Output is a reference to a list of relative abundance. The abundance adds up
  # to 1
  my ($distrib, $param, $diversity) = @_;
  # First calculate rank-abundance values
  my $rel_ab;
  if ($distrib eq 'uniform') {
    # no parameter
    my $val = 1 / $diversity;
    for (my $index = 0 ; $index < $diversity ; $index++) {
      $$rel_ab[$index] = $val;
    }
  } elsif ($distrib eq 'linear') {
    # no parameter
    my $slope = 1 / $diversity;
    my $total = 0;
    for (my $index = 0 ; $index < $diversity ; $index++) {
      $$rel_ab[$index] = 1 - $slope * $index;
      $total += $$rel_ab[$index];
    }
    for (my $index = 0 ; $index < $diversity ; $index++) {
      $$rel_ab[$index] /= $total;
    }
  } elsif ($distrib eq 'powerlaw') {
    # 1 parameter
    die "Error: The powerlaw model requires an input parameter (-p option)\n"
      if not defined $param;
    my $total = 0;
    for (my $index = 0 ; $index < $diversity ; $index++) {
      $$rel_ab[$index] = ($index+1)**-$param;
      $total += $$rel_ab[$index];
    }
    for (my $index = 0 ; $index < $diversity ; $index++) {
      $$rel_ab[$index] /= $total;
    }
  } elsif ($distrib eq 'logarithmic') {
    # 1 parameter
    die "Error: The logarithmic model requires an input parameter (-p option)\n"
      if not defined $param;
    my $total = 0;
    for (my $index = 0 ; $index < $diversity ; $index++) {
      $$rel_ab[$index] = log($index+2)**-$param;
      $total += $$rel_ab[$index];
    }
    for (my $index = 0 ; $index < $diversity ; $index++) {
      $$rel_ab[$index] /= $total;
    }
  } elsif ($distrib eq 'exponential') {
    # 1 parameter
    die "Error: The exponential model requires an input parameter (-p option)\n"
      if not defined $param;
    my $total = 0;
    for (my $index = 0 ; $index < $diversity ; $index++) {
      $$rel_ab[$index] = exp(-($index+1)*$param);
      $total += $$rel_ab[$index];
    }
    for (my $index = 0 ; $index < $diversity ; $index++) {
      $$rel_ab[$index] /= $total;
    }
  } else {
    die "Error: $distrib is not a valid rank-abundance distribution\n";
  }
  return $rel_ab;
}


sub community_calculate_amplicon_abundance {
  my ($r_spp_abs, $r_spp_ids, $seq_ids) = @_;
  # Convert abundance of species into abundance of their amplicons because there
  # can be multiple amplicon per species and the amplicons have a different ID
  # from the species. The r_spp_ids and r_spp_abs arrays are the ID and abundance
  # of the species, sorted by decreasing abundance.

  # Give amplicons from the same species the same sampling probability
  my $sum = 0;
  for (my $i  = 0; $i < scalar @$r_spp_ids; $i++) {
    my $species_ab    = $$r_spp_abs[$i];
    my $species_id    = $$r_spp_ids[$i];
    my @amplicon_ids  = keys %{$seq_ids->{$species_id}};
    my $nof_amplicons = scalar @amplicon_ids;
    my @amplicon_abs  = ($species_ab) x $nof_amplicons;
    splice @$r_spp_abs, $i, 1, @amplicon_abs;
    splice @$r_spp_ids, $i, 1, @amplicon_ids;
    $sum += $species_ab * $nof_amplicons;
    $i += $nof_amplicons - 1;
  }

  # Normalize the relative abundance
  if ($sum != 1) {
    for my $i (0 .. scalar @$r_spp_abs - 1) {
      $$r_spp_abs[$i] /= $sum;
    }
  }

  return $r_spp_abs, $r_spp_ids;
}


sub next_single_read {
  # Generate a single shotgun or amplicon read
  my ($self) = @_;
  my $oids           = $self->{c_structs}->[$self->{cur_lib}-1]->{ids};
  my $mid            = $self->{multiplex_ids}->[$self->{cur_lib}-1];
  my $lib_num        = $self->{num_libraries} > 1 ? $self->{cur_lib} : undef;
  my $max_nof_tries  = $self->{forward_reverse} ? 1 : 10;

  # Choose a random genome or amplicon
  my $genome = $self->rand_seq($self->{positions}, $oids);
  my $nof_tries = 0;
  my $shotgun_seq;
  do {
    # Error if we have exceeded the maximum number attempts
    $nof_tries++;
    if ($nof_tries > $max_nof_tries) {
      my $message = "Error: Could not take a random shotgun read without ".
        "forbidden characters from reference sequence ".$genome->id;
      $message .= " ($max_nof_tries attempts made)" if ($max_nof_tries > 1);
      $message .= ".\n";
      die $message;
    }
    # Take a random orientation if needed
    my $orientation = ($self->{unidirectional} != 0) ? 1 : rand_seq_orientation();
    # Choose a read size according to the specified distribution
    my $length = rand_seq_length($self->{read_length}, $self->{read_model},
      $self->{read_delta});
    # Shorten read length if too long
    $length = $genome->length if $length > $genome->length;
    # Read position on genome or amplicon
    my ($start, $end) = rand_seq_pos($genome, $length, $self->{forward_reverse},
      $mid);
    # Chimerize the template sequence if needed
    $genome = $self->rand_seq_chimera($genome, $self->{chimera_perc}, $start,
      $end, $self->{positions}, $oids) if $self->{chimera_perc};
    # New sequence object
    $shotgun_seq = new_subseq($self->{cur_read}, $genome, $self->{unidirectional},
      $orientation, $start, $end, $mid, undef, $lib_num, $self->{desc_track},
      $self->{qual_levels});
    # Simulate sequence aberrations and sequencing error if needed
    $shotgun_seq = $self->rand_seq_errors($shotgun_seq)
      if ($self->{homopolymer_dist} || $self->{mutation_freq});
  } while (
    $self->{exclude_chars} && not is_valid($shotgun_seq, $self->{exclude_chars})
  );
  return $shotgun_seq;
}


sub next_mate_pair {
  # Generate a shotgun mate pair
  my ($self) = @_;
  my $oids           = $self->{c_structs}->[$self->{cur_lib}-1]->{ids};
  my $mid            = $self->{multiplex_ids}->[$self->{cur_lib}-1];
  my $lib_num        = $self->{num_libraries} > 1 ? $self->{cur_lib} : undef;
  my $pair_num       = int( $self->{cur_read} / 2 + 0.5 );
  my $max_nof_tries  = $self->{forward_reverse} ? 1 : 10;

  # Choose a random genome
  my $genome = $self->rand_seq($self->{positions}, $oids);

  my $nof_tries      = 0;
  my ($shotgun_seq_1, $shotgun_seq_2);
  while (1) {
    # Error if we have exceeded the maximum number of attempts
    $nof_tries++;
    if ($nof_tries > $max_nof_tries) {
      my $message = "Error: Could not take a pair of random shotgun read ".
        "without forbidden characters from reference sequence ".$genome->id;
      $message .= " ($max_nof_tries attempts made)" if ($max_nof_tries > 1);
      $message .= ".\n";
      die $message;
    }
    # Take a random orientation if needed
    my $orientation = ($self->{unidirectional} != 0) ? 1 : rand_seq_orientation();
    # Choose a mate pair length according to the specified distribution
    my $mate_length = rand_seq_length($self->{mate_length}, $self->{mate_model},
      $self->{mate_delta});
    # Shorten mate length if too long
    $mate_length = $genome->length if $mate_length > $genome->length;
    # Mate position on genome or amplicon
    my ($mate_start, $mate_end) = rand_seq_pos($genome, $mate_length,
      $self->{forward_reverse}, $mid);
    # Chimerize the template sequence if needed
    $genome = $self->rand_seq_chimera($genome, $self->{chimera_perc},
      $mate_start, $mate_end, $self->{positions}, $oids) if $self->{chimera_perc};
    # First mate read
    my $read_length = rand_seq_length($self->{read_length}, $self->{read_model},
      $self->{read_delta});
    my $seq_start   = $mate_start;
    my $seq_end     = $mate_start + $read_length - 1;
    $shotgun_seq_1  = new_subseq($pair_num, $genome, $self->{unidirectional},
      $orientation, $seq_start, $seq_end, $mid, '1', $lib_num, $self->{desc_track},
      $self->{qual_levels});
    $shotgun_seq_1 = $self->rand_seq_errors($shotgun_seq_1)
      if ($self->{homopolymer_dist} || $self->{mutation_freq});
    if ($self->{exclude_chars} && not is_valid($shotgun_seq_1, $self->{exclude_chars})) {
      next;
    }
    # Second mate read
    $read_length   = rand_seq_length($self->{read_length}, $self->{read_model},
      $self->{read_delta});
    $seq_start     = $mate_end - $read_length + 1;
    $seq_end       = $mate_end;
    $shotgun_seq_2 = new_subseq($pair_num, $genome, $self->{unidirectional},
      $orientation, $seq_start, $seq_end, $mid, '2', $lib_num, $self->{desc_track},
      $self->{qual_levels});
    $shotgun_seq_2 = $self->rand_seq_errors($shotgun_seq_2)
      if ($self->{homopolymer_dist} || $self->{mutation_freq});
    if ($self->{exclude_chars} && not is_valid($shotgun_seq_2, $self->{exclude_chars})) {
      next;
    }
    # Both shotgun reads were valid
    last;
  }
  return $shotgun_seq_1, $shotgun_seq_2;
}


sub is_valid {
  # Return 1 if the sequence object is valid (is not empty and does not have any
  # of the specified forbidden characters), 0 otherwise. Specify the forbidden
  # characters as a single string, e.g. 'N-' to prevent any reads to have 'N' or
  # '-'. The search is case-insensitive.
  my ($seq, $fchars) = @_;
  if ( (not defined $fchars) || ($fchars eq '') ) {
    # no forbidden chars
    return 1;
  }
  if ((not defined $seq) || ($seq->seq eq '') || ($seq->seq =~ /[$fchars]/gi)) {
    # invalid sequence
    return 0;
  }
  # sequence passed all checks
  return 1;
}


sub proba_create {
  my ($self, $c_struct, $size_dep, $copy_bias) = @_;
  # 1/ Calculate size-dependent, copy number-dependent probabilities
  my $probas = $self->proba_bias_dependency($c_struct, $self->{database}->{db},
    $size_dep, $copy_bias);

  # 2/ Generate proba starting position
  my $positions = $self->proba_cumul($probas);
  return $positions;
}


sub proba_bias_dependency {
  # Affect probability of picking a species by considering genome length or gene
  # copy number bias
  my ($self, $c_struct, $seq_db, $size_dep, $copy_bias) = @_;

  # Calculate probability
  my @probas;
  my $totproba = 0;
  my $diversity = scalar @{$c_struct->{'ids'}};
  for my $i (0 .. scalar $diversity - 1) {
    my $proba = $c_struct->{'abs'}[$i];

    if ( defined $self->{forward_reverse} ) {
      # Gene copy number bias
      if ($copy_bias) {
        my $refseq_id = $self->database_get_parent_id($c_struct->{'ids'}[$i]);
        my $nof_amplicons = scalar @{ $self->database_get_children_seq($refseq_id) };
        $proba *= $nof_amplicons;
      }
    } else {
      # Genome length bias
      if ($size_dep) {
        my $id  = $c_struct->{'ids'}[$i];
        my $seq = $self->database_get_seq($id);
        my $len = $seq->length;
        $proba /= $len;
      }
    }

    push @probas, $proba;
    $totproba += $proba;
  }

  # Normalize if necessary
  if ($totproba != 1) {
    for my $i (0 .. scalar $diversity - 1) {
      $probas[$i] /= $totproba;
    }
  }

  return \@probas;
}


sub proba_cumul {
  # Put the probas end to end on a line and generate their start position on the
  # line (cumulative distribution). This will help with picking genomes or 
  # nucleotides at random using the rand_weighted() subroutine.
  my ($self, $probas) = @_;
  my $sum = 0;
  return [ 0, map { $sum += $_ } @$probas ];
}


sub rand_weighted {
  # Pick a random number based on the given cumulative probabilities.
  # Cumulative weights can be obtained from the proba_cumul() subroutine.
  my ($cum_probas, $pick, $index) = (shift, rand, -1);
  map { $pick >= $_ ? $index++ : return $index } @$cum_probas;
}


sub rand_seq {
  # Choose a sequence object randomly using a probability distribution
  my ($self, $positions, $oids) = @_;
  return $self->database_get_seq( $$oids[rand_weighted($positions)] ); 
}


sub rand_seq_chimera {
  my ($self, $sequence, $chimera_perc, $start, $end, $positions, $oids) = @_;
  # Produce an amplicon that is a chimera of two sequences, starting with the
  # input sequence until a random position between the given start and end, and
  # ending with the end of another sequence taken at random from the database
  # and going at least as far as the specified end
  my $chimera;
  # Sanity check
  if ( (scalar @$oids < 2) && ($chimera_perc > 0) ) {
    die "Error: Not enough sequences to produce chimeras\n";
  }
  # Fate now decides to produce a chimera or not
  if ( rand(100) <= $chimera_perc ) { 
    my $t1_seq = $sequence; # first template sequence
    my $t2_seq;             # second template sequence
    do {
      $t2_seq = $self->rand_seq($positions, $oids);
    } while ($t2_seq->id eq $t1_seq->id);

    my $t1_start = 1;
    my $t2_end   = $t2_seq->length;

    my $diff = $end - $t2_end;
    $start = $diff if $diff > 0;

    my $t1_end   = $start + int( rand($end-$start) ); # start <= t1_end < end
    my $t2_start = $t1_end - $diff + 1;

    # Join chimera fragments
    $chimera      = $t1_seq->trunc($t1_start, $t1_end);
    $chimera->seq( $chimera->seq . $t2_seq->subseq($t2_start, $t2_end) );
    $chimera->id( $chimera->id . ',' . $t2_seq->id );
    $chimera->{_amplicon} = $t1_seq->{_amplicon}.','.$t2_seq->{_amplicon};

  } else {
    # No chimera needed
    $chimera = $sequence;
  }
  return $chimera;
}


sub rand_seq_orientation {
  # Return a random read orientation: 1 for uncomplemented, or -1 for complemented
  return int(rand()+0.5) ? 1 : -1;
}


sub rand_seq_errors {
  # Introduce sequencing errors (point mutations, homopolymers) in a sequence
  # based on error models
  my ($self, $seq) = @_;
  my $seq_str = $seq->seq();
  my $error_specs = {}; # Error specifications

  # First, specify errors in homopolymeric stretches
  $error_specs = $self->rand_homopolymer_errors($seq_str, $error_specs)
    if $self->{homopolymer_dist};

  # Then, specify point sequencing errors: substitutions, insertions, deletions
  $error_specs = $self->rand_point_errors($seq_str, $error_specs) 
    if $self->{mutation_freq};

  # Finally, actually implement the errors as per the specifications
  $seq->errors($error_specs) if (scalar keys %$error_specs > 0);

  return $seq;
}


sub rand_homopolymer_errors {
  # Specify sequencing errors in a sequence's homopolymeric stretches
  my ($self, $seq_str, $error_specs) = @_;
  while ( $seq_str =~ m/(.)(\1+)/g ) {

    # Found a homopolymer
    my $res = $1;                       # residue in homopolymer
    my $len = length($2) + 1;           # length of the homopolymer
    my $pos = pos($seq_str) - $len + 1; # start of the homopolymer (residue no.)

    # Apply homopolymer model based on normal distribution N(mean, standard deviation)
    #   Balzer:    N(n, 0.03494 + n * 0.06856)  Balzer et al. 2010
    #   Richter:   N(n, 0.15 * sqrt(n))         Richter et al. 2008
    #   Margulies: N(n, 0.15 * n)               Margulies et al. 2005
    my ($stddev, $new_len, $diff) = (0, 0, 0);
    if ( $self->{homopolymer_dist} eq 'balzer' ) {
      $stddev = 0.03494 + $len * 0.06856;
    } elsif ($self->{homopolymer_dist} eq 'richter') {
      $stddev = 0.15 + sqrt($len);
    } elsif ($self->{homopolymer_dist} eq 'margulies') {
      $stddev = 0.15 + $len;
    } else {
      die "Error: Unknown homopolymer distribution '".$self->{homopolymer_dist}."'\n";
    }
    $new_len = int( $len + $stddev * randn() + 0.5 );
    $new_len = 0 if $new_len < 0;
    # We're done if no error was introduced
    $diff = $new_len - $len;
    next unless $diff;
    # Otherwise, track the error generated
    if ($diff > 0) { # Homopolymer extension
      my $ext = $res x $diff;
      if (exists $$error_specs{$pos}{'+'}) {
        $$error_specs{$pos}{'+'} .= $ext;
      } else {
        $$error_specs{$pos}{'+'} = $ext;
      }
    } elsif ($diff < 0) { # Homopolymer shrinkage
      for my $offset ( 0 .. abs($diff)-1 ) {
        $$error_specs{$pos+$offset}{'-'} = '';
      }
    }

  }

  return $error_specs;
}


sub rand_point_errors {
  # Do some random point sequencing errors on a sequence based on a model
  my ($self, $seq_str, $error_specs) = @_;

  my $seq_len = length($seq_str);

  # Number of mutations to make in this sequence is assumed to follow a Normal
  # distribution N( mutation_freq, 0.3 * mutation_freq )
  my $read_mutation_freq = $self->{mutation_freq} + 0.3 * $self->{mutation_freq} * randn();
  my $nof_mutations = int( $seq_len*$read_mutation_freq/100 + 0.5 );
 
  # Exit without doing anything if there are no mutations to do
  return $error_specs if $nof_mutations == 0;

  # Mutation cumulative density functions (cdf)
  my $mut_pdf = [ 0 .. $seq_len - 1 ]; # probability density function
  if ($self->{mutation_model} eq 'uniform') {
    my $proba = 1 / $seq_len;
    $mut_pdf = [ map { $proba } (@$mut_pdf) ];
  } elsif ($self->{mutation_model} eq 'linear') {
    my $start = (2 * $self->{mutation_freq} - $self->{mutation_end}) / ($seq_len * $self->{mutation_freq});
    if ($start < 0) {
      die "Error: A 3' end mutation frequency of ".$self->{mutation_end}.
      " % is not possible in combination with an average mutation frequency".
      " of ".$self->{mutation_freq}." %\n";
    }
    my $slope = 2 * ($self->{mutation_end} - $self->{mutation_freq}) / (($seq_len-1) * $seq_len * $self->{mutation_freq});
    $mut_pdf = [ map { $start + $_ * $slope } (@$mut_pdf) ];
  } else {
    die "Error: '".$self->{mutation_model}."' is not a supported error distribution\n";
  }
  my $mut_cdf = $self->proba_cumul($mut_pdf);

  # Make as many mutations in read as needed based on model
  my $subst_frac = $self->{mutation_ratio} / ($self->{mutation_ratio} + 1);
  for ( 1 .. $nof_mutations ) {

    # Position to mutate
    my $idx = rand_weighted($mut_cdf);

    # Do a substitution or indel
    if ( rand() <= $subst_frac ) {
      # Substitute at given position by a random replacement nucleotide
      $$error_specs{$idx+1}{'%'} = rand_nuc( substr($seq_str, $idx, 1) );

    } else {
      # Equiprobably insert or delete
      if ( rand() < 0.5 ) {
        # Insertion after given position
        my $add = rand_nuc();
        if (exists $$error_specs{$idx+1}{'+'}) {
          $$error_specs{$idx+1}{'+'} .= $add;
        } else {
          $$error_specs{$idx+1}{'+'} = $add;
        }
      } else {
        # Make a deletion at given position
        next if length($seq_str) == 1; # skip this deletion to avoid a 0 length
        $$error_specs{$idx+1}{'-'} = '';
      }
    }

  }

  return $error_specs;
}

sub rand_nuc {
  # Pick a nucleotide at random. An optional nucleotide to exclude can be given.
  my $not_nuc = shift;
  my @cdf;
  my @nucs;
  if (defined $not_nuc) {
    # Exclude a specific nucleotide from the search
    my %nucs_hash = ( 'A' => undef,
                      'C' => undef,
                      'G' => undef,
                      'T' => undef  );
    delete $nucs_hash{uc($not_nuc)};
    @nucs = keys %nucs_hash;
    @cdf  = (0, 0.33333333, 0.66666666, 1);
  } else {
    # All nucleotides are possible
    @nucs = ('A', 'C', 'G', 'T');
    @cdf  = (0, 0.25, 0.5, 0.75, 1);
  }
  return $nucs[rand_weighted(\@cdf)]; 
}



sub rand_seq_length {
  # Choose the sequence length following a given probability distribution
  my($avg, $model, $stddev) = @_;
  my $length;
  if (not $model) {
    # No specified distribution: all the sequences have the length of the average
    $length = $avg;
  } else {
    if ($model eq 'uniform') {
      # Uniform distribution: decimal number uniformly distributed in [min, max)
      my ($min, $max) = ($avg - $stddev, $avg + $stddev);
      $length = $min + int( ($max - $min + 1) * rand() );
    } elsif ($model eq 'normal') {
      # Gaussian distribution: decimal number normally distribution in N(avg,stddev)
      $length = $avg + $stddev * randn();
      $length = int( $length + 0.5 );
    } else {
      die "Error: '$model' is not a supported read or insert length distribution\n";
    }
  }
  $length = 1 if ($length < 1);
  return $length;
}


sub rand_seq_pos {
  # Pick the coordinates (start and end) of an amplicon or random shotgun read.
  # Coordinate system: the first base is 1 and the number is inclusive, ie 1-2
  # are the first two bases of the sequence
  my ($seq_obj, $read_length, $amplicon, $mid) = @_;
  # Read length includes the MID
  my $length;
  if (defined $mid) {
    $length = $read_length - length($mid);
  } else {
    $length = $read_length;
  }
  # Pick starting position
  my $start;
  if (defined $amplicon) {
    # Amplicon always start at first position of amplicon
    $start = 1;
  } else {
    # Shotgun reads start at a random position in genome
    $start = int( rand($seq_obj->length - $length + 1) ) + 1;
  }
  # End position
  my $end = $start + $length - 1;
  return $start, $end;
}


sub randn {
  # Normally distributed random value (mean 0 and standard deviation 1) using
  # the Box-Mueller transformation method, adapted from the Perl Cookbook
  my ($g1, $g2, $w);
  do {
    $g1 = 2 * rand() - 1; # uniformly distributed
    $g2 = 2 * rand() - 1;
    $w = $g1**2 + $g2**2; # variance
  } while ( $w >= 1 );
  $w = sqrt( (-2 * log($w)) / $w ); # weight
  $g1 *= $w; # gaussian-distributed
  if ( wantarray ) {
    $g2 *= $w;
    return ($g1, $g2);
  } else {
    return $g1;
  }
}


sub randig {
   # Random value sampled from the inverse gaussian (a.k.a. Wald) distribution,
   # using the method at http://en.wikipedia.org/wiki/Inverse_Gaussian_distribution
   my ($mu, $lambda) = @_;
   my $y = randn()**2;
   my $x = $mu + ($mu**2 * $y)/(2 * $lambda) - $mu / (2 * $lambda)
           * sqrt(4 * $mu * $lambda * $y + $mu**2 * $y**2);
   if ( rand() <= $mu / ($mu + $x) ) {
      $y = $x;
   } else {
      $y = $mu**2 / $x;
   }
   return $y;
}


sub randomize {
  # Randomize an array using the Fisher-Yates shuffle described in the Perl
  # cookbook.
  my ($array) = @_;
  my $i;
  for ($i = @$array; --$i; ) {
   my $j = int rand($i+1);
   next if $i == $j;
   @$array[$i,$j] = @$array[$j,$i];
  }
  return $array;
}


sub database_create {
  # Read and import sequences
  # Parameters:
  #   * FASTA file containing the sequences
  #   * Sequencing unidirectionally? 0: no, 1: yes forward, -1: yes reverse
  #   * Amplicon PCR primers (optional): Should be provided in a FASTA file and
  #     use the IUPAC convention. If a primer sequence is given, any sequence
  #     that does not contain the primer (or its reverse complement for the
  #     reverse primer) is skipped, while any sequence that match is trimmed so
  #     that it is flush with the primer sequence.
  #   * Abundance file (optional): To avoid registering sequences in the database
  #     unless they are needed
  #   * Delete chars (optional): Characters to delete form the sequences.
  my ($self, $fasta_file, $unidirectional, $forward_reverse_primers,
    $abundance_file, $delete_chars) = @_;
  # Input filehandle
  my $in;
  if ($fasta_file eq '-') {
    $in = Bio::SeqIO->newFh(
      -fh     => \*STDIN,
      -format => 'fasta',
    );
  } else {
    $in = Bio::SeqIO->newFh(
      -file   => $fasta_file,
      -format => 'fasta',
    );
  }
  # Regular expression to catch amplicons present in the database
  my ($forward_regexp, $reverse_regexp);
  if (defined $forward_reverse_primers) {
    # Read primers from FASTA file
    my $primer_in = Bio::SeqIO->newFh(
      -file   => $forward_reverse_primers,
      -format => 'fasta',
    );
    my $first_primer = <$primer_in>;
    if (not defined $first_primer) {
      die "Error: The file '$forward_reverse_primers' does not contain any primers\n";
    } else {
      # Force the alphabet since degenerate primers can look like protein sequences
      $first_primer->alphabet('dna');
    }
    my $second_primer = <$primer_in>;
    if (defined $second_primer) {
      $second_primer->alphabet('dna');
    }
    undef $primer_in;
    # Regexp for forward primer and optional reverse complement of reverse primer
    $forward_regexp = iupac_to_regexp($first_primer->seq);
    if (defined $second_primer) {
      $second_primer = $second_primer->revcom;
      $reverse_regexp = iupac_to_regexp($second_primer->seq);
    }
  }
  # Get list of all IDs with a manually-specified abundance
  my %ids_to_keep;
  if ($abundance_file) {
    my ($ids) = community_read_abundances($abundance_file);
    for my $comm_num (0 .. $#$ids) {
      for my $gen_num ( 0 .. scalar @{$$ids[$comm_num]} - 1 ) {
        my $id = $$ids[$comm_num][$gen_num];
        $ids_to_keep{$id} = undef;
      }
    }
  }
  # Process database sequences
  my %seq_db;        # hash of BioPerl sequence objects (all amplicons)
  my %seq_ids;       # hash of reference sequence IDs and IDs of their amplicons
  while ( my $ref_seq = <$in> ) {
    # Skip unwanted sequences
    my $ref_seq_id = $ref_seq->id;
    next if (scalar keys %ids_to_keep > 0) && (not exists $ids_to_keep{$ref_seq_id});
    # If we are sequencing from the reverse strand, reverse complement now
    if ($unidirectional == -1) {
      $ref_seq = $ref_seq->revcom;
    }
    # Extract amplicons if needed
    my $amp_seqs;
    if (defined $forward_regexp) {
      $amp_seqs = $self->database_extract_amplicons($ref_seq, $forward_regexp,
        $reverse_regexp, \%ids_to_keep);
      next if scalar @$amp_seqs == 0;
    } else {
      $amp_seqs = [$ref_seq];
    }

    for my $amp_seq (@$amp_seqs) {
      # Remove forbidden chars
      if ( (defined $delete_chars) && (not $delete_chars eq '') ) {
        my $clean_seq = $amp_seq->seq;
        $clean_seq =~ s/[$delete_chars]//gi;
        $amp_seq->seq($clean_seq);
      }
      # Save amplicon sequence and identify them by their unique object reference
      $seq_db{$amp_seq} = $amp_seq;
      $seq_ids{$ref_seq_id}{$amp_seq} = undef;
    }

  }
  undef $in; # close the filehandle (maybe?!)

  # Sanity check
  if (scalar keys %seq_ids == 0) {
    die "Error: No genome sequences could be used. If you specified a file of".
      " abundances for the genome sequences, make sure that their ID match the".
      " ID in the FASTA file. If you specified amplicon primers, verify that ".
      "they match some genome sequences.\n";
  }

  my $database = { 'db' => \%seq_db, 'ids' => \%seq_ids };
  return $database;
}


sub database_extract_amplicons {
  my ($self, $seq, $forward_regexp, $reverse_regexp, $ids_to_keep) = @_;
  # A database sequence can have several amplicons, e.g. a genome can have 
  # several 16S genes. Extract all amplicons from a sequence. Only the shortest
  # amplicons are returned, which means that for a sequence that matches the
  # primer twice, there are 3 primer combinations, but only two amplicons are
  # returned, not three.

  my $seqstr = $seq->seq;
  my $seqid  = $seq->id;
  my @amplicons;
  if ( (defined $forward_regexp) && (not defined $reverse_regexp) ) {
    while ( $seqstr =~ m/($forward_regexp)/g ) {
      my $start    = pos($seqstr) - length($1) + 1;
      my $end      = $seq->length;
      my $amplicon = $seq->trunc($start, $end);
      $amplicon->{_amplicon} = $start.'-'.$end;
      push @amplicons, $amplicon;
    }
  } elsif ( (defined $forward_regexp) && (defined $reverse_regexp) ) {
    while ( $seqstr =~ m/($forward_regexp.*?$reverse_regexp)/g ) {
      my $end      = pos($seqstr);
      my $start    = $end - length($1) + 1;
      my $amplicon = $seq->trunc($start, $end);
      $amplicon->{_amplicon} = $start.'-'.$end;
      push @amplicons, $amplicon;
    }
  } else {
    die "Error: Need to provide at least a forward primer\n";
  }
  # Complain if primers did not match explicitly specified database sequence
  if ( (scalar keys %{$ids_to_keep} > 0) &&
       (exists $$ids_to_keep{$seqid}   ) &&
       (scalar @amplicons == 0         ) ) {
    die "Error: Requested sequence $seqid does not match the specified forward primer.\n";
  }

  return \@amplicons;
}


sub database_get_seq {
  # Retrieve a sequence object from the database based on its object ID
  my ($self, $oid)  = @_;
  my $db = $self->{database}->{db};
  my $seq_obj;
  if (not exists $$db{$oid}) {
    warn "Warning: Could not find sequence with object ID '$oid' in the database\n";
  }
  $seq_obj = $$db{$oid};
  return $seq_obj;
}


sub database_get_children_seq {
  # Retrieve all the sequences object made from a reference sequence based on the
  # ID of the reference sequence
  my ($self, $refseqid)  = @_;
  my @children;
  for my $child_oid ( keys %{$self->{database}->{ids}->{$refseqid}} ) {
    push @children, $self->database_get_seq($child_oid);
  }
  return \@children;
}

sub database_get_parent_id {
  # Based on a sequence object ID, retrieve the ID of the reference sequence it
  # came from
  my ($self, $oid) = @_;
  my $seq_id = $self->database_get_seq($oid)->id;
###  $seq_id =~ s/_amplicon.*$//;
  return $seq_id;
}


sub iupac_to_regexp {
  # Create a regular expression to match a nucleotide sequence that contain
  # degeneracies (in IUPAC standard)
  my ($seq) = @_;
  # Basic IUPAC code
  #my %iupac = (
  #  'A' => ['A'],
  #  'C' => ['C'],
  #  'G' => ['G'],
  #  'T' => ['T'],
  #  'U' => ['U'],
  #  'R' => ['G', 'A'],
  #  'Y' => ['T', 'C'],
  #  'K' => ['G', 'T'],
  #  'M' => ['A', 'C'],
  #  'S' => ['G', 'C'],
  #  'W' => ['A', 'T'],
  #  'B' => ['G', 'T', 'C'],
  #  'D' => ['G', 'A', 'T'],
  #  'H' => ['A', 'C', 'T'],
  #  'V' => ['G', 'C', 'A'],
  #  'N' => ['A', 'G', 'C', 'T'],
  #);
  # IUPAC code
  #   + degenerate primer residues matching ambiguous template residues
  #   + degenerate primer residues matching uracil U
  my %iupac = (
    'A' => ['A'],
    'C' => ['C'],
    'G' => ['G'],
    'T' => ['T'],
    'U' => ['U'],
    'R' => ['G', 'A', 'R'],
    'Y' => ['T', 'U', 'C', 'Y'],
    'K' => ['G', 'T', 'U', 'K'],
    'M' => ['A', 'C', 'M'],
    'S' => ['G', 'C', 'S'],
    'W' => ['A', 'T', 'U', 'W'],
    'B' => ['G', 'T', 'U', 'C', 'Y', 'K', 'S', 'B'],
    'D' => ['G', 'A', 'T', 'U', 'R', 'K', 'W', 'D'],
    'H' => ['A', 'C', 'T', 'U', 'Y', 'M', 'W', 'H'],
    'V' => ['G', 'C', 'A', 'R', 'M', 'S', 'V'],
    'N' => ['A', 'G', 'C', 'T', 'U', 'R', 'Y', 'K', 'M', 'S', 'W', 'B', 'D', 'H', 'V', 'N'],
  );
  # Regular expression to catch this sequence
  my $regexp;
  for my $pos (0 .. length($seq)-1) {
    my $res = substr $seq, $pos, 1;
    my $iupacs = $iupac{$res};
    if (not defined $iupacs) {
      die "Error: Primer sequence '$seq' is not a valid IUPAC sequence. ".
        "Offending character is '$res'.\n";
    }
    if (scalar @$iupacs > 1) {
      $regexp .= '['.join('',@$iupacs).']';
    } else {
      $regexp .= $$iupacs[0];
    }
  }
  $regexp = qr/$regexp/i;
  return $regexp;
}


sub lib_coverage {
  # Calculate number of sequences needed to reach a given coverage. If the
  # number of sequences is provided, calculate the coverage
  my ($self, $c_struct) = @_;
  my $coverage    = $self->{coverage_fold};
  my $nof_seqs    = $self->{total_reads};
  my $read_length = $self->{read_length};
  # 1/ Calculate library length and size
  my $ref_ids    = $c_struct->{'ids'};
  my $diversity  = scalar @$ref_ids;
  my $lib_length = 0;
  for my $ref_id (@$ref_ids) {
    my $seqobj = $self->database_get_seq($ref_id);
    my $seqlen = $seqobj->length;
    $lib_length += $seqlen;
  }
  # 2/ Calculate number of sequences to generate based on desired coverage. If
  # both number of reads and coverage fold were given, coverage has precedence.
  if ($coverage) {
    $nof_seqs = ($coverage * $lib_length) / $read_length;
    if ( int($nof_seqs) < $nof_seqs ){
      $nof_seqs = int($nof_seqs + 1); # ceiling
    }
  }
  $coverage = ($nof_seqs * $read_length) / $lib_length;
  # 3/ Sanity check

  ####
  # Warn only if diversity was explicitely specified on the command line
  ####

  if ( $nof_seqs < $diversity) {
    warn "Warning: The number of reads to produce is lower than the required ".
      "diversity. Increase the coverage or number of reads to achieve this ".
      "diversity.\n";
    $self->{diversity}->[$self->{cur_lib}-1] = $nof_seqs;
  }
  return $nof_seqs, $coverage;
}


sub new_subseq {
  # Create a new sequence object as a subsequence of another one and name it so
  # we can trace back where it came from
  my ($fragnum, $seq_obj, $unidirectional, $orientation, $start, $end, $mid,
    $mate_number, $lib_number, $tracking, $qual_levels) = @_;
  # If the length is too short for this read, no choice but to decrease it.
  $start = 1 if $start < 1;
  $end   = $seq_obj->length if $end > $seq_obj->length;

  # Build the sequence ID
  my $name_sep  = '_';
  my $field_sep = ' ';
  my $mate_sep  = '/'; # mate pair indicator, by convention
  my $newid = $fragnum;
  if (defined $lib_number) {
    $newid = $lib_number.$name_sep.$newid;
  }
  if (defined $mate_number) {
    $newid .= $mate_sep.$mate_number;
  }

    # Create a new simulated read object
  my $newseq = Bio::Seq::SimulatedRead->new(
     -id          => $newid,
     -reference   => $seq_obj,
     -start       => $start,
     -end         => $end,
     -strand      => $orientation,
     -mid         => $mid,
     -track       => $tracking,
     -qual_levels => $qual_levels,
  );
  
  # Record location of amplicon on reference sequence in the sequence description
  my $amplicon_desc = $seq_obj->{_amplicon};
  if (defined $amplicon_desc) {
    $amplicon_desc = 'amplicon='.$amplicon_desc;
    my $desc = $newseq->desc;
    $desc =~ s/(reference=\S+)/$1 $amplicon_desc/;
    $newseq->desc($desc);
  }

  # Database genomes were already reverse complemented if reverse sequencing was requested
  if ($unidirectional == -1) {
    if ($orientation == -1) {
      $orientation = '+1';
    } else {
      $orientation = '-1';
    }
    $newseq->strand($orientation);
    my $new_desc = $newseq->desc;
    $newseq->desc( $new_desc =~ s/strand=(\S+)/strand=$orientation/g );
  }

  return $newseq;
}


sub two_array_sort {
  # Sort 2 arrays by taking the numeric sort of the first one and keeping the 
  # element of the second one match those of the first one
  my ($l1, $l2) = @_;
  my @ids = map { [ $$l1[$_], $$l2[$_] ] } (0..$#$l1);
  @ids = sort { $a->[0] <=> $b->[0] } @ids;
  my @k1;
  my @k2;
  for (my $i = 0; $i < scalar @ids; $i++) {
    $k1[$i] = $ids[$i][0];
    $k2[$i] = $ids[$i][1];
  }
  return \@k1, \@k2;
}


1;

