use inc::Module::Install;

name     'Grinder';
all_from 'lib/Grinder.pm';

build_requires 'Test::More';

# Regular Perl module dependencies for GAAS
requires 'Getopt::Euclid'          => '0.2.4';
requires 'Math::Random::MT::Perl'  => '1.06';
requires 'Bio::SeqIO'              =>  0;
#required 'Bio::Seq::SimulatedRead' =>  0; # required but packaged with Grinder since it is so recent

install_script 'script/Grinder';
install_script 'utils/average_genome_size';
install_script 'utils/change_paired_read_orientation';

auto_install;

WriteAll;

