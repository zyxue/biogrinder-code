use inc::Module::Install;


# Package information

name     'Grinder';
all_from 'lib/Grinder.pm';
license  'gpl3'; # Module::Install 1.04 does not parse the GPL version number

resources
   homepage   'http://sourceforge.net/projects/biogrinder/';
   bugtracker 'http://sourceforge.net/tracker/?group_id=244196&atid=1124737';
   repository 'git://biogrinder.git.sourceforge.net/gitroot/biogrinder/biogrinder';


# Dependencies for everyone

build_requires  'Test::More';
requires        'Getopt::Euclid'             => '0.3.4';
requires        'Math::Random::MT'           => '1.16';
requires        'Bio::SeqIO'                 =>  0;
#requires       'Bio::Root::Root'            =>  0;
# Modules required, but packaged with Grinder since they are too recent to have had a release yet
#requires       'Bio::Seq::SimulatedRead'    =>  0;
#requires       'Bio::Tools::AmpliconSearch' =>  0;
#requires       'Bio::SeqFeature::SubSeq'    =>  0;


# Dependencies for authors

author_requires 'Module::Install';
author_requires 'Module::Install::AuthorRequires';
author_requires 'Module::Install::AutoLicense';
author_requires 'Module::Install::ReadmeFromPod' => 0.14;
author_requires 'Module::Install::AutoManifest';
author_requires 'Statistics::R'  => '0.21';


# Bundle dependencies
# This system does not support Build.PL based dependencies

#perl_version( '5.005' );
#auto_bundle_deps();


# Install dependencies

auto_install;


# Extra scripts to install

install_script 'script/grinder';
install_script 'utils/average_genome_size';
install_script 'utils/change_paired_read_orientation';


# Generate MANIFEST file

auto_manifest();


# Generate Makefile and META.yml files

WriteAll;


# Generate the LICENSE file

auto_license();


# Generate the README and manpage files from the POD docs

auto_doc();



#--------- UTILS --------------------------------------------------------------#


sub auto_doc {
   print "*** Building doc...\n";
   # Generate script/grinder.pod
   my @args = ($^X, '-Ilib', 'script/grinder', '--podfile');
   system(@args) == 0 or die "System call to '@args' failed:\n$?\n";
   my $grinder = 'script/grinder.pod';
   my $script1 = 'utils/average_genome_size';
   my $script2 = 'utils/change_paired_read_orientation';
   my $man_dir = 'man';
   if (not -d $man_dir) {
      mkdir $man_dir or die "Could not write folder $man_dir:\n$!\n";
   }
   my $clean   = 1;
   readme_from $grinder, $clean, 'txt', 'README';
   readme_from $grinder, $clean, 'htm', 'README.htm';
   readme_from $grinder, $clean, 'man', "$man_dir/grinder.1";
   readme_from $script1, $clean, 'man', "$man_dir/average_genome_size.1";
   readme_from $script2, $clean, 'man', "$man_dir/change_paired_read_orientation.1";
   return 1;
}
