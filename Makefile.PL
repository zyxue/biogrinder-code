use inc::Module::Install;


# Package information

name     'Grinder';
all_from 'lib/Grinder.pm';
license  'gpl3'; # As of version 1.01, Module::Install does not parse the GPL version number

resources
   homepage   'http://sourceforge.net/projects/biogrinder/';
   bugtracker 'http://sourceforge.net/tracker/?group_id=244196&atid=1124737';
   repository 'git://biogrinder.git.sourceforge.net/gitroot/biogrinder/biogrinder';


# Dependencies for everyone

build_requires  'Test::More';

requires        'Getopt::Euclid'          => '0.2.4';
requires        'Math::Random::MT'        =>  0;
requires        'Math::Random::MT::Perl'  => '1.06';
requires        'Bio::SeqIO'              =>  0;
#requires       'Bio::Seq::SimulatedRead' =>  0; # required but packaged here since it is so recent


# Dependencies for authors

author_requires 'Module::Install';
author_requires 'Module::Install::AuthorRequires';
author_requires 'Module::Install::AutoLicense';
author_requires 'Module::Install::ReadmeFromPod';
author_requires 'Module::Install::AutoManifest';
author_requires 'Statistics::R'  => '0.21';


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


# Generate the README files from the module's POD

auto_readme();




#--------- UTILS --------------------------------------------------------------#


sub auto_readme {
   #print "*** Building README files...\n";
   my $source = 'lib/Grinder.pm';
   readme_from $source;
   readme_html_from( $source, 'README.htm' );
   return 1;
}

sub readme_html_from {
   my ($in, $out) = @_;
   `pod2html --infile=$in --outfile=$out`;
   warn "Warning: Could not generate $out.\n$!\n" if $? == -1;
   rm_files(['pod2htmd.tmp', 'pod2htmi.tmp']);
   return $?; # exit status
}

sub rm_files {
   my ($files) = @_;
   for my $file (@$files) {
      if (-e $file) {
         unlink $file or warn "Warning: Could not remove file '$file'.\n$!\n";
      }
   }
   return 1;
}
