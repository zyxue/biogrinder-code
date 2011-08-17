use inc::Module::Install;
use File::Spec::Functions;


# Package information

name     'Grinder';
all_from 'lib/Grinder.pm';

resources
   homepage   'http://sourceforge.net/projects/biogrinder/';
   bugtracker 'http://sourceforge.net/tracker/?group_id=244196&atid=1124737';
   repository 'git://biogrinder.git.sourceforge.net/gitroot/biogrinder/biogrinder';


# Build the README from the module's POD using Perl's pod2text and pod2html
generate_readme_text( catfile('lib','Grinder.pm'), 'README'     );
generate_readme_html( catfile('lib','Grinder.pm'), 'README.htm' );


# Dependencies

build_requires 'Test::More';

# Regular Perl module dependencies for Grinder
requires 'Getopt::Euclid'          => '0.2.4';
requires 'Math::Random::MT::Perl'  => '1.06';
requires 'Bio::SeqIO'              =>  0;
#required 'Bio::Seq::SimulatedRead' =>  0; # required but packaged here since it is so recent


# Extra scripts to install

install_script 'script/grinder';
install_script 'utils/average_genome_size';
install_script 'utils/change_paired_read_orientation';

auto_install;

WriteAll;




sub generate_readme_text {
   my ($in, $out) = @_;
   `pod2text $in $out`;
   warn "Warning: Could not generate $out.\n$!\n" if $? == -1;
   return $?; # exit status
}

sub generate_readme_html {
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
