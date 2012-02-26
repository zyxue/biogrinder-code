package Grinder::AmpliconSearch;

use strict;
use warnings;
use Bio::SeqIO;

use base qw(Bio::Root::Root); # using throw() and _rearrange() methods

=head1

Search for amplicons in a given sequence. This is akin to an in silico, or simulated
PCR reaction. Given a primer or a set of degenerate primers, this module makes a regexp and
extract any amplicon in the given reference sequence

=cut

sub new {
   my ($class, @args) = @_;
   my $self = $class->SUPER::new(@args);
   my ($template, $primer_file, $forward_primer, $reverse_primer) =
      $self->_rearrange([qw(TEMPLATE PRIMER_FILE FORWARD_PRIMER REVERSE_PRIMER)],
      @args);

   if (defined $primer_file) {
      $self->_set_primers($primer_file);
   } else {
      $self->_set_forward_primer($forward_primer) if defined $forward_primer;
      $self->_set_reverse_primer($reverse_primer) if defined $reverse_primer;
   }

   $self->_set_template($template) if defined $template;

   return $self;
}


sub _initialize {
   ### convert forward primer to regexp
   #$self->_set_forward_regexp( _iupac_to_regexp($forward_primer->seq) );
   #if (defined $reverse_primer) {
   #   $primer = $reverse_primer->revcom;
   #   $self->_set_reverse_regexp( _iupac_to_regexp($reverse_primer->seq) );
   #}
   ### take optional reverse primer, reverse complement it and convert to regexp
   ### start regexp iterator
}


sub get_template {
   my ($self) = @_;
   return $self->{'template'};
}


sub _set_template {
   my ($self, $val) = @_;
   $self->{'template'} = $val;
   return $self->get_template;
}


sub _set_primers {
   my ($self, $primer_file) = @_;
   # Read primer file and convert primers into regular expressions to catch
   # amplicons present in the database

   if (not defined $primer_file) {
      $self->throw("Need to provide an input file\n");
   }

   # Read primers from sequence file
   my $in = Bio::SeqIO->newFh( -file => $primer_file );

   # Mandatory first primer
   my $primer = <$in>;
   if (not defined $primer) {
      $self->throw("The file '$primer_file' contains no primers\n");
   }
   $primer->alphabet('dna'); # Force the alphabet since degenerate primers can look like protein sequences
   $self->_set_forward_primer($primer);


   $primer = undef;

   # Optional reverse primers
   $primer = <$in>;
   if (defined $primer) {
      $primer->alphabet('dna');
      $self->_set_reverse_primer($primer);
   }

   return 1;
}


sub get_forward_primer {
   my ($self) = @_;
   return $self->{'forward_primer'};
}


sub _set_forward_primer {
   my ($self, $val) = @_;
   $self->{'forward_primer'} = $val;
   return $self->get_forward_primer;
}


sub get_reverse_primer {
   my ($self) = @_;
   return $self->{'reverse_primer'};
}


sub _set_reverse_primer {
   my ($self, $val) = @_;
   $self->{'reverse_primer'} = $val;
   return $self->get_reverse_primer;
}


sub get_forward_regexp {
   my ($self) = @_;
   return $self->{'forward_regexp'};
}


sub _set_forward_regexp {
   my ($self, $val) = @_;
   $self->{'forward_regexp'} = $val;
   return $self->get_forward_regexp;
}


sub get_reverse_regexp {
   my ($self) = @_;
   return $self->{'reverse_regexp'};
}


sub _set_reverse_regexp {
   my ($self, $val) = @_;
   $self->{'reverse_regexp'} = $val;
   return $self->get_reverse_regexp;
}


sub next_amplicon {
   my ($self) = @_;
   $self->_initialize if not defined $self->get_forward_regexp();

   #### retrieve next amplicon or return undef
}


#### Use Bio::Tools::SeqPattern
sub _iupac_to_regexp {
  my ($val) = @);
  return $val;
}
####

1;
