#!/usr/bin/env python

"""
Move files create by Grinder to a location where it is going to be recognized by
Galaxy as multiple output files with the right format. See
http://wiki.g2.bx.psu.edu/Admin/Tools/Multiple Output Files
Example: python grinder_move_outputs output_dir output_id
Author: Florent Angly
"""

import sys, os, re

assert sys.version_info[:2] >= ( 2, 4 )

def stop_err( msg ):
    sys.stderr.write( "%s\n" % msg )
    sys.exit()

def __main__():
    # Get output dir and ID
    args = sys.argv
    output_dir = args[1]
    output_id  = args[2]

    # Move Grinder files to the proper output
    # Grinder filenames look like this
    #   grinder-ranks.txt
    #   grinder-reads.fa
    #   grinder-reads.qual
    #   grinder-1-ranks.txt
    #   grinder-1-reads.fa
    #   grinder-1-reads.qual
    #   grinder-2-ranks.txt
    #   grinder-2-reads.fa
    #   grinder-2-reads.qual

    p = re.compile(output_id)
    q = re.compile('-(\d+)-')
    r = re.compile('-(\w+)$')
    

    for fname in os.listdir(output_dir):

        # Skip files that do not start with the output_id
        source = os.path.join( output_dir, fname )
        basename, extension = os.path.splitext(fname)
        if not p.match(fname):
           continue

        # Assign the dataset format
        if extension == '.txt': 
           format = 'text'
        elif extension == '.fq':
           format = 'fastqsanger'
        elif extension == '.fastq':
           format = 'fastqsanger'
        elif extension == '.fa':
           format = 'fasta'
        elif extension == '.fna':
           format = 'fasta'
        elif extension == '.faa':
           format = 'fasta'
        elif extension == '.fasta':
           format = 'fasta'
        elif extension == '.qual':
           format = 'qual'
        else:
           stop_err( 'Error: File %s had the unknown extension %s' % ( fname, extension ) )
        
        # Assign the dataset name
        name = ''
        match = q.search(basename)
        if match != None:
          lib_num = match.group(1)
          name = 'lib%s-' % lib_num

        match = r.search(basename)
        if match == None:
          stop_err( 'Error: File with basename %s did not have a recognized name' % (basename) )
        
        lib_type = match.group(1)
        if format == 'qual':
          lib_type = 'qual'

        name = name + lib_type        

        # Move the dataset to the proper place
        #db_ref = ''
        #destination = os.path.join( output_dir, 'primary_%s_%s_visible_%s_%s' % (output_id, name, format, db_ref) )
        destination = os.path.join( output_dir, 'primary_%s_%s_visible_%s' % (output_id, name, format) )

        print "moving %s to %s" % (source, destination)

        try:
          os.rename(source, destination)
        except Exception, e:
          stop_err( 'Error: ' + str( e ) )

if __name__ == "__main__": __main__()
