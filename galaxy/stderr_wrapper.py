#!/usr/bin/env python

"""
Wrapper that execute a program and its arguments but reports standard error
messages only if the program exit status was not 0
Example: ./stderr_wrapper.py myprog arg1 -f arg2
Author: Florent Angly
"""

import sys, subprocess

assert sys.version_info[:2] >= ( 2, 4 )

def stop_err( msg ):
    sys.stderr.write( "%s\n" % msg )
    sys.exit()

def __main__():
    # Get command-line arguments
    args = sys.argv
    # Remove name of calling program, i.e. ./stderr_wrapper.py
    args.pop(0)
    # If there are no arguments left, we're done
    if len(args) == 0:
        return
   
    # If one needs to silence stdout 
    #args.append( ">" )
    #args.append( "/dev/null" )

    cmdline = " ".join(args)
    print cmdline
    try:
        # Run program
        proc = subprocess.Popen( args=cmdline, shell=True, stderr=subprocess.PIPE )
        returncode = proc.wait()
        # Capture stderr, allowing for case where it's very large
        stderr = ''
        buffsize = 1048576
        try:
            while True:
                stderr += proc.stderr.read( buffsize )
                if not stderr or len( stderr ) % buffsize != 0:
                    break
        except OverflowError:
            pass
        # Running Grinder failed: write error message to stderr
        if returncode != 0:
            raise Exception, stderr
    except Exception, e: 
        # Running Grinder failed: write error message to stderr
        stop_err( 'Error:\n' + str( e ) )


if __name__ == "__main__": __main__()
