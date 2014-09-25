#!/usr/bin/env python
"""

SYNOPSIS
========

    pflexible [-h] [-v,--verbose] [--version]

DESCRIPTION
===========

    pflexible: A python module for working with FLEXPART Output.


EXAMPLES
========

    #TODO:

    A lot! This is just a starting point. See the doc strings
    for information about the various functions.


AUTHOR
======

    JFB: John F Burkhart <jfburkhart@gmail.com>

CONTRIBUTORS
============

    HSO: Harald Sodemann
    SEC: Sabine Eckhardt
    AST: Andreas Stohl

    Many functions are adaptations of Fortran / NCARG programs (AST)
    or Matlab functions (HSO/SEC).

LICENSE
=======

    This script follows creative commons usage.

VERSION
=======

    ID: $Id$: $Rev$ 
"""
# builtin imports
# import pdb
import sys
import os
import optparse
import time
import traceback


# local imports
import mapping as mp

__version__ = '0.9.6'
__path__ = os.path.abspath(os.curdir)



#### Below for use from command line  ##########################################

def main ():

    global options, args
    here = os.path.curdir
    dataDir = options.pathname
    H = read_header(dataDir, **{'readp':1})  # readheader and release pointspiu
    G = readgridV8(H, **{'unit':'time'})
    D = get_slabs(H, G)

    plt = mp.plot_map(D[0], H)
    plot_name = 'test'
    plt.title(plot_name, fontsize=10)
    plt.savefig(os.path.join(here, plot_name + '.png'))
    print plot_name


if __name__ == '__main__':
    print 'testing'
    try:
        start_time = time.time()
        parser = optparse.OptionParser(
            formatter=optparse.TitledHelpFormatter(),
            usage=globals()['__doc__'],
            version='$Id: py.tpl 327 2008-10-05 23:38:32Z root $')
        parser.add_option ('-v', '--verbose', action='store_true',
                           default=False, help='verbose output')
        parser.add_option ('-T', '--latex', action='store_true',
                           default=False, help='latex usage flag')
        parser.add_option ('-d', '--directory', dest='pathname',
                           default=False, help='pathname of directory containg FLEXPART Output')
        (options, args) = parser.parse_args()
        # if len(args) < 1:
        #    parser.error ('missing argument')
        if options.verbose: print time.asctime()

        #### LateX #####
        if options.latex:
            mpl.rc('font', **{'family':'sans-serif',
                              'sans-serif':['Helvetica']}
                   )

            mpl.rc('text', usetex=True)

        exit_code = main()
        if exit_code is None:
            exit_code = 0
        if options.verbose: print time.asctime()
        if options.verbose: print 'TOTAL TIME IN MINUTES:',
        if options.verbose: print (time.time() - start_time) / 60.0
        sys.exit(exit_code)
    except KeyboardInterrupt, e:  # Ctrl-C
        raise e
    except SystemExit, e:  # sys.exit()
        raise e
    except Exception, e:
        print 'ERROR, UNEXPECTED EXCEPTION'
        print str(e)
        traceback.print_exc()
        os._exit(1)

# vim:set sr et ts=4 sw=4 ft=python fenc=utf-8: // See Vim, :help 'modeline'

