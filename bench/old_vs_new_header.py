"""
Small benchmark that compares the resource consumption (CPU and memory)
of using the old (coming from pflexible) and new Header data structure.

For running this, just do:

$ python bench/old_vs_new_header.py

without parameters.

This will automatically download a test data file and will run the
benchamrks on it.

Note: this only runs on Linux.

:Author: Francesc Alted

"""

import sys
import os.path
import subprocess

import reflexible as rf
import reflexible.conv2netcdf4 as conv
from reflexible import fprof

fprof.enable()

wdir = os.path.dirname(__file__)


def check_test_data(dirname):
    datadir = os.path.join(wdir, dirname)
    if not os.path.isdir(datadir):
        # Download the test data from John's site
        tarball = "flexpart_V8data.tgz"
        ptarball = os.path.join(wdir, tarball)
        if not os.path.isfile(ptarball):
            retcode = subprocess.call(
                "wget -O %s http://folk.uio.no/johnbur/sharing/%s" % (ptarball, tarball),
                shell=True)
            if retcode < 0:
                print >>sys.stderr, "wget was terminated by signal", -retcode
        retcode = subprocess.call("tar xvfz %s -C %s" % (ptarball, wdir), shell=True)
        if retcode < 0:
            print >>sys.stderr, "tar was terminated by signal", -retcode
    return datadir


test_data = check_test_data("test_data")


ncfile = os.path.join(wdir, "test_data.nc")
if not os.path.isfile(ncfile):
    rf.create_ncfile(test_data, nested=False, wetdep=True, drydep=True, outfile=ncfile)

with fprof.lmprof(), fprof.ctime("Read the original Header"):
    H = conv.Header(test_data)

with fprof.lmprof(), fprof.ctime("Read the new Header"):
    Hnc = rf.Header(ncfile)

with fprof.lmprof(), fprof.ctime("Read all the concentrations with old Header"):
    H.fill_backward(nspec=range(H.nspec))
    for nspec, pointspec in zip(range(H.nspec), range(H.numpointspec)):
        print "C[(%d, %d)]:" % (nspec, pointspec), H.C[(nspec, pointspec)]

with fprof.lmprof(), fprof.ctime("Read all the concentrations with new Header"):
    for nspec, pointspec in zip(range(Hnc.nspec), range(Hnc.numpointspec)):
        print "C[(%d, %d)]:" % (nspec, pointspec), Hnc.C[(nspec, pointspec)]
