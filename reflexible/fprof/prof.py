"""
Utilities for profiling code.

Here you can find some small number of context managers that can be
used to annotate the code and get some nice indicators of where the
resources are consumed.

In the future that could grow into its own package.
"""

from __future__ import print_function

import contextlib
import time
import cProfile
import pstats
import StringIO


# Global variable for disabling all the profiles
disabled = True


def disable():
    """Disable the profiling code"""
    global disabled
    disabled = True


def enable():
    """Enable the profiling code"""
    global disabled
    disabled = False


@contextlib.contextmanager
def ctime(explain=''):
    "Counts the time spent in some context"
    if disabled:
        yield
        return
    t = time.time()
    yield
    print("ctime output for: ******* {} *******".format(explain))
    print("%.3f sec" % (time.time() - t))


@contextlib.contextmanager
def cprof(explain='', nlines=20):
    """Does a profile of some context using cProfile

    `explain` allows user to provide some annotations.
    """
    if disabled:
        yield
        return
    pr = cProfile.Profile()
    pr.enable()
    yield
    pr.disable()
    s = StringIO.StringIO()
    #sortby = ('time', 'calls')
    sortby = ('cumulative',)
    ps = pstats.Stats(pr, stream=s).sort_stats(*sortby)
    print("cProf output for: ******* {} *******".format(explain))
    ps.print_stats(nlines)
    print(s.getvalue())


def show_stats(explain, tref, before):
    "Show the used memory (only works for Linux > 2.6)."
    # Build the command to obtain memory info
    import os
    import subprocess
    global vmrss_, vmsize_, vmdata_, vmstk_
    cmd = "cat /proc/%s/status" % os.getpid()
    sout = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout
    for line in sout:
        if line.startswith("VmSize:"):
            vmsize = int(line.split()[1])
        elif line.startswith("VmRSS:"):
            vmrss = int(line.split()[1])
        elif line.startswith("VmData:"):
            vmdata = int(line.split()[1])
        elif line.startswith("VmStk:"):
            vmstk = int(line.split()[1])
    sout.close()
    if before:
        vmrss_ = vmrss
        vmsize_ = vmsize
        vmdata_ = vmdata
        vmstk_ = vmstk
        return time.time()
    print("Memory usage: ******* {}".format(explain))
    print("VmSize: %7s kB\tVmRSS: %7s kB"
          "\tdelta(VmSize): %7s kB\tdelta(VmRSS): %7s kB" % (
        vmsize, vmsize - vmsize_, vmrss, vmrss - vmrss_))
    print("VmData: %7s kB\tVmStk: %7s kB"
          "\tdelta(VmData): %7s kB\tdelta(VmStk): %7s kB" % (
        vmdata, vmdata - vmdata_, vmstk, vmstk - vmstk_))
    tnow = time.time()
    #print("WallClock time:", round(tnow - tref, 3))
    return tnow


@contextlib.contextmanager
def lmprof(explain=''):
    """Counts the memory consumed in some context

    `explain` allows user to provide some annotations.
    """
    if disabled:
        yield
        return
    t = show_stats(explain + " (before)", time.time(), before=True)
    yield
    show_stats(explain, t, before=False)
