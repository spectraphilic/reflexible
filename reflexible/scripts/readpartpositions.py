"""
Python Python script that reads raw datafiles of FLEXPART particle positions.

This algorithm strives to use as less memory as possible.  For this, the
output array is created initially with some small size and enlarged as
necessary.

NOTE: This code reads directly from un UNFORMATTED SEQUENTIAL data Fortran
file so care has been taken to skip the record length at the beginning and
the end of every record.  See:
http://stackoverflow.com/questions/8751185/fortran-unformatted-file-format


Author: Francesc Alted
Date: 2015-01-30
"""

from __future__ import print_function

from time import time
import numpy as np


CHUNKSIZE = 50 * 1000


def readpartpositions(filename, nspec, xlon0, ylat0, dx, dy, dataframe=False):
    """Read the particle positions in filename.

    Parameters
    ----------
    filename : string
      the file name of the particle raw files
    nspec : int
      number of species
    xlon0 : float
      the longitude
    xlat0 : float
      the latitude
    dx : float
      the cell size (x)
    dy : float
      the cell size (u)
    dataframe : bool
      Return a pandas DataFrame

    Returns
    -------
    structured_numpy_array OR pandas DataFrame
    """

    xmass_dtype = [('xmass_%d' % (i+1), 'f4') for i in range(nspec)]
    out_fields = [
        ('npoint', 'i4'), ('xtra1', 'f4'), ('ytra1', 'f4'), ('ztra1', 'f4'),
        ('itramem', 'i4'), ('topo', 'f4'), ('pvi', 'f4'), ('qvi', 'f4'),
        ('rhoi', 'f4'), ('hmixi', 'f4'), ('tri', 'f4'), ('tti', 'f4')] + xmass_dtype
    raw_fields = [('begin_recsize', 'i4')] + out_fields + [('end_recsize', 'i4')]
    rectype = np.dtype(out_fields)
    raw_rectype = np.dtype(raw_fields)
    recsize = raw_rectype.itemsize

    out = np.empty(CHUNKSIZE, dtype=rectype)
    chunk = np.empty(CHUNKSIZE, dtype=raw_rectype)
    chunk.flags.writeable = True

    numparticlecount = 0
    with open(filename, "rb", buffering=1) as f:
        # The timein value is at the beginning of the file
        reclen = np.ndarray(shape=1, buffer=f.read(4), dtype="i4")[0]
        assert reclen == 4
        itimein = np.ndarray(shape=1, buffer=f.read(4), dtype="i4")
        reclen = np.ndarray(shape=1, buffer=f.read(4), dtype="i4")[0]
        assert reclen == 4
        i = 0
        while True:
            # Try to read a complete chunk
            data = f.read(CHUNKSIZE * recsize)
            read_records = int(len(data) / recsize)  # the actual number of records read
            chunk = chunk[:read_records]
            chunk.data[:] = data
            # Recompute xtra1 and ytra1 fields for the recently added chunk
            chunk['xtra1'][:] = (chunk['xtra1'] - xlon0) / dx
            chunk['ytra1'][:] = (chunk['ytra1'] - ylat0) / dy
            # Add the chunk to the out array
            out[i:i+read_records] = chunk
            i += read_records
            if read_records < CHUNKSIZE:
                # We reached the end of the file
                break
            else:
                # Enlarge out by CHUNKSIZE more elements
                out = np.resize(out, out.size + CHUNKSIZE)

        # Truncate at the max length (last row is always a sentinel, so remove it)
        out = np.resize(out, i - 1)

    if dataframe:
        import pandas as pd
        return pd.DataFrame(out)
    return out


def main():
    import argparse
    import reflexible.conv2netcdf4 as conv

    parser = argparse.ArgumentParser()
    parser.add_argument("--partdir", "-d",
                        help="Directory where the particle metadata is (including the 'header' and 'dates' files)")
    parser.add_argument("--datafile", "-f",
                        help="The datafile to read (typically inside the partdir)")
    parser.add_argument("--dataframe", "-df",  action='store_true',
                        help="Return a pandas dataframe instead of a structured array (default)")
    parser.add_argument("--timming", "-t",  action='store_true',
                        help="Print timmings for reading the different parts")
    args = parser.parse_args()

    t0 = time()
    H = conv.Header(args.partdir)
    if args.timming:
        print("Time for reading the header: %.3fs" % (time()-t0,))
    print("nspec, outlon0, outlat0, dxout, dyout:", H.nspec, H.outlon0, H.outlat0, H.dxout, H.dyout)

    # Structured arrays
    t0 = time()
    readout = readpartpositions(args.datafile, H.nspec, H.outlon0, H.outlat0, H.dxout, H.dyout, args.dataframe)
    if args.timming:
        print("Time for reading the data: %.3fs" % (time()-t0,))
    print("records read:", len(readout))
    print("3 first records:\n", repr(readout[:3]))
    print("3 last records:\n", repr(readout[-3:]))


if __name__ == "__main__":
    main()
