"""
read_part_positions() reads raw datafiles of FLEXPART particle positions.

This function strives to use as less memory as possible; for this, a
bcolz ctable container is used for holding the data.  Besides to be compressed
in-memory, its chunked nature makes a natural fit for data that needs to
be appended because it does not need expensive memory resize operations.

NOTE: This code reads directly from un UNFORMATTED SEQUENTIAL data Fortran
file so care has been taken to skip the record length at the beginning and
the end of every record.  See:
http://stackoverflow.com/questions/8751185/fortran-unformatted-file-format

"""

import sys
from time import time
import numpy as np
import bcolz

# The amount of records that are read from disk in a single shot.  Keep this as
# small as possible.
CHUNKSIZE = 10 * 1000


def get_quantized_ctable(dtype, cparams, quantize=None, expectedlen=None):
    """Return a ctable with the quantize filter enabled for floating point cols."""
    columns, names = [], []
    for fname, ftype in dtype.descr:
        names.append(fname)
        if 'f' in ftype:
            cparams2 = bcolz.cparams(clevel=cparams.clevel, cname=cparams.cname, quantize=quantize)
            columns.append(bcolz.zeros(0, dtype=ftype, cparams=cparams2, expectedlen=expectedlen))
        else:
            columns.append(bcolz.zeros(0, dtype=ftype, cparams=cparams, expectedlen=expectedlen))
    return bcolz.ctable(columns=columns, names=names)


def read_part_positions(filename, nspec, ctable=True, clevel=5, cname="lz4", quantize=None):
    """Read the particle positions in `filename`.

    Parameters
    ----------
    filename : string
        The file name of the particle raw data
    nspec : int
        number of species in particle raw data
    ctable : bool
        Return a bcolz ctable container.  If not, a numpy structured array is returned instead.
    clevel : int
        Compression level for the ctable container
    cname : string
        Codec name for the ctable container.  Can be 'blosclz', 'lz4', 'zlib' or 'zstd'.
    quantize : int
        Quantize data to improve (lossy) compression.  Data is quantized using
        np.around(scale*data)/scale, where scale is 2**bits, and bits is
        determined from the quantize value.  For example, if quantize=1, bits
        will be 4.  0 means that the quantization is disabled.

    Returns
    -------
    ctable object OR structured_numpy_array

    Returning a ctable is preferred because it is used internally so it does not require to be
    converted to other formats, so it is faster and uses less memory.

    Note: Passing a `quantize` param > 0 can increase the compression ratio of the ctable
    container, but it may also slow down the reading speed significantly.
    """

    xmass_dtype = [('xmass_%d' % (i+1), 'f4') for i in range(nspec)]
    # note age is calculated from itramem by adding itimein
    out_fields = [
        ('npoint', 'i4'), ('xtra1', 'f4'), ('ytra1', 'f4'), ('ztra1', 'f4'),
        ('itramem', 'i4'), ('topo', 'f4'), ('pvi', 'f4'), ('qvi', 'f4'),
        ('rhoi', 'f4'), ('hmixi', 'f4'), ('tri', 'f4'), ('tti', 'f4')] + xmass_dtype
    raw_fields = [('begin_recsize', 'i4')] + out_fields + [('end_recsize', 'i4')]
    raw_rectype = np.dtype(raw_fields)
    recsize = raw_rectype.itemsize

    cparams = bcolz.cparams(clevel=clevel, cname=cname)
    if quantize is not None and quantize > 0:
        out = get_quantized_ctable(raw_rectype, cparams=cparams, quantize=quantize, expectedlen=int(1e6))
    else:
        out = bcolz.zeros(0, dtype=raw_rectype, cparams=cparams, expectedlen=int(1e6))

    with open(filename, "rb", buffering=1) as f:
        # The timein value is at the beginning of the file
        reclen = np.ndarray(shape=(1,), buffer=f.read(4), dtype="i4")[0]
        assert reclen == 4
        itimein = np.ndarray(shape=(1,), buffer=f.read(4), dtype="i4")
        reclen = np.ndarray(shape=(1,), buffer=f.read(4), dtype="i4")[0]
        assert reclen == 4
        nrec = 0
        while True:
            # Try to read a complete chunk
            data = f.read(CHUNKSIZE * recsize)
            read_records = int(len(data) / recsize)  # the actual number of records read
            chunk = np.ndarray(shape=(read_records,), buffer=data, dtype=raw_rectype)
            # Add the chunk to the out array
            out.append(chunk[:read_records])
            nrec += read_records
            if read_records < CHUNKSIZE:
                # We reached the end of the file
                break

    # Truncate at the max length (last row is always a sentinel, so remove it)
    out.trim(1)
    # Remove the first and last columns
    out.delcol("begin_recsize")
    out.delcol("end_recsize")

    if ctable:
        return out
    else:
        return out[:]


if __name__ == "__main__":
    # ctable
    t0 = time()
    readout = read_part_positions(sys.argv[1], 2)
    print("Time for reading a bcolz array: %.3fs" % (time()-t0,))
    print("records read:", len(readout))
    print("3 first records:\n", repr(readout[:3]))
    # print(" last records:\n", repr(readout[-3:]))
    print("Compression ratio: %.2f" % (readout.nbytes / readout.cbytes))

    # ctable and quantize > 0
    t0 = time()
    readout = read_part_positions(sys.argv[1], 2, clevel=1, cname="zstd", quantize=2)
    print("Time for reading a bcolz array: %.3fs" % (time()-t0,))
    print("records read:", len(readout))
    print("3 first records:\n", repr(readout[:3]))
    # print(" last records:\n", repr(readout[-3:]))
    print("Compression ratio: %.2f" % (readout.nbytes / readout.cbytes))

    # Structured arrays
    t0 = time()
    readout = read_part_positions(sys.argv[1], 2, ctable=False)
    print("Time for reading a structured array: %.3fs" % (time()-t0,))
    print("records read:", len(readout))
    # print("3 first records:\n", repr(readout[:3]))
    # print(" last records:\n", repr(readout[-3:]))
