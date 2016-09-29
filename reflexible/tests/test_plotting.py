#!/usr/bin/env python
from __future__ import absolute_import
from __future__ import print_function

import itertools
import pytest

import reflexible as rf
import reflexible.plotting as pf


# tuples locating test data, nested(True) and global(False)
test_datasets = [('Fwd1_V10.0', True), ('Fwd1_V10.0', False)]


# Small test on color maps
def test_flexpart_colormap():
    from matplotlib.colors import ListedColormap
    assert isinstance(pf._gen_flexpart_colormap(), ListedColormap)


class Dataset:
    def __init__(self, fp_dataset):
        self.fp_name = fp_dataset[0]
        self.nested = fp_dataset[1]
        self.fp_path = rf.datasets[self.fp_name]

    def setup(self):
        self.H = rf.Header(self.fp_path, nested=self.nested,
                           absolute_path=False)
        self.nc_path = self.H.ncfile
        return self.H, self.fp_path, self.nc_path

    def cleanup(self):
        pass


class TestPlotting:

    @pytest.fixture(autouse=True, params=test_datasets)
    def setup(self, request):
        dataset = Dataset(request.param)
        self.H, self.fp_path, self.nc_path = dataset.setup()
        request.addfinalizer(dataset.cleanup)

    def test_quickplot(self, tmpdir):
        self.H.fill_grids()
        keys = [(s,k) for s, k in itertools.product(
            range(self.H.nspec), range(self.H.pointspec))]
        for k in keys:
            tc = self.H.C[k].total_column
            dr = [tc.min(), tc.max()]
            if self.H.nested:
                pf.plot_sensitivity(self.H, tc, data_range=dr,
                                    map_region='svalbard')
                # Create a temporary file for the PNG output
                p = tmpdir.join('stads_{0}-{1}_nested.png'.format(*k))
                png_file = str(p.realpath())
                pf.plt.savefig(png_file)
                assert p.size() > 0   # check that PNG file has some content
            else:
                pf.plot_sensitivity(self.H, tc, data_range=dr,
                                    map_region='north_atlantic')
                # Create a temporary file for the PNG output
                p = tmpdir.join('stads_{0}-{1}.png'.format(*k))
                png_file = str(p.realpath())
                pf.plt.savefig(png_file)
                assert p.size() > 0  # check that PNG file has some content
