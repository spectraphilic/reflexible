#!/usr/bin/env python
from __future__ import absolute_import
from __future__ import print_function

import os
import itertools
import pytest
import matplotlib.pyplot as plt

import reflexible as rf


# tuples locating test data, nested(True) and global(False)
test_datasets = [('Fwd1_V9.02', True), ('Fwd1_V9.02', False)]
# TODO: support the next use cases too
#test_datasets = [('Fwd1_V10.1', True), ('Fwd1_V10.1', False)]


# Small test on color maps
def test_flexpart_colormap():
    from matplotlib.colors import ListedColormap
    assert isinstance(rf.plotting._gen_flexpart_colormap(), ListedColormap)


class Dataset:
    def __init__(self, fp_dataset):
        self.fp_name = fp_dataset[0]
        self.nested = fp_dataset[1]
        self.fp_path = rf.datasets[self.fp_name]
        self.fp_pathnames = os.path.join(self.fp_path, "pathnames")

    def setup(self):
        self.fprun = rf.Flexpart(self.fp_pathnames, nested=self.nested)
        return self.fprun.Header

    def cleanup(self):
        pass


class TestPlotting:

    @pytest.fixture(autouse=True, params=test_datasets)
    def setup(self, request):
        dataset = Dataset(request.param)
        self.H = dataset.setup()
        request.addfinalizer(dataset.cleanup)

    def test_quickplot(self, tmpdir):
        self.H.fill_grids()
        keys = [(s,k) for s, k in itertools.product(
            range(self.H.nspec), range(self.H.pointspec))]
        for k in keys:
            tc = self.H.C[k].total_column
            dr = [tc.min(), tc.max()]
            if self.H.nested:
                rf.plot_sensitivity(self.H, tc, data_range=dr,
                                    map_region='svalbard')
                # Create a temporary file for the PNG output
                p = tmpdir.join('stads_{0}-{1}_nested.png'.format(*k))
                png_file = str(p.realpath())
                plt.savefig(png_file)
                assert p.size() > 0   # check that PNG file has some content
            else:
                rf.plot_sensitivity(self.H, tc, data_range=dr,
                                    map_region='north_atlantic')
                # Create a temporary file for the PNG output
                p = tmpdir.join('stads_{0}-{1}.png'.format(*k))
                png_file = str(p.realpath())
                plt.savefig(png_file)
                assert p.size() > 0  # check that PNG file has some content
