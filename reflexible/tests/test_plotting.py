#!/usr/bin/env python
from __future__ import absolute_import
from __future__ import print_function

import pytest

import reflexible.plotting as pf
import reflexible as rf

# tuples locating test data, nested(True) and global(False)
test_datasets = [('Fwd1_V10.0', False), ('Fwd1_V10.0', True)]

# TODO: finish code for plotting something related with the datasets above


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
        pass
        # dataset = Dataset(request.param)
        # self.H, self.fp_path, self.nc_path = dataset.setup()
        # request.addfinalizer(dataset.cleanup)

    def test_flexpart_colormap(self):
        from matplotlib.colors import ListedColormap
        assert isinstance(pf._gen_flexpart_colormap(), ListedColormap)
