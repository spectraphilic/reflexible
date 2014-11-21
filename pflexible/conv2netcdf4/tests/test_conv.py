import os, os.path

from unittest import TestCase

import pflexible as pf
import pflexible.conv2netcdf4 as conv


class Test_Fwd(TestCase):

    H = conv.Header(pf.Fwd1_data)
    fdkeys = [
            'dry', 'grid', 'gridfile', 'itime', 'max', 'min', 'rel_i',
            'shape', 'slabs', 'spec_i', 'species', 'timestamp', 'wet']

    def test_nc_create(self):
        pf.create_ncfile(pf.Fwd1_data, nested=False, outfile="/tmp/Fwd1.nc")
        assert os.path.exists("/tmp/Fwd1.nc")

    def test_read_grid(self):
        FD = conv.read_grid(self.H, time_ret=0, nspec_ret=0)
        fdkeys = sorted(FD.keys())
        assert fdkeys == ['grid_dates', 'options', (0, '20070121100000')]
        fdkeys_ = sorted(FD[(0, '20070121100000')].keys())
        assert fdkeys_ == self.fdkeys

    def test_H_read_grid(self):
        self.H.read_grid(time_ret=0, nspec_ret=0)
        fdkeys = sorted(self.H.FD.keys())
        assert fdkeys == ['grid_dates', 'options', (0, '20070121100000')]
        fdkeys_ = sorted(self.H.FD[(0, '20070121100000')].keys())
        assert fdkeys_ == self.fdkeys

    def test_read_trajectories(self):
        T = conv.read_trajectories(self.H)
        tkeys = sorted(T.keys())
        assert tkeys == [
            'RELEASE_TEST1', 'Trajectories', 'date', 'info', 'labels',
            'version']

class Test_Bwd(TestCase):

    H = conv.Header(pf.Bwd1_data)
    fdkeys = [
            'dry', 'grid', 'gridfile', 'itime', 'max', 'min', 'rel_i',
            'shape', 'slabs', 'spec_i', 'species', 'timestamp', 'wet']

    def test_fill_backward(self):
        self.H.fill_backward(nspec=(0,))
        ckeys = self.H.C[(0,0)].keys()
        assert sorted(ckeys) == self.fdkeys
        fdkeys = self.H.FD[(0, self.H.available_dates[0])].keys()
        assert sorted(fdkeys) == self.fdkeys

    def test_read_grid(self):
        FD = conv.read_grid(self.H, time_ret=0, nspec_ret=0)
        fdkeys = sorted(FD.keys())
        assert fdkeys == ['grid_dates', 'options', (0, self.H.available_dates[0])]
        fdkeys_ = sorted(FD[(0, self.H.available_dates[0])].keys())
        assert fdkeys_ == self.fdkeys

    def test_H_read_grid(self):
        self.H.read_grid(time_ret=0, nspec_ret=0)
        fdkeys = sorted(self.H.FD.keys())
        assert fdkeys == ['grid_dates', 'options', (0, self.H.available_dates[0])]
        fdkeys_ = sorted(self.H.FD[(0, self.H.available_dates[0])].keys())
        assert fdkeys_ == self.fdkeys