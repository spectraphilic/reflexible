import os, os.path

from unittest import TestCase

import reflexible as pf
import reflexible.conv2netcdf4 as conv

fd_keys = [
    'dry', 'grid', 'gridfile', 'itime', 'max', 'min', 'rel_i',
    'shape', 'slabs', 'spec_i', 'species', 'timestamp', 'wet']


class Fwd_API(object):

    def test_nc_create(self):
        pf.create_ncfile(self.dataset, nested=False, outfile="/tmp/Fwd1.nc")
        assert os.path.exists("/tmp/Fwd1.nc")

    def test_read_grid(self):
        FD = conv.read_grid(self.H, time_ret=0, nspec_ret=0)
        fdkeys = sorted(FD.keys())
        assert fdkeys == ['grid_dates', 'options', (0, '20070121100000')]
        fdkeys_ = sorted(FD[(0, '20070121100000')].keys())
        assert fdkeys_ == fd_keys

    def test_H_read_grid(self):
        self.H.read_grid(time_ret=0, nspec_ret=0)
        fdkeys = sorted(self.H.FD.keys())
        assert fdkeys == ['grid_dates', 'options', (0, '20070121100000')]
        fdkeys_ = sorted(self.H.FD[(0, '20070121100000')].keys())
        assert fdkeys_ == fd_keys

    def test_read_trajectories(self):
        T = conv.read_trajectories(self.H)
        tkeys = sorted(T.keys())
        assert tkeys == [
            'RELEASE_TEST1', 'Trajectories', 'date', 'info', 'labels',
            'version']

class Test_Fwd1_V9_02(Fwd_API, TestCase):
    dataset = pf.datasets['Fwd1_V9.02']
    H = conv.Header(dataset)


class Bwd_API(object):

    def test_fill_backward(self):
        self.H.fill_backward(nspec=(0,))
        ckeys = self.H.C[(0,0)].keys()
        assert sorted(ckeys) == fd_keys
        fdkeys = self.H.FD[(0, self.H.available_dates[0])].keys()
        assert sorted(fdkeys) == fd_keys

    def test_read_grid(self):
        FD = conv.read_grid(self.H, time_ret=0, nspec_ret=0)
        fdkeys = sorted(FD.keys())
        assert fdkeys == ['grid_dates', 'options', (0, self.H.available_dates[0])]
        fdkeys_ = sorted(FD[(0, self.H.available_dates[0])].keys())
        assert fdkeys_ == fd_keys

    def test_H_read_grid(self):
        self.H.read_grid(time_ret=0, nspec_ret=0)
        fdkeys = sorted(self.H.FD.keys())
        assert fdkeys == ['grid_dates', 'options', (0, self.H.available_dates[0])]
        fdkeys_ = sorted(self.H.FD[(0, self.H.available_dates[0])].keys())
        assert fdkeys_ == fd_keys

class Test_Bwd_V9_02(Bwd_API, TestCase):
    dataset = pf.datasets['Bwd1_V9.02']
    H = conv.Header(dataset)

class Test_Bwd_V9_2beta(Bwd_API, TestCase):
    dataset = pf.datasets['Bwd2_V9.2beta']
    H = conv.Header(dataset)
