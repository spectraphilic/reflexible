import os.path

import pytest

import reflexible as rf
import reflexible.conv2netcdf4 as conv

fd_keys = [
    'dry', 'grid', 'gridfile', 'itime', 'max', 'min', 'rel_i',
    'shape', 'slabs', 'spec_i', 'species', 'timestamp', 'wet']


class Dataset:
    def __init__(self, fp_name):
        self.fp_name = fp_name
        self.fp_path = rf.datasets[fp_name]
        self.fp_pathnames = os.path.join(self.fp_path, "pathnames")
        self.fp_options, self.fp_output = conv.get_fpdirs(self.fp_pathnames)

    def setup(self, tmpdir):
        self.tmpdir = tmpdir  # bring the fixture to the Dataset instance
        self.H = conv.Header(self.fp_output)
        self.nc_path = tmpdir.join("%s.nc" % self.fp_name).strpath
        return self.H, self.fp_pathnames, self.nc_path

    def cleanup(self):
        self.tmpdir.remove(self.nc_path)


class TestFwdAPI:
    @pytest.fixture(autouse=True, params=['Fwd1_V9.02', 'Fwd2_V9.02'])
    def setup(self, request, tmpdir):
        dataset = Dataset(request.param)
        self.H, self.fp_pathnames, self.nc_path = dataset.setup(tmpdir)
        request.addfinalizer(dataset.cleanup)

    def test_nc_create(self):
        rf.create_ncfile(self.fp_pathnames, nested=False, outfile=self.nc_path)
        assert os.path.exists(self.nc_path)

    def test_read_grid(self):
        FD = conv.read_grid(self.H, time_ret=0, nspec_ret=0)
        fdkeys = sorted([k for k in FD.keys() if type(k) is str])
        assert fdkeys == ['grid_dates', 'options']
        fdkeys_ = sorted(FD[(0, '20070121100000')].keys())
        assert fdkeys_ == fd_keys

    def test_H_read_grid(self):
        self.H.read_grid(time_ret=0, nspec_ret=0)
        fdkeys = sorted([k for k in self.H.FD.keys() if type(k) is str])
        assert fdkeys == ['grid_dates', 'options']
        fdkeys_ = sorted(self.H.FD[(0, '20070121100000')].keys())
        assert fdkeys_ == fd_keys

    def test_read_trajectories(self):
        T = conv.read_trajectories(self.H)
        tkeys = sorted(T.keys())
        assert tkeys == [
            'RELEASE_TEST1', 'Trajectories', 'date', 'info', 'labels',
            'version']


class TestBwdAPI:
    @pytest.fixture(autouse=True, params=['Bwd1_V9.02', 'Bwd2_V9.2beta'])
    def setup(self, request, tmpdir):
        dataset = Dataset(request.param)
        self.H, self.fp_path, self.nc_path = dataset.setup(tmpdir)
        request.addfinalizer(dataset.cleanup)

    def test_fill_backward(self):
        self.H.fill_backward(nspec=(0,))
        ckeys = self.H.C[(0, 0)].keys()
        assert sorted(ckeys) == fd_keys
        fdkeys = self.H.FD[(0, self.H.available_dates[0])].keys()
        assert sorted(fdkeys) == fd_keys

    def test_read_grid(self):
        FD = conv.read_grid(self.H, time_ret=0, nspec_ret=0)
        fdkeys = sorted([k for k in FD.keys() if type(k) is str])
        assert fdkeys == ['grid_dates', 'options']
        fdkeys_ = sorted(FD[(0, self.H.available_dates[0])].keys())
        assert fdkeys_ == fd_keys

    def test_H_read_grid(self):
        self.H.read_grid(time_ret=0, nspec_ret=0)
        fdkeys = sorted([k for k in self.H.FD.keys() if type(k) is str])
        assert fdkeys == ['grid_dates', 'options']
        fdkeys_ = sorted(self.H.FD[(0, self.H.available_dates[0])].keys())
        assert fdkeys_ == fd_keys
