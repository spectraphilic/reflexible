import os, os.path

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
        self.H = conv.Header(self.fp_path)

    def cleanup(self):
        self.tmpdir.remove(self.nc_path)


@pytest.fixture(scope="module", params=['Fwd1_V9.02'])
def dataset_fwd(request):
    return Dataset(request.param)


class TestFwdAPI:
    @pytest.fixture(autouse=True)
    def setup(self, request, dataset_fwd, tmpdir):
        self.dataset = dataset = dataset_fwd
        dataset.tmpdir = tmpdir   # bring the fixture to the Dataset instance
        self.H = dataset.H
        self.dataset.nc_path = tmpdir.join("%s.nc" % dataset.fp_name).strpath
        self.fp_path, self.nc_path = dataset.fp_path, dataset.nc_path
        request.addfinalizer(dataset.cleanup)

    def test_nc_create(self):
        rf.create_ncfile(self.fp_path, nested=False, outfile=self.nc_path)
        assert os.path.exists(self.nc_path)

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


@pytest.fixture(scope="module", params=['Bwd1_V9.02', 'Bwd2_V9.2beta'])
def dataset_bwd(request):
    return Dataset(request.param)


class TestBwdAPI:
    @pytest.fixture(autouse=True)
    def setup(self, request, dataset_bwd, tmpdir):
        self.dataset = dataset = dataset_bwd
        dataset.tmpdir = tmpdir   # bring the fixture to the Dataset instance
        self.H = dataset.H
        self.dataset.nc_path = tmpdir.join("%s.nc" % dataset.fp_name).strpath
        self.fp_path, self.nc_path = dataset.fp_path, dataset.nc_path
        request.addfinalizer(dataset.cleanup)

    def test_fill_backward(self):
        self.H.fill_backward(nspec=(0,))
        ckeys = self.H.C[(0, 0)].keys()
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
