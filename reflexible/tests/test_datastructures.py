import os
import pytest
import numpy as np

import reflexible as rf

fd_keys = [
    'dry', 'grid', 'gridfile', 'itime', 'max', 'min', 'rel_i',
    'shape', 'slabs', 'spec_i', 'species', 'timestamp', 'wet']

# tuples locating test data, nested(True) and global(False)
test_datasets = [('Fwd1_V10.1', False),
                 ('Fwd1_V9.02', False), ('Fwd1_V9.02', True)]


def monotonically_increasing(l):
    return all(x < y for x, y in zip(l, l[1:]))


class Dataset:
    def __init__(self, fp_dataset):
        self.fp_name = fp_dataset[0]
        self.nested = fp_dataset[1]
        self.fp_path = rf.datasets[self.fp_name]
        self.fp_pathnames = os.path.join(self.fp_path, "pathnames")

    def setup(self):
        self.fprun = rf.Flexpart(self.fp_pathnames, nested=self.nested)
        return self.fp_name, self.fprun

    def cleanup(self):
        pass


class TestHeader:
    @pytest.fixture(autouse=True, params=test_datasets)
    def setup(self, request):
        dataset = Dataset(request.param)
        self.fp_name, self.fprun = dataset.setup()
        self.H = self.fprun.Header
        request.addfinalizer(dataset.cleanup)

    def test_ncattrs(self):
        header_attrs = ('outlon0', 'outlat0', 'dxout',
            'dyout', 'ibdate', 'ibtime', 'iedate', 'ietime',
            'loutstep', 'loutaver', 'loutsample',
            'lsubgrid', 'lconvection', 'ind_source',
            'ind_receptor', 'ldirect', 'iout')
        for attr in header_attrs:
            assert getattr(self.H, attr) == getattr(self.H.nc, attr)

    def test_direction(self):
        assert self.H.direction in ['forward', 'backward']

    def test_nspec(self):
        assert self.H.nspec == self.H.nc.dims['numspec']

    def test_species(self):
        """This test needs better consideration."""
        assert self.H.species[0] in ['I-131', 'TRACER', 'AEROTRACER']

    def test_dimensions(self):
        header_dims = ('numpoint', 'pointspec', 'nageclass')
        for attr in header_dims:
            assert getattr(self.H, attr) == self.H.nc.dims[attr]

    def test_grid_dimensions(self):
        grid_dims = (('numxgrid', 'longitude'),
            ('numygrid', 'latitude'), ('numzgrid', 'height'))
        for hdim, ncdim in grid_dims:
            assert getattr(self.H, hdim) == self.H.nc.dims[ncdim]

    def test_lon(self):
        assert monotonically_increasing(self.H.longitude)

    def test_lat(self):
        assert monotonically_increasing(self.H.latitude)

    def test_available_dates(self):
        assert isinstance(self.H.available_dates[0], str)

    def test_releasetimes(self):
        assert hasattr(self.H, 'releasetimes')

    def test_variables(self):
        ncvars = ('ORO', 'zpoint1', 'zpoint2')
        for ncvar in ncvars:
            assert hasattr(self.H, ncvar)

    def test_options(self):
        assert self.H.nested == self.H.options['nested']

    def test_area(self):
        H = self.H
        assert isinstance(H.area, np.ndarray)
        assert H.area.shape == H.ORO.shape


class TestTrajectory:
    @pytest.fixture(autouse=True, params=test_datasets)
    def setup(self, request):
        dataset = Dataset(request.param)
        self.fp_name, self.fprun = dataset.setup()
        self.H = self.fprun.Header
        request.addfinalizer(dataset.cleanup)

    def test_read_trajectories(self):
        if self.fp_name == 'Fwd1_V10.1':
            # Fwd1_V10.1 does not have a trajectories file
            return
        T = rf.read_trajectories(self.H)
        assert isinstance(T.Trajectories, np.ndarray)


class TestCommand:
    @pytest.fixture(autouse=True, params=test_datasets)
    def setup(self, request):
        self.C = rf.Command()
        #dataset = Dataset(request.param)
        #self.H, self.fp_path, self.nc_path = dataset.setup()
        #request.addfinalizer(dataset.cleanup)

    def test_attributes(self):
        for k in self.C._OPTIONS.keys():
            assert hasattr(self.C, k.lower())

    def test_options(self):
        C = rf.Command(**{'LINIT_COND':1})
        assert C.linit_cond == 1
    # def test_read_grid(self):
    #     FD = conv.read_grid(self.H, time_ret=0, nspec_ret=0)
    #     fdkeys = sorted(FD.keys())
    #     assert fdkeys == ['grid_dates', 'options', (0, '20070121100000')]
    #     fdkeys_ = sorted(FD[(0, '20070121100000')].keys())
    #     assert fdkeys_ == fd_keys

    # def test_H_read_grid(self):
    #     self.H.read_grid(time_ret=0, nspec_ret=0)
    #     fdkeys = sorted(self.H.FD.keys())
    #     assert fdkeys == ['grid_dates', 'options', (0, '20070121100000')]
    #     fdkeys_ = sorted(self.H.FD[(0, '20070121100000')].keys())
    #     assert fdkeys_ == fd_keys

    # def test_read_trajectories(self):
    #     T = conv.read_trajectories(self.H)
    #     tkeys = sorted(T.keys())
    #     assert tkeys == [
    #         'RELEASE_TEST1', 'Trajectories', 'date', 'info', 'labels',
    #         'version']


# class TestBwdAPI:
#     @pytest.fixture(autouse=True, params=['Bwd1_V9.02', 'Bwd2_V9.2beta'])
#     def setup(self, request, tmpdir):
#         dataset = Dataset(request.param)
#         self.H, self.fp_path, self.nc_path = dataset.setup(tmpdir)
#         request.addfinalizer(dataset.cleanup)

#     def test_fill_backward(self):
#         self.H.fill_backward(nspec=(0,))
#         ckeys = self.H.C[(0, 0)].keys()
#         assert sorted(ckeys) == fd_keys
#         fdkeys = self.H.FD[(0, self.H.available_dates[0])].keys()
#         assert sorted(fdkeys) == fd_keys

#     def test_read_grid(self):
#         FD = conv.read_grid(self.H, time_ret=0, nspec_ret=0)
#         fdkeys = sorted(FD.keys())
#         assert fdkeys == ['grid_dates', 'options', (0, self.H.available_dates[0])]
#         fdkeys_ = sorted(FD[(0, self.H.available_dates[0])].keys())
#         assert fdkeys_ == fd_keys

#     def test_H_read_grid(self):
#         self.H.read_grid(time_ret=0, nspec_ret=0)
#         fdkeys = sorted(self.H.FD.keys())
#         assert fdkeys == ['grid_dates', 'options', (0, self.H.available_dates[0])]
#         fdkeys_ = sorted(self.H.FD[(0, self.H.available_dates[0])].keys())
#         assert fdkeys_ == fd_keys
