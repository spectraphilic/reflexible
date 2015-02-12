import pytest
import numpy as np

import reflexible as rf
from reflexible.conv2netcdf4 import Header as OldHeader


output_list = ['Fwd1_V9.02', 'Fwd2_V9.02', 'Bwd1_V9.02', 'Bwd2_V9.2beta']


class Dataset:
    def __init__(self, fp_name):
        self.fp_name = fp_name
        self.fp_path = rf.datasets[fp_name]

    def setup(self, tmpdir, nested=False, wetdep=True, drydep=True):
        self.tmpdir = tmpdir   # bring the fixture to the Dataset instance
        self.nc_path = tmpdir.join("%s.nc" % self.fp_name).strpath
        rf.create_ncfile(self.fp_path, nested, wetdep, drydep, outfile=self.nc_path)
        self.oldH = OldHeader(self.fp_path, nested=False)
        self.oldH.fill_backward(nspec=(0,))
        self.wetdep = wetdep
        self.drydep = drydep
        return self.nc_path, self.oldH

    def cleanup(self):
        self.tmpdir.remove(self.nc_path)


class TestHeader:
    @pytest.fixture(autouse=True, params=output_list)
    def setup(self, request, tmpdir):
        dataset = Dataset(request.param)
        self.nc_path, self.oldH = dataset.setup(tmpdir)
        self.H = rf.Header(self.nc_path)
        request.addfinalizer(dataset.cleanup)

    # Test header attributes
    def test_attrs(self):
        req_attrs = (
            # XXX Still small differences (~1 hour) in release times
            # 'releasestart', 'releaseend', 'releasetimes',
            'available_dates', 'outlon0', 'outlat0', 'dxout', 'dyout')
        for attr in req_attrs:
            assert getattr(self.H, attr) == getattr(self.oldH, attr)

    # Test FD structures
    def test_FD(self):
        for date in self.H.available_dates:
            np.testing.assert_array_almost_equal(
                self.H.FD[(0, date)].grid, self.oldH.FD[(0, date)].grid)

    # Test concentrations
    def test_concentrations(self):
        np.testing.assert_array_almost_equal(
            self.H.C[(0, 0)].grid, self.oldH.C[(0, 0)].grid)

    # Test slabs for concentrations
    # XXX The slabs in new Header seems to be transposed wrt to old Header
    def _test_concentration_slabs(self):
        for k in self.H.C[(0, 0)].slabs:
            np.testing.assert_array_almost_equal(
                self.H.C[(0, 0)].slabs[k], self.oldH.C[(0, 0)].slabs[k])
