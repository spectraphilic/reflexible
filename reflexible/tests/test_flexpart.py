"""Tests to the Flexpart interface."""

import os
import pytest

import reflexible as rf

output_list = ['Fwd1_V9.02', 'Fwd2_V9.02', 'Bwd1_V9.02', 'Bwd2_V9.2beta',
               'Fwd1_V10.1']

class Dataset:
    def __init__(self, fp_name):
        self.fp_name = fp_name
        self.fp_path = rf.datasets[fp_name]
        self.fp_pathnames = os.path.join(self.fp_path, "pathnames")


class TestHeader:
    @pytest.fixture(autouse=True, params=output_list)
    def setup(self, request):
        self.dataset = dataset = Dataset(request.param)
        self.fprun = rf.Flexpart(dataset.fp_pathnames)

    def test_releases(self):
        # Only test a few FP runs here
        if self.dataset.fp_name == 'Fwd1_V10.1':
            releases = self.fprun.Releases[:]
            assert len(releases) == 1
            assert tuple(releases[0]) == (
                20050501, 0, 20050502, 0, 10.5, 11.0, 59.75, 60.25, 0.0,
                100.0, 1, 1.0, 10000, b'OsloRelease')
        elif self.dataset.fp_name == 'Fwd1_V9.02':
            releases = self.fprun.Releases
            assert len(releases) == 1
            assert releases['release_point_names'][0] == b'RELEASE_TEST1'

    def test_command(self):
        command = self.fprun.Command
        # Only test a few FP runs here
        if self.dataset.fp_name == 'Fwd1_V10.1':
            assert command == {
                'IND_RECEPTOR': 1, 'IETIME': 0, 'LSYNCTIME': 900, 'LDIRECT': 1,
                'IEDATE': 20050510, 'SURF_ONLY': 0, 'LOUTAVER': 86400,
                'IPOUT': 0, 'CBLFLAG': 0, 'IFLUX': 0, 'LCONVECTION': 1,
                'NESTED_OUTPUT': 0, 'LOUTSAMPLE': 3600, 'LSUBGRID': 1,
                'CTL': -5, 'LOUTSTEP': 86400, 'LAGESPECTRA': 1,
                'LINIT_COND': 1, 'IOUT': 1, 'IBDATE': 20050501,
                'IOUTPUTFOREACHRELEASE': 0, 'ITSPLIT': 9999999, 'MQUASILAG': 0,
                'MDOMAINFILL': 0, 'IBTIME': 0, 'IND_SOURCE': 1, 'IFINE': 4,
                'IPIN': 0}
        elif self.dataset.fp_name == 'Fwd1_V9.02':
            assert command == {
                'SIM_DIR': 1, 'CTL': 3.0, 'IPIN': 0, 'IFLUX': 0,
                'NESTED_OUTPUT': 1, 'MDOMAINFILL': 0, 'IND_RECEPTOR': 1,
                'LAGESPECTRA': 0, 'SIM_START': ['20070121', '090000'],
                'IOUT': 5, 'AVG_CNC_INT': 3600, 'IPOUT': 0, 'MQUASILAG': 0,
                'OUTPUTFOREACHRELEASE': 0, 'IFINE': 4, 'SYNC': 300,
                'LSUBGRID': 1, 'SIM_END': ['20070122', '180000'],
                'LCONVECTION': 0, 'AVG_CNC_TAVG': 3600, 'LINIT_COND': 2,
                'T_PARTSPLIT': 999999999, 'IND_SOURCE': 1,
                'CNC_SAMP_TIME': 300}

    def test_species(self):
        species = self.fprun.Species
        assert species['dryvel'] == [-9.99]
