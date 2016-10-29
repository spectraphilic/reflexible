"""Tests to the Flexpart interface."""

import os
import pytest

import reflexible as rf

# Only tests for FP v10 here because this is what we will be using the most
output_list = ['Fwd1_V10.1']


class Dataset:
    def __init__(self, fp_name):
        self.fp_name = fp_name
        self.fp_path = rf.datasets[fp_name]
        self.fp_pathnames = os.path.join(self.fp_path, "pathnames")


class TestHeader:
    @pytest.fixture(autouse=True, params=output_list)
    def setup(self, request):
        dataset = Dataset(request.param)
        self.fprun = rf.Flexpart(dataset.fp_pathnames)

    def test_releases(self):
        releases = self.fprun.Releases[:]
        assert len(releases) == 1
        assert tuple(releases[0]) == (
            20050501, 0, 20050502, 0, 10.5, 11.0, 59.75, 60.25, 0.0,
            100.0, 1, 1.0, 10000, b'OsloRelease')

    def test_command(self):
        command = self.fprun.Command
        assert command == {
            'IND_RECEPTOR': 1, 'IETIME': 0, 'LSYNCTIME': 900, 'LDIRECT': 1,
            'IEDATE': 20050510, 'SURF_ONLY': 0, 'LOUTAVER': 86400, 'IPOUT': 0,
            'CBLFLAG': 0, 'IFLUX': 0, 'LCONVECTION': 1, 'NESTED_OUTPUT': 0,
            'LOUTSAMPLE': 3600, 'LSUBGRID': 1, 'CTL': -5, 'LOUTSTEP': 86400,
            'LAGESPECTRA': 1, 'LINIT_COND': 1, 'IOUT': 1, 'IBDATE': 20050501,
            'IOUTPUTFOREACHRELEASE': 0, 'ITSPLIT': 9999999, 'MQUASILAG': 0,
            'MDOMAINFILL': 0, 'IBTIME': 0, 'IND_SOURCE': 1, 'IFINE': 4,
            'IPIN': 0}

    def test_species(self):
        species = self.fprun.Species
        assert species['dryvel'] == [-9.99]
