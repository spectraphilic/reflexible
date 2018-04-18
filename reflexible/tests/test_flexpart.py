"""Tests to the Flexpart interface."""

import os
import pytest

import reflexible as rf

output_list = ['Fwd1_V9.02', 'Fwd2_V9.02', 'Bwd1_V9.02', 'Bwd2_V9.2beta',
               'Fwd1_V10.1', 'HelloWorld_V9.02', 'Only_Outputs_V9.02']

class Dataset:
    def __init__(self, fp_name):
        self.fp_name = fp_name
        self.fp_path = rf.datasets[fp_name]
        if fp_name != 'Only_Outputs_V9.02':
            self.fp_pathnames = os.path.join(self.fp_path, "pathnames")
        else:
            self.fp_pathnames = self.fp_path


class TestHeader:
    @pytest.fixture(autouse=True, params=output_list)
    def setup(self, request):
        self.dataset = dataset = Dataset(request.param)
        self.fprun = rf.Flexpart(dataset.fp_pathnames)

    def test_Header(self):
        # Only test a few FP runs here
        if self.dataset.fp_name in ('Fwd1_V10.1', 'Fwd1_V9.02'):
            H = self.fprun.Header
            assert len(H.C.keys()) >= 1
            assert len(H.FD.keys()) >= 1
            assert len(H.C[0,0].data_cube.shape) == 5

    def test_releases(self):
        # Only test a few FP runs here
        if self.dataset.fp_name == 'Fwd1_V10.1':
            releases = self.fprun.Releases['releases']
            assert len(releases) == 1
            assert tuple(releases[0]) == (
                20050501, 0, 20050502, 0, 10.5, 11.0, 59.75, 60.25, 0.0,
                100.0, 1, 1.0, 10000, b'OsloRelease')
        elif self.dataset.fp_name == 'Fwd1_V9.02':
            releases = self.fprun.Releases['release_point_names']
            assert len(releases) == 1
            assert releases[0] == b'RELEASE_TEST1'

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
        if self.dataset.fp_name == 'Fwd1_V10.1':
            # species id = 15
            assert species == {
                'decay': [691200.0],
                # TODO: Search string in reflexible.scripts.create_ncfile.read_species
                # needs to be updated to handle below/in-cloud scavenging in V10
                #'weta': [1.0E-04],
                #'wetb': [0.80],
                'reldiff': [-9.9],
                'henry': [-9.9],
                'f0': [-9.9],
                'dquer': [6.0E-7],
                'dsigma': [3.0E-1],
                'dryvel': [-9.99],
                'weightmolar': [-9.99],
                'ohreact': [-9.9E-09],
                'spec_ass': [-9],
                'kao': [-99.99]
            }
        elif self.dataset.fp_name == 'Fwd1_V9.02':
            # species id = 1
            assert species == {
                'decay': [-999.9],
                'weta': [-9.9E-09],
                'wetb': [-9.9E-09],
                'reldiff': [-9.9],
                'henry': [-9.9],
                'f0': [-9.9],
                'dquer': [-9.9],
                'dsigma': [-9.9],
                'dryvel': [-9.99],
                'weightmolar': [350.00],
                'ohreact': [-9.9E-09],
                'spec_ass': [-9],
                'kao': [-99.99]
            }
        elif self.dataset.fp_name == 'Fwd2_V9.02':
            # species id = 12
            assert species == {
                'decay': [-999.9],
                'weta': [5.0E-06],
                'wetb': [0.62],
                'reldiff': [-9.9],
                'henry': [1.0E-09],
                'f0': [1.0E-09],
                'dquer': [4.0E-7],
                'dsigma': [3.0E-1],
                'dryvel': [-9.99],
                'weightmolar': [-9.99],
                'ohreact': [-9.9E-09],
                'spec_ass': [-9],
                'kao': [-99.99]
            }
