from unittest import TestCase

import pflexible as pf


class Test_API(TestCase):

    H = pf.Header(pf.Fwd1_data)
    fdkeys = [
            'dry', 'grid', 'gridfile', 'itime', 'max', 'min', 'rel_i',
            'shape', 'slabs', 'spec_i', 'species', 'timestamp', 'wet']

    def test_fill_backward(self):
        self.H.fill_backward(nspec=(0,))
        ckeys = self.H.C[(0,2)].keys()
        assert sorted(ckeys) == self.fdkeys
        fdkeys = self.H.FD[(0, '20070121220000')].keys()
        assert sorted(fdkeys) == self.fdkeys

    def test_read_grid(self):
        FD = pf.read_grid(self.H, time_ret=0, nspec_ret=0)
        fdkeys = sorted(FD.keys())
        assert fdkeys == ['grid_dates', 'options', (0, '20070121100000')]
        fdkeys_ = sorted(FD[(0, '20070121100000')].keys())
        assert fdkeys_ == [
            'dry', 'grid', 'gridfile', 'itime', 'max', 'min', 'rel_i',
            'shape', 'spec_i', 'species', 'timestamp', 'wet']

    def test_read_trajectories(self):
        T = pf.read_trajectories(self.H)
        tkeys = sorted(T.keys())
        assert tkeys == [
            'RELEASE_TEST1', 'Trajectories', 'date', 'info', 'labels',
            'version']
