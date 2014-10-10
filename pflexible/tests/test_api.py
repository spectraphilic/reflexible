from unittest import TestCase

import pflexible as pf


class Test_API(TestCase):

    def test_(self):
        H = pf.Header(pf.tests.Fwd1_data)
        assert sorted(H.keys()) == sorted([
            'zpoint2', 'zpoint1', 'maxpoint', 'compoint', 'nxmax',
            'maxageclass', 'jjjjmmdd', 'releasestart', 'loutsample',
            'latitude', 'nspec', 'numzgrid', 'nz_list', 'dxout', 'hhmmss',
            'maxspec', '_3', '_2', '_1', '_0', 'output_unit', 'Area',
            'available_dates_dt', 'lconvection', 'ibtime', 'plot_unit',
            'ind_receptor', 'xpoint', 'fp_version', 'ireleasestart',
            'version', 'dyout', 'numpointspec', 'nymax', 'mpart', 'oro',
            'path', 'releaseend', 'numxgrid', 'ireleaseend', 'drydep',
            'last_date', 'ibdate', 'pathname', 'loutstep', 'ind_source',
            'numpoint', 'options', 'outheight', 'nested', 'Heightnn', 'unit',
            'area', 'simulationstart', 'first_date', 'lage', 'nageclass',
            'direction', 'decayconstant', 'numygrid', 'longitude', 'xp2',
            'xp1', 'nzmax', 'outlon0', 'wetdep', 'lsubgrid', 'species',
            'npart', 'nx', 'ny', 'nz', 'ypoint', 'outlat0', 'junk',
            'releasetimes', 'ageclasses', 'layerthickness', 'loutaver',
            'kindz', 'flexpart', 'numageclasses', 'xmass', 'alt_unit',
            'yp2', 'yp1', 'available_dates'])
