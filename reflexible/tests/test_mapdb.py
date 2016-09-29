import pytest
import os
from ..mapping import map_regions

region_list = ['default', 'polarcat']

user_mapdb = """
northern_hemisphere:
    alias: default
    descr: Northern Hemisphere (default)
    map_par:
        projection: cyl
        llcrnrlat: 0.5
        urcrnrlat: 90
        llcrnrlon: -180.7
        urcrnrlon: 180
        resolution: c
        # anchor: W
    fig_par:
        figsize: [8, 3]  # w,h tuple
        axlocs: [0.1, 0.1, .7, .8]

polarcat:
    descr: Whole North Atlantic and Norwegian, Greenland Seas
    map_par:
        llcrnrlat: 90.
        llcrnrlon: 95.
        urcrnrlat: 81.
        urcrnrlon: 40.
        area_thresh: 1000.
        resolution: 'l'
        projection: 'npstere'
        lat_1: 5.
        lon_0: -20.
        rsphere: [6378137.00, 6356752.3142]
        boundinglat: 45.
        anchor: W
    fig_par:
        figsize: [7, 7]
        axlocs: [0.1, 0.05, 0.7, .85]
        # axlocs = [0.1, 0.05, 0.9, 0.8]
"""


class TestMapping:
    @pytest.fixture(autouse=True, params=region_list)
    def setup(self, request):
        self.region = request.param

    # Test the system mapping DB
    def test_system_mapdb(self):
        self.map_par, self.fig_par = map_regions(self.region)
        if self.region == "default":
            assert self.map_par['llcrnrlat'] == 0
            assert self.map_par['llcrnrlon'] == -180
        if self.region == "polarcat":
            assert self.map_par['llcrnrlat'] == 35.
            assert self.map_par['llcrnrlon'] == -95.

    # Test a user-provided mapping DB
    def test_user_mapdb(self, tmpdir):
        # Create a temporay file as the user DB
        p = tmpdir.join("userdb.yml")
        p.write(user_mapdb)
        userdb_file = str(p.realpath())
        os.environ['REFLEXIBLE_MAPDB'] = userdb_file

        self.map_par, self.fig_par = map_regions(self.region)
        if self.region == "default":
            assert self.map_par['llcrnrlat'] == 0.5
            assert self.map_par['llcrnrlon'] == -180.7
        if self.region == "polarcat":
            assert self.map_par['llcrnrlat'] == 90.
            assert self.map_par['llcrnrlon'] == 95.

    # Test passing map_par as param
    def test_map_par(self):
        map_par = {'llcrnrlat': 1, 'llcrnrlon': 2}
        self.map_par, self.fig_par = map_regions(self.region, map_par=map_par)
        assert self.map_par['llcrnrlat'] == 1
        assert self.map_par['llcrnrlon'] == 2

    # Test passing fig_par as param
    def test_fig_par(self):
        fig_par = {'figsize': [7, 3], 'axlocs': [0.1, 0.2, .7, .9]}
        self.map_par, self.fig_par = map_regions(self.region, fig_par=fig_par)
        assert self.fig_par['figsize'] == [7, 3]
        assert self.fig_par['axlocs'] == [0.1, 0.2, .7, .9]
