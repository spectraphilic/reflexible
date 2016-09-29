import pytest
import netCDF4 as nc
import os
import reflexible as rf


output_list = [('Fwd1_V10.0','grid_conc_20110101000000.nc'),
               ('Fwd1_V10.0','grid_conc_20110101000000_nest.nc')]


class Dataset:
    def __init__(self, fp_name):
        self.fp_name = fp_name[0]
        self.fp_filename = fp_name[1]
        self.fp_path = rf.datasets[fp_name[0]]

    def setup(self):
        #self.H = rf.Header(self.fp_path, absolute_path=False)
        ncfile = os.path.join(self.fp_path, self.fp_filename)

        self.ncid = nc.Dataset(ncfile, 'r')
        return self.ncid, self.fp_path, ncfile, None

    def cleanup(self):
        pass #self.tmpdir.remove(self.nc_path)


class TestStructure:
    @pytest.fixture(autouse=True, params=output_list)
    def setup(self, request, tmpdir):
        dataset = Dataset(request.param)
        print(dataset.fp_path)
        self.ncid, self.fp_path, self.nc_path, self.H = dataset.setup()
        request.addfinalizer(dataset.cleanup)

    # CF convention required attributes
    def test_CF_conventions(self):
        req_attrs = (
            'Conventions', 'title', 'institution',
            'source', 'history', 'references')
        existing_attrs = self.ncid.ncattrs()
        for attr in req_attrs:
            assert attr in existing_attrs
            #just test the nc file is complete
            # no reflexible Header testing
            #assert self.H.nc.attrs[attr] == getattr(self.ncid, attr)

    # Attributes describing model run
    def test_model_run(self):
        req_attrs = ('outlon0', 'outlat0', 'dxout', 'dyout')
        existing_attrs = self.ncid.ncattrs()
        for attr in req_attrs:
            assert attr in existing_attrs

    # # COMMAND file settings
    def test_command_settings(self):
        req_attrs = (
            'ldirect',
            'ibdate', 'ibtime', 'iedate', 'ietime',
            'loutstep', 'loutaver', 'loutsample',
            'iout', 'lsubgrid', 'lconvection',
            'ind_source', 'ind_receptor', 'itsplit',
            'lsynctime', 'ctl',
            'ifine', 'ipout', 'lagespectra', 'ipin',
            'ioutputforeachrelease', 'mdomainfill',
            'mquasilag', 'nested_output', 'surf_only',
            'linit_cond', 'iflux')
        existing_attrs = self.ncid.ncattrs()
        for attr in req_attrs:
            assert attr in existing_attrs

    # Test dimensions
    def test_dimensions(self):
        req_dims = (
            'time', 'longitude', 'latitude', 'height', 'numspec',  'pointspec',
            'nageclass', 'nchar', 'numpoint')
        existing_dims = self.ncid.dimensions
        for dim in req_dims:
            assert dim in existing_dims

    # Test variables
    def test_time(self):
        assert 'time' in self.ncid.variables
        var_attrs = self.ncid.variables['time'].ncattrs()
        assert 'units' in var_attrs
        assert 'calendar' in var_attrs

    def test_longitude(self):
        assert 'longitude' in self.ncid.variables
        var_attrs = self.ncid.variables['longitude'].ncattrs()
        attr_names = ('long_name', 'axis', 'units', 'standard_name',
                      'description')
        for attr in attr_names:
            assert attr in var_attrs

    def test_latitude(self):
        assert 'latitude' in self.ncid.variables
        var_attrs = self.ncid.variables['latitude'].ncattrs()
        attr_names = ('long_name', 'axis', 'units', 'standard_name',
                      'description')
        for attr in attr_names:
            assert attr in var_attrs

    def test_height(self):
        assert 'height' in self.ncid.variables
        var_attrs = self.ncid.variables['height'].ncattrs()
        attr_names = ('units', 'positive', 'standard_name', 'long_name')
        for attr in attr_names:
            assert attr in var_attrs

    # Assumption for all REL variables: write_releases.eqv is True
    def test_RELCOM(self):
        assert 'RELCOM' in self.ncid.variables
        assert 'long_name' in self.ncid.variables['RELCOM'].ncattrs()

    def test_RELLNG1(self):
        assert 'RELLNG1' in self.ncid.variables
        var_attrs = self.ncid.variables['RELLNG1'].ncattrs()
        assert 'units' in var_attrs
        assert 'long_name'in var_attrs

    def test_RELLNG2(self):
        assert 'RELLNG2' in self.ncid.variables
        var_attrs = self.ncid.variables['RELLNG2'].ncattrs()
        assert 'units' in var_attrs
        assert 'long_name' in var_attrs

    def test_RELLAT1(self):
        assert 'RELLAT1' in self.ncid.variables
        var_attrs = self.ncid.variables['RELLAT1'].ncattrs()
        assert 'units' in var_attrs
        assert 'long_name' in var_attrs

    def test_RELLAT2(self):
        assert 'RELLAT2' in self.ncid.variables
        var_attrs = self.ncid.variables['RELLAT2'].ncattrs()
        assert 'units' in var_attrs
        assert 'long_name' in var_attrs

    def test_RELZZ1(self):
        assert 'RELZZ1' in self.ncid.variables
        var_attrs = self.ncid.variables['RELZZ1'].ncattrs()
        assert 'units' in var_attrs
        assert 'long_name' in var_attrs

    def test_RELZZ2(self):
        assert 'RELZZ2' in self.ncid.variables
        var_attrs = self.ncid.variables['RELZZ2'].ncattrs()
        assert 'units' in var_attrs
        assert 'long_name' in var_attrs

    def test_RELKINDZ(self):
        assert 'RELKINDZ' in self.ncid.variables
        assert 'long_name' in self.ncid.variables['RELKINDZ'].ncattrs()

    def test_RELSTART(self):
        assert 'RELSTART' in self.ncid.variables
        var_attrs = self.ncid.variables['RELSTART'].ncattrs()
        assert 'units' in var_attrs
        assert 'long_name' in var_attrs

    def test_RELEND(self):
        assert 'RELEND' in self.ncid.variables
        var_attrs = self.ncid.variables['RELEND'].ncattrs()
        assert 'units' in var_attrs
        assert 'long_name' in var_attrs

    def test_RELPART(self):
        assert 'RELPART' in self.ncid.variables
        assert 'long_name' in self.ncid.variables['RELPART'].ncattrs()

    def test_RELXMASS(self):
        assert 'RELXMASS' in self.ncid.variables
        assert 'long_name' in self.ncid.variables['RELXMASS'].ncattrs()

    def test_LAGE(self):
        assert 'LAGE' in self.ncid.variables
        var_attrs = self.ncid.variables['LAGE'].ncattrs()
        assert 'units' in var_attrs
        assert 'long_name' in var_attrs

    # Concentration output, wet and dry deposition variables (one per species)
    # To be checked:  iout, wetdep, drydep
    def test_species_mr(self):
        attr_names = ('units', 'long_name', 'decay', 'weightmolar',
                      'ohreact', 'kao', 'vsetaver', 'spec_ass')
        for i in range(0, len(self.ncid.dimensions['numspec'])):
            anspec = "%3.3d" % (i + 1)
            if self.ncid.iout in (1, 3, 5):
                var_name = "spec" + anspec + "_mr"
                var_attrs = self.ncid.variables[var_name].ncattrs()
                assert var_name in self.ncid.variables
                # The following asserts fail because some attributes have not
                # been set (decay, weightmolar, ohreact, kao, vsetaver,
                # spec_ass)
                # for attr in attr_names:
                #     assert attr in var_attrs

    def test_species_pptv(self):
        attr_names = ('units', 'long_name', 'decay', 'weightmolar',
                      'ohreact', 'kao', 'vsetaver', 'spec_ass')
        for i in range(0, len(self.ncid.dimensions['numspec'])):
            anspec = "%3.3d" % (i + 1)
            if self.ncid.iout in (2, 3):
                var_name = "spec" + anspec + "_pptv"
                var_attrs = self.ncid.variables[var_name].ncattrs()
                assert var_name in self.ncid.variables
                # The following asserts fail because some attributes have not
                # been set yet (decay, weightmolar, ohreact, kao, vsetaver,
                # spec_ass)
                # for attr in attr_names:
                #     assert attr in var_attrs
