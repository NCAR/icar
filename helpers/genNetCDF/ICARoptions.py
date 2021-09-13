
class ICARoptions:
    def __init__(self):
        print("ICAR OPTIONS")

        self.f = open('icar_options.nm', 'w')
        self.gen(self.model_version)
        self.gen(self.output_list)
        self.gen(self.physics)
        self.gen(self.files_list)
        self.gen(self.parameters)
        self.gen(self.var_list)
        self.gen(self.z_info)
        self.close()

    def gen(self, nml):
        f = self.f
        f.write("&"+nml['name'])

        i = 0
        for name, val in nml.items():
            if name == 'name':
                continue
            if i != 0:
                f.write(',')
            else:
                f.write('\n')
                i += 1
            f.write(str(name)+'='+str(val))
        f.write('\n/\n')

    def close(self):
        self.f.close()


    # namelist options to write
    model_version = {
        'name': 'model_version',
        'version': 2.0,
        'comment': '"Unit Test Data"'
    }

    output_list = {
        'name': 'output_list',
        'names': '"u","v","ta2m","hus2m", "precipitation", "swe"',
        'outputinterval': 3600,
        'output_file': '"out_"'
    }

    physics = {
        'name': 'physics',
        'pbl': 0,  'lsm': 0,
        'water': 0, 'mp': 1,
        'rad': 0, 'conv': 0,
        'adv': 1, 'wind': 1
    }

    files_list = {
        'name': 'files_list',
        'init_conditions_file': '"topography.nc"'
        # 'boundary_files':'TBD',
        # 'forcing_file_list': '"forcing.nc"'
    }

    z_info = {
        'name': 'z_info',
        'dz_levels': str([50., 75., 125., 200., 300., 400.] + [500.] * 35)[1:-1]
    }


    var_list = {
        'name': 'var_list',
        'pvar': 'pressure',
        'tvar': 'temperature',
        'qvvar': 'water_vapor',
        'hgtvar': 'height',
        'zvar': 'z',
        'latvar': 'lat',
        'lonvar': 'lon',
        'uvar':'u',
        'vvar':'v',
        'wvar':'w',
    }

    files_list = {
        'name': 'files_list',
        'init_conditions_file': 'topography.nc',
        'boundary_files': 'forcing.nc'
    }

    parameters = {
        'name': 'parameters',
        'forcing_start_date': '2020-12-01 00:00:00',
        # 'start_date': '2020-12-01 00:01:00',
        'end_date': '2020-12-01 00:01:00',
        'calendar': 'standard',
        'inputinterval': '3600',
        'dx': '4000.0',
        'readdz': 'true',
        'nz': '15',
        'z_is_geopotential': 'True',
        'z_is_on_interface': 'True',
        'time_varying_z': 'False',
        'use_agl_height': 'False',
        'smooth_wind_distance': '72000'
    }
