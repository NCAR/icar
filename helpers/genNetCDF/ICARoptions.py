import sys

# class generates ICAR namelist file
class ICARoptions:
    def __init__(self):
        self.f = open('icar_options.nm', 'w')
        self.gen(self.model_version)
        self.gen(self.output_list)
        self.gen(self.physics)
        self.gen(self.files_list)
        self.gen(self.parameters)
        self.clean_var_list()
        self.gen(self.var_list)
        self.gen(self.z_info)
        self.close()
        print("generated ICAR options")


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

    def clean_var_list(self):
        for name, val in self.var_list.items():
            if name == 'name':
                continue
            self.var_list[name] = '"' + val + '"'
            # print(name, val)

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
        'water': 0, 'mp': 2,
        'rad': 0, 'conv': 0,
        'adv': 1, 'wind': 0
    }

    files_list = {
        'name': 'files_list',
        'init_conditions_file': '"init.nc"',
        # 'boundary_files':'TBD',
        'forcing_file_list': '"forcing.nc"'
    }

    z_info = {
        'name': 'z_info',
        'dz_levels': str([50., 75., 125., 200., 300., 400.] + [500.] * 35)[1:-1]
    }


    var_list = {
        # forcing variables
        'name': 'var_list',
        'uvar':'u',
        'vvar':'v',
        'pvar': 'pressure',
        'tvar': 'theta',
        'qvvar': 'qv', # water_vapor
        'hgtvar': 'height',
        'zvar': 'z',
        'latvar': 'lat_m',
        'lonvar': 'lon_m',
        # init conditions variables
        'lat_hi': 'lat_hi',
        'lon_hi': 'lon_hi',
        'hgt_hi': 'hgt_hi', # surface elevation
        'time_var':'time'
    }
        # 'wvar':'w'


    files_list = {
        'name': 'files_list',
        'init_conditions_file': '"init.nc"',
        'boundary_files': '"forcing.nc"'
    }

    parameters = {
        'name': 'parameters',
        'forcing_start_date': '"2020-12-01 00:00:00"',
        'end_date': '"2020-12-02 00:00:00"',
        'calendar': '"standard"',
        'dx': '1000.0',
        'qv_is_relative_humidity':'true',
        't_is_potential':'true',
        'readdz': 'true',
        'nz': '15',
        'z_is_geopotential': 'False',
        'z_is_on_interface': 'False',
        'time_varying_z': 'False',
        'use_agl_height': 'False',
        'smooth_wind_distance': '2000'
    }
