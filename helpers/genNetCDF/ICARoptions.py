# class generates ICAR namelist file
class ICARoptions:
    def __init__(self,
                 filename = 'icar_options.nml',
                 # model namelist
                 model_version = 2.0,
                 model_comment = 'Unit Test Data',
                 # output namelist
                 output_vars = ['u','v','precipitation','swe'],
                 output_interval = 3600,
                 output_file = 'icar_out_',
                 restart_interval = 3600,
                 restart_file = 'icar_rst_',
                 # physics namelist
                 phys_opt_pbl = 0,
                 phys_opt_lsm = 0,
                 phys_opt_water = 0,
                 phys_opt_mp = 2,
                 phys_opt_rad = 0,
                 phys_opt_conv = 0,
                 phys_opt_adv = 1,
                 phys_opt_wind = 0,
                 # files namelist
                 init_conditions_file = 'init.nc',
                 boundary_files = 'forcing.nc',
                 forcing_file_list = [],
                 # z_info namelist
                 dz_levels = [50., 75., 125., 200., 300., 400.] + [500.] * 50,
                 # forcing variables namelist
                 forc_u_var = 'u',
                 forc_v_var = 'v',
                 forc_p_var = 'pressure',
                 forc_t_var = 'theta',
                 forc_qv_var = 'qv',
                 forc_hgt_var = 'height',
                 forc_z_var = 'z',
                 forc_lat_var = 'lat_m',
                 forc_lon_var = 'lon_m',
                 forc_lat_hi_var = 'lat_hi',
                 forc_lon_hi_var = 'lon_hi',
                 forc_hgt_hi_var = 'hgt_hi',
                 forc_time_var = 'time',
                 # parameters namelist
                 start_date = '2020-12-01 00:00:00',
                 end_date = '2020-12-02 00:00:00',
                 calendar = 'standard',
                 input_interval = '3600',
                 dx = '4000.0',
                 qv_is_relative_humidity ='true',
                 readdz = 'true',
                 nz = '15',
                 z_is_geopotential = 'False',
                 z_is_on_interface = 'False',
                 t_is_potential = 'True',
                 time_varying_z = 'False',
                 use_agl_height = 'False',
                 smooth_wind_distance = '72000',
                 # parcels namelist
                 total_parcels = 0):

        # Open file, create namelist objects, then write
        f = open(filename, 'w')
        self.model_version = ModelVersion(filename=f,
                                          version=model_version,
                                          comment=model_comment)

        self.output_list = OutputList(filename=f,
                                      names=output_vars,
                                      outputinterval=output_interval,
                                      output_file=output_file,
                                      restartinterval=restart_interval,
                                      restart_file=restart_file)

        self.physics_list = PhysicsList(filename=f,
                                        pbl=phys_opt_pbl,
                                        lsm=phys_opt_lsm,
                                        water=phys_opt_water,
                                        mp=phys_opt_mp,
                                        rad=phys_opt_rad,
                                        conv=phys_opt_conv,
                                        adv=phys_opt_adv,
                                        wind=phys_opt_wind)

        self.files_list = FilesList(filename=f,
                                    init_conditions_file=init_conditions_file,
                                    boundary_files=boundary_files,
                                    forcing_file_list=forcing_file_list)

        self.z_info_list = ZInfoList(filename=f,
                                     dz_levels=dz_levels)

        self.forcing_var_list = ForcingVarList(filename=f,
                                               uvar=forc_u_var,
                                               vvar=forc_v_var,
                                               pvar=forc_p_var,
                                               tvar=forc_t_var,
                                               qvvar=forc_qv_var,
                                               hgtvar=forc_hgt_var,
                                               zvar=forc_z_var,
                                               latvar=forc_lat_var,
                                               lonvar=forc_lon_var,
                                               lat_hi=forc_lat_hi_var,
                                               lon_hi=forc_lon_hi_var,
                                               hgt_hi=forc_hgt_hi_var,
                                               time_var=forc_time_var)

        self.parameters_list = ParametersList(filename=f,
                                              forcing_start_date=start_date,
                                              end_date=end_date,
                                              calendar=calendar,
                                              inputinterval=input_interval,
                                              dx=dx,
                                              qv_is_relative_humidity =\
                                              qv_is_relative_humidity,
                                              readdz=readdz,
                                              nz=nz,
                                              z_is_geopotential =\
                                              z_is_geopotential,
                                              z_is_on_interface =\
                                              z_is_on_interface,
                                              t_is_potential =\
                                              t_is_potential,
                                              time_varying_z =\
                                              time_varying_z,
                                              use_agl_height =\
                                              use_agl_height,
                                              smooth_wind_distance =\
                                              smooth_wind_distance)

        self.parcels_list = ParcelsList(filename=f,
                                        total_parcels=total_parcels)

        self.generate_all_namelists()
        f.close()


    def generate_all_namelists(self):
        self.model_version.gen()
        self.output_list.gen()
        self.physics_list.gen()
        self.files_list.gen()
        self.z_info_list.gen()
        self.forcing_var_list.gen()
        self.parameters_list.gen()
        self.parcels_list.gen()


class Namelist:
    def __init__(self, kargs):
        self.filename = kargs['filename']
        del kargs['filename']
        self.nml = {}
        for name, val in kargs.items():
            self.nml[name] = val
        self.remove_empty_values()
    def remove_empty_values(self):
        delete = []
        for name, val in self.nml.items():
            if val in ['', ""]:
                delete.append(name)
        for name in delete:
                del self.nml[name]
    def gen(self):
        f = self.filename
        f.write("&"+self.nml['name'])
        del self.nml['name']
        i = 0
        for name, val in self.nml.items():
            if i != 0:
                f.write(', \n')
            else:
                f.write('\n')
                i += 1
            f.write(str(name)+'='+str(val))
        f.write('\n/\n\n')


class ModelVersion(Namelist):
    def __init__(self, **kargs):
        Namelist.__init__(self, kargs)
        self.nml['name'] = 'model_version'

    def gen(self):
        for name, val in self.nml.items():
            if name == 'comment':
                self.nml[name] = '"' + val + '"'
        super().gen()


class OutputList(Namelist):
    def __init__(self, **kargs):
        Namelist.__init__(self, kargs)
        self.nml['name'] = 'output_list'
    def gen(self):
        for name, val in self.nml.items():
            if name == 'names':
                self.nml[name] = '"' + '","'.join(val) + '"'
            if name in ['output_file', 'restart_file']:
                self.nml[name] = '"' + val + '"'
        super().gen()


class PhysicsList(Namelist):
    def __init__(self, **kargs):
        Namelist.__init__(self, kargs)
        self.nml['name'] = 'physics'


class FilesList(Namelist):
    def __init__(self, **kargs):
        Namelist.__init__(self, kargs)
        self.nml['name'] = 'files_list'
    def gen(self):
        for name, val in self.nml.items():
            if name == 'forcing_file_list':
                self.nml[name] = '"' + '","'.join(val) + '"'
            elif name != 'name':
                self.nml[name] = '"' + val + '"'
        super().gen()


class ZInfoList(Namelist):
    def __init__(self, **kargs):
        Namelist.__init__(self, kargs)
        self.nml['name'] = 'z_info'
    def gen(self):
        for name, val in self.nml.items():
            if name == 'dz_levels':
                self.nml[name] = str(val)[1:-1]
        super().gen()


class ForcingVarList(Namelist):
    def __init__(self, **kargs):
        Namelist.__init__(self, kargs)
        self.nml['name'] = 'var_list'

    def gen(self):
        for name, val in self.nml.items():
            if name == 'name':
                continue
            self.nml[name] = '"' + val + '"'
        super().gen()


class ParametersList(Namelist):
    def __init__(self, **kargs):
        Namelist.__init__(self, kargs)
        self.nml['name'] = 'parameters'

    def gen(self):
        for name, val in self.nml.items():
            if name in ['forcing_start_date', 'end_date', 'calendar']:
                self.nml[name] = '"' + val + '"'
        super().gen()


class ParcelsList(Namelist):
    def __init__(self, **kargs):
        Namelist.__init__(self, kargs)
        self.nml['name'] = 'parcel_parameters'
