#!/usr/bin/env python
import os,traceback,sys
import config
import io_routines
import output
import convert

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Requirements: (on NCAR's Cheyenne)
#   - module load ncl
#   - the config.py file requires a link to mygis.py, found in icar/helpers/lib. Either add this
#     to your path with sys.append({yourpath}/icar/helpers/lib) (in config.py) or ln -s in the folder dir.
#   - since large amounts of data are generated, it makes sense to copy the entire /erai dir to scratch
#     before running. (note previous link to mygis.py in this case!)
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def main(info):

    for i in range(info.ntimes):
        raw_data=io_routines.load_data(info.times[i],info)
        processed_data=convert.era2icar(raw_data)
        output.write_file(info.times[i],info,processed_data)


if __name__ == '__main__':
    try:
        info=config.parse()
        config.update_info(info)

        exit_code = main(info)
        if exit_code is None:
            exit_code = 0
        sys.exit(exit_code)
    except KeyboardInterrupt as e: # Ctrl-C
        raise e
    except SystemExit as e: # sys.exit()
        raise e
    except Exception as e:
        print('ERROR, UNEXPECTED EXCEPTION')
        print(str(e))
        traceback.print_exc()
        os._exit(1)
