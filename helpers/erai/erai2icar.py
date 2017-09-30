#!/usr/bin/env python
import os,traceback,sys
import config
import io_routines
import output
import convert

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
    