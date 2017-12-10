#!/usr/bin/env python
import os,traceback,sys
import config
import io_routines
import output
import convert

def main(info):
    
    for k in info.keys():
        if k!="times" and k!="lat_data" and k!="lon_data":
            print(k,info[k])
    
    print(info.times[0],info.times[-1])

    curyear=info.times[0].year
    lastyear=info.times[0].year-1
    for curtime in info.times:
        if curtime.year>lastyear:
            raw_data=io_routines.load_data(curtime,info)
            processed_data=convert.cmip2icar(raw_data)
            output.write_file(curtime,info,processed_data)
            lastyear=curtime.year
            
            
    

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
    