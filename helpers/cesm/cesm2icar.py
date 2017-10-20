#!/usr/bin/env python
import os,traceback,sys
import datetime

import config
import io_routines
import output
import convert

def main(info):

    for k in info.keys():
        if k!="times" and k!="lat_data" and k!="lon_data":
            print(k,info[k])

    print(info.times[0],info.times[-1])

    curtime=info.times[0]
    firsttime=curtime

    timesteps_per_year=365*4 # no leap calendar, 4 steps per day
    starttime=0
    endtime=timesteps_per_year
    for i in range(info.nyears):
        print(curtime)
        raw_data=io_routines.load_data(firsttime,info,starttime,endtime)
        processed_data=convert.cesm2icar(raw_data)
        output.write_file(curtime,info,processed_data)

        curtime=datetime.datetime(curtime.year+1,curtime.month,curtime.day)
        starttime=endtime
        endtime=endtime+timesteps_per_year

        # curpos+=raw_data.atm.ntimes
        # curtime=info.times[curpos]
        # curtime=info.times[-1]



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
