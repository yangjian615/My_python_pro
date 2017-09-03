#!/usr/bin/env python

import numpy as np
import math
import datetime as dt

def var_dict(line):
    splitline = line.split('=')
    for x,y in enumerate(splitline):
        splitline[x] = y.strip()
    return splitline

def str2flt(yyyymmddThhmmssZ):
    '''From yyyymmddThhmmssZ to seconds since 19700101'''
    yyyymmddhhmmss = (yyyymmddThhmmssZ[:-1]).replace('T', '')
    notz = 'UTC' #%Z Time zone name UTC, EST, CST; Should be naive but no...
    dttime = dt.datetime.strptime(yyyymmddhhmmss+notz, '%Y-%m-%d%H:%M:%S.%f%Z')
    timestamp = (dttime - dt.datetime(1970, 1, 1)) / dt.timedelta(seconds=1)
    dayfloat = timestamp
    return dayfloat

def read_cef(filename):
    # Get the dataset and file information
    with open(filename) as f:
        d = {}
        data = 0
        meta = 0
        ndata = 0
        var_order = [] # because the dictionary won't preserve the order
        for line in f:
            if data == 0 and line[0:16] != 'DATA_UNTIL = EOF':
                if line[:12] == 'END_VARIABLE':
                    meta = 0
                if line[:14] == 'START_VARIABLE' or meta == 1:
                    meta = 1
                    var_data = var_dict(line)
                    if var_data[0] =='START_VARIABLE':
                        varname = var_data[1]
                        d[varname] = {}   
                        var_order.append(varname)
                    else:
                        d[varname][var_data[0]] = var_data[1]
            elif data == 0 and line[:16] == 'DATA_UNTIL = EOF': 
                data = 1
            elif data == 1 and line[:8] != '!RECORDS':
                ndata += 1
        d['npts'] = ndata
    
    # Construct the data array
    build = []
    full_list = []
    for var in var_order:
        n = int(d[var]['SIZES'])
        t = d[var]['VALUE_TYPE'].lower()
        if t == 'iso_time':
            t = 'float'
        if n > 1:
            vect = d[var]['LABEL_1'].split(',')
            for i,j in enumerate(vect):
                vect[i] = (j.strip())[1:-1]
                build.append((vect[i], t))
                full_list.append(var)
        else:
            build.append((d[var]['LABLAXIS'][1:-1], t))
            full_list.append(var)
    
    data_array = np.zeros(ndata, dtype = build)
    
    # Read in the data
    with open(filename) as f:
        data = 0
        for i,line in enumerate(f):
            if data == 0 and line[:16] == 'DATA_UNTIL = EOF':
                data = 1
                j = i+1
            elif data == 1 and line[:8] != '!RECORDS':
                # take dollar sign from end of line
                lineparts = line[:line.find('$')].split(',')
                # this assumes that the time is first
                timestamp = str2flt(lineparts[0].strip())
                linelist = [timestamp]
                x = 1
                for dpt in full_list[1:]:
                    # As there are multiple identical varnames in full_list, 
                    # this assumes that a data point in the metadata is a scalar.
                    if 'DATA' in d[dpt]:
                        linelist.append(float(d[dpt]['DATA']))
                        if float(d[dpt]['SIZES']) > 1:
                            print('Non-scalar metadata data point')
                    else:
                        linelist.append(float(lineparts[x]))
                        x += 1
                data_array[i-j] = tuple(linelist)
    
    return data_array