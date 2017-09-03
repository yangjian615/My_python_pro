#!/usr/bin/env python

import datetime as dt
import glob
import numpy as np
from collections import Counter # for size_mode()

def mon_no(month):
    '''To return the Mon as a string number''' #This could return a string, e.g. '01', '02', etc?
    months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    return str(months.index(month) + 1)

def flt2date(decday):
    '''From decimal day to datetime.datetime in UTC'''
    # Convert back to seconds
    seconds = decday * (24*60*60)
    return dt.datetime.utcfromtimestamp(seconds)

def days_in_month(any_day):
    '''Given a decimal day, return how many days in that month'''
    # Convert to date object
    any_day = flt2date(any_day)
    next_month = any_day.replace(day=28) + dt.timedelta(days=4)  # this will never fail
    ldom = next_month - dt.timedelta(days=next_month.day) #Last Day Of Month
    return ldom.day

def str2flt(yyyymmddhhmmss):
    '''From yyyymmddhhmmss to decimal day'''
    notz = 'UTC' #%Z	Time zone name UTC, EST, CST; Should be naive but no...
    dttime = dt.datetime.strptime(yyyymmddhhmmss+notz, '%Y%m%d%H%M%S%Z')
    #Convert it to an integer for the recarray
    #dt.datetime.timestamp(dttime) is unreliable as I can't make it UTC
    #timestamp = dttime.replace(tzinfo=timezone.utc).timestamp() timezone??
    #or by calculating the timestamp directly:
    timestamp = (dttime - dt.datetime(1970, 1, 1)) / dt.timedelta(seconds=1)
    dayfloat = timestamp / (24*60*60)
    return dayfloat

def date2flt(dttime):
    '''From datetime.datetime to decimal day''' # subsection of str2flt()
    timestamp = (dttime - dt.datetime(1970, 1, 1)) / dt.timedelta(seconds=1)
    dayfloat = timestamp / (24*60*60)
    return dayfloat

def days_in_year(year):
    '''Given a year, return how many days in that year'''
    # Convert to date object
    nyd = str2flt(year+'0101000000')
    next_nyd = str2flt(str(int(year) + 1) +'0101000000')
    ndays = next_nyd - nyd
    return ndays

def lineofinventory(lineoffile):
    '''
    Takes a line (string) from the inventory csv as input:
    '"C1_CP_FGM_SPIN","2001-02-28 10:36:32.348","2001-03-01 00:00:00.0","12016","10"'
    Output is a tuple of creation day as float, filesize in bytes and plot time as float
    '''

    # Split by comma
    lineoffile = lineoffile.split(',')
    # For first date&time, take off leading and trailing double quote and convert to date using given format
    stad = dt.datetime.strptime(lineoffile[1][1:-1], '%Y-%m-%d %H:%M:%S.%f')
    # Convert date to flt (date2flt gives decimal day)
    stad = date2flt(stad)
    # Repeat for end of period
    endd = dt.datetime.strptime(lineoffile[2][1:-1], '%Y-%m-%d %H:%M:%S.%f')
    endd = date2flt(endd)
    # Number of data points in that line, without double quotes
    npts = int(lineoffile[3][1:-1])
    # Return tuple of results
    return (stad, endd, npts)

def read_inventory(sc, path, y, m):
    '''
    Given a spacecraft ID (e.g. 'C1'), a path to the inventory lists, a year and a month,
    return an ndarray containing the times () and points in each line.
    '''
    # inventory filename example: C1_CG_FGM_SPIN_200101Inventory.csv
    # C1_CP_EFW_L3_E3D_INERT_200101Inventory.csv
    inv_req = sc+'*'+y+m+'Inventory.csv'
    # Find all matching files
    inventory = glob.glob(path+inv_req)
    if len(inventory) == 0:
        print("Are all the directories correct? read_inventory can't find any inventory files.")
    if len(inventory) > 1:
        print('More than one inventory file?!')
        
    # Open file and count lines
    with open(inventory[0]) as f:
        size=len([0 for _ in f])
    
    if size > 0:    
    # Set up an np.ndarray
        inv_file_data = np.zeros(size-1, dtype=[('time1',float), 
                                           ('time2', float),
                                           ('npts', int)])
        
    # Open the file again and read it with lineofinventory
        with open(inventory[0]) as f:
            for c, line in enumerate(f):
                if line[1] == 'C': #ignore the header
                    dataline = (lineofinventory(line))
                    inv_file_data[c-1] = dataline #ignoring the header
    else:
        inv_file_data = np.zeros(1, dtype=[('time1',float), 
                                           ('time2', float),
                                           ('npts', int)])

    return inv_file_data

def query_inventory(inv_data, plot_st):
    '''
    Check inventory file for lines enclosing the plot time and 
    return number of seconds in that hour with data.
    
    20160609 Revised to not care about the order of the lines, 
    after advice from Stratos.
    '''
    
    # Plots cover 1 hour 
    plot_end = plot_st + (1/24)
    
    # All lines with start time before plot end
    inv_indices = [x for x in range(0,len(inv_data)) if inv_data['time1'][x] <= plot_end]
    # Now narrow that down to lines with end time after plot start
    inv_indices = [x for x in inv_indices if inv_data['time2'][x] >= plot_st]

    # Relevant Inventory Lines:
    ril = inv_data[inv_indices]
    #print(ril)
    
    # Calculate time in inventory lines if they contain points
    total_overlap = 0.0
    for p in ril:
        # From http://stackoverflow.com/questions/325933/determine-whether-two-date-ranges-overlap
        # by Vitalii Fedorenko
        overlap = max(0, min(p['time2'], plot_end) - max(p['time1'], plot_st))
        if overlap > 0 and p['npts'] > 0:
            total_overlap += overlap
            
    return total_overlap

def lineofplotlist(lineoffile):
    '''
    Takes a line of the plot list; modification date (system time, so LT), size of file, and filename
    containing start and end of plot (1 hour):
    'Sun Nov  1 09:25:57 2015 7366 C1_CG_FGM_BMAG_CAA__20010301_000000_20010301_010000_V01.png'
    Returns modification date in seconds, size of file and decimal day of start of plot.
    '''
    
    # Split the whole line by spaces
    lparts = lineoffile.split()

    # Creation date and time    
    # Convert and replace the 'Mmm' to 'n' month number (string)
    lparts[1] = mon_no(lparts[1])
    # Join them back together, convert to datetime object (which won't work in a list)
    # using the given format
    creation_time = dt.datetime.strptime(' '.join(lparts[1:5]), '%m %d %H:%M:%S %Y')
    # Then convert the datetime object to decimal day since 1/1/1970 LT (not UTC since this is system time)
    #creation = dt.datetime.timestamp(creation_time)
    creation = date2flt(creation_time)
    
    # Filesize in bytes
    filesize = lparts[5]
    
    # Split up the plot filename by underscores
    fnparts = lparts[-1].split('_')
    # Start time, taken from end of filename since parameter may contain '_', converted to decimal day
    start = str2flt(fnparts[-5]+fnparts[-4])
    
    # Return a tuple of the results
    return (creation, filesize, start)

def read_png_lists(sc, path, y, m):
    '''
    Given a spacecraft ID (e.g. 'C1'), a path to the png lists, a year and a month,
    return a dictionary where each key is SCParameter_yyyymm
    '''
    
    # Filename example: C1FGM_BGSE_PHI200106, C1EFW_L3_E3D_INERT_X200101
    pnglist_req = sc+'*'+y+m+'.txt'
    # Find all matching files; these will include the full path
    plotlistlist = glob.glob(path+pnglist_req)
    # Extract all the parameter names into a list, this should retain its order. 
    files = [x[x.find(sc):] for x in plotlistlist]
    parameters = [x[x.find('_')+1:x.find(y+m)-1] for x in files]
    
    # Set up an empty dictionary
    scym_pngs = {}

    # This assumes that the parameters are singular
    for i in range(0,len(parameters)):
        
        # Set up the key for the dictionary
        scymkey = sc + parameters[i] + '_' + y + m

        # Open file and check size
        with open(plotlistlist[i]) as f:
            size=len([0 for _ in f])
        
        # Set up receiving ndarray
        pngs = np.zeros(size, dtype=[('creation',float), 
                                     ('filesize', int), 
                                     ('plot_time', float)])
        
        # Open and read file again, send each line to Pregen1hPlot 
        # and save the returned tuple in the array
        with open(plotlistlist[i]) as f:
            for c, line in enumerate(f):
                dataline = (lineofplotlist(line))
                pngs[c] = dataline
                
        # Once one file (parameter) has been read in, store 
        # one month together, 1 s/c, 1 month, all parameters
        scym_pngs[scymkey] = pngs
    
    # Return whole dictionary
    return scym_pngs

def size_mode(d_sizes):
    '''
    Given a list of integers, find the mode. If more than one mode, return the integers 
    with the same count as the most common. If no mode, return all integers.
    '''
    # Count up the different sizes
    d_size_counts = Counter(d_sizes) #<class 'collections.Counter'>
    #e.g. Counter({35137: 14, 96384: 1, 43906: 1, 106550: 1, 105157: 1, 132854: 1, 132839: 1, 68484: 1, 102534: 1, 44717: 1, 65438: 1})#print(type(d_size_counts)) 

    # How many has the most common got?
    modeinfo = d_size_counts.most_common(1) #<class 'list'>
    # e.g., modeinfo = [(7046, 21)]
    modecount = modeinfo[0][1]

    # Count the counts of the different sizes
    mode_check = Counter(d_size_counts.values()) #<class 'collections.Counter'>
    #e.g. Counter({1: 10, 14: 1})
    if mode_check[modecount] > 1:
        #if mode_check[modecount] != 24: # This doesn't appear to be needed...
            #print('Non-singular mode!', mode_check[modecount])
        # This will return only the modes that have a count equal to the most common
        # So if there are 'twins', all twins will be returned
        # If all the values are different, count=1 and it will return all
        return [key for key in d_size_counts if d_size_counts[key] == modecount]
    else:
        return [modeinfo[0][0]]



def Tally(all_citizens, key):
    # these list comprehensions should return lists of indices where text matches
    humans = [k for k,x in enumerate(all_citizens[key]['status']) if x == b'Human']
    empties = [k for k,x in enumerate(all_citizens[key]['status']) if x == b'Empty']
    cylons = [k for k,x in enumerate(all_citizens[key]['status']) if x == b'Cylon']
    suspects = [k for k,x in enumerate(all_citizens[key]['status']) if x == b'Suspect']
    twins = [k for k,x in enumerate(all_citizens[key]['status']) if x == b'Twin']
    others = [k for k,x in enumerate(all_citizens[key]['status']) if x != b'Suspect' \
              and x != b'Empty' and x != b'Cylon' and x != b'Human' and x != b'Twin']

    return [len(humans), len(empties), len(cylons), len(twins), len(suspects), len(others)]

