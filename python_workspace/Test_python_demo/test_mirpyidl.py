#!/usr/bin/env pyton
# coding=utf-8



# region TEST
print  '''
    '==================================================================================='
     # 测试python 调用 IDl
     #
     #
    '==================================================================================='
    '''

import os
#获取当前工作目录
os.getcwd()
#更改当前工作目录
#'/Applications/exelis/idl82/bin/bin.darwin.x86_64/libXm.3.0.2.dylib'
os.chdir('/Applications/exelis/idl82/bin/bin.darwin.x86_64/')
os.getcwd()




# Import mirpyidl.
import mirpyidl as idl
import numpy as np

# Execute a command in IDL.
# This will print 'Hello mom!' from IDL.
idl.execute("PRINT, 'Hello mom!'")

# Create a numpy array in python.
py_array = np.random.normal(size=1000)

# Copy this data to IDL.
idl.setVariable('idl_array', py_array)


# Calculate the standard devation and the mean in IDL.
idl.execute('idl_stddev = STDDEV(idl_array)')
idl.execute('idl_mean = MEAN(idl_array)')

# Copy the results back to python.
py_stddev = idl.getVariable('idl_stddev')
py_mean = idl.getVariable('idl_mean')

# Print out the results.
print('Mean: {}, StdDev: {}'.format(py_mean, py_stddev))


# endregion