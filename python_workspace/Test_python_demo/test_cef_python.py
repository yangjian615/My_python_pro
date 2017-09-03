# coding=utf-8


import sys
from copy import deepcopy
sys.path.append('/Users/yangjian/CEFLIB/PYTHON')  #载入ceflib.so
from ceflib import *

file="/Users/yangjian/Documents/MYgithub/My_python_pro/python_workspace" \
     "/Test_python_demo/test_data/C1_CP_FGM_SPIN__20030303_000000_20030304_000000_V060926.cef"


def run():
    read(file)
    print '''CEF varnames'''
    print '''============================'''
    print "\n".join(varnames())
    close()

    read(file)
    print '''READ metadata'''
    print '''============================'''
    read_metadata(file)
    close()

    read(file)
    print '''READ time_tags'''
    print '''============================'''
    print records()
    print var('time_tags')
    time = deepcopy(var('time_tags'))
    close()

    print '''Converts internal time-tags to ISOTIME strings'''
    print '''======================================='''

    print milli_to_isotime(time[0], 0)

if __name__ == '__main__':
    print '======Debug========='
    run()











