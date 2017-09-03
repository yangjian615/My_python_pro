#!/usr/bin/env pyton
# coding=utf-8


def show ():
    # coding=utf-8

    print  '''
        '==================================================================================='
        pycharm console 2016.1 与 ipython5.0 冲突
        用/Users/yangjian/Documents/MYgithub/My_python_pro/pydev_ipython_console_011_change.py
        替代/Applications/PyCharm.app/Contents/helpers/pydev/_pydev_bundle/pydev_ipython_console_011.py

        '==================================================================================='


        '==================================================================================='
        # 1.3.4.2.4. Well-known (& more obscure) file formats
        #
        # HDF5: h5py, PyTables
        # NetCDF: scipy.io.netcdf_file, netcdf4-python, ...
        # Matlab: scipy.io.loadmat, scipy.io.savemat
        # MatrixMarket: scipy.io.mmread, scipy.io.mmwrite
        # IDL: scipy.io.readsav
        '==================================================================================='

        '==================================================================================='
        # 中文字符前面加入
        # coding=utf-8
        '==================================================================================='

        '==================================================================================='
        #PYTTHON 路径添加的三种方式
        # (1) sys.path.append()函数添加
        import sys
        sys.path
        sys.path.append('/Library/Python/2.7/site-packages')

        # (2) 修改环境变量
        # .bashrc文件添加
        # export PYTHONPATH=$PYTHONPATH:/Library/Python/2.7/site-packages

        # (3)在site-packages路径下添加一个路径配置文件,文件的扩展名为.pth,内容为要添加的路径即可
        # 如mypython_path.pth，必须以.pth为后缀
        '==================================================================================='

        #(4)__name__ 作用可以调试代码当模块被其他调用不执行; 调试的时候又可以执行
        '==================================================================================='

        '==================================================================================='
        #(5) __init__.py 的用法:
        可以为空是指包的标识

        将包下的module全部导入
         __all__ = ['*']
        from package import *


        '==================================================================================='

        '==================================================================================='

        '''



if __name__ == '__main__':
    print '======Test========='
    show()
