#!/usr/bin/env python

def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    import os.path as op

    auto_name = op.split(op.dirname(__file__))[-1]
    config = Configuration(auto_name, parent_package, top_path)

    config.add_data_files('matlab_eig.m')

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
