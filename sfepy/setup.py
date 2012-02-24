import sys
print sys.path

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    import os.path as op

    auto_name = op.split(op.dirname(__file__))[-1]
    config = Configuration('sfepy', parent_package, top_path)

    subdirs = [
        'applications',
        'base',
        'fem',
        'mesh',
        'homogenization',
        'interactive',
        'linalg',
        'mechanics',
        'optimize',
        'physics',
        'postprocess',
        'solvers',
        'terms'
    ]
    for subdir in subdirs:
        config.add_subpackage(subdir)
    import numpy.distutils.misc_util
    print numpy.distutils.misc_util.__file__
    config.make_config_py()
    print config
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
