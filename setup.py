import setuptools
from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration

sources = [
    'src/basis.f90',
    'src/evaluate.f90',
    'src/jacobian.f90',
    'src/knots.f90',
    'src/paramuni.f90',
    ]

config = Configuration('MBI')
config.add_extension('MBIlib', sources=sources)#, extra_compile_args=['--fcompiler=gnu95'])

kwds = {'install_requires': ['numpy', 'scipy'],
        'version': '0.1',
        'zip_safe': False,
        'license': 'Apache (v2.0)',
        'packages': ['MBI'],
        'package_data': { 'MBI': ['examples/*.py'] },
        'include_package_data': True,
        #'script_args': ['build', '--fcompiler=gnu95', 'install']
        }
kwds.update(config.todict())

setup(**kwds)

