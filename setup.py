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

config = Configuration(name='MBI')
config.add_extension('MBIlib', sources=sources)

kwds = {'install_requires':['numpy','scipy'],
        'version': '0.1',
        'zip_safe': False,
        'license': 'LGPL',
        }
kwds.update(config.todict())

setup(**kwds)
