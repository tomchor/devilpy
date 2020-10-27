__version__ = '0.1'
print('Initiating devilpy version %s'%__version__, end='...')
#from .core import *
#from .routines import *
from . import utils, routines, plot, plot2, physics
from . import vector, numerical, stats, spectral#, nml
from . import io, input
print('done')
