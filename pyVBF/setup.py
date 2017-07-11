# pyVBF install script

# Check for numpy, future versions won't strictly require it.
try:
  import numpy
except:
  print("This package currently requires numpy")
  from sys import exit
  exit(1)

from os.path import expandvars
from distutils.core import setup
from distutils.extension import Extension

# Set this to False to use the included .cpp file, rather than
# generating your own version with cython.  Future versions 
# will not use cython by default.
USE_CYTHON = True

ext = '.pyx' if USE_CYTHON else '.cpp'

extensions = [
  Extension("pyVBF", ["pyVBF"+ext],
    language = "c++",
    include_dirs = [expandvars('${VBFPATH}/include'), numpy.get_include()],
    libraries = ['VBF', 'bz2', 'z'],
    library_dirs = [expandvars('${VBFPATH}/lib/')],
    runtime_library_dirs = [expandvars('${VBFPATH}/lib/')],)
]

if USE_CYTHON:
  from Cython.Build import cythonize

setup(name = "VBF Reader",
      version="0.1",
      description="extension module for python to read VERITAS VBF files",
      author="Jonathan Eisch",
      author_email="jeisch@iastate.edu",
      include_dirs = [numpy.get_include()],
      ext_modules = cythonize(extensions) if USE_CYTHON else extensions,
)


