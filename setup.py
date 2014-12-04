from numpy.distutils.core import setup, Extension
from Cython.Build import cythonize

libinterferometry = cythonize("dishes/interferometry/libinterferometry.pyx")[0]
libimaging = cythonize("dishes/imaging/libimaging.pyx")[0]

setup(ext_modules=[libinterferometry, libimaging])
