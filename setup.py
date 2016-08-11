from numpy.distutils.core import setup, Extension
from Cython.Build import cythonize

libinterferometry = cythonize([Extension('dishes.interferometry.libinterferometry',["dishes/interferometry/libinterferometry.pyx"],libraries=["m"],extra_compile_args=['-ffast-math'])])[0]
#libinterferometry = cythonize("dishes/interferometry/libinterferometry.pyx")[0]
libimaging = cythonize("dishes/imaging/libimaging.pyx")[0]

setup(ext_modules=[libinterferometry, libimaging])
