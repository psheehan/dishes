from numpy.distutils.core import setup, Extension
from Cython.Build import cythonize

# Set up the extension modules.

libinterferometry = cythonize([\
        Extension('dishes.interferometry.libinterferometry',\
            ["dishes/interferometry/libinterferometry.pyx"],\
            libraries=["m"], extra_compile_args=['-ffast-math'])])[0]

libimaging = cythonize([Extension('dishes.imaging.libimaging',\
        ["dishes/imaging/libimaging.pyx"], libraries=[], \
        extra_compile_args=[])])[0]

setup(name="dishes", version="1.0.0", \
        packages=[\
        "dishes",\
        "dishes.constants", \
        "dishes.imaging",\
        "dishes.interferometry", \
        "dishes.spectroscopy",\
        "dishes.table"], \
        ext_modules=[libinterferometry, libimaging])
