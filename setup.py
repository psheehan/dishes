from setuptools import setup
from setuptools.extension import Extension
from Cython.Build import cythonize
import numpy as np

# Set up the extension modules.

libinterferometry = cythonize([\
        Extension('dishes.interferometry.libinterferometry',\
            ["dishes/interferometry/libinterferometry.pyx"],\
            libraries=["m"], extra_compile_args=['-ffast-math'], \
            include_dirs=[np.get_include()], \
            define_macros=[('NPY_NO_DEPRECATED_API', 0)])])[0]

libimaging = cythonize([Extension('dishes.imaging.libimaging',\
        ["dishes/imaging/libimaging.pyx"], libraries=[], \
        extra_compile_args=[], include_dirs=[np.get_include()])])[0]

# Now define the setup for the package.

setup(name="dishes", \
        version="0.0.1", \
        author="Patrick Sheehan", \
        author_email="psheehan@nrao.edu", \
        description="Tools for working with radio interferometry visibilities and images", \
        long_description=open("README.md","r").read(), \
        long_description_content_type="text/markdown", \
        url="https://github.com/psheehan/dishes", \
        packages=[\
        "dishes",\
        "dishes.imaging",\
        "dishes.interferometry", \
        "dishes.spectroscopy"], \
        package_data={\
        'dishes.imaging': ['*.pyx'], \
        'dishes.interferometry': ['*.pyx']}, \
        ext_modules=[libinterferometry, libimaging], \
        install_requires=['numpy','scipy','matplotlib',\
        'h5py','Cython','astropy','dynesty'])
