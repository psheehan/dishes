from setuptools import setup
from setuptools.extension import Extension
from Cython.Build import cythonize
import numpy as np

# Now define the setup for the package.

setup(name="dishes", \
        version="0.0.1", \
        author="Patrick Sheehan", \
        author_email="psheehan@nrao.edu", \
        description="Tools for working with radio interferometry visibilities and images", \
        long_description=open("README.md","r").read(), \
        long_description_content_type="text/markdown", \
        url="https://github.com/psheehan/dishes", \
        packages=["dishes","dishes.imaging","dishes.interferometry",\
        "dishes.spectroscopy"], \
        install_requires=['numpy','scipy','matplotlib',\
        'h5py','Cython','astropy','dynesty','numba'])
