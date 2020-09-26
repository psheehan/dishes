============
Installation
============

Installing the code with Anaconda
"""""""""""""""""""""""""""""""""

Anaconda is probably the easiest way to install dishes because the dependency GALARIO is easily conda-installable.

1. Download the code from GitHub: https://github.com/psheehan/dishes.git

2. In a terminal, in the directory where the code was downloaded to:
   ::

       conda build dishes -c conda-forge
       conda install dishes -c conda-forge --use-local

Installing the code with pip
""""""""""""""""""""""""""""

1. In a terminal, run:
   ::

       pip install dishes

2. Install GALARIO. Unfortunately, GALARIO is not pip-installable, so you will need to follow the instructions `here <https://mtazzari.github.io/galario/>`_.

Installing the code manually
""""""""""""""""""""""""""""

1. Download the code from this webpage. Git clone is recommended if you would like to be able to pull updates:
   ::

       git clone https://github.com/psheehan/dishes.git

2. Install the Python dependencies:

   * numpy  
   * scipy  
   * matplotlib  
   * emcee  
   * corner  
   * hyperion  
   * h5py  
   * mpi4py  
   * galario  
   * Cython  
   * astropy
   * schwimmbad  
   * dynesty

3. In a terminal, go to the directory where the code was downloaded, and into the code directory. Run:
   ::

        python setup.py install
   
   or

   ::
   
        pip install -e .
   
