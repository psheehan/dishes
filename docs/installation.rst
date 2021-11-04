============
Installation
============

Installing the code with Anaconda
"""""""""""""""""""""""""""""""""

Anaconda is probably the easiest way to install dishes because the dependency GALARIO is easily conda-installable. dishes is not yet available for a standard conda install, but you can first install galario with Anaconda and then pip install dishes:

   ::

       conda install galario -c conda-forge
       pip install dishes

Installing the code with pip
""""""""""""""""""""""""""""

1. In a terminal, run:
   ::

       pip install dishes

2. Install GALARIO. Unfortunately, GALARIO is not pip-installable by default, so you will need to follow the instructions `here <https://mtazzari.github.io/galario/>`_. Alternatively, a pip-installable version of GALARIO is in the works and can be installed like so:
   ::

       pip install git+https://github.com/psheehan/galario.git@add_unstructured

but note that this is a fork of GALARIO and is not yet completely merged.

Installing the code manually
""""""""""""""""""""""""""""

If you would like the most cutting-edge version of dishes, with updates that go beyond the pre-packaged versions, you can download and compile the code yourself following these instructions:

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
   
   or

   ::

       conda build dishes -c conda-forge
       conda install dishes -c conda-forge --use-local

   depending on what method you prefer best.

