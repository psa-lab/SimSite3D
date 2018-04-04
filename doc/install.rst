.. _install_guide:

**************************
Steps to Install SimSite3D
**************************

.. _virtualenv:

Setup a Virtual Environment for SimSite3D
=========================================

.. note:: You may skip this step if you have a recent version of Python 2
          (2.6+) and you want to install the SimSite3D Python dependencies
          either in your home directory or in the system (e.g. /usr/local).

.. note:: you only need to setup the tools virtualenv and virtualenvwrapper 
          once per machine if you install it in a standard location.  
          Of course, you may have multiple virtual environments, but don't 
          need to replicate the two tools used to work with virtual 
          environments.

Choose a directory to store the virtual environments.  This should be a
directory that is group read/writable and is relatively low latency or local
to the machine.
On the PSA network mounts, this location is currently::

  /soft/linux64/.virtualenvs

.. _Python:

Install a recent version of Python
----------------------------------

Go to python.org and get the latest 2.* version and install it.  Currently the
installed version is at::

  /soft/linux64/Python-2.7.2

This version was installed as an alternate install and should probably only
be used for virtualenvs.

Install a recent version of Boost
---------------------------------

Go to boost.org and get the lastest stable Boost version and install it.
Currently the installed version is at::

  /soft/linux64/boost_1_48_0

.. _note: You will need to move your user-config.jam out of the way if you
          have one, otherwise there will be build issues.  Don't forget to
          move it back once you have built the libraries.

Build the boost::python libary using latest version of Python::

  ./bootstrap.sh --with-libraries=python --with-python-root=/soft/linux64/Python-2.7.2 --prefix=/soft/linux64

Edit the tools/build/v2/site-config.jamfile (in the installed boost directory)
so that we use this version of boost::python and the right python library::

  project site-config ;

  using gcc ;
  # We are using an alternate location install of Python 2.7
  using python : 2.7 : /soft/linux64/Python-2.7.2 ;
  # We had to build boost using the alternate install of Python for ....
  lib boost_python : : <file>/soft/linux64/boost_1_48_0/stage/lib/libboost_python.so ;

  # These are dependencies for SimSite3D -- if someone felt the need they could
  # only link them where needed.
  lib blas : : <file>/usr/lib64/libblas.so.3 ;
  lib lapack : : <file>/usr/lib64/liblapack.so.3 ;
  lib popt : : <file>/usr/lib64/libpopt.so ;


Edit the SimSite3D/python/boost-build.jam to point to this version of boost
build; currently this may be done as::

  boost-build /soft/linux64/boost_1_48_0/tools/build/v2 ;

Install virtualenv
------------------

Next we want to install virtualenv
Get virtualenv from the website::

  cd /soft/linux64/src
  wget http://pypi.python.org/packages/source/v/virtualenv/virtualenv-X.Y.Z.tar.gz

Unzip it in a src directory::

  tar -xvzf virtualenv-X.Y.Z.tar.gz

Build and install it in a local directory::

  cd /soft/linux64/src/virtualenv-X.Y.Z
  /soft/linux64/Python-X.Y.Z/bin/pythonX.Y setup.py install --prefix=/soft/linux64/Python-X.Y.Z

Setup the virtualenv using python and virtualenv.py
---------------------------------------------------

.. note:: If you don't want to use bash, you must follow this section as
          virtualenvwrapper is written in bash

Initialize a virtual environment using::
  
  python virtualenv.py --no-site-packages --distribute <name or path of VE>

Edit the activate script to set the necessary environment variables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Set the environment variables in $VIRTUAL_ENV/bin/activate and/or
$VIRTUAL_ENV/bin/activate.csh

Set the ASCBASE_INSTALL_DIR variable::

  setenv ASCBASE_INSTALL_DIR=</path/to/SimSite3D>

.. _activate_virtualenv:

Activate the virtual environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Activate your desired virtual environment (virtualenv) by changing to that
directory and sourcing the correct activate file (for your shell).  
If you are using csh/tcsh, you may use::

  cd <path to my virtualenv>
  source bin/activate.csh

.. _virtualenvwrapper:

Setup the virtualenv using virtualenvwrapper
--------------------------------------------
.. note:: virtualenvwrapper is quite nice, but it isn't necessary and it
          does require that you use the bash shell

.. note:: you only need to setup virtualenv and virtualenvwrapper once per 
          machine if you install it in a standard location.

I chose to do it in this manner.
Edit your environment and login files to set the following environment
variables (using bash)::

  # the following are needed for virtualenv and virtualenvwrappers
  export PYTHONPATH=$PYTHONPATH:/soft/linux64/Python-X.Y.Z/lib/pythonX.Y/site-packages
  export WORKON_HOME=/soft/linux64/Python-X.Y.Z/.virtualenvs
  source /soft/linux64/Python-X.Y.Z/bin/virtualenvwrapper.sh

Download, build, and install virtualenvwrapper::

  cd /soft/linux64/src/
  wget http://pypi.python.org/packages/source/v/virtualenvwrapper/virtualenvwrapper-2.11.tar.gz 
  tar -xvzf virtualenvwrapper-2.11.tgz
  cd /soft/linux64/src/virtualenvwrapper-2.11
  /soft/linux64/Python-X.Y.Z/bin/pythonX.Y setup.py install --prefix=/soft/linux64/Python-X.Y.Z

.. See the notes from the previous section and set above variables in
   the postactivate script in the virtualenvwrapper directory.

..
   note:: The easiest method is to copy ~/virtualenv_utils/postactivate to 
   to $VIRTUAL_ENV/bin/postactivate.  
   Then, make the edits required for your setup.

Initialize a virtual environment using::

  mkvirtualenv --no-site-packages --distribute <name of new VE>
  
Activate the drugsite-prod virtualenv::
  
  workon <name of new VE>

Installing SimSite3D
====================
   
Before cloning the SimSite3D repository, we will need a recent version of
mercurial.  This can be done by installing mercurial in the virtualenv
using pip::

  pip install mercurial

.. _clone_repo:

Cloning the SimSite3D Mercurial Repository on BitBucket
-------------------------------------------------------

.. _note:: Mercurial is a distributed version control system (dvcs).  This 
           means that, among other things, *every* Mercurial clone is just
           as much a full repository as another (except possibly for those
           clones that only store necessary data to fully rebuild a repository
           -- that is those repositories that one does not use to make
           changes to except via Mercurial).  In other words, the repo on 
           BitBucket is only special as long as we want it to be.

Clone SimSite3D from BitBucket::

  hg clone ssh://hg@bitbucket.org/jeff_vanvoorst/simsite3d <local dir name>

Change directory to the root directory of the repository::

  cd <local dir name>

Install SimSite3D's Requirements
--------------------------------

SimSite3D's code can be divided among three categories: pure C++ code, 
C++ code for boost::python extension modules, and Python code.

  * The core of SimSite3D is the pure C++ code
  * The core code has several dependencies that depend on what you seek to do
    * The molecular surfaces for binding sites require one of the following:
      * Generate surfaces yourself:
        * Generating surfaces using your tool of choice
        * Converting them to the format written by MSMS
        * Passing the surfaces to gen_points via command line flags
      * Having an MSMS binary at the path $ASCBASE_INSTALL_DIR/bin/linux_msms
      * Providing the location of an MSMS binary to gen_points via a 
        command line flag
    * ArtSurf uses numerical linear algebra that should be as robust as 
      reasonably possible
      * We chose to use LAPACK, and require:
        * A Fortran77 LAPACK library
        * A BLAS library that can be loaded by the LAPACK library
  * Boost::Python extension modules exist for some of the higher level C++   
    classes, and they require:
    * Boost headers
    * a recent version of Boost::Python
    * Python development tools (header files and libraries)
  * There are a number of utility functions, etc. that do not require 
    high performance, but rather require less time coding; these tools are
    written in Python and will require:
    * numpy
    * matplotlib
    * sphinx
    * mercurial

