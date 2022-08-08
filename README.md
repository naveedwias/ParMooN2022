ParMooN - Parallel Mathematics and object oriented Numerics
==========

Introduction
------------

This is the introduction.


Quick start with ParMooN
------------

You are eager to get started and have no no time to read any verbose documentation?

This is the minimal out-of-source build instruction.
Write the following into a shell 

    mkdir parmoon
    cd parmoon
    mkdir code
    git clone https://lab.wias-berlin.de/parmoon/parmoon.git code
    mkdir build
    cd build/
    cmake ../code
    make check doc -j

The '-j' for using all available cores is optional.

Requirements
----------
- PETSc requires gfortran
- Doxygen shall be installed in order to build documentation correctly

Recommended start with ParMooN
----------

ParMooN depends on and uses several external libraries. These partly depend on
each other in somewhat difficult ways. In order to facilitate the process as
much as possible you can use PETSc to handle external libraries completely. This
is the recommended way even though it is also possible to use libraries
installed on your system. The following guidelines create a certain directory
structure which can of course be changed. First download and configure PETSc as
given below

    cd <path to software directory>
    mkdir petsc
    git clone -b release https://gitlab.com/petsc/petsc.git petsc
    cd petsc
    ./configure --with-cc=gcc --with-cxx=g++ --with-fc=gfortran --with-debugging=0 COPTFLAGS='-O3 -march=native -mtune=native' CXXOPTFLAGS='-O3 -march=native -mtune=native' FOPTFLAGS='-O3 -march=native -mtune=native' --download-fblaslapack --download-openmpi --download-scalapack --download-mumps --download-metis --download-suitesparse --download-parmetis --download-hdf5 --download-zlib

This will take time and requires internet access because many external libraries
are downloaded and compiled. Then you have to compile PETSc itself, the exact
command is written at the end of the configure output. Something like

    make PETSC_DIR=<path to software directory>/petsc PETSC_ARCH=arch-linux2-c-opt all

This again will take time. Finally you can (but need not to) check the
installation, again the exact command is given at the end of the previous output

    make PETSC_DIR=<path to software directory>/petsc PETSC_ARCH=arch-linux2-c-opt check

You will be able to compile executable files for a some applications, such as
time-dependent or stationary convection--diffusion and Navier--Stokes equations.
Additionally, ParMooN can be used as a library, which is installed at a 
user-chosen directory. ParMooN can be set up using the PETSC_DIR and PETSC_ARCH
arguments with a specified directory for the documentation and installation as
follows:

    cd <path to software directory>
    mkdir parmoon
    cd parmoon
    git clone https://lab.wias-berlin.de/parmoon/parmoon.git code
    mkdir documentation
    mkdir build
    mkdir install_dir
    cd build
    cmake -C ../initial_cache.cmake ../code

The last command requires a text file 'initial_cache.cmake' in your `parmoon`
directory. A template for such a file is as follows:

    # decide if ParMooN is built to run sequentially or in parallel
    set(PARMOON_PARALLEL_TYPE "MPI"
        CACHE STRING "The used parallelism. Choose between MPI or SEQUENTIAL")
    # tell ParMooN where you installed PETSc
    set(PETSC_DIR "<path to software directory>/petsc"
        CACHE PATH "path to where you downloaded PETSc")
    set(PETSC_ARCH "arch-linux2-c-opt"
        CACHE STRING "name of subdirectory within PETSC_DIR where you installed PETSc")
    # in case you want to use ParMooN as a library this is where ParMooN will be installed
    set(CMAKE_INSTALL_PREFIX "<path to software directory>/parmoon/install_dir"
        CACHE PATH "path to where you want this build to install")
    # ParMooN has a doxygen-generated documentation, this is where it is built
    set(PARMOON_DOCUMENTATION_DIRECTORY "<path to software directory>/parmoon/documentation"
        CACHE PATH "path to where the doxygen documentation should be created")
    # determine level of optimizations and debug information
    set(CMAKE_BUILD_TYPE "Release"
        CACHE STRING "choose between 'Debug', 'Release', 'MinSizeRel', 'RelWithDebInfo'")
    # choose a C compiler
    set(CMAKE_C_COMPILER gcc
        CACHE FILEPATH "C compiler, typically 'gcc' or 'clang'")
    # choose a C++ compiler
    set(CMAKE_CXX_COMPILER g++
        CACHE FILEPATH "C++ compiler, typically 'g++' or 'clang++'")

Make sure to properly adjust all paths according to your setup. Alternatively,
you can pass all these variables to `cmake`, each as a separate parameter,
e.g., `cmake -D PARMOON_PARALLEL_TYPE=MPI -D PETSC_DIR=...`.

**Note**:
For MacOS users you might need `PETSC_ARCH=arch-darwin-c-opt`, on Linux also 
`PETSC_ARCH=arch-linux-c-opt` (the '2' is missing) occurs. Make sure to use the
one that PETSc tells you. In all cases, it is the name of a directory in
`<path to software directory>/petsc`

This will configure ParMooN. Then in the build directory you can compile.
If you specify `CMAKE_BUILD_TYPE=DEBUG` instead of `RELEASE` then a debug build
is done.
In addition, if you choose `PARMOON_PARALLEL_TYPE=MPI` instead of `SEQUENTIAL`
a parallel version is compiled.

As a start it makes sense to compile and run all the unit and regression tests:

    make check

Use the option '-j4' to use 4 parallel builds to speed up this process. Next

    make doc

will build the Doxygen documentation (this may take a while), open the file 
<path to software directory>/parmoon/documentation/html/index.html in a web
browser to see it. There is more information there for you to read. Finally,
to build and install the libraries, type

    make install


Writing H5 files in the parallel case
----------

Depending on the filesystem and various other factors, parallel writing of H5
files may take a very long time, stall entirely, or even crash. This is
especially common when writing to network storage through NFS. There are
a handful of options for dealing with this issue:

- Write somewhere else. This may not be feasible, particularly when
distributing computations over multiple physical devices or using limited user
accounts on HPC servers.
- Fiddle with NFS settings. Also not feasible without privileged access.
- Change the I/O driver used by MPI. For example: when using OpenMPI, the
default I/O driver is OMPIO, which will tend to (near-)deadlock itself in
`fcntl` as each process tries to acquire a write lock on the file. ROMIO,
some version of which should be available in a default OpenMPI install, handles
parallel writing differently; try passing `-mca io romio321` (or similar,
depending on your version) to `mpiexec`/`mpirun`.
- Set the option `output_collective_h5_files` to `false`. This forces ParMooN
to collect all data at the root node and write it from there, for improved
stability at the cost of performance. For large computations this aggregation
step is not usually the bottleneck.

The safest fallback is disabling collective writing altogether by setting
`output_collective_h5_files` to `false`. Depending on your setup, using ROMIO
may work for a few hundred writes and then unpredictably fail with a
nonspecific I/O error, for instance.

Writing compressed H5 files
----------

Some computations may generate very large amounts of data. Although XDMF/HDF5
output is considerably more space efficient than pure XML, it may nevertheless
become appropriate to compress the output by setting the option
`output_compress_h5_files` to `true`.

This works only if your HDF5 library was built with zlib. If activating
compressed output results in strange HDF5 errors, you may need to rebuild HDF5
with zlib. If you're using PETSC, the best way to make sure of this is to run
its configuration script with the `--download-zlib` argument before building
everything.

Building with the right HDF5 library
----------

ParMooN can use the HDF5 library whenever solutions are written into files which
then can be read by visualization programs such as ParView. HDF5 distinguishes
between a parallel and a serial version and it is important to link to the right
one, depending on which configuration you chose for ParMooN, see the cmake flag
PARMOON_PARALLEL_TYPE. If the wrong one is found you will see a message
indicating that a "non-compatible hdf5 library [was] found". Then you have to
rerun cmake with a parameter HDF5_ROOT and give the path to the correct version
of HDF5, e.g.
'cmake -DHDF5_ROOT=/usr/lib/x86_64-linux-gnu/hdf5/openmpi/ ../code'. Currently,
we have no way of automating this easily, so it has to be done by the user.


Using installed ParMooN libraries in another project
----------

Assume you have installed ParMooN as described above and want to use it within
another cmake project/application. The directory is set up as follows

     cd <path to application directory>
     mkdir application_code
     mkdir build

For this demonstration, assume a C++-file (main.cpp) containing
'int main(int, char* argv[])' and a CMakeLists.txt is in the application_code
directory. Then the CMakeLists.txt file couldbe

    cmake_minimum_required(VERSION 3.9)

    find_package(parmoon)
    add_executable(main-example "${CMAKE_CURRENT_SOURCE_DIR}/main.cpp" )
    target_link_libraries(main-example PUBLIC parmoon::parmoon_2d_SEQUENTIAL)

This uses the sequential 2d library of ParMooN, you can use the 3d library as
well. The parallel version (parmoon::parmoon_3d_MPI) can be compiled with the
appropriate cmake-flag, `-D PARMOON_PARALLEL_TYPE=MPI`, when building the
ParMooN core libraries. Next, you have to tell cmake where ParMooN was 
installed:

    cd build
    cmake ../application_code -Dparmoon_DIR=<path to software directory>/parmoon/install_dir/lib/cmake/parmoon



Note that you have to provide the specified subdirectory ('lib/cmake/parmoon')
of the install directory where you installed ParMooN. 
Check the directory where the ParMooN library has been installed, it might depend on the system (e.g. `lib`, `lib64`, etc.)
Finally, you can build
your application, yielding the executable 'main-example':

    make main-example
