# libOL version 1.40
Quick and easy spatial localization with _octree_

# Overview
the **libOL** first store a mesh made of vertices, edges and triangles in an octree structure with a very small memory footprint.
Subsequently, you can perform geometrical queries very quickly on this mesh:
- retrieve the closest entity from a given set of coordinates
- build the list of mesh entities than are include in a given bounding box

# Build for *Linux* or *macOS*
Simply follow these steps:
- unarchive the ZIP file
- `cd libOL-master`
- `mkdir build`
- `cd build`
- `cmake -DCMAKE_INSTALL_PREFIX=$HOME/local ../`
- `make`
- `make install`

# Build for *Windows*
- You first need to install [CMake](https://cmake.org/files/v3.7/cmake-3.7.2-win64-x64.msi). Do not forget to choose "add cmake to the path for all users", from the install panel.
- Then you need a valid C compiler like the free [Visual Studio Express 2015](https://www.visualstudio.com/vs/visual-studio-express/)
- unarchive the ZIP file
- open the windows shell
- `cd libOL-master`
- `mkdir build`
- `cd build`
- `cmake -DCMAKE_INSTALL_PREFIX=%HOMEPATH%\local ..\`
- `cmake --build . --config Release --target INSTALL`

Optionally, you may download some sample meshes to run the examples:
- you need to install the [libMeshb](https://github.com/LoicMarechal/libMeshb) from GitHub
- manually download files from the *Git LFS* repository: [sample files](sample_meshes/)
- move them into /opt/libOL/sample_meshes/
- uncompress them with `lzip -d *.meshb.lz`
- you may now enter /opt/libOL/examples directory and run the various examples

# Usage
It is made of a single *ANSI C* file and a header file to be compiled and linked alongside the calling program.  
It may be used in C, C++, Fortran 77 and 90 programs.  
Tested on *Linux*, *macOS*, and *Windows 7->10*.
