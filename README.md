# libOL version 1.40
Quick and easy spatial localization with _octree_

# Overview
the libOL first store a mesh made of vertices, edges and triangles in an octree structure with a very small memory footprint.
Subsequently, you can perform geometrical queries very quickly on this mesh:
- retrieve the closest entity from a given set of coordinates
- build the list of mesh entities than are include in a given bounding box

# Build
Simply follow these steps:
- unarchive the ZIP file
- `cd libOL-master`
- `cmake .`
- `make`
- `make install`

Optionally, you may download some sample meshes to run the examples:
- you need to install the [libMeshb](https://github.com/LoicMarechal/libMeshb) from GitHub
- manually download files from the *Git LFS* repository: [sample files](sample_meshes/)
- move them into /opt/libOL/sample_meshes/
- uncompress them with `lzip -d *.meshb.lz`
- you may now enter /opt/libOL/examples directory and run the various examples

# Usage
It is made of a single ANSI C file and a header file to be compiled and linked alongside the calling program.  
It may be used in C, C++, Fortran 77 and 90 programs.  
Tested on Linux, Mac OS X, and Windows 7-10.
