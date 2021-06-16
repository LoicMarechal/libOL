[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## libOL version 1.77
Quick and easy spatial localization with _octree_

## Overview
the **libOL** first store a mesh made of vertices, edges and triangles in an octree structure with a very small memory footprint.
Subsequently, you can perform geometrical queries very quickly on this mesh:
- retrieve the closest entity from a given set of coordinates
- build the list of mesh entities than are include in a given bounding box
- project a vertex on any kind of geometrical entity
- launch a ray and get the first intersected entity
- all query operations can be performed in parallel as the library is thread safe

## Build

### Prerequisites for *Linux* or *macOS*
- Install [CMake](https://cmake.org/files/v3.7/cmake-3.7.2-win64-x64.msi)
- A valid C99 compiler
- Open a shell window

### Prerequisites for *Windows*
- You first need to install [CMake](https://cmake.org/files/v3.7/cmake-3.7.2-win64-x64.msi). Do not forget to choose "add cmake to the path for all users", from the install panel.
- Then you need a valid C compiler like the free [Visual Studio Community 2019](https://www.visualstudio.com/vs/visual-studio-express/)
- Open the x64 Native Tools Command Prompt for VS (or x86 if you need to build a 32-bit version)

### Build commands for all platforms
- unarchive the ZIP file
- `cd libOL-master`
- `mkdir build`
- `cd build`
- `cmake ..`
- `cmake --build . --target install`

### Optional steps
Optionally, you may download some sample meshes to run the examples:
- you need to install the [libMeshb](https://github.com/LoicMarechal/libMeshb) from GitHub
- manually download files from the *Git LFS* repository: [sample files](sample_meshes/)
- move them into /opt/libOL/sample_meshes/
- uncompress them with `lzip -d *.meshb.lz`
- you may now enter /opt/libOL/examples directory and run the various examples

When speed is is more critical than memory, you can compile the library with the `-DWITH_FAST_MODE` option in order to speed-up queries by 35%, at the cost of 2.5 times the memory footprint.

## Usage
It is made of a single *ANSI C* file and a header file to be compiled and linked alongside the calling program.  
It may be used in C, C++, Fortran 77 and 90 programs.  
Tested on *Linux*, *macOS*, and *Windows 7->10*.

Here is a basic example that builds an octree from a surface mesh made of vertices (NmbVer vertices stored in VerTab[]) and triangles (NmbTri elements stored in TriTab[]), then searches for the closest triangle from a set of coordinates:

```C++
   int64_t OctIdx;
   int     TriIdx;
   double  dis, crd[3] = { 0.2, 2.5, -3.1 };

   // Build an octree from a surface mesh
   OctIdx = LolNewOctree(  NmbVer, VerTab[1], VerTab[2],
                           0, NULL, NULL,
                           NmbTri, TriTab[1], TriTab[2],
                           0, NULL, NULL,
                           0, NULL, NULL,
                           0, NULL, NULL,
                           0, NULL, NULL,
                           0, NULL, NULL,
                           1, 1);

   // Find the closest triangle from a given point
   TriIdx = LolGetNearest(OctIdx, LolTypTri, crd, &dis, 0., NULL, NULL, 0);

   // Print the triangle's index and the distance from the source point
   printf("The closest triangle from (%g %g %g) is: %d, distance = %g\n",
          crd[0], crd[1], crd[2], TriIdx, dis);

   // Free the octree
   LolFreeOctree(OctIdx);
```
