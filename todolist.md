### HIGH PRIORITY

### STANDARD PRIORITY
- reserve keywords and data structures for hybrid and 2nd order elements
- add a procedure to detect volume subdomains
- build the octree in parallel by inserting batches of entities whose bounding box do not overlap
- add test scenarios: surface proximity, level set, distributed surface volume collision, all to all mapping

### DONE
- generate local diameter size map
- add volume elements kinds
- handle 32 or 64-bit integers
- handle 32 or 64-bit floating points
- add a mode flag to handle user's tables that either range from (0 -> n-1) or (1 -> n)
- make the GetNearest() thread safe: allocate thread local stacks
- translate the documentation to english
- optional mode for faster queries (+35%) at the cost of 2.5 the memory footprint
- add a NewOctreeFromSTL() procedure for simpler and faster octree build and queries
- added an example to illustrate the use of filtering procedures
- added an example that computes the distances between a volume mesh's vertices and a surface mesh
