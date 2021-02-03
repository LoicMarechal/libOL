### STANDARD PRIORITY
- reserve keywords and data structures for hybrid and 2nd order elements
- optionally include the LPlib and build the octree in parallel
- add a procedure to detect volume subdomains

### DONE
- generate local diameter size map
- add volume elements kinds
- handle 32 or 64-bit integers
- handle 32 or 64-bit floating points
- add a mode flag to handle user's tables that either range from (0 -> n-1) or (1 -> n)
- make the GetNearest() thread safe: allocate thread local stacks
