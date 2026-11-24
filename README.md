# Faster inertia for triangle meshes.

a C99 implementation of the algorithm described in the paper "Faster inertia for triangle meshes. A moment's appreciation for PGA". 

## Use

Copy `inertia.h` and `inertia.c` into your project. `inertia.c` is C99.

```c
#include "inertia.h"

inertia_result r;
inertia(vtx, idx, ntri, &r);    // vtx: double[][3], idx: uint32_t[][3]
// r.vol, r.com[3], r.I[3][3]   // inertia tensor about the center of mass
```

The mesh must be closed, with triangles counter-clockwise seen from outside.
Density is 1, so the mass is the volume. `inertia()` accumulates moments
relative to the first vertex, so results stay accurate when the mesh is far
from the origin. `inertia_unshifted()` accumulates them about the origin,
slightly faster, for meshes that are already near the origin. Both work in 
`double`; change the `inertia_real` typedef in `inertia.h` for another type.
