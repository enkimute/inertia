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

## Performance

From the paper. The methods compared:

- **edge, shifted** and **edge, unshifted** are ours.
- **Kallay**: M. Kallay, [Computing the Moment of Inertia of a Solid Defined
  by a Triangle Mesh](https://doi.org/10.1080/2151237X.2006.10129220),
  Journal of Graphics Tools 11(2), 2006. Signed origin tetrahedra with a
  four-point rule per face, his published code, unmodified.
- **DiCarlo**: A. DiCarlo, A. Paoluzzi, [Fast Computation of Inertia through
  Affinely Extended Euler Tensor](https://doi.org/10.1016/j.cad.2006.03.005),
  Computer-Aided Design 38(11), 2006. Mass, first and second moments in one
  4x4 tensor per tetrahedron. Our implementation.
- **Mirtich**: B. Mirtich, [Fast and Accurate Computation of Polyhedral Mass
  Properties](https://doi.org/10.1080/10867651.1996.10487458), Journal of
  Graphics Tools 1(2), 1996. His published `volInt.c`.
- **Eberly**: D. Eberly, [Polyhedral Mass Properties
  (Revisited)](https://www.geometrictools.com/Documentation/PolyhedralMassProperties.pdf),
  2002. Mirtich's method specialised to triangle meshes; The Geometric Tools implementation.

Floating-point operations per triangle for volume, center of mass and
inertia:

| method | mul | add/sub | total |
|---|---|---|---|
| edge, unshifted | 36 | 39 | 75 |
| edge, shifted | 36 | 48 | 84 |
| Kallay (unshifted) | 45 | 39 | 84 |
| Eberly | 52 | 61 | 113 |
| DiCarlo | 73 | 51 | 124 |

Time per triangle in nanoseconds, double precision. On 5k and 5M triangles.

| | edge, shifted | edge, unshifted | Kallay | DiCarlo | Mirtich | Eberly |
|---|---|---|---|---|---|---|
| i9-13980HX, 5k | 3.81 | 3.28 | 4.04 | 6.28 | 33.8 | 5.93 |
| i9-13980HX, 5M | 5.06 | 4.36 | 5.01 | 7.29 | 39.8 | 7.23 |
| Apple M2, 5k | 3.02 | 2.54 | 3.34 | 4.56 | 32.8 | 5.83 |
| Apple M2, 5M | 3.70 | 3.22 | 3.97 | 5.28 | 34.3 | 6.57 |

Relative error of the inertia tensor with the mesh moved away from the origin:

| offset | edge, shifted | edge, unshifted | Kallay | DiCarlo | Mirtich | Eberly |
|---|---|---|---|---|---|---|
| 0 | 1.4e-14 | 2.5e-14 | 2.6e-14 | 9.0e-15 | 1.0e-14 | 6.3e-15 |
| 10 | 1.6e-14 | 2.6e-10 | 1.6e-10 | 6.5e-11 | 5.7e-11 | 1.1e-11 |
| 100 | 1.3e-14 | 1.4e-6 | 7.8e-5 | 2.9e-5 | 2.9e-7 | 3.3e-8 |
| 1000 | 1.0e-14 | 2.1 | 1.8 | 4.6e-2 | 4.1e-3 | 3.3e-5 |
