// Volume, center of mass and inertia tensor of a closed, consistently
// oriented triangle mesh (0th, 1st and 2nd moments, unit density).
// C99, also callable from C++.
#ifndef INERTIA_H
#define INERTIA_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

// Change to float for single precision.
typedef double inertia_real;

typedef struct {
  inertia_real vol;      // volume (= mass at unit density)
  inertia_real com[3];   // center of mass
  inertia_real I[3][3];  // inertia tensor about the center of mass
} inertia_result;

// vtx = [[x,y,z],...]   idx = [[a,b,c],...]   ntri = number of triangles.
// Triangles are counter-clockwise seen from outside.

// Moments are accumulated relative to the first vertex, so accuracy does not
// degrade when the mesh is far from the origin.
void inertia(const inertia_real (*vtx)[3], const uint32_t (*idx)[3], size_t ntri,
             inertia_result *out);

// Moments are accumulated about the origin,
// slightly faster, for meshes that are already near the origin.
void inertia_unshifted(const inertia_real (*vtx)[3], const uint32_t (*idx)[3], size_t ntri,
                       inertia_result *out);

#ifdef __cplusplus
}
#endif

#endif
