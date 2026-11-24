#include "inertia.h"

typedef inertia_real real;

// Calculate volume, com and inertia. (shifted version)
// Moments are accumulated relative to the first vertex, so accuracy does not
// degrade when the mesh is far from the origin.
void inertia(const real (*vtx)[3], const uint32_t (*idx)[3], size_t ntri,
             inertia_result *out) {

  // triangles * ( 36 mul, 48 add )

  real s00 = 0, s11 = 0, s22 = 0, s01 = 0, s02 = 0, s12 = 0;
  real m0 = 0, m1 = 0, m2 = 0;
  real vol = 0;

  // Use first vertex as offset for numerical stability.
  const real *p = ntri ? vtx[idx[0][0]] : (const real[3]){0, 0, 0};
  const real p0 = p[0], p1 = p[1], p2 = p[2];

  for (size_t i = 0; i < ntri; ++i) {

    // grab vertex coords for triangle i, relative to the reference point.
    const real *a = vtx[idx[i][0]], *b = vtx[idx[i][1]], *c = vtx[idx[i][2]];
    const real a0 = a[0] - p0, a1 = a[1] - p1, a2 = a[2] - p2;
    const real b0 = b[0] - p0, b1 = b[1] - p1, b2 = b[2] - p2;
    const real c0 = c[0] - p0, c1 = c[1] - p1, c2 = c[2] - p2;

    // calculate volume contribution. (scale 6)
    // 9 mul, 6 add.
    const real tvol = (c2 * ((a0 * b1) - (a1 * b0)))
                    + (c1 * ((b0 * a2) - (a0 * b2)))
                    + (c0 * ((b2 * a1) - (b1 * a2)));
    vol += tvol;

    // centers of edges
    // 9 add
    const real A0 = a0 + b0, A1 = a1 + b1, A2 = a2 + b2;
    const real B0 = b0 + c0, B1 = b1 + c1, B2 = b2 + c2;
    const real C0 = c0 + a0, C1 = c1 + a1, C2 = c2 + a2;

    // accumulate com. (scale 6*4*vol)
    // 3 mul, 6 add.
    m0 += tvol * (A0 + c0);
    m1 += tvol * (A1 + c1);
    m2 += tvol * (A2 + c2);

    // accumulate second moments (scale 6*4*5*vol)
    // 24 muls, 18 adds.
    s00 += tvol * (A0 * A0 + B0 * B0 + C0 * C0);
    s11 += tvol * (A1 * A1 + B1 * B1 + C1 * C1);
    s22 += tvol * (A2 * A2 + B2 * B2 + C2 * C2);
    s01 += tvol * (A0 * A1 + B0 * B1 + C0 * C1);
    s02 += tvol * (A0 * A2 + B0 * B2 + C0 * C2);
    s12 += tvol * (A1 * A2 + B1 * B2 + C1 * C2);
  }

  // Fix vol and com scale.
  vol = vol / 6;
  const real iv = 1 / (24 * vol), iv2 = (real)(1.0 / 120);
  m0 *= iv; m1 *= iv; m2 *= iv;

  // scale inertia and parallel transport
  const real c00 = vol * m0 * m0, c11 = vol * m1 * m1, c22 = vol * m2 * m2;
  const real c01 = vol * m0 * m1, c02 = vol * m0 * m2, c12 = vol * m1 * m2;

  // move com back.
  out->vol = vol;
  out->com[0] = m0 + p0; out->com[1] = m1 + p1; out->com[2] = m2 + p2;

  // second moments about the com, calculate full inertia tensor.
  const real x00 = iv2 * s00 - c00, x11 = iv2 * s11 - c11, x22 = iv2 * s22 - c22;
  const real x01 = iv2 * s01 - c01, x02 = iv2 * s02 - c02, x12 = iv2 * s12 - c12;
  out->I[0][0] = x11 + x22;
  out->I[1][1] = x00 + x22;
  out->I[2][2] = x00 + x11;
  out->I[1][0] = out->I[0][1] = -x01;
  out->I[2][0] = out->I[0][2] = -x02;
  out->I[2][1] = out->I[1][2] = -x12;
}

// Calculate volume, com and inertia. (unshifted version)
// Moments are accumulated about the origin,
// slightly faster, for meshes that are already near the origin.
void inertia_unshifted(const real (*vtx)[3], const uint32_t (*idx)[3], size_t ntri,
                       inertia_result *out) {

  // triangles * ( 36 mul, 39 add )

  real s00 = 0, s11 = 0, s22 = 0, s01 = 0, s02 = 0, s12 = 0;
  real m0 = 0, m1 = 0, m2 = 0;
  real vol = 0;

  for (size_t i = 0; i < ntri; ++i) {

    // grab vertex coords for triangle i.
    const real *a = vtx[idx[i][0]], *b = vtx[idx[i][1]], *c = vtx[idx[i][2]];
    const real a0 = a[0], a1 = a[1], a2 = a[2];
    const real b0 = b[0], b1 = b[1], b2 = b[2];
    const real c0 = c[0], c1 = c[1], c2 = c[2];

    // calculate volume contribution. (scale 6)
    // 9 mul, 6 add.
    const real tvol = (c2 * ((a0 * b1) - (a1 * b0)))
                    + (c1 * ((b0 * a2) - (a0 * b2)))
                    + (c0 * ((b2 * a1) - (b1 * a2)));
    vol += tvol;

    // centers of edges
    // 9 add
    const real A0 = a0 + b0, A1 = a1 + b1, A2 = a2 + b2;
    const real B0 = b0 + c0, B1 = b1 + c1, B2 = b2 + c2;
    const real C0 = c0 + a0, C1 = c1 + a1, C2 = c2 + a2;

    // accumulate com. (scale 6*4*vol)
    // 3 mul, 6 add.
    m0 += tvol * (A0 + c0);
    m1 += tvol * (A1 + c1);
    m2 += tvol * (A2 + c2);

    // accumulate second moments (scale 6*4*5*vol)
    // 24 muls, 18 adds.
    s00 += tvol * (A0 * A0 + B0 * B0 + C0 * C0);
    s11 += tvol * (A1 * A1 + B1 * B1 + C1 * C1);
    s22 += tvol * (A2 * A2 + B2 * B2 + C2 * C2);
    s01 += tvol * (A0 * A1 + B0 * B1 + C0 * C1);
    s02 += tvol * (A0 * A2 + B0 * B2 + C0 * C2);
    s12 += tvol * (A1 * A2 + B1 * B2 + C1 * C2);
  }

  // Fix vol and com scale.
  vol = vol / 6;
  const real iv = 1 / (24 * vol), iv2 = (real)(1.0 / 120);
  m0 *= iv; m1 *= iv; m2 *= iv;

  // scale inertia and parallel transport
  const real c00 = vol * m0 * m0, c11 = vol * m1 * m1, c22 = vol * m2 * m2;
  const real c01 = vol * m0 * m1, c02 = vol * m0 * m2, c12 = vol * m1 * m2;

  out->vol = vol;
  out->com[0] = m0; out->com[1] = m1; out->com[2] = m2;

  // second moments about the com, calculate full inertia tensor.
  const real x00 = iv2 * s00 - c00, x11 = iv2 * s11 - c11, x22 = iv2 * s22 - c22;
  const real x01 = iv2 * s01 - c01, x02 = iv2 * s02 - c02, x12 = iv2 * s12 - c12;
  out->I[0][0] = x11 + x22;
  out->I[1][1] = x00 + x22;
  out->I[2][2] = x00 + x11;
  out->I[1][0] = out->I[0][1] = -x01;
  out->I[2][0] = out->I[0][2] = -x02;
  out->I[2][1] = out->I[1][2] = -x12;
}
