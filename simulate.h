
#ifndef _SIMULATE_
#define _SIMULATE_

#include <eigen3/Eigen/Core>
#include <ctime>
#include <cmath>

typedef Eigen::Vector2f vec2;
typedef Eigen::Matrix2f mat2;
#define EPSILON 0.00001f



float sign(const float x) {
  return (x > 0) - (x < 0);
}

float clamp(const float v, const float min, const float max) {
  return (v < min) ? min : ((v > max) ? max : v);
}

inline int min(const int a, const int b) {
  return (a < b) ? (a) : (b);
}

inline int max(const int a, const int b) {
  return (a > b) ? (a) : (b);
}

// returns a "vector" in the zhat direction
float cross(const vec2 a, const vec2 b) {
  return a[0]*b[1] - a[1]*b[0];
}

// scalar argument is assumed to be a vector in the zhat direction
vec2 cross(const float a, const vec2 b) {
  return a * vec2(-b[1], b[0]);
}

// takes an angle, t, in radians
// returns a 2x2 CCW rotation matrix
mat2 RotationMatrix(const float t) {
  float s = sin(t);
  float c = cos(t);
  mat2 rot;
  rot << c, -s,
         s, c;
  return rot;
}

// reference: https://en.wikipedia.org/wiki/Centroid#Of_a_polygon
vec2 CenterOfMass(const vec2* verts, const int nverts, float& outArea) {
  float area=0, x=0, y=0;
  float a;

  for (int i=0; i < nverts-1; i++) {
    a = (verts[i][0] * verts[i+1][1]) - (verts[i+1][0] * verts[i][1]);
    x += a * (verts[i][0] + verts[i+1][0]);
    y += a * (verts[i][1] + verts[i+1][1]);

    area += a;
  }

  // final iteration: wrap around
  a = (verts[nverts-1][0] * verts[0][1]) - (verts[0][0] * verts[nverts-1][1]);
  x += a * (verts[nverts-1][0] + verts[0][0]);
  y += a * (verts[nverts-1][1] + verts[0][1]);
  area += a;

  area /= 2;
  x /= (6 * area);
  y /= (6 * area);

  outArea = area;
  return vec2(x, y);
}

// returns the moment of inertia of a polygon about the origin (defined implicitly by verts)
// computed as the sum of the moments of the subtriangles (with one vertex at the origin)
float MomentOfInertia(const vec2* verts, const int nverts, const float mass) {
  float area=0, I=0;
  float da, dI;
  vec2 p, q;

  for (int i=0; i < nverts-1; i++) {
    p = verts[i]; q = verts[i+1];
    da = p[0]*q[1] - p[1]*q[0];
    dI = p.dot(p) + p.dot(q) + q.dot(q);

    area += da;
    I += dI * da;
  }

  // final iteration: wrap around
  p = verts[nverts-1]; q = verts[0];
  da = p[0]*q[1] - p[1]*q[0];
  dI = p.dot(p) + p.dot(q) + q.dot(q);
  area += da;
  I += dI * da;

  return mass/area * I/6; // da is double-area, but cancels with /2area
}

#endif //_SIMULATE_