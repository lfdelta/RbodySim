
#ifndef _SIMULATE_
#define _SIMULATE_

#include <eigen3/Eigen/Core>
#include <ctime>
#include <cmath>

typedef Eigen::Vector2f vec2;
typedef Eigen::Matrix2f mat2;
#define EPSILON 0.00000001f



float sign(float x) {
  return (x > 0) - (x < 0);
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

#endif //_SIMULATE_