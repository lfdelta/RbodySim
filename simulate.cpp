#include "simulate.h"
#include <iostream>
#include <string>

int main(int argc, char* argv[]) {
  float t;
  t = (argc > 1) ? std::stof(argv[1]) : 0;
  mat2 rot = RotationMatrix(t);

  vec2 verts[3];
  verts[0] = vec2(0,0); verts[1] = vec2(0,1); verts[2] = vec2(1,0);
  Rigidbody rbody(verts, 3, vec2(0,0), t);

  std::cout << "AABB {" << rbody.AABB[0][0] << ", " << rbody.AABB[0][1] << "} {" << rbody.AABB[1][0] << ", " << rbody.AABB[1][1] << "}\n";
  //std::cout << rot * vec2(1,0) << "\n";

  return 0;
}