#include <iostream>
#include <string>

#include "simulate.h"
#include "rigidbody.h"

void printProperties(const Rigidbody* rbody, const std::string title="") {
  if (title[0])
    std::cout << title << "\n";

  vec2 worldvs[rbody->nverts], worldnormals[rbody->nverts];
  rbody->worldCoords(worldvs, worldnormals);
  float m;
  vec2 CoM = CenterOfMass(rbody->vertices, rbody->nverts, m);

  std::cout << "AABB {" << rbody->AABB[0][0] << ", " << rbody->AABB[0][1] << "} {" << rbody->AABB[1][0] << ", " << rbody->AABB[1][1] << "}\n";
  std::cout << "World:";
  for (int i=0; i < rbody->nverts; i++)
    std::cout << " {" << worldvs[i][0] << ", " << worldvs[i][1] << "}";
  std::cout << "\n";
  std::cout << "Nhats:";
  for (int i=0; i < rbody->nverts; i++)
    std::cout << " {" << worldnormals[i][0] << ", " << worldnormals[i][1] << "}";
  std::cout << "\n";
  std::cout << "Local CoM: {" << CoM[0] << ", " << CoM[1] << "}\n";
  std::cout << "Pos: {" << rbody->position[0] << ", " << rbody->position[1] << "}\n";
  std::cout << "Mass: " << rbody->mass << "\n";
  std::cout << "\n";
}


// takes args of form {t} or { x, y, [t] }
int main(int argc, char* argv[]) {
  vec2 x;
  float t;
  t = (argc == 2) ? std::stof(argv[1]) : (argc >= 4 ? std::stof(argv[3]) : 0);
  x = (argc >= 3) ? vec2(std::stof(argv[1]), std::stof(argv[2])) : vec2(0,0);

  // generate test objects
  vec2 triverts[3] = {vec2(0,0), vec2(2,0), vec2(0,2)};
  vec2 squareverts[4] = {vec2(0,0), vec2(1,0), vec2(1,1), vec2(0,1)};
  vec2 diamondverts[4] = {vec2(0, -1), vec2(1,0), vec2(0, 1), vec2(-1, 0)};
  Rigidbody triangle(triverts, 3, x, t);
  Rigidbody square(squareverts, 4, x, t);
  Rigidbody diamond(diamondverts, 4, x, t);

  printProperties(&triangle, "Right Tri");
  printProperties(&square, "Square");
  printProperties(&diamond, "Diamond");

  return 0;
}