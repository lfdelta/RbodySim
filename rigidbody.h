#ifndef _RIGIDBODY_
#define _RIGIDBODY_

#include "simulate.h"



// a class representing a convex 2D polygon
class Rigidbody {
public:
  vec2* vertices; // local space, counterclockwise
  int nverts;
  vec2 position; // world space
  float rotation; // radians
  float mass;
  vec2* normals; // normal vector to the edge; normal[i] = vertices[i+1] - vertices[i]
  mat2 rotmat; // rotation matrix; should be updated as needed with updateRotMat()
  vec2 AABB[2]; // local axis-aligned bounding box; lower-left and upper-right

  Rigidbody(const vec2* verts, const int sz, const vec2& pos=vec2(0,0), const float rot=0, const float m=0)
  :vertices(new vec2[sz]), normals(new vec2[sz]), nverts(sz), position(pos), rotation(rot), mass(m) {
    initRbody(verts);
  }

  ~Rigidbody() {
    delete[] vertices;
    delete[] normals;
  }

  // used for instance construction
  void initRbody(const vec2* verts) {
    float area;
    vec2 CoM;
    CoM = CenterOfMass(verts, nverts, area);

    // shift the object's local origin to its center of mass
    for (int i=0; i < nverts; i++)
      vertices[i] = verts[i] - CoM;

    // compute normal vectors
    for (int i=0; i < nverts; i++) {
      vec2 v = vertices[(i+1)%nverts] - vertices[i];
      if (abs(v[0]) < EPSILON) {
        normals[i] = vec2(sign(v[1]), 0);
      } else {
        float m = v[1]/v[0];
        float ny = -sign(v[0]) * sqrt(1 / (1 + m*m));
        float nx = -m * ny;
        normals[i] = vec2(nx, ny);
      }
    }

    if (mass == 0)
      mass = area;

    updateRotMat();
    updateAABB();
  }

  void updateRotMat() {
    rotmat = RotationMatrix(rotation);
  }

  // stores the AABB as the rotated but local coordinates
  void updateAABB() {
    float minX=0, maxX=0, minY=0, maxY=0;
    for(int i=0; i<nverts; i++) {
      vec2 vert = rotmat * vertices[i];
      float x = vert[0];
      float y = vert[1];
      if (x < minX) {minX = x;}
      if (x > maxX) {maxX = x;}
      if (y < minY) {minY = y;}
      if (y > maxY) {maxY = y;}
    }
    AABB[0] = vec2(minX, minY);
    AABB[1] = vec2(maxX, maxY);
  }

  // calculates the vertices and normals in world coordinates,
  //   and stores the values in the out-parameters,
  //   which must be of length nverts
  void worldCoords(vec2* outVerts, vec2* outNormals) const {
    for (int i=0; i < nverts; i++) {
      outVerts[i] = position + (rotmat * vertices[i]);
      outNormals[i] = rotmat * normals[i];
    }
  }
};


// class StaticBody: public Rigidbody {
// public:

//   StaticBody(const vec2* verts, const int sz, const vec2& pos=vec2(0,0), const float rot=0)
//   :Rigidbody(verts, sz, pos, rot) {
//     bakeWorldCoords();
//   }

//   // stores the AABB and vertices in world coordinates, for collision-checking efficiency
//   void bakeWorldCoords() {
//     for (int i=0; i < nverts; i++) {
//       vertices[i] = position + (rotmat * vertices[i]);
//       normals[i] = rotmat * normals[i];
//     }

//     float minX=vertices[0][0], minY=vertices[0][1];
//     float maxX=minX, maxY=minY;
//     for(int i=1; i<nverts; i++) {
//       vec2 vert = vertices[i];
//       float x = vert[0];
//       float y = vert[1];
//       if (x < minX) {minX = x;}
//       if (x > maxX) {maxX = x;}
//       if (y < minY) {minY = y;}
//       if (y > maxY) {maxY = y;}
//     }
//     AABB[0] = vec2(minX, minY);
//     AABB[1] = vec2(maxX, maxY);
//   }
// };



// a Dynamic Rigid Body, with simulated physics
class Dynamic : public Rigidbody {
public:

  Dynamic(const vec2* verts, const int& sz, const vec2& pos=vec2(0,0), const float& rot=0, const float& m=0, const vec2& vel=vec2(0,0), const float& rotvel=0)
  :Rigidbody(verts, sz, pos, rot, m), velocity(vel), force(vec2(0,0)), rotspeed(rotvel), torque(0) {}

  vec2 velocity; // units per second
  vec2 force; // units per second squared; flushed with each timestep
  float rotspeed; // rad/s
  float torque; // rad/s/s; flushed with each timestep
};

#endif //_RIGIDBODY_