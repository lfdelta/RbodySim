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
  float mass, invmass;
  float inertia, invinertia; // moment of inertia
  float radius; // distance from CoM to farthest vertex
  vec2* normals; // normal vector to the edge; normal[i] between vertices[i] and vertices[i+1]
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
    // accumulate radius
    float radsq = 0;
    for (int i=0; i < nverts; i++) {
      vertices[i] = verts[i] - CoM;
      radsq = fmax(radsq, vertices[i].dot(vertices[i]));
    }
    radius = sqrt(radsq);

    if (mass == 0)
      mass = area;

    inertia = MomentOfInertia(vertices, nverts, mass);
    invmass = 1/mass;
    invinertia = 1/inertia;

    // compute normal vectors
    for (int i=0; i < nverts; i++) {
      vec2 v = vertices[(i+1)%nverts] - vertices[i];
      vec2 n(v[1], -v[0]);
      normals[i] = n.normalized();
    }

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
  virtual void worldCoords(vec2* outVerts, vec2* outNormals) const = 0;

  // takes a vertex index and a time to integrate along
  // returns the vertex's world position, linearly interpolated
  virtual vec2 worldVertLerp(int ind, float dt) = 0;

  // takes a normal vector index and a time to integrate along
  // returns the normal's world orientation, linearly interpolated
  virtual vec2 worldNormLerp(int ind, float dt) = 0;

  // calculates all vertices' world positions, lerped by timestep t,
  //   and stores the values in the out-parameter,
  //   which must be of length nverts
  virtual void worldCoordsLerp(vec2* outVerts, float dt) = 0;
};



// a Dynamic Rigid Body, with simulated physics
class Dynamic : public Rigidbody {
public:

  Dynamic(const vec2* verts, const int& sz, const vec2& pos=vec2(0,0), const float& rot=0, const float& m=0, const vec2& vel=vec2(0,0), const float& rotvel=0)
  :Rigidbody(verts, sz, pos, rot, m), velocity(vel), force(vec2(0,0)), rotspeed(rotvel), torque(0) {}

  vec2 velocity; // units per second
  vec2 force; // units per second squared; flushed with each timestep
  float rotspeed; // rad/s
  float torque; // rad/s/s; flushed with each timestep

  // calculates the vertices and normals in world coordinates,
  //   and stores the values in the out-parameters,
  //   which must be of length nverts
  void worldCoords(vec2* outVerts, vec2* outNormals) const {
    for (int i=0; i < nverts; i++) {
      outVerts[i] = position + (rotmat * vertices[i]);
      outNormals[i] = rotmat * normals[i];
    }
  }

  // takes a vertex index and a time to integrate along
  // returns the vertex's world position, lerped by the timestep t
  vec2 worldVertLerp(int ind, float dt) {
    mat2 rot = RotationMatrix(rotation + dt * rotspeed);
    vec2 pos = position + velocity*dt;
    return pos + (rot * vertices[ind]);
  }

  // takes a normal vector index and a time to integrate along
  // returns the normal's world orientation, linearly interpolated
  vec2 worldNormLerp(int ind, float dt) {
    mat2 rot = RotationMatrix(rotation + dt * rotspeed);
    return rot * normals[ind];
  }

  // calculates all vertices' world positions, lerped by timestep t,
  //   and stores the values in the out-parameter,
  //   which must be of length nverts
  void worldCoordsLerp(vec2* outVerts, float dt) {
    for (int i=0; i < nverts; i++)
      outVerts[i] = worldVertLerp(i, dt);
  }
};


// a Static Rigid Body, which is fully constrained and does not respond to physics
class Static : public Rigidbody {
public:
  vec2* bakedVerts;
  vec2* bakedNorms;

  Static(const vec2* verts, const int sz, const vec2& pos=vec2(0,0), const float rot=0, const float m=0)
  :Rigidbody(verts, sz, pos, rot, m), bakedVerts(new vec2[sz]), bakedNorms(new vec2[sz]) {
    // bake in final world coordinates
    for (int i=0; i < nverts; i++) {
      bakedVerts[i] = position + (rotmat * vertices[i]);
      bakedNorms[i] = rotmat * normals[i];
    }

    // bake in final AABB
    float minX=0, maxX=0, minY=0, maxY=0;
    for(int i=0; i<nverts; i++) {
      vec2 vert = bakedVerts[i];
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

  ~Static() {
    delete[] bakedVerts;
    delete[] bakedNorms;
  }

  void worldCoords(vec2* outVerts, vec2* outNormals) const {
    for (int i=0; i < nverts; i++) {
      outVerts[i] = bakedVerts[i];
      outNormals[i] = bakedNorms[i];
    }
  }

  vec2 worldVertLerp(int ind, float dt) {
    return bakedVerts[ind];
  }

  vec2 worldNormLerp(int ind, float dt) {
    return bakedNorms[ind];
  }

  void worldCoordsLerp(vec2* outVerts, float dt) {
    for (int i=0; i < nverts; i++)
      outVerts[i] = bakedVerts[i];
  }
};

#endif //_RIGIDBODY_