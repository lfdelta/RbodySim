#include <eigen3/Eigen/Core>
#include <cmath>

typedef Eigen::Vector2f vec2;
typedef Eigen::Matrix2f mat2;



// takes an angle, t, in radians
// returns a 2x2 CCW rotation matrix
mat2 RotationMatrix(float t) {
  float s = sin(t);
  float c = cos(t);
  mat2 rot;
  rot << c, -s,
         s, c;
  return rot;
}



// a class representing a convex 2D polygon
class Rigidbody {
public:
  Rigidbody(const vec2* verts, const int& sz, const vec2& pos, const float& rot, const float& m)
  :vertices(new vec2[sz]), normals(new vec2[sz]), nverts(sz), position(pos), rotation(rot), mass(m) {
    initRbody(verts);
  }

  Rigidbody(const vec2* verts, const int& sz, const vec2& pos, const float& rot)
  :vertices(new vec2[sz]), normals(new vec2[sz]), nverts(sz), position(pos), rotation(rot), mass(1) {
    initRbody(verts);
  }

  Rigidbody(const vec2* verts, const int& sz)
  :vertices(new vec2[sz]), normals(new vec2[sz]), nverts(sz), position(vec2(0,0)), rotation(0), mass(1) {
    initRbody(verts);
  }

  // used for instance construction
  void initRbody(const vec2* verts) {
    for (int i=0; i < nverts; i++)
      vertices[i] = verts[i];
    resetAABB();
    updateRotMat();
  }

  ~Rigidbody() {
    delete[] vertices;
    delete[] normals;
  }

  void resetAABB() {
    float minX=0, maxX=0, minY=0, maxY=0;
    for(int i=0; i<nverts; i++) {
      vec2 vert = RotationMatrix(rotation) * vertices[i];
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

  void updateRotMat() {
    rotmat = RotationMatrix(rotation);
  }


  // calculates the vertices and normals in world coordinates,
  //   and stores the values in the out-parameters,
  //   which must be of length nverts
  void worldCoords(vec2* worldVerts, vec2* worldNormals) {
    for (int i=0; i < nverts; i++) {
      worldVerts[i] = position + (rotmat * vertices[i]);
      worldNormals[i] = rotmat * normals[i];
    }
  }


  vec2* vertices; // local space, counterclockwise
  int nverts;
  vec2 position; // world space
  float rotation; // radians
  float mass;
  vec2* normals; // normal vector for each edge
  mat2 rotmat; // rotation matrix; should be updated as needed with updateRotMat()
  vec2 AABB[2]; // local axis-aligned bounding box; upper-left and lower-right
};



// a Dynamic Rigid Body, with simulated physics
class Dynamic : public Rigidbody {
public:

  /* should have a constructor with initial velocity and rotspeed */

  vec2 velocity; // units per second
  vec2 force; // units per second squared; flushed with each timestep
  float rotspeed; // rad/s
  float torque; // rad/s/s; flushed with each timestep
};
