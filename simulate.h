#include <eigen3/Eigen/Core>
#include <ctime>
#include <cmath>
#include <iostream>

typedef Eigen::Vector2f vec2;
typedef Eigen::Matrix2f mat2;



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



// https://en.wikipedia.org/wiki/Centroid#Of_a_polygon
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



// a class representing a convex 2D polygon
class Rigidbody {
public:
  vec2* vertices; // local space, counterclockwise
  int nverts;
  vec2 position; // world space
  float rotation; // radians
  float mass;
  vec2* normals; // normal vector for each edge
  mat2 rotmat; // rotation matrix; should be updated as needed with updateRotMat()
  vec2 AABB[2]; // local axis-aligned bounding box; upper-left and lower-right

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

    if (mass == 0)
      mass = area;

    resetAABB();
    updateRotMat();
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
  void worldCoords(vec2* outVerts, vec2* outNormals) const {
    for (int i=0; i < nverts; i++) {
      outVerts[i] = position + (rotmat * vertices[i]);
      outNormals[i] = rotmat * normals[i];
    }
  }
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
};



enum SimType {
  Sim_Timer,      // only run the raw simulation, and compute the time required
  Sim_FullOutput  // print the system state for each timestep
};

class Simulator {
public:
  Dynamic* physObjs;
  int nphys;
  Rigidbody* staticObjs;
  int nstatics;
  double tstep; // seconds
  int nsteps;
  vec2 globalForce; // units/second/second
  SimType mode;

  Simulator(const Dynamic* phys, const int nphysObjs, const Rigidbody* statics, const int nstatObjs,
            const double t=0.01, const int n=1, const vec2& g=vec2(0,0))
  :tstep(t), nsteps(n), globalForce(g) {
    loadObjects(phys, nphys, statics, nstatics);
  }

  ~Simulator() {
    if (physObjs) {
      free(physObjs);
      free(staticObjs);
    }
  }

  // copy object arrays into the instance members
  void loadObjects(const Dynamic* phys, const int nphysObjs, const Rigidbody* statics, const int nstatObjs) {
    nphys = nphysObjs;
    nstatics = nstatObjs;

    if (physObjs) {
      free(physObjs);
      free(staticObjs);
    }

    physObjs = (Dynamic*)malloc(nphys * sizeof(Dynamic));
    staticObjs = (Rigidbody*)malloc(nstatics * sizeof(Rigidbody));

    for (int i=0; i < nphys; i++)
      physObjs[i] = phys[i];
    for (int i=0; i < nstatics; i++)
      staticObjs[i] = statics[i];
  }

  void simulationSettings(const double t, const int n, const vec2& g) {
    tstep = t;
    nsteps = n;
    globalForce = g;
  }

  // void loadSystem(filestream) //http://www.cplusplus.com/doc/tutorial/files/  //https://www.learncpp.com/cpp-tutorial/186-basic-file-io/

  // performs a complete simulation, running for nsteps
  void runFullSimulation() {
    clock_t timer;

    switch (mode) {
      case Sim_Timer:
      default:
        timer = clock();
        for (int i=0; i < nsteps; i++)
          simulateTimestep();

        timer = clock() - timer;
        fprintf(stdout, "%f\n", timer/(double)CLOCKS_PER_SEC);
        break;

      case Sim_FullOutput:
        outputSystemData();
        for (int i=0; i < nsteps; i++) {
          simulateTimestep();
          outputStepData();
        }
        break;
    }
  }

  // prints
  void outputSystemData() {}

  // prints CSV-formatted information about the current timestep to stdout
  void outputStepData() {}

  // uses the instance properties to apply physics, handle collisions, and update the system
  virtual void simulateTimestep();
};