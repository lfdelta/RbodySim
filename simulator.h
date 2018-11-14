#ifndef _SIMULATORS_
#define _SIMULATORS_

#include "simulate.h"
#include "rigidbody.h"


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
            const double t=0.01, const int n=1, const vec2& g=vec2(0,0), const SimType type=Sim_Timer)
  :tstep(t), nsteps(n), globalForce(g), mode(type) {
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



bool AABBoverlap(const Rigidbody a, const Rigidbody b) {
  vec2 apos = a.position, bpos = b.position;

  float aL = apos[0]+a.AABB[0][0], aR = apos[0]+a.AABB[1][0], bL = bpos[0]+b.AABB[0][0], bR = bpos[0]+b.AABB[1][0];
  if (aL > bR || bL > aR) // if x doesn't overlap
    return false;

  float aT = apos[1]+a.AABB[1][1], aB = apos[1]+a.AABB[0][1], bT = bpos[1]+b.AABB[1][1], bB = bpos[1]+b.AABB[0][1];
  if (aB > bT || bB > aT) // if y doesn't overlap
    return false;

  return true; // x and y both overlap
}


// https://www.gamedev.net/forums/topic/563009-sat-translation-vector/
bool MinimumTranslationVector(const float amin, const float amax, const float bmin, const float bmax, const vec2 axis, vec2& MTDaxis, float& sqDepth) {
  return false;
}


/////////////////////////////// http://www.dyn4j.org/2010/01/sat/#sat-mtv

// https://developer.mozilla.org/en-US/docs/Games/Techniques/2D_collision_detection#Separating_Axis_Theorem
// returns true if the two convex polygons are colliding (a separating axis does NOT exist)
// "recurse" is used to check the other permutation without copying code
bool SeparatingAxisOverlap(const Rigidbody a, const Rigidbody b, const bool recurse=true) {
  vec2 nhat, origin, axis, r;
  float amin, amax, bmin, bmax; // signed distances from origin, along the defined axis

  vec2 aVerts[a.nverts], aNorms[a.nverts], bVerts[b.nverts], bNorms[b.nverts];
  a.worldCoords(aVerts, aNorms);
  b.worldCoords(bVerts, bNorms);

  for (int i=0; i < a.nverts; i++) {
    origin = aVerts[i];
    axis = aVerts[(i+1)%a.nverts] - origin;
    nhat = aNorms[i];

    amin = amax = 0;
    for (int j=0; j < a.nverts; j++) {
      r = aVerts[j] - origin;
      float d = axis.dot(r - nhat * (nhat.dot(r)));
      if (d < amin)
        amin = d;
      if (d > amax)
        amax = d;
    }

    r = bVerts[0] - origin;
    float d = axis.dot(r - nhat * (nhat.dot(r)));
    bmin = bmax = d;
    for (int j=1; j < b.nverts; j++) {
      r = bVerts[j] - origin;
      d = axis.dot(r - nhat * (nhat.dot(r)));
      if (d < bmin)
        bmin = d;
      if (d > bmax)
        bmax = d;
    }

    if (bmax < amin || amax < bmin)
      return false; // nhat defines a separating axis
  }

  if (recurse)
    return SeparatingAxisOverlap(b, a, false);
  else
    return true; // there is no separating axis
}



class EulerPairwise : public Simulator {
  int physLoops; // number of physics loops per timestep

  void dynamicCollision(Dynamic a, Dynamic b) {
    // if SAT collision, move/apply forces to both objects
  }

  void staticCollision(Dynamic d, Rigidbody s) {
    // if SAT collision, move/apply forces to d
  }

  void simulateTimestep() {
    // speculatively apply physics to each dynamic object, ignoring potential collisions
    for (int i=0; i < nphys; i++) {
      Dynamic thisphys = physObjs[i];
      thisphys.force += globalForce * tstep;
      thisphys.velocity += thisphys.force * tstep;
      thisphys.rotspeed += thisphys.torque * tstep;
      thisphys.position += thisphys.velocity * tstep;
      thisphys.rotation += thisphys.rotspeed * tstep;
      thisphys.force = vec2(0,0); // flush the accumulated forces
      thisphys.torque = 0;
      thisphys.updateAABB();
    }

    // iterate through objects to check collisions (pairwise)
    for (int t=0; t < physLoops; t++) {

      for (int i=0; i < nphys-1; i++) {
        Dynamic thisphys = physObjs[i];
        for (int j=i+1; j < nphys; j++) {
          Dynamic other = physObjs[j];
          if (AABBoverlap(thisphys, other))
            dynamicCollision(thisphys, other);
        }

        for (int j=0; j < nstatics; j++) {
          Rigidbody other = staticObjs[j];
          if (AABBoverlap(thisphys, other))
            staticCollision(thisphys, other);
        }
      } // end collision
    } // end physloop
  }

};

#endif //_SIMULATORS_