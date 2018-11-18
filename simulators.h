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
  Dynamic** physObjs;
  int nphys;
  Rigidbody** staticObjs;
  int nstatics;
  double tstep; // seconds
  int nsteps;
  vec2 globalForce; // units/second/second
  SimType mode;

  Simulator(Dynamic** phys, const int nphysObjs, Rigidbody** statics, const int nstatObjs,
            const double t=0.01, const int n=1, const vec2& g=vec2(0,0), const SimType type=Sim_Timer)
  :tstep(t), nsteps(n), globalForce(g), mode(type) {
    loadObjects(phys, nphysObjs, statics, nstatObjs);
  }

  ~Simulator() {
    if (physObjs) {
      for (int i=0; i < nphys; i++)
        free(physObjs[i]);
      for (int i=0; i < nstatics; i++)
        free(staticObjs[i]);

      free(physObjs);
      free(staticObjs);
    }
  }

  // copy object arrays into the instance members
  void loadObjects(Dynamic** phys, const int nphysObjs, Rigidbody** statics, const int nstatObjs) {
    if (physObjs) {
      for (int i=0; i < nphys; i++)
        free(physObjs[i]);
      for (int i=0; i < nstatics; i++)
        free(staticObjs[i]);

      free(physObjs);
      free(staticObjs);
    }

    nphys = nphysObjs;
    nstatics = nstatObjs;
    physObjs = (Dynamic**)malloc(nphys * sizeof(Dynamic*));
    staticObjs = (Rigidbody**)malloc(nstatics * sizeof(Rigidbody*));

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
  void runFullSimulation(FILE* outSIM=stdout, FILE* outSYS=stdout, FILE* outSTA=stdout, FILE* outDYN=stdout) {
    clock_t timer;

    switch (mode) {
      case Sim_Timer:
      default:
        timer = clock();
        for (int i=0; i < nsteps; i++)
          simulateTimestep();

        timer = clock() - timer;
        fprintf(outSIM, "%f\n", timer/(double)CLOCKS_PER_SEC);
        break;

      case Sim_FullOutput:
        outputSystemData(outSYS, outSTA, outDYN);
        outputStepData(outSIM);
        for (int i=0; i < nsteps; i++) {
          simulateTimestep();
          outputStepData(outSIM);
        }
        break;
    }
  }

  // prints CSV-formatted information about the simulation variables and static objects
  void outputSystemData(FILE* outSYS=stdout, FILE* outSTA=stdout, FILE* outDYN=stdout) {
    fprintf(outSYS, "tstep,nsteps,g_x,g_y\n");
    fprintf(outSYS, "%f,%d,%f,%f\n", tstep, nsteps, globalForce[0], globalForce[1]);

    for (int i=0; i < nphys-1; i++) {
      Dynamic* rbody = physObjs[i];
      fprintf(outDYN, "%f,%f,", rbody->mass, rbody->inertia);
    }
    fprintf(outDYN, "%f,%f\n\n", physObjs[nphys-1]->mass, physObjs[nphys-1]->inertia);

    // print world-space vertex xs, then ys
    for (int i=0; i < nstatics; i++) {
      Rigidbody* rbody = staticObjs[i];

      vec2 rVerts[rbody->nverts], rNorms[rbody->nverts];
      rbody->worldCoords(rVerts, rNorms);

      for (int j=0; j < rbody->nverts-1; j++)
        fprintf(outSTA, "%f,", rVerts[j][0]);
      fprintf(outSTA, "%f\n", rVerts[rbody->nverts-1][0]);

      for (int j=0; j < rbody->nverts-1; j++)
        fprintf(outSTA, "%f,", rVerts[j][1]);
      fprintf(outSTA, "%f\n", rVerts[rbody->nverts-1][1]);
    }
  }

  // prints CSV-formatted information about the current timestep to stdout
  void outputStepData(FILE* outf=stdout) {
    for (int i=0; i < nphys-1; i++) {
      fprintf(outf, "%f,%f,%f,%f,%f,%f,", physObjs[i]->position[0], physObjs[i]->position[1],
                                          physObjs[i]->velocity[0], physObjs[i]->velocity[1],
                                          physObjs[i]->rotation, physObjs[i]->rotspeed);
    }
    fprintf(outf, "%f,%f,%f,%f,%f,%f\n", physObjs[nphys-1]->position[0], physObjs[nphys-1]->position[1],
                                         physObjs[nphys-1]->velocity[0], physObjs[nphys-1]->velocity[1],
                                         physObjs[nphys-1]->rotation, physObjs[nphys-1]->rotspeed);
  }

  // uses the instance properties to apply physics, handle collisions, and update the system
  virtual void simulateTimestep() = 0;
};



bool AABBoverlap(const Rigidbody* a, const Rigidbody* b) {
  vec2 apos = a->position, bpos = b->position;

  float aL = apos[0]+a->AABB[0][0], aR = apos[0]+a->AABB[1][0], bL = bpos[0]+b->AABB[0][0], bR = bpos[0]+b->AABB[1][0];
  if (aL > bR || bL > aR) // if x doesn't overlap
    return false;

  float aT = apos[1]+a->AABB[1][1], aB = apos[1]+a->AABB[0][1], bT = bpos[1]+b->AABB[1][1], bB = bpos[1]+b->AABB[0][1];
  if (aB > bT || bB > aT) // if y doesn't overlap
    return false;

  return true; // x and y both overlap
}


// helper function for SeparatingAxisOverlap; MTVdist should be initialized to a negative value
// projects each vertex from A and B onto each of A's normal axes, to look for overlap
// if no separation is found, returns true, with out-parameters to accumulate minimum translation vector
bool MinimumTranslationVector(const vec2* aVerts, const vec2* aNorms, const int aSize, const vec2* bVerts, const int bSize,
                              vec2& MTVdir, float& MTVdist) {
  vec2 origin, axis, r; // origin vertex, associated normal, distance to current vertex of interest
  float amin, amax, bmin, bmax; // signed distances from origin vertex, along the current axis

  // for each normal axis of A
  for (int i=0; i < aSize; i++) {
    origin = aVerts[i];
    axis = aNorms[i];

    // project each of A's vertices onto the normal axis
    amin = amax = 0;
    for (int j=0; j < aSize; j++) {
      r = aVerts[j] - origin;
      float d = axis.dot(r);
      if (d < amin)
        amin = d;
      if (d > amax)
        amax = d;
    }

    // project each of B's vertices onto the normal axis
    r = bVerts[0] - origin;
    float d = axis.dot(r);
    bmin = bmax = d;
    for (int j=1; j < bSize; j++) {
      r = bVerts[j] - origin;
      d = axis.dot(r);
      if (d < bmin)
        bmin = d;
      if (d > bmax)
        bmax = d;
    }

    // minimum overlap in the ranges defined by A and B's projections
    float depth = std::fmin(bmax-amin, amax-bmin);

    // if A and B's projections do not overlap, then the normal defines a separating axis
    if (depth < 0)
      return false;

    // accumulate minimum translation vector
    if (depth < MTVdist || MTVdist < 0) {
      MTVdir = axis;
      MTVdist = depth;
    }
  }

  return true; // none of the normals is a separating axis
}


// returns true if the two convex polygons are colliding (a separating axis does NOT exist)
// if A and B are overlapping, outMTV will store the minimum translating vector required to separate them
// reference: http://www.dyn4j.org/2010/01/sat/#sat-mtv
bool SeparatingAxisOverlap(const Rigidbody* a, const Rigidbody* b, vec2& outMTVdir, float& outMTVdist) {
  vec2 aVerts[a->nverts], aNorms[a->nverts], bVerts[b->nverts], bNorms[b->nverts];
  a->worldCoords(aVerts, aNorms);
  b->worldCoords(bVerts, bNorms);

  // these functions accumulate and store MTV info
  if (MinimumTranslationVector(aVerts, aNorms, a->nverts, bVerts, b->nverts, outMTVdir, outMTVdist)) {
    if (MinimumTranslationVector(bVerts, bNorms, b->nverts, aVerts, a->nverts, outMTVdir, outMTVdist)) {
      return true; // no separating axis was found
    }
  }

  return false; // a separating axis was found
}



class EulerPairwise : public Simulator {
  public:
  int physLoops; // number of physics loops per timestep

  EulerPairwise(Dynamic** phys, const int nphysObjs, Rigidbody** statics, const int nstatObjs, const int loops=5,
               const double t=0.01, const int n=1, const vec2& g=vec2(0,0), const SimType type=Sim_Timer)
  :Simulator(phys, nphysObjs, statics, nstatObjs, t, n, g, type), physLoops(loops) {}

  // if SAT collision, move d and apply forces
  void staticCollision(Dynamic* d, Rigidbody* s) {
    vec2 MTVdir; // normal vector to the plane of impact
    float MTVdist = -1;

    if (SeparatingAxisOverlap(d, s, MTVdir, MTVdist)) {
      d->position -= MTVdir * MTVdist; // remove d from its overlap with s
      vec2 perpvel = MTVdir * MTVdir.dot(d->velocity);
      d->velocity -= 2 * perpvel; // reflects velocity about the impact normal (elastically)
    }
  }

  // if SAT collision, change position and velocity of both objects
  void dynamicCollision(Dynamic* a, Dynamic* b) {
    vec2 MTVdir; // normal vector to the plane of impact
    float MTVdist = -1;

    if (SeparatingAxisOverlap(a, b, MTVdir, MTVdist)) {
      // move A and B in proportion to their relative momenta
      float av = a->mass * a->velocity.norm();
      float bv = b->mass * b->velocity.norm();
      float aRatio = av / (av + bv);
      float bRatio = 1 - aRatio;

      vec2 MTV = MTVdir * MTVdist;
      a->position -= aRatio * MTV;
      b->position += bRatio * MTV;

      vec2 aPerpvel = MTVdir * MTVdir.dot(a->velocity);
      vec2 bPerpvel = MTVdir * MTVdir.dot(b->velocity);
      a->velocity -= 2 * aPerpvel;
      b->velocity -= 2 * bPerpvel;
    }
  }

  void simulateTimestep() {
    // speculatively apply physics to each dynamic object, ignoring potential collisions
    for (int i=0; i < nphys; i++) {
      Dynamic* thisphys = physObjs[i];
      thisphys->force += globalForce * tstep;
      thisphys->velocity += thisphys->force * tstep;
      thisphys->rotspeed += thisphys->torque * tstep;
      thisphys->position += thisphys->velocity * tstep;
      thisphys->rotation += thisphys->rotspeed * tstep;
      thisphys->force = vec2(0,0); // flush the accumulated forces
      thisphys->torque = 0;
      thisphys->updateAABB();
    }

    // iterate through objects to check collisions (pairwise)
    for (int t=0; t < physLoops; t++) {

      if (nphys == 1) {
        Dynamic* thisphys = physObjs[0];
        for (int j=0; j < nstatics; j++) {
          Rigidbody* other = staticObjs[j];
          if (AABBoverlap(thisphys, other))
            staticCollision(thisphys, other);
        }
      } else { // nphys > 1
        for (int i=0; i < nphys-1; i++) {
          Dynamic* thisphys = physObjs[i];
          for (int j=i+1; j < nphys; j++) { // all dynamic-dynamic
            Dynamic* other = physObjs[j];
            if (AABBoverlap(thisphys, other))
              dynamicCollision(thisphys, other);
          }

          for (int j=0; j < nstatics; j++) { // all but one dynamic-static
            Rigidbody* other = staticObjs[j];
            if (AABBoverlap(thisphys, other))
              staticCollision(thisphys, other);
          }
        }

        Dynamic* thisphys = physObjs[nphys-1];
        for (int j=0; j < nstatics; j++) { // final dynamic-static
            Rigidbody* other = staticObjs[j];
            if (AABBoverlap(thisphys, other))
              staticCollision(thisphys, other);
        }
      }
    } // end physloop
  }

};

#endif //_SIMULATORS_