#ifndef _CONTINUOUS_
#define _CONTINUOUS_

#include "simulators.h"

// Minkowski sum vertex
// used to determine rigidbody overlap and contact points
struct MinkVert {
  vec2 pos; // Minkowski-space location
  short Aind; // index of corresponding vertex in A
  short Bind;

  MinkVert() {
    pos = vec2(0,0);
    Aind = Bind = -1;
  }

  MinkVert(vec2 Ap, short Ai, vec2 Bp, short Bi) {
    pos = Ap - Bp;
    Aind = Ai;
    Bind = Bi;
  }

  bool operator==(MinkVert& other) {
    return Aind == other.Aind && Bind == other.Bind;
  }
};

struct ContactMesh {
  vec2 pos[2];
  vec2 nhat[2];
  char sz;

  ContactMesh() {
    sz = 0;
  }

  ContactMesh(vec2 p, vec2 n) {
    pos[0] = p;
    nhat[0] = n;
    sz = 1;
  }

  ContactMesh(vec2 p1, vec2 p2, vec2 n1, vec2 n2) {
    pos[0] = p1;
    pos[1] = p2;
    nhat[0] = n1;
    nhat[1] = n2;
    sz = 2;
  }

  ContactMesh(vec2 p1, vec2 p2, vec2 n1, vec2 n2, char ct) {
    pos[0] = p1;
    pos[1] = p2;
    nhat[0] = n1;
    nhat[1] = n2;
    sz = ct;
  }

  ContactMesh operator-() {
    return ContactMesh(pos[0], pos[1], -nhat[0], -nhat[1], sz);
  }
};

// outputs the lerped coordinates such that u*A + v*B is
//   the point between A and B which is nearest the origin
// if u or v is negative, the point is not "between" A and B
void edgeUV(const vec2 A, const vec2 B, float& u, float& v) {
  vec2 d = B - A;
  float normz = 1/d.dot(d);
  u = normz * B.dot(d);
  v = -normz * A.dot(d);
}

// returns the point which is farthest in the direction of proj
//   (whose projection along the vector is greatest)
MinkVert SupportPoint(MinkVert* verts, int nverts, vec2 proj) {
  MinkVert maxvert = verts[0];
  float maxdist = proj.dot(maxvert.pos);

  for (int i=1; i < nverts; i++) {
    float d = proj.dot(verts[i].pos);
    if (d > maxdist) {
      maxdist = d;
      maxvert = verts[i];
    }
  }

  return maxvert;
}

// returns the index of the vertex which is farthest in the direction
//   of proj (whose projection along the vector is greatest)
int SupportPoint(vec2* verts, int nverts, vec2 proj) {
  int maxvert = 0;
  float maxdist = proj.dot(verts[0]);

  for (int i=1; i < nverts; i++) {
    float d = proj.dot(verts[i]);
    if (d > maxdist) {
      maxdist = d;
      maxvert = i;
    }
  }

  return maxvert;
}


// represents a point, segment, or triangle
// used by GJK to determine collisions
class Simplex {
public:
  MinkVert verts[3];
  char nverts;

  Simplex(const MinkVert v) {
    initialize(v);
  }

  Simplex() {
    nverts = 0;
  }


  void initialize(const MinkVert v) {
    verts[0] = v;
    nverts = 1;
  }

  // adds a new simplex point at the end
  void insertPoint(MinkVert v) {
    if (nverts >= 3)
      return;
    verts[nverts] = v;
    nverts++;
  }

  // removes the given simplex point, shifting the others down
  void removePoint(char index) {
    if (nverts <= 0)
      return;
    for (int i=index; i < nverts-1; i++)
      verts[i] = verts[i+1];
    nverts--;
  }

  // keeps only the provided simplex vertex
  void reduceTo(char index) {
    verts[0] = verts[index];
    nverts = 1;
  }

  bool contains(MinkVert& elem) {
    for (int i=0; i < nverts; i++) {
      if (verts[i] == elem)
        return true;
    }
    return false;
  }

  // returns the point enclosed by S which is closest to the origin
  // reduces the simplex to only contributing points before returning
  // reference: Erin Catto, "Computing Distance", GDC 2010
  vec2 pointNearestOrigin() {
    vec2 O(0,0);
    vec2 a, d, A, B, C;
    float normz;
    float u, v, w;
    float uAB, vAB, uBC, vBC, uCA, vCA;

    switch(nverts) {
      case 1:
        return verts[0].pos;

      case 2:
        A = verts[0].pos; B = verts[1].pos;
        edgeUV(A, B, u, v);
        if (u <= 0) {
          removePoint(0);
          return B;
        }
        if (v <= 0) {
          removePoint(1);
          return A;
        }
        return u*A + v*B; // lerp

      case 3:
        A = verts[0].pos; B = verts[1].pos; C = verts[2].pos;
        a = B-A;
        d = C-A;
        normz = 1 / (a[0]*d[1] - a[1]*d[0]); // total double-area
        u = normz * (B[0]*C[1] - B[1]*C[0]); // fractional areas
        v = normz * (C[0]*A[1] - C[1]*A[0]);
        w = normz * (A[0]*B[1] - A[1]*B[0]);

        // test if origin enclosed
        if (u > 0 && v > 0 && w > 0)
          return O;

        // compute and test edge regions
        edgeUV(A, B, uAB, vAB);
        if (uAB > 0 && vAB > 0 && w <= 0) {
          removePoint(2);
          return uAB*A + vAB*B;
        }
        edgeUV(B, C, uBC, vBC);
        if (uBC > 0 && vBC > 0 && u <= 0) {
          removePoint(0);
          return uBC*B + vBC*C;
        }
        edgeUV(C, A, uCA, vCA);
        if (uCA > 0 && vCA > 0 && v <= 0) {
          removePoint(1);
          return uCA*C + vCA*A;
        }

        // test vertex regions
        if (uCA <= 0 && vAB <= 0) {
          reduceTo(0);
          return A;
        }
        if (uAB <= 0 && vBC <= 0) {
          reduceTo(1);
          return B;
        }
        if (uBC <= 0 && vCA <= 0) {
          reduceTo(2);
          return C;
        }

        fprintf(stderr, "ERROR: nearestPoint (triangle) failed to return\n");
        return O;

      default:
        fprintf(stderr, "ERROR: nearestPoint given simplex of size %d\n", nverts);
        return O;
    }
  }
};


// takes the nearest features on two rigidbodies, and performs root-finding on the
//   point where a separating axis is crossed in order to perform conservative advancement
class RootFinder {
private:
  enum AxisType { // stores which features define the separating edge
    AxisA,
    AxisB,
    AxisVert,
    AxisDouble
  };

public:
  Rigidbody *A, *B;
  int Aind, Bind; // indices to the points on A and B which we are separating
  int Aind2, Bind2; // used for edge-edge resolution
  int nhat; // used for edge-vertex resolution (refers to index of edge normal)
  vec2 axis; // used for vertex-vertex resolution (no edge normal to reference)
  AxisType type;
  float t1, t2, t;
  float tolerance;

  RootFinder() {}

  // takes tstart in [0, 1) as the fraction of the timestep to begin searching from
  // returns the time, [tstart, tstep], at which the bodies cross the separating axis within distance tol
  RootFinder(Rigidbody* rbA, Rigidbody* rbB, Simplex* S, float tstart, float tstep, float tol)
  : A(rbA), B(rbB), t1(tstart), tolerance(tol) {
    MinkVert a, b;

    if (S->nverts == 1) {
      a = S->verts[0];
      type = AxisVert;

      // Aind, Bind are support points (vertices with deepest penetration) at end time
      axis = B->worldVertLerp(a.Bind, t1) - A->worldVertLerp(a.Aind, t1);
      vec2 Afin[A->nverts], Bfin[B->nverts];
      A->worldCoordsLerp(Afin, tstep);
      B->worldCoordsLerp(Bfin, tstep);
      Aind = SupportPoint(Afin, A->nverts, axis);
      Bind = SupportPoint(Bfin, B->nverts, -axis);

    } else if (S->nverts == 2) {
      a = S->verts[0];
      b = S->verts[1];
      Aind = a.Aind;
      Bind = b.Bind;

      // take the normal vector from the separating edge
      if (a.Aind == b.Aind) {
        type = AxisB;
        if ((a.Bind == 0 && b.Bind != 1) || (b.Bind == 0 && a.Bind != 1)) // handle wraparound
          nhat = max(a.Bind, b.Bind);
        else
          nhat = min(a.Bind, b.Bind);
      } else if (a.Bind == b.Bind) {
        type = AxisA;
        if ((a.Aind == 0 && b.Aind != 1) || (b.Aind == 0 && a.Aind != 1)) // handle wraparound
          nhat = max(a.Aind, b.Aind);
        else
          nhat = min(a.Aind, b.Aind);
      } else {
        type = AxisDouble;
        // the two edges must be parallel, so we arbitrarily choose edge A normal
        if ((a.Aind == 0 && b.Aind != 1) || (b.Aind == 0 && a.Aind != 1)) // handle wraparound
          nhat = max(a.Aind, b.Aind);
        else
          nhat = min(a.Aind, b.Aind);
        Aind2 = b.Aind; // track these to potentially use for contact generation
        Bind2 = a.Bind;
      }

    } else { // overlapping -> no separating axis
      fprintf(stderr, "ERROR: RootFinder give bad simplex size %d\n", S->nverts);
      Aind = Bind = nhat = -1;
    }

    t = t1;
    t2 = tstep;
  }

  // returns the time, [t1, t2), where the reference features cross their separating axis
  // returns t < 0 if the objects are not penetrating at the end of the timestep
  float FindRoot() {
    float d;
    vec2 n;

    switch(type) {
      case AxisA:
      case AxisDouble:
        // if the root is not bracketed, exit
        n = A->worldNormLerp(nhat, t2);
        d = n.dot(B->worldVertLerp(Bind, t2)) - n.dot(A->worldVertLerp(Aind, t2));
        if (d > 0)
          return -1;
        n = A->worldNormLerp(nhat, t1);
        d = n.dot(B->worldVertLerp(Bind, t1)) - n.dot(A->worldVertLerp(Aind, t1));
        if (d < tolerance)
          return t1;

        while (true) {
          n = A->worldNormLerp(nhat, t);
          d = n.dot(B->worldVertLerp(Bind, t)) - n.dot(A->worldVertLerp(Aind, t));

          if (abs(d) < tolerance)
            return t;

          // bisection root-finding
          if (d < 0)
            t2 = t;
          else
            t1 = t;
          t = 0.5 * (t1+t2);
        }

      case AxisB:
        // if the root is not bracketed, exit
        n = B->worldNormLerp(nhat, t2);
        d = n.dot(A->worldVertLerp(Aind, t2)) - n.dot(B->worldVertLerp(Bind, t2));
        if (d > 0)
          return -1;
        n = B->worldNormLerp(nhat, t1);
        d = n.dot(A->worldVertLerp(Aind, t1)) - n.dot(B->worldVertLerp(Bind, t1));
        if (d < tolerance)
          return t1;

        while (true) {
          n = B->worldNormLerp(nhat, t);
          d = n.dot(A->worldVertLerp(Aind, t)) - n.dot(B->worldVertLerp(Bind, t));

          if (abs(d) < tolerance)
            return t;

          // bisection root-finding
          if (d < 0)
            t2 = t;
          else
            t1 = t;
          t = 0.5 * (t1+t2);
        }
        break;

      case AxisVert:
        // if the root is not bracketed, exit
        d = axis.dot(B->worldVertLerp(Bind, t2)) - axis.dot(A->worldVertLerp(Aind, t2));
        if (d > 0)
          return -1;
        d = axis.dot(B->worldVertLerp(Bind, t1)) - axis.dot(A->worldVertLerp(Aind, t1));
        if (d < tolerance)
          return t1;

        while (true) {
          d = axis.dot(B->worldVertLerp(Bind, t)) - axis.dot(A->worldVertLerp(Aind, t));

          if (abs(d) < tolerance)
            return t;

          // bisection root-finding
          if (d < 0)
            t2 = t;
          else
            t1 = t;
          t = 0.5 * (t1+t2);
        }
        break;

      default:
        fprintf(stderr, "ERROR: FindRoot given bad type enum %d\n", type);
        return -1;
    }
  }

  // returns a struct containing the world coordinates of the nearest features of the two rigidbodies
  // nhat will always point from A to B
  ContactMesh GetContacts(float dt) {
    switch(type) {
      case AxisA:
      case AxisDouble:
        return ContactMesh(B->worldVertLerp(Bind, dt), A->worldNormLerp(nhat, dt));
      case AxisB:
        return ContactMesh(A->worldVertLerp(Aind, dt), -B->worldNormLerp(nhat, dt));
      case AxisVert: // we assume Aind position is approximately Bind position
        return ContactMesh(A->worldVertLerp(Aind, dt), axis);

      default:
        return ContactMesh();
    }
  }
};


// returns the "time", between 0 and 1, at which two spheres A and B (with radius R) collide,
//   assuming that they move at constant velocity between their initial and final positions
// returns t < 0 if there is no collision between time 0 and time 1
// reference: https://www.toptal.com/game/video-game-physics-part-ii-collision-detection-for-solid-objects
//   (and a bunch of algebra)
float SphereCast(vec2 A1, vec2 A2, float Ra, vec2 B1, vec2 B2, float Rb) {
  // set up helper variables
  float Rsq = Ra+Rb; Rsq *= Rsq;
  vec2 BmA = B1 - A1;
  float BmAsq = BmA.dot(BmA);
  vec2 D = A2 - A1 - B2 + B1;
  float Dsq = D.dot(D);
  float t;

  // apply the quadratic formula
  float b = 2*BmA.dot(D); // actually -b
  float discr = b*b - 4*Dsq*(BmAsq - Rsq);
  if (discr < 0)
    return -1;
  if (discr < EPSILON)
    t = b/(2*Dsq);
  else // take the first root (collision entry; second is collision exit)
    t = (b - sqrt(discr)) / (2*Dsq);
  return (t <= 1) ? t : -1;
}

// wrapper for SphereCast
float SphereOverlap(Dynamic* A, Dynamic* B, float tstep) {
  vec2 A1 = A->position;
  vec2 A2 = A1 + tstep * A->velocity;
  vec2 B1 = B->position;
  vec2 B2 = B1 + tstep * B->velocity;
  return SphereCast(A1, A2, A->radius, B1, B2, B->radius);
}

float SphereOverlap(Dynamic* A, Static* B, float tstep) {
  vec2 A1 = A->position;
  vec2 A2 = A1 + tstep * A->velocity;
  vec2 B1 = B->position;
  return SphereCast(A1, A2, A->radius, B1, B1, B->radius);
}



// returns the nearest square distance between minkowski difference and origin
// out-parameter "nearest" is the position itself
// in/out-parameter simplex gets updated; stores points contributing to nearest
// reference: so many, but mostly Erin Catto GDC 2013
float GJK(MinkVert* minkverts, int nverts, Simplex* S, vec2& nearest, float sqtolerance) {
  while (true) {
    // return if the minkowski difference is within threshold of the origin
    nearest = S->pointNearestOrigin();
    float sqdist = nearest.dot(nearest);
      if (sqdist < sqtolerance) {
      nearest = vec2(0,0);
      return 0;
    }

    // return if this iteration didn't bring us any closer to the origin
    MinkVert support = SupportPoint(minkverts, nverts, -nearest);
    if (S->contains(support))
      return sqdist;
    S->insertPoint(support);
  }
}



// performs root-finding on the distance between A and B, approximated as constant-velocity movers
// takes bracketing start and end times to search for collision, in the range of [0, tstep)
// returns the time [0,tstep] at which collision first occurs
// returns a negative value if the objects are not colliding at t=tstep
// out-parameter cmesh is a list of world-space contact points between the two objects
float TimeOfImpact(Rigidbody* A, Rigidbody* B, float t1, float tstep, ContactMesh& cmesh, float threshold=0.00001) {
  float sqtolerance = threshold*threshold;
  float t = t1;

  // set up helper objects
  Simplex* S = new Simplex();
  int nMVs = A->nverts * B->nverts;
  MinkVert* minkDiff = (MinkVert*)malloc(nMVs * sizeof(MinkVert*));
  RootFinder advancement;

  // stop if we detect non-collision, or if we advance beyond the end of the current timestep
  while (t >= 0 && t < tstep) {
    // compute the minkowski difference at the current time
    // critical hypothesis: since GJK relies on support points, computing the convex hull on
    //   the minkowski difference is not necessary for algorithmic accuracy, just for speed
    for (int i=0; i < A->nverts; i++) {
      for (int j=0; j < B->nverts; j++)
        minkDiff[B->nverts*i + j] = MinkVert(A->worldVertLerp(i, t), i, B->worldVertLerp(j, t), j);
    }

    // initialize the simplex, and run GJK to compute nearest points
    S->initialize(minkDiff[0]);
    vec2 nearest;
    float sqdist = GJK(minkDiff, nMVs, S, nearest, sqtolerance);
    if (sqdist < sqtolerance)
      break; // current distance (and time) is sufficiently close

    // perform conservative advancement
    advancement = RootFinder(A, B, S, t, tstep, threshold);
    t = advancement.FindRoot();
  }

  delete S;
  free(minkDiff);

  // return -1 if we advanced beyond the end of the timestep
  if (t > tstep)
    return -1;

  cmesh = advancement.GetContacts(t); // store the contact mesh
  return t;
}


class ContinuousSim : public Simulator {
  public:
  int physLoops;

  ContinuousSim(Dynamic** phys, const int nphysObjs, Static** statics, const int nstatObjs, const int loops=5,
                const double t=0.01, const int n=1, const vec2& g=vec2(0,0), const float e=1, const SimType type=Sim_Timer)
  :Simulator(phys, nphysObjs, statics, nstatObjs, t, n, g, e, type), physLoops(loops) {}


  void staticCollision(Dynamic* d, Static* s, ContactMesh c) {
    return;
  }

  void dynamicCollision(Dynamic* a, Dynamic* b, ContactMesh c) {
    return;
  }

  // update the position of every dynamic object
  void integratePositions(float dt) {
    for (int i=0; i < nphys; i++) {
      Dynamic* thisphys = physObjs[i];
      thisphys->position += thisphys->velocity * dt;
      if (thisphys->rotspeed != 0) {
        thisphys->rotation += thisphys->rotspeed * dt;
        thisphys->updateAABB();
      }
    }
  }

  void simulateTimestep() {
    float remaining = tstep; // amount of time left to integrate in this timestep
    bool collision = false; // was a collision detected?
    float tcol; // time of earliest collision

    // apply forces to each object
    for (int i=0; i < nphys; i++) {
      Dynamic* thisphys = physObjs[i];
      thisphys->force += globalForce;
      thisphys->velocity += thisphys->force * tstep * thisphys->invmass;
      thisphys->rotspeed += thisphys->torque * tstep * thisphys->invinertia;

      thisphys->force = vec2(0,0); // flush the accumulated forces
      thisphys->torque = 0;
    }

    if (nphys == 1) {
      Dynamic* thisphys = physObjs[0];
      Static* B; // first object collided with
      ContactMesh cmesh; // contact mesh for earliest collision

      while (remaining > 0) {
        collision = false;
        tcol = remaining+1;

        // accumulate earliest impact
        for (int j=0; j < nstatics; j++) {
          Static* other = staticObjs[j];
          float t = SphereOverlap(thisphys, other, remaining); // test possible collision
          if (t >= 0) {
            ContactMesh tmpmesh;
            t = TimeOfImpact(thisphys, other, t, remaining, tmpmesh); // test actual collision
            if (t >= 0 && t < tcol) {
              collision = true;
              cmesh = tmpmesh;
              tcol = t;
              B = other;
            }
          }
        }

        if (collision) {
          integratePositions(tcol);
          staticCollision(thisphys, B, cmesh);
          remaining -= tcol;
        } else {
          integratePositions(remaining);
          remaining = -1;
        }
      }
    } else { // nphys > 1

    }
    //while (remaining > 0) {
      // iterate through objects to detect the first collision (accumulate objects and time of impact)
      // handle that collision
      // integrate every object to that time tC
      // remaining -= tC


      // if (nphys == 1) {
      //   Dynamic* thisphys = physObjs[0];
      //   for (int j=0; j < nstatics; j++) {
      //     Static* other = staticObjs[j];
      //     if (SphereOverlap(thisphys, other, remaining))
      //       staticCollision(thisphys, other);
      //   }
      // } else { // nphys > 1
      //   for (int i=0; i < nphys-1; i++) {
      //     Dynamic* thisphys = physObjs[i];
      //     for (int j=i+1; j < nphys; j++) { // all dynamic-dynamic
      //       Dynamic* other = physObjs[j];
      //       if (SphereOverlap(thisphys, other, remaining))
      //         dynamicCollision(thisphys, other);
      //     }

      //     for (int j=0; j < nstatics; j++) { // all but one dynamic-static
      //       Static* other = staticObjs[j];
      //       if (SphereOverlap(thisphys, other, remaining))
      //         staticCollision(thisphys, other);
      //     }
      //   }

      //   Dynamic* thisphys = physObjs[nphys-1];
      //   for (int j=0; j < nstatics; j++) { // final dynamic-static
      //     Static* other = staticObjs[j];
      //     if (SphereOverlap(thisphys, other, remaining))
      //       staticCollision(thisphys, other);
      //   }
      // }
    //} // end physloop
  }

};

#endif //_CONTINUOUS_