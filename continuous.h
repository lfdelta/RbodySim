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

// outputs the lerped coordinates such that u*A + v*B is
//   the point between A and B which is nearest the origin
// if u or v is negative, the point is not "between" A and B
void edgeUV(const vec2 A, const vec2 B, float& u, float& v) {
  vec2 d = B - A;
  float normz = 1/d.dot(d);
  u = normz * B.dot(d);
  v = -normz * A.dot(d);
}

// represents a point, segment, or triangle
// used by GJK to determine collisions
class Simplex {
public:
  MinkVert verts[3];
  char nverts;

  Simplex(const MinkVert v) {
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
  // which features define the separating edge?
  enum AxisType {
    AxisA,
    AxisB,
    AxisBoth
  };

public:
  Rigidbody *A, *B;
  int Aind, Bind; // indices to the points on A and B which we are separating
  int edge[2]; // indices to the edge defining separating axis
  AxisType type;
  float t1, t2, t, dt;

  RootFinder(Rigidbody* rbA, Rigidbody* rbB, Simplex* S, float tstart, float tend, float tstep)
  : A(rbA), B(rbB), t1(tstart), t2(tend), dt(tstep) {
    MinkVert a, b;
    vec2 axis;
    switch (S->nverts) {
      case 1:
        a = S->verts[0];
        Aind = a.Aind;
        Bind = b.Bind;
        edge[1] = Bind;
        edge[0] = Aind;
        type = AxisBoth;
        break;

      case 2:
        a = S->verts[0];
        b = S->verts[1];
        Aind = a.Aind;
        Bind = a.Bind;
        if (a.Aind == b.Aind) {
          edge[1] = a.Bind;
          edge[0] = b.Bind;
          type = AxisB;
          axis = B->vertices[edge[1]] - B->vertices[edge[0]]; // parallel to edge
          axis = vec2(-axis[1], axis[0]); // perpendicular to edge
          if (axis.dot(B->vertices[Bind]) > axis.dot(A->vertices[Aind])) {
            edge[0] = a.Bind; // ensure axis points from A to B
            edge[1] = b.Bind;
          }
        } else if (a.Bind == b.Bind) {
          edge[1] = a.Aind;
          edge[0] = b.Aind;
          type = AxisA;
          axis = A->vertices[edge[1]] - A->vertices[edge[0]];
        }

        break;

      case 3:
      default:
        Aind = Bind = -1;
        //nhat = vec2(0,0);
        break; // overlapping -> no separating axis
    }

    t = t1;
  }
};


// returns the "time", between 0 and 1, at which two spheres A and B (with radius R)
//   collide, assuming that they move at constant velocity between their initial
//   and final positions
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


// returns the nearest square distance between minkowski difference and origin
// out-parameter "nearest" is the position itself
// in/out-parameter simplex gets updated; stores points contributing to nearest
// reference: so many, but mostly Erin Catto GDC 2013
float GJK(MinkVert* minkverts, int nverts, Simplex* S, vec2& nearest, float sqtolerance) {
  while (true) {
    nearest = S->pointNearestOrigin();
    float sqdist = nearest.dot(nearest);
      if (sqdist < sqtolerance) {
      nearest = vec2(0,0);
      return 0;
    }

    MinkVert support = SupportPoint(minkverts, nverts, -nearest);
    if (S->contains(support))
      return sqdist;
    S->insertPoint(support);
  }
}


// uses the minkowksi vertices to determine which features on the corresponding
//   rigidbodies are nearest one another, and returns the normal vector
//   of a separating axis, pointing from A to B
// vec2 ComputeSeparatingAxis(Rigidbody* A, Rigidbody* B, Simplex* S) {
  // int Ainds[2], Binds[2];
  // MinkVert a, b;
  // switch (S->nverts) {
  //   case 1:
  //     a = S->verts[0];
  //     return A->vertices[a.Aind] - B->vertices[a.Bind];
  //   case 2:
  //     a = S->verts[0];
  //     b = S->verts[1];
  //     if (a.Aind == b.Aind)
  //       return
  //     return;
  //   case 3:
  //   default:
  //     return vec2(0,0); // overlapping -> no separating axis
  // }
// }


// performs root-finding on the distance between A and B, approximated as constant-velocity movers
// takes bracketing start and end times to search for collision, in the range of [0, 1)
// returns the time [0,1] at which collision first occurs
float TimeOfImpactGJK(Rigidbody* A, Rigidbody* B, float t1, float t2=1, float threshold=0.00001) {
  t1 = fmin(0, t1 - 0.1); // backtrack to slightly before the earliest possible penetration point
  t2 = fmin(1, t2); // only advance to the end of the timestep
  float t = 0.5 * (t1+t2); // current probe time
  float sqtolerance = threshold*threshold;

  // compute the minkowski difference
  // critical hypothesis: since GJK relies on support points, computing the
  //   convex hull on the minkowski difference is not necessary for algorithmic accuracy


  //// ***** need to compute minkDiff using current positions based upon partial time integration
  int nMVs = A->nverts * B->nverts;
  MinkVert* minkDiff = (MinkVert*)malloc(nMVs * sizeof(MinkVert*));
  for (int i=0; i < A->nverts; i++) {
    for (int j=0; j < B->nverts; j++)
      minkDiff[B->nverts*i + j] = MinkVert(A->vertices[i], i, B->vertices[j], j);
  }

  // initialize the simplex, and run GJK to compute nearest points
  Simplex* S = new Simplex(minkDiff[0]);
  vec2 nearest;
  float sqdist = GJK(minkDiff, nMVs, S, nearest, sqtolerance);
  float lastSqdist = sqdist;

  if (sqdist == 0) {
    t2 = t;
    t = 0.5 * (t1+t2);
    /// ................................... reverse time and repeat GJK
  } else if (sqdist < sqtolerance) {
    // ...................... compute nearest features (and collision point)
    return t; // IMPACT !
  } else {
    t1 = t;
    t = 0.5 * (t1+t2);
    // ............................. advance time and repeat
  }

  // extract nearest features, and root-find separating axis crossing

  free(minkDiff);


  // bilateral advancement:
  // if (p > oldP)
  //   return -1 ??
  // if (p == O)
  //   t2 = t
  //   t = (t1 + t2)/2
  //   repeat GJK
  // if (p dot p < threshold)
  //   return t
  // else
  //   extract nearest features from simplex -> separating axis
  //   conservative advancement until axis crossed
  //   repeat GJK
  return 0;
}


class ContinuousSim : public Simulator {
  public:
  ContinuousSim(Dynamic** phys, const int nphysObjs, Static** statics, const int nstatObjs,
                const double t=0.01, const int n=1, const vec2& g=vec2(0,0), const float e=1, const SimType type=Sim_Timer)
  :Simulator(phys, nphysObjs, statics, nstatObjs, t, n, g, e, type) {}

  void simulateTimestep() {
    return;
  }
};

#endif //_CONTINUOUS_