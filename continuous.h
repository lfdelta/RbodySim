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

  ContactMesh(vec2 p1, vec2 n1, vec2 p2, vec2 n2) {
    pos[0] = p1;
    pos[1] = p2;
    nhat[0] = n1;
    nhat[1] = n2;
    sz = 2;
  }

  ContactMesh(vec2 p1, vec2 n1, vec2 p2, vec2 n2, char ct) {
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
// if dup is true, out-parameter duplicate stores a second vertex which is equally far
MinkVert SupportPoint(MinkVert* verts, int nverts, vec2 proj, bool& dup, MinkVert* duplicate) {
  MinkVert maxvert = verts[0];
  float maxdist = proj.dot(maxvert.pos);
  dup = false;

  for (int i=1; i < nverts; i++) {
    float d = proj.dot(verts[i].pos);
    if (d > maxdist) {
      maxdist = d;
      maxvert = verts[i];
      dup = false;
    } else if (d == maxdist) {
      *duplicate = verts[i];
      dup = true;
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

// calculates the vertices which are farthest in the direction of proj
// returns the number of support points found (1 or 2)
// out-parameter supports must be length 2, and stores the support vertices
int SupportPoints(vec2* verts, int nverts, vec2 proj, vec2* supports) {
  float maxdist = proj.dot(verts[0]);
  supports[0] = verts[0];
  int sz = 1;

  for (int i=1; i < nverts; i++) {
    float d = proj.dot(verts[i]);
    if (d > maxdist) {
      maxdist = d;
      supports[0] = verts[i];
      sz = 1;
    } else if (d == maxdist) {
      supports[1] = verts[i];
      sz = 2;
    }
  }

  return sz;
}


// uses separating axis theorem (SAT) to generate nearest points on two
//   just-barely-overlapping meshes
// assumes that this is NOT a vertex-vertex collision; those should be handled separately
ContactMesh EdgeContactPoints(vec2* Averts, vec2* Anorms, int Asz, vec2* Bverts, vec2* Bnorms, int Bsz) {
  vec2 MTVdirA, MTVdirB, MTVdir;
  float MTVdistA = -1;
  float MTVdistB = -1;
  float MTVdist;

  if (MinimumTranslationVector(Averts, Anorms, Asz, Bverts, Bsz, MTVdirA, MTVdistA)) {
    if (MinimumTranslationVector(Bverts, Bnorms, Bsz, Averts, Asz, MTVdirB, MTVdistB)) {
      if (MTVdistA < MTVdistB) {
        MTVdist = MTVdistA;
        MTVdir = -MTVdirA;
      } else {
        MTVdist = MTVdistB;
        MTVdir = MTVdirB;
      }
      vec2 Asupp[2], Bsupp[2];
      int Anum = SupportPoints(Averts, Asz, MTVdir, Asupp);
      int Bnum = SupportPoints(Bverts, Bsz, -MTVdir, Bsupp);
      if (Anum == 1) {
        return ContactMesh(Asupp[0], MTVdir); // B edge -> A vertex
      } else if (Bnum == 1) {
        return ContactMesh(Bsupp[0], MTVdir); // A edge -> B vertex

      } else {
        vec2 edge = Asupp[1] - Asupp[0]; // A edge -> B edge (should be roughly collinear)
        vec2 edgeverts[4] = {Asupp[0], Asupp[1], Bsupp[0], Bsupp[1]};
        int rightmost = SupportPoint(edgeverts, 4, edge);
        int leftmost = SupportPoint(edgeverts, 4, -edge);
        vec2 innerpoints[2], normals[2];
        int j = 0;
        for (int i=0; i < 2; i++) {
          if (i != leftmost && i != rightmost) { // find the two innermost (collinear) points
            innerpoints[j] = edgeverts[i];
            normals[j++] = MTVdir; // A vertex -> B edge
          }
        }
        for (int i=2; i < 4; i++) {
          if (i != leftmost && i != rightmost) {
            innerpoints[j] = edgeverts[i];
            normals[j++] = MTVdir; // B vertex -> A edge
          }
        }
        return ContactMesh(innerpoints[0], normals[0], innerpoints[1], normals[1]);
      }
    }
  }

  return ContactMesh(); // no collision
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
    AxisA = 0,
    AxisB = 1,
    AxisVert = 2,
    AxisDouble = 3
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

          if (d < 0 && abs(d) < tolerance)
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

          if (d < 0 && abs(d) < tolerance)
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

          if (d < 0 && abs(d) < tolerance)
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
    vec2 Avs[A->nverts], Ans[A->nverts];
    vec2 Bvs[B->nverts], Bns[B->nverts];

    switch(type) {
      case AxisVert: // we assume Aind position is approximately Bind position
        return ContactMesh(A->worldVertLerp(Aind, dt), axis.normalized());
      case AxisA:
      case AxisB:
      case AxisDouble:
        A->worldCoordsLerp(Avs, dt); A->worldNormalsLerp(Ans, dt);
        B->worldCoordsLerp(Bvs, dt); B->worldNormalsLerp(Bns, dt);
        return EdgeContactPoints(Avs, Ans, A->nverts, Bvs, Bns, B->nverts);

      default:
        fprintf(stderr, "ERROR: GetContacts given invalid type %d\n", type);
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
  return (t <= 1) ? fmax(0, t) : -1;
}

// wrapper for SphereCast
float SphereOverlap(Dynamic* A, Dynamic* B, float tstep) {
  vec2 A1 = A->position;
  vec2 A2 = A1 + tstep * A->velocity;
  vec2 B1 = B->position;
  vec2 B2 = B1 + tstep * B->velocity;
  return tstep * SphereCast(A1, A2, A->radius, B1, B2, B->radius);
}

float SphereOverlap(Dynamic* A, Static* B, float tstep) {
  vec2 A1 = A->position;
  vec2 A2 = A1 + tstep * A->velocity;
  vec2 B1 = B->position;
  return tstep * SphereCast(A1, A2, A->radius, B1, B1, B->radius);
}


// uses cross product to determine if a->b->c is counterclockwise
// returns positive if CCW, 0 if collinear, negative if CW
float Handedness(vec2 a, vec2 b, vec2 c) {
  return ((b[0]-a[0])*(c[1]-b[1]) - (b[1]-a[1])*(c[0]-b[0]));
}

// uses "gift wrapping" algorithm to compute the convex hull of a point cloud
// out-parameter hull stores the culled vertices, in counterclockwise order
// returns the number of vertices in the hull
int ConvexHull(MinkVert* cloud, int nverts, MinkVert* hull) {
  int startind = 0;
  float startx = cloud[0].pos[0];
  float starty = cloud[0].pos[1];

  // look for the leftmost point; if multiple, take bottommost
  for (int i=0; i < nverts; i++) {
    vec2 pos = cloud[i].pos;
    if (pos[0] < startx || (pos[0] == startx && pos[1] < starty)) {
      startind = i;
      startx = pos[0];
      starty = pos[1];
    }
  }

  int nextind = startind; // the one we're about to add
  int lastind; // the one we just added

  hull[0] = cloud[startind];
  int sz = 1;

  // keep adding hull points until we've wrapped back to the start
  while (true) {
    lastind = nextind;
    nextind = 0; // start searching through the point cloud again

    // find the next hull point, given the previous one
    for (int i=1; i < nverts; i++) {
      float h = Handedness(cloud[lastind].pos, cloud[i].pos, cloud[nextind].pos);
      if (h > 0) {
        nextind = i; // this point is "more clockwise" than the accumulator
      } else if (h == 0) {
        vec2 p1 = cloud[i].pos - cloud[lastind].pos;
        vec2 p2 = cloud[nextind].pos - cloud[lastind].pos;
        if (p1.dot(p1) > p2.dot(p2))
          nextind = i; // collinear -> take the further one (outermost)
      }
    }

    if (nextind == startind)
      break;
    else
      hull[sz++] = cloud[nextind]; // add to hull
  }

  return sz;
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
    if (sqdist < sqtolerance)
      return sqdist;

    // return if this iteration didn't bring us any closer to the origin
    MinkVert alt;
    bool dup;
    MinkVert support = SupportPoint(minkverts, nverts, -nearest, dup, &alt);
    if (S->contains(support)) {
      if (!dup || S->contains(alt))
        return sqdist;
      else
        S->insertPoint(alt);
    } else {
      S->insertPoint(support);
    }
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
  MinkVert* minkDiff = (MinkVert*)malloc(nMVs * sizeof(MinkVert));
  RootFinder advancement;

  //fprintf(stderr, "ToI called with t1=%f, tstep=%f\n", t, tstep);

  // stop if we detect non-collision, or if we advance beyond the end of the current timestep
  while (t >= 0 && t < tstep) {
    // compute the minkowski difference at the current time
    // critical hypothesis: since GJK relies on support points, computing the convex hull on
    //   the minkowski difference is not necessary for algorithmic accuracy, just for speed
    for (int i=0; i < A->nverts; i++) {
      for (int j=0; j < B->nverts; j++)
        minkDiff[B->nverts*i + j] = MinkVert(A->worldVertLerp(i, t), i, B->worldVertLerp(j, t), j);
    }

    // MinkVert minkHull[nMVs];
    // int hullsz = ConvexHull(minkDiff, nMVs, minkHull);
    MinkVert* minkHull = minkDiff;
    int hullsz = nMVs;

    // initialize the simplex, and run GJK to compute nearest points
    S->initialize(minkHull[0]);
    vec2 nearest;
    float sqdist = GJK(minkHull, hullsz, S, nearest, sqtolerance);
    fprintf(stderr, "sqdist: %f\n", sqdist);

    advancement = RootFinder(A, B, S, t, tstep, threshold);

    // if (sqdist < sqtolerance) {
    //   cmesh = advancement.GetContacts(t); // current distance (and time) is sufficiently close
    //   break;
    // } else {
    // // perform conservative advancement
    // t = advancement.FindRoot();
    // cmesh = advancement.GetContacts(t); // store the contact mesh
    // fprintf(stderr, "root found at %f\n", t);
    // }

    // wait for them to just barely penetrate
    if (sqdist == 0) {
      cmesh = advancement.GetContacts(t);
      break;
    } else {
      // perform conservative advancement
      t = advancement.FindRoot();
      cmesh = advancement.GetContacts(t); // store the contact mesh
      fprintf(stderr, "root found at %f\n", t);
    }
  }

  delete S;
  free(minkDiff);

  // return -1 if we advanced beyond the end of the timestep
  if (t > tstep)
    return -1;

  return t;
}

float cross(const vec2 a, const vec2 b) {
  return a[0]*b[1] - a[1]*b[0];
}


class ContinuousSim : public Simulator {
  public:
  int physLoops;

  ContinuousSim(Dynamic** phys, const int nphysObjs, Static** statics, const int nstatObjs, const int loops=5,
                const double t=0.01, const int n=1, const vec2& g=vec2(0,0), const float e=1, const SimType type=Sim_Timer)
  :Simulator(phys, nphysObjs, statics, nstatObjs, t, n, g, e, type), physLoops(loops) {}


  // reference: Realtime Physics (John Dingliana, presented by Michael Manzke)
  // reference: Erin Catto, GDC 2014
  void staticCollision(Dynamic* d, Static* s, ContactMesh c) {
    // vec2 dp = vec2(0,0);
    // float dw = 0;

    // for (int i=0; i < c.sz; i++) {
    //   //   fprintf(stderr, "nhat {%0.1f, %0.1f}\n", nhat[0], nhat[1]);
    //   //   fprintf(stderr, "velocity: {%f, %f}\n", d->velocity[0], d->velocity[1]);
    //   //   fprintf(stderr, "perpvel: %f\n", perpvel);
    //   vec2 nhat = -c.nhat[i]; // points from static to dynamic
    //   vec2 r = c.pos[i] - d->position;
    //   vec2 wxr = d->rotspeed * vec2(-r[1], r[0]); // omega cross r
    //   vec2 vrel = nhat * nhat.dot(d->velocity + wxr);
    //   vec2 changeP = (1 + elasticity) * d->mass * vrel;
    //   dp += changeP;
    //   dw += cross(r, changeP) * d->invinertia;
    // }

    // if (c.sz == 2)
    //   dp *= 0.5;

    // d->velocity -= dp * d->invmass;
    // d->rotspeed -= dw;

    vec2 dp = vec2(0,0);
    float dalpha = 0;

    for (int i=0; i < c.sz; i++) {
      vec2 nhat = -c.nhat[i]; // points from static to dynamic
      fprintf(stderr, "nhat {%0.1f, %0.1f}\n", nhat[0], nhat[1]);
      vec2 r = c.pos[i] - d->position;
      vec2 wxr = d->rotspeed * vec2(-r[1], r[0]); // omega cross r
      float vrel = nhat.dot(d->velocity + wxr);
      float rxn = cross(r, nhat);
      float j = -(1 + elasticity) * vrel / (d->invmass + d->invinertia*rxn*rxn);

      vec2 J = nhat * j;
      dp += J;
      dalpha += cross(r, J);
    }

    // if (c.sz == 2)
    //   dp *= 0.5;

    d->velocity += dp * d->invmass;
    d->rotspeed += dalpha * d->invinertia;
  }

  void dynamicCollision(Dynamic* A, Dynamic* B, ContactMesh c) {
    // vec2 dpA = vec2(0,0), dpB = vec2(0,0);
    // float dwA = 0, dwB = 0;

    // for (int i=0; i < c.sz; i++) {
    //   vec2 nhat = c.nhat[i]; // points from A to B
    //   vec2 Ar = c.pos[i] - A->position;
    //   vec2 Awxr = A->rotspeed * vec2(-Ar[1], Ar[0]); // omega cross r
    //   vec2 Br = c.pos[i] - B->position;
    //   vec2 Bwxr = B->rotspeed * vec2(-Br[1], Br[0]);

    //   vec2 vrel = nhat * nhat.dot(A->velocity + Awxr - B->velocity - Bwxr);
    //   vec2 changePA = (1 + elasticity) * A->mass * vrel;
    //   dpA += changePA;
    //   dw += cross(r, changeP) * d->invinertia;
    // }

    // if (c.sz == 2)
    //   dp *= 0.5;

    // d->velocity -= dp * d->invmass;
    // d->rotspeed -= dw;
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
              fprintf(stderr, "collision at time %f\n", t);
              collision = true;
              cmesh = tmpmesh;
              tcol = t;
              B = other;
            }
          }
        }

        if (collision) {
          fprintf(stderr, "COLLISION!!!!!!\n");
          integratePositions(tcol);
          staticCollision(thisphys, B, cmesh);
          remaining -= tcol;
          integratePositions(remaining); // hopefully a temporary measure
          remaining = -1;
        } else {
          integratePositions(remaining);
          remaining = -1;
        }
      }
    } else { // nphys > 1

    }
  }

};

#endif //_CONTINUOUS_