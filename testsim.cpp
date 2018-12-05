#include <string>
#include <unistd.h>
#include "stdio.h"

#include "simulate.h"
#include "rigidbody.h"
#include "simulators.h"
#include "continuous.h"

int main(int argc, char* argv[]) {
  int scene = 0;
  int nsteps = 1000;
  int physloops = 5;
  bool continuous = false;
  SimType simt = Sim_Timer;
  float dt = 0.01;
  float elasticity = 1.0;
  char* fname;
  FILE *outSIM, *outSYS, *outSTA, *outDYN;

  int c;
  extern char *optarg;
  extern int optind, opterr, optopt;

  while ((c = getopt(argc, argv, "s:n:t:f:p:e:c")) != -1) {
    switch(c) {
    case 'c':
      continuous = true;
      break;
    case 'n':
      nsteps = std::stoi(optarg);
      break;
    case 't':
      dt = std::stof(optarg);
      break;
    case 's':
      scene = std::stoi(optarg);
      break;
    case 'p':
      physloops = std::stoi(optarg);
      break;
    case 'e':
      elasticity = std::stof(optarg);
      break;
    case 'f':
      simt = Sim_FullOutput;

      fname = (char*)malloc(strlen(optarg) + 4+1);
      strcpy(fname, optarg);
      strcat(fname, ".sim");
      outSIM = fopen(fname, "w");

      strcpy(fname, optarg);
      strcat(fname, ".sys");
      outSYS = fopen(fname, "w");

      strcpy(fname, optarg);
      strcat(fname, ".sta");
      outSTA = fopen(fname, "w");

      strcpy(fname, optarg);
      strcat(fname, ".dyn");
      outDYN = fopen(fname, "w");
      free(fname);
      break;
    case ':':
    case '?':
    default:
      fprintf(stderr, "Usage: [-c] [-s test_scenario] [-n nsteps] [-t timestep] [-f outfile_prefix] [-p physloops] [-e coeff_of_restitution]\n");
      return 1;
    }
  }

  int nphys=0, nstat=0;
  vec2 grav(0, -10);
  Dynamic* physs[11];
  Static* stats[10];
  const float pi = 3.141592653589793f;
  const float halfpi = 1.5707963267948966f;

  // reusable shapes
  vec2 boxverts[] = {vec2(0,0), vec2(1,0), vec2(1,1), vec2(0,1)};
  vec2 floorverts[] = {vec2(0,0), vec2(20, 0), vec2(20,1), vec2(0,1)};
  vec2 thinfloorverts[] = {vec2(0,0), vec2(20, 0), vec2(20,0.1), vec2(0,0.1)};
  vec2 hexadecverts[16]; float hexstep = pi/8;
  for (int i=0; i<16; i++) {
    float theta = i * hexstep;
    hexadecverts[i] = vec2(cos(theta), sin(theta));
  }

  if (scene <= 1) { // box falling onto floor
    Dynamic* box = new Dynamic(boxverts, 4, vec2(0, 5));
    Static* floor = new Static(floorverts, 4, vec2(0, -0.5));

    nphys = 1; nstat = 1;
    physs[0] = box;
    stats[0] = floor;

  } else if (scene == 2) { // two boxes hitting each other; two walls
    Dynamic* boxL = new Dynamic(boxverts, 4, vec2(-5, 0), 0,0, vec2(3,  0));
    Dynamic* boxR = new Dynamic(boxverts, 4, vec2(5,  0), 0,0, vec2(-3, 0));
    Static* wallL = new Static(floorverts, 4, vec2(-10.5, 0), halfpi);
    Static* wallR = new Static(floorverts, 4, vec2(10.5, 0), halfpi);
    Static* wallB = new Static(floorverts, 4, vec2(0, -10.5));

    nphys = 2; nstat = 3;
    physs[0] = boxL; physs[1] = boxR;
    stats[0] = wallL; stats[1] = wallR; stats[2] = wallB;

  } else if (scene == 3) { // tower of bouncing boxes
    int nboxes = 3;
    for (int i=0; i < nboxes; i++)
      physs[i] = new Dynamic(boxverts, 4, vec2(0, 1 + 2*i));
    Static* floor = new Static(floorverts, 4, vec2(0, -0.5));

    nphys = nboxes; nstat = 1;
    stats[0] = floor;

  } else if (scene == 4) { // billiard squares, to test basic conservation of momentum
    float mH1 = 1;
    float mH2 = 1000;
    float mD = 1;
    Dynamic* boxHoriz1 = new Dynamic(boxverts, 4, vec2(0,0), 0,mH1, vec2(1, 0));
    Dynamic* boxDiag1 = new Dynamic(boxverts, 4, vec2(0,3), 0,mD, vec2(1, -1));

    Dynamic* boxHoriz2 = new Dynamic(boxverts, 4, vec2(0,5), 0,mH2, vec2(1, 0));
    Dynamic* boxDiag2 = new Dynamic(boxverts, 4, vec2(0,8), 0,mD, vec2(1, -1));

    nphys = 4; nstat = 0;
    physs[0] = boxHoriz1; physs[1] = boxDiag1; physs[2] = boxHoriz2; physs[3] = boxDiag2;
    grav = vec2(0,0);

  } else if (scene == 5) { // box approaching a corner
    Dynamic* box = new Dynamic(boxverts, 4, vec2(0, 0), 0,0, vec2(1,-1));
    Static* wallR = new Static(floorverts, 4, vec2(10.5, 0), halfpi);
    Static* wallB = new Static(floorverts, 4, vec2(0, -10.5));

    nphys = 1; nstat = 2;
    physs[0] = box;
    stats[0] = wallR; stats[1] = wallB;
    grav = vec2(0,0);

  } else if (scene == 6)  { // fast-moving particle (probable tunneling)
    Dynamic* box = new Dynamic(boxverts, 4, vec2(0, 0), 0,0, vec2(100,20), pi);
    Static* wallR = new Static(thinfloorverts, 4, vec2(10.05, 0), halfpi);
    Static* wallL = new Static(thinfloorverts, 4, vec2(-10.05, 0), halfpi);
    Static* wallT = new Static(thinfloorverts, 4, vec2(0, 10.05));
    Static* wallB = new Static(thinfloorverts, 4, vec2(0, -10.05));

    nphys = 1; nstat = 4;
    physs[0] = box;
    stats[0] = wallR; stats[1] = wallL; stats[2] = wallT; stats[3] = wallB;
    grav = vec2(0,0);

  } else if (scene == 7) { // tower of boxes at rest under gravity
    int nboxes = 10;
    for (int i=0; i < nboxes; i++)
      physs[i] = new Dynamic(boxverts, 4, vec2(0, 0.5+i));
    Static* floor = new Static(floorverts, 4, vec2(0, -0.5));

    nphys = nboxes; nstat = 1;
    stats[0] = floor;

  } else if (scene == 8) { // cue ball breaking a tower -> particles in a box
    float gridsz = 2.0;
    int k=0;
    for (int i=0; i < 4; i++) {
      float y = gridsz * i;
      for (int j=0; j <= i; j++) {
        float x = (2*j - i) * gridsz;
        physs[k++] = new Dynamic(hexadecverts, 16, vec2(x,y));
      }
    }
    Dynamic* cue = new Dynamic(hexadecverts, 16, vec2(0, -5), 0,0, vec2(0,5));

    Static* wallR = new Static(floorverts, 4, vec2(10.5, 0), halfpi);
    Static* wallL = new Static(floorverts, 4, vec2(-10.5, 0), halfpi);
    Static* wallT = new Static(floorverts, 4, vec2(0, 10.5));
    Static* wallB = new Static(floorverts, 4, vec2(0, -10.5));

    nphys = 11; nstat = 4;
    physs[10] = cue;
    stats[0] = wallR; stats[1] = wallL; stats[2] = wallT; stats[3] = wallB;
    grav = vec2(0,0);

  } else {
    return 0;
  }

  // run with full output
  Simulator* sim;
  if (continuous)
    sim = new ContinuousSim(physs, nphys, stats, nstat, dt, nsteps, grav, elasticity, simt);
  else
    sim = new EulerPairwise(physs, nphys, stats, nstat, physloops, dt, nsteps, grav, elasticity, simt);

  if (outSIM) {
    sim->runFullSimulation(outSIM, outSYS, outSTA, outDYN);
    fclose(outSIM);
    fclose(outSYS);
    fclose(outSTA);
    fclose(outDYN);
  } else {
    sim->runFullSimulation();
  }

  return 0;
}