#include <string>
#include <unistd.h>
#include "stdio.h"

#include "simulate.h"
#include "rigidbody.h"
#include "simulators.h"

int main(int argc, char* argv[]) {
  int scene = 0;
  int nsteps = 100;
  int physloops = 5;
  float dt = 0.1;
  char* fname;
  FILE *outSIM, *outSYS, *outSTA, *outDYN;

  int c;
  extern char *optarg;
  extern int optind, opterr, optopt;

  while ((c = getopt(argc, argv, "s:n:t:f:p:")) != -1) {
    switch(c) {
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
    case 'f':
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
      fprintf(stderr, "Usage: [-n nsteps] [-t timestep] [-s scenario] [-f outfile_prefix]\n");
      return 1;
    }
  }

  int nphys=0, nstat=0;
  vec2 grav(0, -10);
  Dynamic* physs[10];
  Rigidbody* stats[10];

  // reusable shapes
  vec2 boxverts[] = {vec2(0,0), vec2(1,0), vec2(1,1), vec2(0,1)};
  vec2 floorverts[] = {vec2(0,0), vec2(20, 0), vec2(20,1), vec2(0,1)};

  if (scene <= 1) { // box falling onto floor
    Dynamic* box = new Dynamic(boxverts, 4, vec2(0, 5));
    Rigidbody* floor = new Rigidbody(floorverts, 4, vec2(0, -0.5));

    nphys = 1; nstat = 1;
    physs[0] = box;
    stats[0] = floor;

  } else if (scene == 2) { // two boxes hitting each other; two walls
    Dynamic* boxL = new Dynamic(boxverts, 4, vec2(-5, 0), 0,0, vec2(1,  0));
    Dynamic* boxR = new Dynamic(boxverts, 4, vec2(5,  0), 0,0, vec2(-1, 0));
    Rigidbody* wallL = new Rigidbody(floorverts, 4, vec2(-10.5, 0), 1.57);
    Rigidbody* wallR = new Rigidbody(floorverts, 4, vec2(10.5, 0), 1.57);
    Rigidbody* wallB = new Rigidbody(floorverts, 4, vec2(0, -10.5));

    nphys = 2; nstat = 3;
    physs[0] = boxL; physs[1] = boxR;
    stats[0] = wallL; stats[1] = wallR; stats[2] = wallB;

  } else if (scene == 3) { // tower of boxes
    for (int i=0; i < 10; i++)
      physs[i] = new Dynamic(boxverts, 4, vec2(0, 1 + 2*i));
    Rigidbody* floor = new Rigidbody(floorverts, 4, vec2(0, -0.5));

    nphys = 10; nstat = 1;
    stats[0] = floor;
  } else {
    return 0;
  }

  EulerPairwise* sim = new EulerPairwise(physs, nphys, stats, nstat, physloops, dt, nsteps, grav, Sim_FullOutput);
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