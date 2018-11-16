#include <string>
#include <unistd.h>
#include "stdio.h"

#include "simulate.h"
#include "rigidbody.h"
#include "simulators.h"

int main(int argc, char* argv[]) {
  int scene = 0;
  int nsteps = 100;
  float dt = 0.1;
  char* fname;
  FILE *outCSV, *outSYS;

  int c;
  extern char *optarg;
  extern int optind, opterr, optopt;
  //int nsteps = (argc >= 2) ? std::stoi(argv[1]) : 100;
  //float dt = (argc >= 3) ? std::stof(argv[2]) : 0.1f;

  while ((c = getopt(argc, argv, "s:n:t:f:")) != -1) {
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
    case 'f':
      fname = (char*)malloc(strlen(optarg) + 4+1);
      strcpy(fname, optarg);
      strcat(fname, ".csv");
      outCSV = fopen(fname, "w");

      strcpy(fname, optarg);
      strcat(fname, ".sys");
      outSYS = fopen(fname, "w");
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
  } else {
    return 0;
  }

  EulerPairwise* sim = new EulerPairwise(physs, nphys, stats, nstat, 5, dt, nsteps, grav, Sim_FullOutput);
  if (outCSV) {
    sim->runFullSimulation(outCSV, outSYS);
    fclose(outCSV);
    fclose(outSYS);
  } else {
    sim->runFullSimulation();
  }

  return 0;
}