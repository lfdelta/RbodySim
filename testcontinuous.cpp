#include "continuous.h"


void TestSphereCast(vec2 A1, vec2 A2, float Ra, vec2 B1, vec2 B2, float Rb, float t) {
  printf("SphereCast: [%0.1f,%0.1f]->[%0.1f,%0.1f](%0.1f) [%0.1f,%0.1f]->[%0.1f,%0.1f](%0.1f)\n",
         A1[0],A1[1],A2[0],A2[1],Ra, B1[0],B1[1],B2[0],B2[1],Rb);
  printf("  Predicted: %f\n", t);
  printf("  Computed:  %f\n", SphereCast(A1, A2, Ra, B1, B2, Rb));
}

int main(int argc, char* argv[]) {
  TestSphereCast(vec2(-10, 0), vec2(0,0), 0, vec2(0, -5), vec2(0,0), 0, 1); // end at origin
  TestSphereCast(vec2(-10, 0), vec2(10,0), 0, vec2(0, -5), vec2(0,5), 0, 0.5); // cross
  TestSphereCast(vec2(0, 0), vec2(10,0), 0, vec2(0, 0), vec2(0,5), 0, 0); // start at origin

  TestSphereCast(vec2(2,2), vec2(2,2), 0, vec2(4,0), vec2(-1,5), 0, 0.4); // one stationary
  TestSphereCast(vec2(0,0), vec2(0,5), 1, vec2(0,5), vec2(0,5), 1, 0.6); // radius testing
  TestSphereCast(vec2(0,0), vec2(0,5), 2, vec2(0,5), vec2(0,5), 1, 0.4); // radius testing
  TestSphereCast(vec2(0,0), vec2(0,5), 1, vec2(0,5), vec2(0,5), 2, 0.4); // radius symmetry

  TestSphereCast(vec2(0, 0), vec2(0,1), 0, vec2(0, 1), vec2(0,2), 0.9, -1); // non-intersection
  TestSphereCast(vec2(-5, 7), vec2(15,7), 10, vec2(50,10), vec2(0,50), 10, -1); // non-intersection

  return 0;
}