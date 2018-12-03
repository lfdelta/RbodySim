#include "continuous.h"


void TestSphereCast(vec2 A1, vec2 A2, float Ra, vec2 B1, vec2 B2, float Rb, float t) {
  printf("SphereCast: [%0.1f,%0.1f]->[%0.1f,%0.1f](%0.1f) [%0.1f,%0.1f]->[%0.1f,%0.1f](%0.1f)\n",
         A1[0],A1[1],A2[0],A2[1],Ra, B1[0],B1[1],B2[0],B2[1],Rb);
  printf("  Predicted: %f\n", t);
  printf("  Computed:  %f\n", SphereCast(A1, A2, Ra, B1, B2, Rb));
}



void TestContacts(vec2* v1, vec2* n1, int s1, vec2* v2, vec2* n2, int s2) {
  ContactMesh cmesh = EdgeContactPoints(v1, n1, s1, v2, n2, s2);

  printf("verts1:");
  for (int i=0; i < s1; i++)
    printf(" {%0.1f, %0.1f}", v1[i][0], v1[i][1]);
  printf("\n");

  printf("verts2:");
  for (int i=0; i < s2; i++)
    printf(" {%0.1f, %0.1f}", v2[i][0], v2[i][1]);
  printf("\n");

  printf("  contacts:");
  for (int i=0; i < cmesh.sz; i++)
    printf(" p{%0.1f, %0.1f}n{%0.1f, %0.1f}", cmesh.pos[i][0], cmesh.pos[i][1], cmesh.nhat[i][0], cmesh.nhat[i][1]);
  printf("\n");
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

  vec2 box[] = {vec2(0,0), vec2(1,0), vec2(1,1), vec2(0,1)};
  vec2 floor[] = {vec2(-10,-1), vec2(10,-1), vec2(10,0.001), vec2(-10,0.001)};
  vec2 sqnorms[] = {vec2(0, -1), vec2(1, 0), vec2(0, 1), vec2(-1, 0)};
  TestContacts(box, sqnorms, 4, floor, sqnorms, 4);

  vec2 tri[] = {vec2(1,1), vec2(-1,1), vec2(0,0)};
  vec2 trinorms[] = {vec2(0,1), vec2(-1,-1)/sqrt(2), vec2(1,-1)/sqrt(2)};
  TestContacts(tri, trinorms, 3, floor, sqnorms, 4);

  vec2 tri2[] = {vec2(1+3,1), vec2(-1+3,1), vec2(0+3,0)};
  TestContacts(tri2, trinorms, 3, floor, sqnorms, 4);

  TestContacts(floor, sqnorms, 4, tri, trinorms, 3);

  return 0;
}