#include "../defs.h"
#include <assert.h>

int main() {

  int numPoints = 108;
  char vtkfilename[STRLEN] = "test.vtk";

  vec rr[numPoints], rw[numPoints];

  printf("Test 1: Checking if readPointCloudFromVTK works\n");

  // initialize rw with something
  for (int i = 0; i < numPoints; i++) {
    rw[i].x = i*0.5;
    rw[i].y = i*0.6;
    rw[i].z = i*0.7;
  }
  // write vtk
  writePointCloudToVTK(vtkfilename, rw, numPoints);
  // read vtk again
  int RR = readPointCloudFromVTK(vtkfilename, rr, numPoints);
  if (RR!=0) {return -1;}

  for (int i = 0; i < numPoints; i++) {
    assert((
    rw[i].x == rr[i].x &&
    rw[i].y == rr[i].y &&
    rw[i].z == rr[i].z
    ));
  }

  printf("Test passed 👌.\n");
  printf("All tests passed 👌🏞️🏖️🚀.\n");

  return 0;
}
