/* Definition of built-in neighbourhood functions for Kohonen self-organizing
   maps.

   Authors: Johannes Kruisselbrink and Ron Wehrens
*/

#include "neighbourhood-functions.h"

#include <math.h>

/*
 * Returns a neighbourhood function pointer of the specified type.
 */
NeighbourhoodFunctionPtr CreateNeighbourhoodFunction(NeighbourhoodFunctionType type) {
  NeighbourhoodFunctionPtr neighbourhoodFunction;
  switch (type) {
    case BUBBLE:
      neighbourhoodFunction = &Bubble;
      break;
    case GAUSSIAN:
      neighbourhoodFunction = &Gaussian;
      break;
    default:
      neighbourhoodFunction = &Bubble;
      break;
  }
  return neighbourhoodFunction;
}

/*
 * Gaussian neighbourhood function.
 */
double Gaussian(double distance, double radius) {
  return exp(-(distance * distance) / (2 * radius * radius));
}

/*
 * Bubble neighbourhood function.
 */
double Bubble(double distance, double radius) {
  if (distance <= radius) {
    return 1.0;
  }
  return 0.0;
}
