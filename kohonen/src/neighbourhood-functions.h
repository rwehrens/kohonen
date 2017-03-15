#ifndef NEIGHBOURHOOD_H
#define NEIGHBOURHOOD_H

typedef enum {
  BUBBLE = 1,
  GAUSSIAN = 2
} NeighbourhoodFunctionType;

typedef double (*NeighbourhoodFunctionPtr)(double, double);

/*
 * Returns a neighbourhood function pointer of the specified type.
 */
NeighbourhoodFunctionPtr CreateNeighbourhoodFunction(NeighbourhoodFunctionType type);

/*
 * Gaussian neighbourhood function.
 */
double Gaussian(double distance, double radius);

/*
 * Bubble neighbourhood function.
 */
double Bubble(double distance, double radius);

#endif
