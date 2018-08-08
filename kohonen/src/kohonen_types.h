typedef enum {
    SUMOFSQUARES = 1,
    EUCLIDEAN = 2,
    MANHATTAN = 3,
    TANIMOTO = 4,
    DTW      = 5
} DistanceType;

typedef double (*DistanceFunctionPtr)(double *, double *, int, int);
typedef DistanceFunctionPtr* DistanceFunctionPtrs;
