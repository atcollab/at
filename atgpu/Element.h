#ifndef AT_GPU_ELEMENT_H
#define AT_GPU_ELEMENT_H

// Include file shared with host and GPU
// Define GPU<->Host memory exchange
// DEFINE_TYPE (do not remove this comment line, used to define type for gpu code compiler, see genheader.py)

typedef double AT_FLOAT;

// GPU thread block size
#define GPU_GRID 128

// Type of lattice element
#define NB_PASSMETHOD_TYPE 4

#define IDENTITY  0
#define DRIFT     1
#define BEND      2
#define MPOLE     3

#define MAX_POLYNOMIAL_ORDER 32

// Lattice element
typedef struct {

  uint32_t  Type;
  uint32_t  NumIntSteps;
  uint32_t  MaxOrder;
  uint32_t  doFringe;
  AT_FLOAT  Length;
  AT_FLOAT  PolynomA[MAX_POLYNOMIAL_ORDER];
  AT_FLOAT  PolynomB[MAX_POLYNOMIAL_ORDER];
  AT_FLOAT  irho;
  AT_FLOAT  L1;
  AT_FLOAT  L2;
  AT_FLOAT  K1;
  AT_FLOAT  K2;
  AT_FLOAT  T1[6];
  AT_FLOAT  T2[6];
  AT_FLOAT  R1[6];
  AT_FLOAT  R2[6];
  AT_FLOAT  EApertures[2];
  AT_FLOAT  RApertures[4];

  uint32_t  FringeBendEntrance; // Method: 1 legacy 2 Soleil 3 ThomX
  AT_FLOAT  EntranceAngle;      // Geometrical edge entrance angle
  AT_FLOAT  FringeInt1;
  AT_FLOAT  tgEntranceAngle;

  uint32_t  FringeBendExit;     // Method: 1 legacy 2 Soleil 3 ThomX
  AT_FLOAT  ExitAngle;          // Geometrical edge exit angle
  AT_FLOAT  FringeInt2;
  AT_FLOAT  tgExitAngle;

  AT_FLOAT  FullGap;

  uint32_t  FringeQuadEntrance;
  uint32_t  FringeQuadExit;
  AT_FLOAT  FringeIntM0[5]; // I0m/K1, I1m/K1, I2m/K1, I3m/K1, Lambda2m/K1
  AT_FLOAT  FringeIntP0[5]; // I0p/K1, I1p/K1, I2p/K1, I3p/K1, Lambda2p/K1

} ELEMENT;

#endif //AT_GPU_ELEMENT_H
