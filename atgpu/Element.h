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

#if defined(__CUDACC__) // NVCC
#define STRUCT_ALIGN(n) __align__(n)
#elif defined(__GNUC__) // GCC
#define STRUCT_ALIGN(n) __attribute__((aligned(n)))
#elif defined(_MSC_VER) // MSVC
#define STRUCT_ALIGN(n) __declspec(align(n))
#else
  #error "Please provide a definition for structure alligmenent for your host compiler!"
#endif

// Lattice element
typedef struct STRUCT_ALIGN(16) {

  uint32_t  Type;
  uint32_t  NumIntSteps;
  AT_FLOAT  SL;
  AT_FLOAT  Length;
  AT_FLOAT  *T1;
  AT_FLOAT  *T2;
  AT_FLOAT  *R1;
  AT_FLOAT  *R2;
  AT_FLOAT  *EApertures;
  AT_FLOAT  *RApertures;

  uint32_t  MaxOrder;
  AT_FLOAT  *PolynomA;
  AT_FLOAT  *PolynomB;
  AT_FLOAT  *NormD;
  AT_FLOAT  *NormK;
  AT_FLOAT  *KickAngle;

  uint32_t  FringeQuadEntrance;
  uint32_t  FringeQuadExit;

  /*
  AT_FLOAT  irho;
  uint32_t  doFringe;
  uint32_t  FringeBendEntrance; // Method: 1 legacy 2 Soleil 3 ThomX
  AT_FLOAT  EntranceAngle;      // Geometrical edge entrance angle
  AT_FLOAT  FringeInt1;
  AT_FLOAT  tgEntranceAngle;

  uint32_t  FringeBendExit;     // Method: 1 legacy 2 Soleil 3 ThomX
  AT_FLOAT  ExitAngle;          // Geometrical edge exit angle
  AT_FLOAT  FringeInt2;
  AT_FLOAT  tgExitAngle;

  AT_FLOAT  FullGap;
  AT_FLOAT  FringeIntM0[5]; // I0m/K1, I1m/K1, I2m/K1, I3m/K1, Lambda2m/K1
  AT_FLOAT  FringeIntP0[5]; // I0p/K1, I1p/K1, I2p/K1, I3p/K1, Lambda2p/K1
  */

} ELEMENT;

#endif //AT_GPU_ELEMENT_H
