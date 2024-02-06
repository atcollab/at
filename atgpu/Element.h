#ifndef AT_GPU_ELEMENT_H
#define AT_GPU_ELEMENT_H

// Include file shared with host and GPU
// Define GPU<->Host memory exchange
// DEFINE_TYPE (do not remove this comment line, used to define type for gpu code compiler, see genheader.py)

typedef double AT_FLOAT;
#define SQR(x) ((x)*(x))
#define PNORM(x) ((AT_FLOAT)1/((AT_FLOAT)1+(x)))

// Type of lattice element
// The define must be the corresponding class name in uppercase
#define NB_PASSMETHOD_TYPE           12
#define IDENTITYPASS                 0
#define DRIFTPASS                    1
#define BNDMPOLESYMPLECTIC4PASS      2
#define BNDMPOLESYMPLECTIC4RADPASS   3
#define STRMPOLESYMPLECTIC4PASS      4
#define STRMPOLESYMPLECTIC4RADPASS   5
#define CAVITYPASS                   6
#define RFCAVITYPASS                 7
#define EXACTDRIFTPASS               8
#define EXACTMULTIPOLEPASS           9
#define EXACTMULTIPOLERADPASS        10
#define THINMPOLEPASS                11

// Structure alignment
#if defined(__CUDACC__) // NVCC
#define STRUCT_ALIGN(n) __align__(n)
#elif defined(__GNUC__) // GCC
#define STRUCT_ALIGN(n) __attribute__((aligned(n)))
#elif defined(_MSC_VER) // MSVC
#define STRUCT_ALIGN(n) __declspec(align(n))
#else
  #error "Please provide a definition for structure alligmenent for your host compiler!"
#endif

// Ring global parameter
typedef struct STRUCT_ALIGN(8) {
  AT_FLOAT Energy;       // Nominal energy
  AT_FLOAT Length;       // Ring length
  uint64_t turnCounter;  // Turn counter
} RING_PARAM;

// Lattice element
typedef struct STRUCT_ALIGN(8) {

    // FringeBend: Method: 1 legacy 2 Soleil 3 ThomX
    uint32_t  FringeBendEntrance;
    uint32_t  FringeBendExit;
    AT_FLOAT  irho;
    AT_FLOAT  EntranceAngle;
    AT_FLOAT  FringeCorrEntranceX;
    AT_FLOAT  FringeCorrEntranceY;
    AT_FLOAT  ExitAngle;
    AT_FLOAT  FringeCorrExitX;
    AT_FLOAT  FringeCorrExitY;

} MPOLEBEND;

typedef struct STRUCT_ALIGN(8) {

  AT_FLOAT bax;
  AT_FLOAT bay;

} MPOLETHIN;

typedef struct STRUCT_ALIGN(8) {

  // StrMPole
  uint32_t  NumIntSteps;
  uint32_t  MaxOrder;
  AT_FLOAT  K;
  AT_FLOAT  *PolynomA;
  AT_FLOAT  *PolynomB;
  uint32_t  FringeQuadEntrance;
  uint32_t  FringeQuadExit;
  AT_FLOAT  CRAD;

  union {
      MPOLEBEND bnd;
      MPOLETHIN thin;
  };

} MPOLE;

typedef struct STRUCT_ALIGN(8) {

  // Cavity
  AT_FLOAT  NV; // Voltage/Energy
  AT_FLOAT  FC; // 2.PI.freq/c0
  AT_FLOAT  HC; // C0*(round(freq*T0)/F - T0), T0 = ringLength/C0
  AT_FLOAT  TimeLag;
  AT_FLOAT  PhaseLag;

} CAVITY;

typedef struct STRUCT_ALIGN(8) {

  uint32_t  Type;
  uint32_t  SubType;
  AT_FLOAT  SL;
  AT_FLOAT  Length;
  AT_FLOAT  *T1;
  AT_FLOAT  *T2;
  AT_FLOAT  *R1;
  AT_FLOAT  *R2;
  AT_FLOAT  *EApertures;
  AT_FLOAT  *RApertures;

  union {
    MPOLE mpole;
    CAVITY cavity;
  };

} ELEMENT;

#endif //AT_GPU_ELEMENT_H
