#ifndef AT_GPU_ELEMENT_H
#define AT_GPU_ELEMENT_H

// Include file shared with host and GPU
// Define GPU<->Host memory exchange
// DEFINE_TYPE (do not remove this comment line, used to define type for gpu code compiler, see genheader.py)

typedef double AT_FLOAT;
#define SQR(x) ((x)*(x))

// GPU thread block size
#define GPU_BLOCK_SIZE 128

// Type of lattice element
// The define must be the corresponding class name in uppercase
#define NB_PASSMETHOD_TYPE           8
#define IDENTITYPASS                 0
#define DRIFTPASS                    1
#define BNDMPOLESYMPLECTIC4PASS      2
#define BNDMPOLESYMPLECTIC4RADPASS   3
#define STRMPOLESYMPLECTIC4PASS      4
#define STRMPOLESYMPLECTIC4RADPASS   5
#define CAVITYPASS                   6
#define RFCAVITYPASS                 7

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

  uint32_t  Type;
  uint32_t  SubType;
  uint32_t  NumIntSteps;
  AT_FLOAT  SL;
  AT_FLOAT  Length;
  AT_FLOAT  *T1;
  AT_FLOAT  *T2;
  AT_FLOAT  *R1;
  AT_FLOAT  *R2;
  AT_FLOAT  *EApertures;
  AT_FLOAT  *RApertures;

  // StrMPole
  uint32_t  MaxOrder;
  AT_FLOAT  K;
  AT_FLOAT  *PolynomA;
  AT_FLOAT  *PolynomB;
  uint32_t  FringeQuadEntrance;
  uint32_t  FringeQuadExit;

  // BndMPole
  AT_FLOAT  irho;
  AT_FLOAT  CRAD;
  uint32_t  FringeBendEntrance; // Method: 1 legacy 2 Soleil 3 ThomX
  AT_FLOAT  EntranceAngle;
  AT_FLOAT  FringeCorrEntranceX;
  AT_FLOAT  FringeCorrEntranceY;
  uint32_t  FringeBendExit; // Method: 1 legacy 2 Soleil 3 ThomX
  AT_FLOAT  ExitAngle;
  AT_FLOAT  FringeCorrExitX;
  AT_FLOAT  FringeCorrExitY;

  // Cavity
  AT_FLOAT  NV; // Voltage/Energy
  AT_FLOAT  FC; // 2.PI.freq/c0
  AT_FLOAT  HC; // C0*(round(freq*T0)/F - T0), T0 = ringLength/C0
  AT_FLOAT  TimeLag;
  AT_FLOAT  PhaseLag;

} ELEMENT;

#endif //AT_GPU_ELEMENT_H
