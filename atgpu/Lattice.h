#ifndef AT_GPU_LATTICE_H
#define AT_GPU_LATTICE_H
#include "PassMethodFactory.h"
#include "AbstractGPU.h"

class AbstractElement;
class SymplecticIntegrator;

class Lattice {

public:

  explicit Lattice(int32_t nbElements,SymplecticIntegrator& integrator,double energy,int gpuId);
  ~Lattice();

  // Add an element in the lattice
  // (The new element comes from the current loaded object in the AbstractInterface)
  void addElement();
  // Return number of element
  uint32_t getNbElement();
  // Load (or reload) lattice in GPU. Lattice structure (polynomial orders, kick angles, misalignment configuration) must remain the same.
  void fillGPUMemory();
  // Generate and compile GPU code
  void generateGPUKernel();
  // Run the simulation
  void run(uint64_t nbTurn,uint64_t nbParticles,AT_FLOAT *rin,AT_FLOAT *rout,uint32_t nbRef,uint32_t *refPts,
           uint32_t *lostAtTurn,uint32_t *lostAtElem,AT_FLOAT *lostAtCoord);
  // Return handle to the GPU context
  GPUContext *getGPUContext();
  // Get ring length
  AT_FLOAT getLength();
  // Set turn counter
  void setTurnCounter(uint64_t count);

private:


  PassMethodFactory factory;                // Pass method code generation
  std::vector<AbstractElement *> elements;  // All elements

  GPUContext *gpu;                          // GPU context
  void *gpuRing;                            // Ring in GPU memory
  RING_PARAM ringParams;                    // General ring parameters

};

#endif //AT_GPU_LATTICE_H
