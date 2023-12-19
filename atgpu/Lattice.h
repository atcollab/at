#ifndef AT_GPU_LATTICE_H
#define AT_GPU_LATTICE_H
#include "PassMethodFactory.h"
#include "AbstractGPU.h"

class AbstractElement;
class SymplecticIntegrator;

class Lattice {

public:

  explicit Lattice(SymplecticIntegrator& integrator,int gpuId);
  ~Lattice();

  // Add an element in the lattice
  // (The new element comes from the current loaded object in the AbstractInterface)
  void addElement();
  // Return number of element
  uint32_t getNbElement();
  // Generate device code
  void generateGPUKernel(std::string& code);
  // Fill GPU memory
  void fillGPUMemory();
  // Run the simulation
  void run(uint64_t nbTurn,uint64_t nbParticles,AT_FLOAT *rin,AT_FLOAT *rout,uint32_t nbRef,uint32_t *refPts);
  // Return handle to the GPU context
  GPUContext *getGPUContext();


private:

  void Transpose64(int32_t w,int32_t hm,void *mem);

  PassMethodFactory factory;                // Pass method code generation
  std::vector<AbstractElement *> elements;  // All elements

  GPUContext *gpu;                          // GPU context
  void *gpuRing;                            // Ring in GPU memory
  uint32_t* lost;                           // Particle lost flag

};

#endif //AT_GPU_LATTICE_H
