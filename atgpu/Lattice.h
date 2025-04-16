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
  // Resize lattice (debugging purpose)
  void resize(int32_t nbElements);
  // Return number of element
  uint32_t getNbElement();
  // Load (or reload) lattice in GPU. Lattice structure (element number, polynomial orders, kick angles and misalignment
  // configuration) must remain the same.
  uint64_t fillGPUMemory();
  // Generate and compile GPU code
  void generateGPUKernel();
  // Run the simulation
  void run(uint64_t nbTurn,uint32_t nbParticles,AT_FLOAT *rin,AT_FLOAT *rout,uint32_t nbRef,uint32_t *refPts,
           uint32_t nbElemOffset,uint32_t *elemOffsets,uint32_t *lostAtTurn,uint32_t *lostAtElem,AT_FLOAT *lostAtCoord,
           bool updateRin);
  // Return handle to the GPU context
  GPUContext *getGPUContext();
  // Get ring length
  AT_FLOAT getLength();
  // Set turn counter
  void setTurnCounter(uint64_t count);
  // Return code
  void getCode(std::string& code);

private:

  // Map buffer address to device space
  void mapBuffers(std::string &code);

  void checkLostParticle(std::string &code);
  void storeParticleCoord(std::string &code);

  PassMethodFactory factory;                // Pass method code generation
  std::vector<AbstractElement *> elements;  // All elements

  GPUContext *gpu;                          // GPU context
  void *gpuRing;                            // Ring in GPU memory
  RING_PARAM ringParams;                    // General ring parameters
  std::string code;                         // Kernel code

};

#endif //AT_GPU_LATTICE_H
