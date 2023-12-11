#ifndef AT_GPU_LATTICE_H
#define AT_GPU_LATTICE_H
#include "PassMethodFactory.h"

class AbstractElement;
class SymplecticIntegrator;

class Lattice {

public:

  explicit Lattice(SymplecticIntegrator& integrator);

  // Add an element in the lattice
  // (The new element comes from the current loaded object in the AbstractInterface)
  void addElement();
  // Generate device code
  void generateGPUKernel(std::string& code);

private:

  PassMethodFactory factory;
  std::vector<AbstractElement *> elements;  // All elements

};

#endif //AT_GPU_LATTICE_H
