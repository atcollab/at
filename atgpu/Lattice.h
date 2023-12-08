#ifndef AT_GPU_LATTICE_H
#define AT_GPU_LATTICE_H
#include "AbstractElement.h"

class AbstractElement;

class Lattice {

public:

  explicit Lattice(size_t nElements);
  // Add an element in the lattice
  // (The new element comes from the current loaded object in the AbstractInterface)
  void addElement();
  // Generate pass methods device code
  void generatePassMethods(std::string& code);
  // Generate device code
  void generateGPUKernel(std::string& code);

private:

  std::vector<AbstractElement *> elements;  // All elements

};

#endif //AT_GPU_LATTICE_H
