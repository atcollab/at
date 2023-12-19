//
// AT GPU debug file
//

#include "Lattice.h"
#include "REPRLoader.h"
#include <vector>
#include <iostream>
#include <string.h>
#include "npy.hpp"
#include <inttypes.h>
#include "Element.h"

using namespace std;

#define RINPTR(p) (((AT_FLOAT *)rin) + (p * 6))

AT_FLOAT *createGrid(AT_FLOAT x1,AT_FLOAT y1,AT_FLOAT x2,AT_FLOAT y2,uint64_t nbX,uint64_t nbY) {

  uint64_t nbParticles = nbX*nbY;
  uint32_t rinSize = nbParticles * 6 * sizeof(AT_FLOAT);
  AT_FLOAT* rin = (AT_FLOAT*)malloc(rinSize);
  memset(rin,0,rinSize);

  int i = 0;

  if(nbParticles > 1) {
    AT_FLOAT grid_size_h = x2 - x1;
    AT_FLOAT grid_size_v = y2 - y1;
    AT_FLOAT grid_step_h = grid_size_h / ((AT_FLOAT)nbX - 1.0);
    AT_FLOAT grid_step_v = grid_size_v / ((AT_FLOAT)nbY - 1.0);
    for(uint32_t y = 0; y < (uint32_t)nbY; y++) {
      for(uint32_t x = 0; x < (uint32_t)nbX; x++) {
        RINPTR(i)[0] = x1 + (AT_FLOAT)x * grid_step_h;
        RINPTR(i)[2] = y1 + (AT_FLOAT)y * grid_step_v;
        i++;
      }
    }
  } else {
    RINPTR(i)[0] = x1;
    RINPTR(i)[2] = y1;
  }

  return rin;

}


int main(int argc,char **arv) {

  /*
  ELEMENT e;
  cout << sizeof(ELEMENT) <<endl;
  cout << sizeof(e) <<endl;
  printf("0x%" PRIx64 "\n", &e.Type);
  printf("0x%" PRIx64 "\n", &e.SubType);
  printf("0x%" PRIx64 "\n", &e.NumIntSteps);
  printf("0x%" PRIx64 "\n", &e.SL);
  printf("0x%" PRIx64 "\n", &e.MaxOrder);
  printf("0x%" PRIx64 "\n", &e.Length);
  exit(0);
  */

  SymplecticIntegrator integrator(4);
  CppInterface *dI = new CppInterface();
  AbstractInterface::setHandler(dI);
  vector<CppObject> elements;

  try {
    REPRLoader *loader = new REPRLoader("/segfs/tmp/pons/lattice/simple_ebs.repr");
    loader->parseREPR(elements);
  } catch (string& errStr) {
    cout << "Parse failed: " << errStr << endl;
    exit(0);
  }

  string code;

  try {

    Lattice *l = new Lattice(integrator,1);
    for(auto & element : elements) {
      dI->setObject(&element);
      l->addElement();
    }
    l->generateGPUKernel(code);

    uint64_t nbTurn = 1;
    uint64_t nbX = 1;
    uint64_t nbY = 1;
    uint64_t nbPart = nbX * nbY;
    uint32_t refs[] = {l->getNbElement()};
    uint32_t nbRef = sizeof(refs)/sizeof(uint32_t);

    uint64_t routSize = nbTurn * nbPart * nbRef * 6 * sizeof(AT_FLOAT);
    AT_FLOAT* rout = (AT_FLOAT*)malloc(routSize);

    AT_FLOAT *rin = createGrid(-0.001,-0.001,0.001,0.001,nbX,nbY);

    string gpuName = l->getGPUContext()->name() + "(" + to_string(l->getGPUContext()->coreNumber()) + " cores)";
    cout << "Running " << to_string(nbPart) << " particles, " << to_string(nbTurn) << " turn(s) on " << gpuName << endl;
    l->run(nbTurn,nbPart,rin,rout,nbRef,refs);

    npy::npy_data_ptr<double> d;
    d.data_ptr = (double*)rout;
    d.shape = { 6,nbPart,nbRef,nbTurn };
    d.fortran_order = true;
    std::replace(gpuName.begin(),gpuName.end(),' ','_');
    npy::write_npy("/home/esrf/pons/at/test/data/part_" + gpuName + ".npy",d);

    free(rout);
    free(rin);
    delete l;

  } catch (string& errStr) {
    string err =  "at_gpupass() failed: " + errStr;
    cout << "Error: " << err << endl;
  }


  return 0;

}