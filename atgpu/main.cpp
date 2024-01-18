//
// AT GPU debug file
//

#include "Lattice.h"
#include "REPRLoader.h"
#include <vector>
#include <iostream>
#include <string.h>
#include "npy.hpp"

using namespace std;

#define RINPTR(p) (((AT_FLOAT *)rin) + (p * 6))
#define ROUTPTR(p,r,t) (rout + ((t)* (6 * nbRef * nbPart) + (r)* (6 * nbPart) + (p) * 6))

AT_FLOAT *createGrid(AT_FLOAT x1,AT_FLOAT y1,AT_FLOAT x2,AT_FLOAT y2,uint32_t nbX,uint32_t nbY) {

  uint32_t nbParticles = nbX*nbY;
  uint32_t rinSize = nbParticles * 6 * sizeof(AT_FLOAT);
  AT_FLOAT* rin = (AT_FLOAT*)malloc(rinSize);
  memset(rin,0,rinSize);

  int i = 0;

  AT_FLOAT grid_size_h = x2 - x1;
  AT_FLOAT grid_size_v = y2 - y1;
  AT_FLOAT grid_step_h = (nbX>1)?grid_size_h / ((AT_FLOAT)nbX - 1.0):0.0;
  AT_FLOAT grid_step_v = (nbY>1)?grid_size_v / ((AT_FLOAT)nbY - 1.0):0.0;
  for(uint32_t y = 0; y < (uint32_t)nbY; y++) {
    for(uint32_t x = 0; x < (uint32_t)nbX; x++) {
      RINPTR(i)[0] = x1 + (AT_FLOAT)x * grid_step_h;
      RINPTR(i)[2] = y1 + (AT_FLOAT)y * grid_step_v;
      i++;
    }
  }

  return rin;

}

void printGPUInfo() {

  try {
    vector<GPU_INFO> infos = AbstractGPU::getInstance()->getDeviceList();
    for(size_t i=0;i<infos.size();i++)
      cout << "Device #" << i << ": " << infos[i].name << " [" << infos[i].version << "] " << infos[i].platform << endl;
  } catch (string& errStr) {
    cout << errStr << endl;
  }

}

int main(int argc,char **arv) {

  printGPUInfo();

  SymplecticIntegrator integrator(4);
  CppInterface *dI = new CppInterface();
  AbstractInterface::setHandler(dI);
  vector<CppObject> elements;
  int DEVID = 1;


  try {
    //REPRLoader *loader = new REPRLoader("Z:/tmp/pons/at/test/lattice/betamodel_radon.repr");
    REPRLoader *loader = new REPRLoader("/segfs/tmp/pons/at/test/lattice/betamodel_radon.repr");
    //REPRLoader *loader = new REPRLoader("/segfs/tmp/pons/lattice/simple_ebs.repr");
    loader->parseREPR(elements);
  } catch (string& errStr) {
    cout << "Parse failed: " << errStr << endl;
    exit(0);
  }

  try {

    Lattice *l = new Lattice(0, integrator, 6e9, DEVID);
    double t0 = AbstractGPU::get_ticks();
    for (auto &element: elements) {
      dI->setObject(&element);
      l->addElement();
    }
    l->generateGPUKernel();
    double t1 = AbstractGPU::get_ticks();
    cout << "Ring build: " << (t1 - t0) * 1000.0 << "ms" << endl;

    t0 = AbstractGPU::get_ticks();
    l->fillGPUMemory();
    t1 = AbstractGPU::get_ticks();
    cout << "GPU lattice loading: " << (t1 - t0) * 1000.0 << "ms" << endl;

    string gpuName = l->getGPUContext()->name();
    cout << "Running test on " << gpuName << " (" << l->getGPUContext()->coreNumber() << " units)"
         << " [" << AbstractGPU::getInstance()->implName() << "]" << endl;

    const int nbTest = 1;
    double testTime[2*nbTest];
    testTime[0] = 0;
    testTime[1] = 0;

    for(int count=0;count<nbTest;count++) {

      uint32_t nbTurn = 100;
      uint32_t nbX = 16*(count+1);
      uint32_t nbY = 16;
      uint32_t nbPart = nbX * nbY;
      uint32_t refs[] = {l->getNbElement()};
      uint32_t nbRef = sizeof(refs) / sizeof(uint32_t);
      uint32_t starts[] = {0,200,300,400,500,600,700,200};
      uint32_t nbStart = sizeof(starts) / sizeof(uint32_t);

      uint64_t routSize = nbTurn * nbPart * nbRef * 6 * sizeof(AT_FLOAT);
      AT_FLOAT *rout = (AT_FLOAT *) malloc(routSize);

      AT_FLOAT *rin = createGrid(-0.001, -0.001, 0.001, 0.001, nbX, nbY);
      uint32_t *lostAtTurn = new uint32_t[nbPart];
      uint32_t *lostAtElem = new uint32_t[nbPart];
      AT_FLOAT *lostAtCoord = new AT_FLOAT[nbPart * 6];

      t0 = AbstractGPU::get_ticks();
      l->run(nbTurn, nbPart, rin, rout, nbRef, refs, nbStart, starts, lostAtTurn, lostAtElem, lostAtCoord);
      t1 = AbstractGPU::get_ticks();

      //int pIdx = 0;
      //cout << "lostAtTurn[" << pIdx << "]" << lostAtTurn[pIdx] << endl;
      //cout << "lostAtElem[" << pIdx << "]" << lostAtElem[pIdx] << endl;
      //cout << "lostAtCoord[" << pIdx << "]" << lostAtCoord[pIdx*6+0] << "," << lostAtCoord[pIdx*6+1] << endl;

      AT_FLOAT *P;
      P = ROUTPTR(0,0,nbTurn-1);
      cout << P[0] << " " << P[1] << " " << P[2] << " " << P[3] << " " << P[4] << " " << P[5] << endl;
      P = ROUTPTR(nbPart-1,0,nbTurn-1);
      cout << P[0] << " " << P[1] << " " << P[2] << " " << P[3] << " " << P[4] << " " << P[5] << endl;

      int nbLost = 0;
      for (uint32_t i = 0; i < nbPart; i++)
        if (lostAtTurn[i] != nbTurn)
          nbLost++;

      cout << "Lost: " << nbLost << "/" << nbPart << ": " << (t1-t0)*1000.0 << "ms" << endl;
      testTime[2*count+0] = t1-t0;
      testTime[2*count+1] = nbPart;

      /*
      npy::npy_data_ptr<double> d;
      d.data_ptr = (double*)rout;
      d.shape = { 6,nbPart,nbRef,nbTurn };
      d.fortran_order = true;
      std::replace(gpuName.begin(),gpuName.end(),' ','_');
      npy::write_npy("/segfs/tmp/pons/at/test/part_" + gpuName + ".npy",d);
      */

      free(rout);
      free(rin);

    }

    /*
    npy::npy_data_ptr<double> d;
    d.data_ptr = (double*)testTime;
    d.shape = { 2,(uint32_t)nbTest };
    d.fortran_order = true;
    std::replace(gpuName.begin(),gpuName.end(),' ','_');
    npy::write_npy("/segfs/tmp/pons/at/test/data/timing_" + gpuName + "_w128_fp32.npy",d);
    */

    delete l;

  } catch (string& errStr) {
    string err =  "at_gpupass() failed: " + errStr;
    cout << "Error: " << err << endl;
  }


  return 0;

}