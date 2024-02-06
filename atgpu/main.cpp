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

AT_FLOAT *createArc(AT_FLOAT radius,AT_FLOAT startAngle,AT_FLOAT endAngle,uint32_t nbParticles) {

  uint32_t rinSize = nbParticles * 6 * sizeof(AT_FLOAT);
  AT_FLOAT* rin = (AT_FLOAT*)malloc(rinSize);
  memset(rin,0,rinSize);

  for(int i=0;i<nbParticles;i++) {
    AT_FLOAT angle = startAngle + (endAngle-startAngle)*(AT_FLOAT)i/(AT_FLOAT)(nbParticles-1);
    RINPTR(i)[0] = radius*cos(angle);
    RINPTR(i)[2] = radius*sin(angle);
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

Lattice *loadLattice(string latticeName,SymplecticIntegrator& integrator,int gpu,double stepSize=0.0) {

  CppInterface *dI = new CppInterface();
  AbstractInterface::setHandler(dI);
  vector<CppObject> elements;

  REPRLoader *loader = new REPRLoader(latticeName);
  loader->parseREPR(elements);
  if( stepSize!=0.0 ) {
    // Adjust NumIntSteps to be as close as possible to stepSize
    for(auto & element : elements) {
      try {
        // For element that has NumIntSteps attribute
        string &stepStr = element.getField("NumIntSteps");
        string &lengthStr = element.getField("Length");
        double length = stod(lengthStr);
        int nstep = (int)((length / stepSize)+0.5);
        if(nstep<1) nstep = 1;
        element.addField("NumIntSteps", to_string(nstep));
      } catch (string&) {
      }
    }
  }
  delete loader;

  Lattice *l = new Lattice(0, integrator, 6e9, gpu);
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

  return l;

}

void integratorTest(int gpu,string latticeName) {

  SymplecticIntegrator integrator(4);

  npy::npy_data<double> refPoints = npy::read_npy<double>("/segfs/tmp/pons/at/test/data/ref_arc_yfr_10000.npy");

  try {

    for(int integ=1;integ<7;integ++) {

      integrator.setType(integ);
      Lattice *l = loadLattice(latticeName,integrator,gpu,0.01);

      uint32_t nbTurn = 1;
      uint32_t nbPart = 256;
      uint32_t refs[] = {l->getNbElement()};
      uint32_t nbRef = sizeof(refs) / sizeof(uint32_t);
      uint64_t routSize = nbTurn * nbPart * nbRef * 6 * sizeof(AT_FLOAT);
      AT_FLOAT *rout = (AT_FLOAT *) malloc(routSize);

      // Choose an arc close to 1mm where unexpected tune drift is observed when step size in too small (EBS lattice)
      AT_FLOAT *rin = createArc(0.001,M_PI/2.0,-M_PI/2.0,nbPart);

      l->run(nbTurn, nbPart, rin, rout, nbRef, refs, 0, nullptr, nullptr, nullptr, nullptr, false);

      double err = 0;
      double max = 0;
      int pMax;
      int cMax;
      for(int p=0;p<256;p++) {
        AT_FLOAT *P = ROUTPTR(p, 0, nbTurn - 1);
        err += SQR(P[0] - refPoints.data[p*6+0]);
        err += SQR(P[1] - refPoints.data[p*6+1]);
        err += SQR(P[2] - refPoints.data[p*6+2]);
        err += SQR(P[3] - refPoints.data[p*6+3]);
        err += SQR(P[4] - refPoints.data[p*6+4]);
        err += SQR(P[5] - refPoints.data[p*6+5]);
        if( abs(P[0] - refPoints.data[p*6+0])>max ) {max = abs(P[0] - refPoints.data[p*6+0]);pMax=p;cMax=0;}
        if( abs(P[1] - refPoints.data[p*6+1])>max ) {max = abs(P[1] - refPoints.data[p*6+1]);pMax=p;cMax=1;}
        if( abs(P[2] - refPoints.data[p*6+2])>max ) {max = abs(P[2] - refPoints.data[p*6+2]);pMax=p;cMax=2;}
        if( abs(P[3] - refPoints.data[p*6+3])>max ) {max = abs(P[3] - refPoints.data[p*6+3]);pMax=p;cMax=3;}
        if( abs(P[4] - refPoints.data[p*6+4])>max ) {max = abs(P[4] - refPoints.data[p*6+4]);pMax=p;cMax=4;}
        if( abs(P[5] - refPoints.data[p*6+5])>max ) {max = abs(P[5] - refPoints.data[p*6+5]);pMax=p;cMax=5;}
      }
      err = sqrt(err) / (6.0*256.0);
      AT_FLOAT *P = ROUTPTR(pMax, 0, nbTurn - 1);
      cout << "[" << integ << "]" << err << " max=" << max << " @" << P[cMax] << " " << abs(max/P[cMax]) << endl;

      /*
      if( integ==4 ) {
        // Save ref
        npy::npy_data_ptr<double> d;
        d.data_ptr = (double *) rout;
        d.shape = {6, nbPart, nbRef, nbTurn};
        d.fortran_order = true;
        npy::write_npy("/segfs/tmp/pons/at/test/data/ref_arc_yfr_10000.npy", d);
      }
      */

      free(rout);
      free(rin);
      delete l;

    }

  } catch (string& errStr) {
    string err =  "Fail: " + errStr;
    cout << "Error: " << err << endl;
  }

}

void performanceTest(int gpu,string latticeName) {

  double t0,t1;
  SymplecticIntegrator integrator(4);

  try {

    Lattice *l = loadLattice(latticeName,integrator,gpu);

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
      uint32_t starts[] = {0,100,200,300,400,500,600,700};
      uint32_t nbStride = sizeof(starts) / sizeof(uint32_t);
      //uint32_t nbStride = 0;

      uint64_t routSize = nbTurn * nbPart * nbRef * 6 * sizeof(AT_FLOAT);
      AT_FLOAT *rout = (AT_FLOAT *) malloc(routSize);

      AT_FLOAT *rin = createGrid(-0.001, -0.001, 0.001, 0.001, nbX, nbY);
      uint32_t *lostAtTurn = new uint32_t[nbPart];
      uint32_t *lostAtElem = new uint32_t[nbPart];
      AT_FLOAT *lostAtCoord = new AT_FLOAT[nbPart * 6];

      t0 = AbstractGPU::get_ticks();
      l->run(nbTurn, nbPart, rin, rout, nbRef, refs, nbStride, starts, lostAtTurn, lostAtElem, lostAtCoord,false);
      t1 = AbstractGPU::get_ticks();

      AT_FLOAT *P = &rin[114*6];
      cout << "[" << 114 << "] " << P[0] << " " << P[1] << " " << P[2] << " " << P[3] << " " << P[4] << " " << P[5] << endl;
      cout << "Lost turn:" << lostAtTurn[114] << endl;
      cout << "Lost elem:" << lostAtElem[114] << endl;
      cout << "Lost coord:" << lostAtCoord[114*6] << "," << lostAtCoord[114*6+2] << endl;

      uint32_t strideSize = nbPart / nbStride;
      for(int stride=0; stride < nbStride; stride++) {
        P = ROUTPTR(stride*strideSize, 0, nbTurn - 1);
        cout << "[" << stride << "] " << P[0] << " " << P[1] << " " << P[2] << " " << P[3] << " " << P[4] << " " << P[5] << endl;
      }

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
    string err =  "Fail: " + errStr;
    cout << "Error: " << err << endl;
  }

}

int main(int argc,char **arv) {

  //printGPUInfo();

  int DEVID = 0;
  //performanceTest(DEVID,"Z:/tmp/pons/at/test/lattice/betamodel_radon.repr");
  //performanceTest(DEVID,"/segfs/tmp/pons/at/test/lattice/betamodel_radon.repr");
  performanceTest(DEVID,"/segfs/tmp/pons/at/test/lattice/betamodel_exact.repr");
  //performanceTest(DEVID,"/segfs/tmp/pons/at/test/lattice/simple_ebs_radon.repr");
  //integratorTest(DEVID,"/segfs/tmp/pons/at/test/lattice/betamodel_radon.repr");

  return 0;

}