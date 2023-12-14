//
// AT GPU debug file
//

#include "Lattice.h"
#include "REPRLoader.h"
#include <vector>
#include <iostream>

using namespace std;

int main(int argc,char **arv) {

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

  //try {

    Lattice *l = new Lattice(integrator,0);
    for(auto & element : elements) {
      dI->setObject(&element);
      l->addElement();
    }
    l->generateGPUKernel(code,true);
    //cout << code << endl;
    delete l;

  //} catch (string& errStr) {
  //  string err =  "at_gpupass() failed: " + errStr;
  //  cout << "Error: " << err << endl;
  //}


  return 0;

}