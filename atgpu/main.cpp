//
// AT GPU test file
//

#include "CudaGPU.h"
#include "iostream"
#include <Python.h>
#include <string.h>

using namespace std;

int main(int argc,char **arv) {

  AbstractGPU *gpu = AbstractGPU::getInstance();

  try {



  } catch (string& errStr) {
    cerr << errStr << endl;
    return -1;
  }



  return 0;

}