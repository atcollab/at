#include <mex.h>
#include <MatlabInterface.h>
#include <AbstractGPU.h>

using namespace std;

// Return GPU Info
const char * gpuInfoNames[] = {"Name","Version","CoreNumber","Platform"};

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  if(nlhs != 1 || nrhs !=0)
    mexErrMsgIdAndTxt("AT:WrongParameter","gpuinfo() must have only one output argument");

  try {
    vector<GPU_INFO> gpuInfos;
    gpuInfos = AbstractGPU::getInstance()->getDeviceList();
    const mwSize dims[] = {(mwSize)gpuInfos.size()};
    plhs[0] = mxCreateCellArray(1,dims);
    for(int i=0;i<gpuInfos.size();i++) {
      mxArray *s = mxCreateStructMatrix(1, 1, 4, gpuInfoNames);
      mxArray *gpuName = mxCreateString(gpuInfos[i].name.c_str());
      mxArray *gpuVersion = mxCreateString(gpuInfos[i].version.c_str());
      mxArray *gpuCore = mxCreateDoubleScalar(gpuInfos[i].mpNumber);
      mxArray *gpuPlatform = mxCreateString(gpuInfos[i].platform.c_str());
      mxSetField(s, 0, gpuInfoNames[0], gpuName);
      mxSetField(s, 0, gpuInfoNames[1], gpuVersion);
      mxSetField(s, 0, gpuInfoNames[2], gpuCore);
      mxSetField(s, 0, gpuInfoNames[3], gpuPlatform);
      mxSetCell(plhs[0],i,s);
    }
  } catch (string& errorStr) {
    mexErrMsgIdAndTxt("AT:Error",errorStr.c_str());
  }

}
