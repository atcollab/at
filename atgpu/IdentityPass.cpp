#include "IdentityPass.h"
#include "AbstractGPU.h"
#include <string.h>
#include <iostream>

using namespace std;

IdentityPass::IdentityPass() noexcept {
  memset(&elemData,0,sizeof(elemData));
}

IdentityPass::~IdentityPass() noexcept {
  delete[] R1;
  delete[] R2;
  delete[] T1;
  delete[] T2;
  delete[] EApertures;
  delete[] RApertures;
}

// Retrieve parameters from upper layer (Python, Matlab)
void IdentityPass::getParameters(AbstractInterface *param, PASSMETHOD_INFO *info) {

  param->getOptional1DArray(&R1,"R1", 36);
  param->getOptional1DArray(&R2,"R2", 36);
  param->getOptional1DArray(&T1,"T1", 6);
  param->getOptional1DArray(&T2,"T2", 6);
  param->getOptional1DArray(&EApertures,"EApertures", 2);
  param->getOptional1DArray(&RApertures,"RApertures", 4);

  elemData.Type = IDENTITY;
  info->used = true;
  info->doR1 |= (R1 != nullptr);
  info->doR2 |= (R2 != nullptr);
  info->doT1 |= (T1 != nullptr);
  info->doT2 |= (T2 != nullptr);
  info->doEAperture |= (EApertures != nullptr);
  info->doRAperture |= (RApertures != nullptr);

}

uint64_t IdentityPass::getMemorySize() {

  uint64_t sum = 0;
  if(R1) sum += 36 * sizeof(AT_FLOAT);
  if(R2) sum += 36 * sizeof(AT_FLOAT);
  if(T1) sum += 6 * sizeof(AT_FLOAT);
  if(T2) sum += 6 * sizeof(AT_FLOAT);
  if(EApertures) sum += 2 * sizeof(AT_FLOAT);
  if(RApertures) sum += 4 * sizeof(AT_FLOAT);
  return sum;

}

// Fill device memory
void IdentityPass::fillGPUMemory(void *deviceMem) {

  AbstractGPU *gpu = AbstractGPU::getInstance();
  AT_FLOAT *dest = (AT_FLOAT *)deviceMem;

  if(R1) {
    elemData.R1 = dest;
    gpu->hostToDevice(dest,R1,6*6*sizeof(AT_FLOAT));
    dest += 6*6;
  }
  if(R2) {
    elemData.R2 = (AT_FLOAT *)dest;
    gpu->hostToDevice(dest,R2,6*6*sizeof(AT_FLOAT));
    dest += 6*6;
  }
  if(T1) {
    elemData.T1 = (AT_FLOAT *)dest;
    gpu->hostToDevice(dest,T1,6*sizeof(AT_FLOAT));
    dest += 6;
  }
  if(T2) {
    elemData.T2 = (AT_FLOAT *)dest;
    gpu->hostToDevice(dest,T2,6*sizeof(AT_FLOAT));
    dest += 6;
  }
  if(EApertures) {
    elemData.EApertures = (AT_FLOAT *)dest;
    gpu->hostToDevice(dest,EApertures,2*sizeof(AT_FLOAT));
    dest += 2;
  }
  if(RApertures) {
    elemData.RApertures = (AT_FLOAT *)dest;
    gpu->hostToDevice(dest,EApertures,4*sizeof(AT_FLOAT));
    dest += 4;
  }

}

void IdentityPass::generateCall(std::string& code) noexcept {
  code.append("      case IDENTITY:\n");
  code.append("        IdentityPass(r6,elemPtr);\n");
  code.append("        break;\n");
}

// Generates GPU code
void IdentityPass::generateGPUKernel(std::string& code,PASSMETHOD_INFO *info) noexcept {

  code.append("__device__ void IdentityPass(AT_FLOAT* r6,ELEMENT* elem) {\n");
  generateEnter(code,info);
  generateApertures(code,info);
  generateExit(code,info);
  code.append("}\n");

}

void IdentityPass::generateEnter(std::string& code, PASSMETHOD_INFO *info) noexcept {

  if( info->doEAperture || info->doRAperture )
    code.append("  bool isLost = false;\n");

  if(info->doT1) generateT(code,"T1");
  if(info->doR1) generateR(code,"R1");

}

void IdentityPass::generateExit(std::string& code, PASSMETHOD_INFO *info) noexcept {

  if(info->doR2) generateR(code,"R2");
  if(info->doT2) generateT(code,"T2");

  if( info->doEAperture || info->doRAperture )
    code.append("  if(isLost) r6[5] = INF;\n");

}

void IdentityPass::generateApertures(std::string& code, PASSMETHOD_INFO *info) noexcept {

  if(info->doEAperture) generateEAperture(code);
  if(info->doRAperture) generateRAperture(code);

}

void IdentityPass::generateEAperture(std::string& code) noexcept {
  code.append("  if(elem->EApertures) {\n");
  code.append("    AT_FLOAT xnorm = r6[0]/elem->EApertures[0];\n");
  code.append("    AT_FLOAT ynorm = r6[2]/elem->EApertures[1];\n");
  code.append("    isLost |= (xnorm*xnorm + ynorm*ynorm) >= 1;\n");
  code.append("  }\n");
}

void IdentityPass::generateRAperture(std::string& code) noexcept {
  code.append("  if(elem->RApertures) {\n"
              "    isLost |= r6[0]<elem->RApertures[0] || r6[0]>elem->RApertures[1] ||\n"
              "              r6[2]<elem->RApertures[2] || r6[2]>elem->RApertures[3];\n"
              "  }");
}

void IdentityPass::generateR(std::string& code,const string& pname) noexcept {
  code.append("  if(elem->" + pname + ") transform66(r6,elem->" + pname + ");\n");
}

void IdentityPass::generateT(std::string& code,const string& pname) noexcept {
  code.append("  if(elem->" + pname + ") translate6(r6,elem->" + pname + ");\n");
}
