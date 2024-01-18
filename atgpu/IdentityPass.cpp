#include "IdentityPass.h"
#include "PassMethodFactory.h"
#include <string.h>

using namespace std;

IdentityPass::IdentityPass() noexcept {
  memset(&elemData,0,sizeof(ELEMENT));
  R1 = nullptr;
  R2 = nullptr;
  T1 = nullptr;
  T2 = nullptr;
  EApertures = nullptr;
  RApertures = nullptr;
}

IdentityPass::~IdentityPass() noexcept {
}

// Retrieve parameters from upper layer (Python, Matlab)
void IdentityPass::getParameters(AbstractInterface *param, PassMethodInfo *info) {

  param->getOptional1DArray(&R1,"R1", 36);
  param->getOptional1DArray(&R2,"R2", 36);
  param->getOptional1DArray(&T1,"T1", 6);
  param->getOptional1DArray(&T2,"T2", 6);
  param->getOptional1DArray(&EApertures,"EApertures", 2);
  param->getOptional1DArray(&RApertures,"RApertures", 4);

  elemData.Type = IDENTITYPASS;
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

AT_FLOAT IdentityPass::getLength() {
  return elemData.Length;
}

uint32_t IdentityPass::getType() {
  return elemData.Type;
}

// Fill device memory
void IdentityPass::fillGPUMemory(uint8_t *startAdd,ELEMENT *elemMem,uint64_t* offset) {

  // Store an offset from the beginning of gpuRing memory in ELEMENT
  // for mapping buffers in GPU memory address space (see Lattice::mapBuffers)

  if(R1) {
    elemData.R1 = (AT_FLOAT *)(*offset);
    memcpy(startAdd+*offset,R1,6*6*sizeof(AT_FLOAT));
    *offset += 6*6*sizeof(AT_FLOAT);
  }
  if(R2) {
    elemData.R2 = (AT_FLOAT *)(*offset);
    memcpy(startAdd+*offset,R2,6*6*sizeof(AT_FLOAT));
    *offset += 6*6*sizeof(AT_FLOAT);
  }
  if(T1) {
    elemData.T1 = (AT_FLOAT *)(*offset);
    memcpy(startAdd+*offset,T1,6*sizeof(AT_FLOAT));
    *offset += 6*sizeof(AT_FLOAT);
  }
  if(T2) {
    elemData.T2 = (AT_FLOAT *)(*offset);
    memcpy(startAdd+*offset,T2,6*sizeof(AT_FLOAT));
    *offset += 6*sizeof(AT_FLOAT);
  }
  if(EApertures) {
    elemData.EApertures = (AT_FLOAT *)(*offset);
    ((AT_FLOAT *)(startAdd+*offset))[0] = 1.0 / EApertures[0];
    ((AT_FLOAT *)(startAdd+*offset))[1] = 1.0 / EApertures[1];
    *offset += 2*sizeof(AT_FLOAT);
  }
  if(RApertures) {
    elemData.RApertures = (AT_FLOAT *)(*offset);
    memcpy(startAdd+*offset,RApertures,4*sizeof(AT_FLOAT));
    *offset += 4*sizeof(AT_FLOAT);
  }

  memcpy(elemMem,&elemData,sizeof(ELEMENT));

}

// Generates GPU code
void IdentityPass::generateCode(std::string& code,PassMethodInfo *info,SymplecticIntegrator &integrator) noexcept {

  generateEnter(code,info);
  generateApertures(code,info);
  generateExit(code,info);

}

void IdentityPass::generateEnter(std::string& code, PassMethodInfo *info) noexcept {

  if( info->doEAperture || info->doRAperture )
    code.append("  bool isLost = false;\n");

  if(info->doT1) generateT(code,"T1");
  if(info->doR1) generateR(code,"R1");

}

void IdentityPass::generateExit(std::string& code, PassMethodInfo *info) noexcept {

  if(info->doR2) generateR(code,"R2");
  if(info->doT2) generateT(code,"T2");

  if( info->doEAperture || info->doRAperture )
    code.append("  if(isLost) r6[5] = INF;\n");

}

void IdentityPass::generateApertures(std::string& code, PassMethodInfo *info) noexcept {

  if(info->doEAperture) generateEAperture(code);
  if(info->doRAperture) generateRAperture(code);

}

void IdentityPass::generateEAperture(std::string& code) noexcept {
  code.append("  if(elem->EApertures) {\n");
  code.append("    AT_FLOAT xnorm = r6[0]*elem->EApertures[0];\n");
  code.append("    AT_FLOAT ynorm = r6[2]*elem->EApertures[1];\n");
  code.append("    isLost |= (xnorm*xnorm + ynorm*ynorm) >= (AT_FLOAT)1.0;\n");
  code.append("  }\n");
}

void IdentityPass::generateRAperture(std::string& code) noexcept {
  code.append("  if(elem->RApertures) {\n"
              "    isLost |= r6[0]<elem->RApertures[0] || r6[0]>elem->RApertures[1] ||\n"
              "              r6[2]<elem->RApertures[2] || r6[2]>elem->RApertures[3];\n"
              "  }\n");
}

void IdentityPass::generateR(std::string& code,const string& pname) noexcept {
  code.append("  if(elem->" + pname + ") transform66(r6,elem->" + pname + ");\n");
}

void IdentityPass::generateT(std::string& code,const string& pname) noexcept {
  code.append("  if(elem->" + pname + ") translate6(r6,elem->" + pname + ");\n");
}

void IdentityPass::generateUtilsFunction(std::string& code, PassMethodInfo *info) noexcept {

  string ftype;
  AbstractGPU::getInstance()->getDeviceFunctionQualifier(ftype);

  // 6D transfrom
  code.append(
          ftype +
          "void translate6(AT_FLOAT* r,AT_FLOAT *t) {\n"
          "  r[0] += t[0];  r[1] += t[1];  r[2] += t[2];\n"
          "  r[3] += t[3];  r[4] += t[4];  r[5] += t[5];\n"
          "}\n"
  );
  code.append(
          ftype +
          "void transform66(AT_FLOAT* r,AT_FLOAT *M) {\n"
          "  int i,j;\n"
          "  AT_FLOAT sum[6];\n"
          "  for(i=0;i<6;i++)\n"
          "  {\n"
          "    sum[i]=0;\n"
          "    for(j=0;j<6;j++)\n"
          "      sum[i]+=M[i+j*6]*r[j];\n"
          "  }\n"
          "  for(i=0;i<6;i++)\n"
          "    r[i]=sum[i];\n"
          "}\n"
  );

}
