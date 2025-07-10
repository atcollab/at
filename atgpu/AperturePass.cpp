#include "AperturePass.h"
#include "PassMethodFactory.h"
#include <string.h>

using namespace std;

AperturePass::AperturePass() noexcept {
  memset(&elemData,0,sizeof(ELEMENT));
  RApertures = nullptr;
}

AperturePass::~AperturePass() noexcept {
}

// Retrieve parameters from upper layer (Python, Matlab)
void AperturePass::getParameters(AbstractInterface *param, PassMethodInfo *info) {

  param->get1DArray(&RApertures,"Limits", 4);

  elemData.Type = APERTUREPASS;
  info->used = true;
  info->doRAperture |= (RApertures != nullptr);

}

uint64_t AperturePass::getMemorySize() {

  uint64_t sum = 0;
  if(RApertures) sum += 4 * sizeof(AT_FLOAT);
  return sum;

}

AT_FLOAT AperturePass::getLength() {
  return elemData.Length;
}

uint32_t AperturePass::getType() {
  return elemData.Type;
}

// Fill device memory
void AperturePass::fillGPUMemory(uint8_t *startAdd,ELEMENT *elemMem,uint64_t* offset) {

  // Store an offset from the beginning of gpuRing memory in ELEMENT
  // for mapping buffers in GPU memory address space (see Lattice::mapBuffers)
  // We keep RAperture to avoid increasing element size
  elemData.RApertures = (AT_FLOAT *)(*offset);
  memcpy(startAdd+*offset,RApertures,4*sizeof(AT_FLOAT));
  *offset += 4*sizeof(AT_FLOAT);

  memcpy(elemMem,&elemData,sizeof(ELEMENT));

}

// Generates GPU code
void AperturePass::generateCode(std::string& code,PassMethodInfo *info,SymplecticIntegrator &integrator) noexcept {

  generateApertures(code,info);

}

void AperturePass::generateApertures(std::string& code, PassMethodInfo *info) noexcept {

  if(info->doRAperture) generateRAperture(code);

}

void AperturePass::generateRAperture(std::string& code) noexcept {
  code.append("  if(elem->RApertures) {\n"
              "    isLost |= r6[0]<elem->RApertures[0] || r6[0]>elem->RApertures[1] ||\n"
              "              r6[2]<elem->RApertures[2] || r6[2]>elem->RApertures[3];\n"
              "  }\n");
}

void AperturePass::generateUtilsFunction(std::string& code, PassMethodInfo *info) noexcept {

}
