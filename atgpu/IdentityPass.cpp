#include "IdentityPass.h"
using namespace std;

IdentityPass::IdentityPass() noexcept {
  R1 = nullptr;
  R2 = nullptr;
  T1 = nullptr;
  T2 = nullptr;
  EApertures = nullptr;
  RApertures = nullptr;
}

static const vector<int64_t> SHAPE6x6 = {6,6};
static const vector<int64_t> SHAPE6 = {6};
static const vector<int64_t> SHAPE4 = {4};
static const vector<int64_t> SHAPE2 = {2};

// Retrieve parameters from upper layer (Python, Matlab)
void IdentityPass::getParameters(AbstractInterface *param, PASSMETHOD_INFO *info) {

  R1 = param->getOptionalDoubleArray("R1",SHAPE6x6);
  R2 = param->getOptionalDoubleArray("R2",SHAPE6x6);
  T1 = param->getOptionalDoubleArray("T1",SHAPE6);
  T2 = param->getOptionalDoubleArray("T2",SHAPE6);
  EApertures = param->getOptionalDoubleArray("EApertures",SHAPE2);
  RApertures = param->getOptionalDoubleArray("RApertures",SHAPE4);
  Length = 0;

  info->used = true;
  info->doR1 |= (R1 != nullptr);
  info->doR2 |= (R2 != nullptr);
  info->doT1 |= (T1 != nullptr);
  info->doT2 |= (T2 != nullptr);
  info->doEAperture |= (EApertures != nullptr);
  info->doRAperture |= (RApertures != nullptr);

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
    code.append("  if(isLost) r6[4] = INF;\n");

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
