#include "PassMethodFactory.h"
#include "DriftPass.h"
#include "ExactDriftPass.h"
#include "ExactMultipoleRadPass.h"
#include "ThinMPolePass.h"
#include "StrMPoleSymplectic4RadPass.h"
#include "BndMPoleSymplectic4RadPass.h"
#include "RFCavityPass.h"
#include "AbstractGPU.h"
#include <algorithm>

using namespace std;

// Element creation
template<typename T> AbstractElement *create() {
  return new T();
}
// Pass method code generation
template<typename T> void generate(string& code,PassMethodInfo *info,SymplecticIntegrator& integrator) {
  T::generateCode(code,info,integrator);
}
template<typename T> void generateUtil(string& code,PassMethodInfo *info) {
  T::generateUtilsFunction(code,info);
}

#define REGISTER_PASS(_NAME_) PassMethodInfo(#_NAME_,create<_NAME_>, generate<_NAME_>, generateUtil<_NAME_>)

PassMethodFactory::PassMethodFactory(SymplecticIntegrator& integrator) noexcept:integrator(integrator) {

  // Pass methods map
  passMethodInfos[IDENTITYPASS] = REGISTER_PASS(IdentityPass);
  passMethodInfos[DRIFTPASS] = REGISTER_PASS(DriftPass);
  passMethodInfos[BNDMPOLESYMPLECTIC4PASS] = REGISTER_PASS(BndMPoleSymplectic4Pass);
  passMethodInfos[BNDMPOLESYMPLECTIC4RADPASS] = REGISTER_PASS(BndMPoleSymplectic4RadPass);
  passMethodInfos[STRMPOLESYMPLECTIC4PASS] = REGISTER_PASS(StrMPoleSymplectic4Pass);
  passMethodInfos[STRMPOLESYMPLECTIC4RADPASS] = REGISTER_PASS(StrMPoleSymplectic4RadPass);
  passMethodInfos[CAVITYPASS] = REGISTER_PASS(CavityPass);
  passMethodInfos[RFCAVITYPASS] = REGISTER_PASS(RFCavityPass);
  passMethodInfos[EXACTDRIFTPASS] = REGISTER_PASS(ExactDriftPass);
  passMethodInfos[EXACTMULTIPOLEPASS] = REGISTER_PASS(ExactMultipolePass);
  passMethodInfos[EXACTMULTIPOLERADPASS] = REGISTER_PASS(ExactMultipoleRadPass);
  passMethodInfos[THINMPOLEPASS] = REGISTER_PASS(ThinMPolePass);

}

PassMethodInfo::PassMethodInfo() {
  name = "";
  create = nullptr;
  generate = nullptr;
  ugenerate = nullptr;
  used = false;
  doR1 = false;
  doR2 = false;
  doT1 = false;
  doT2 = false;
  doEAperture = false;
  doRAperture = false;
  doQuadEnter = false;
  doQuadExit = false;
}

PassMethodInfo::PassMethodInfo(const string& name,Constructor c,Generator g,UGenerator ug) : PassMethodInfo() {
  this->name = name;
  create = c;
  generate = g;
  ugenerate = ug;
}

AbstractElement *PassMethodFactory::createElement(std::string& passMethod) {

  AbstractElement *elem;
  AbstractInterface *I = AbstractInterface::getInstance();

  // Create an element and update corresponding pass method infos
  bool found = false;
  size_t i = 0;
  while( !found && i<NB_PASSMETHOD_TYPE ) {
    found = passMethod == passMethodInfos[i].name;
    if(!found) i++;
  }
  if( found  ) {
    elem = passMethodInfos[i].create();
    elem->getParameters(I,&passMethodInfos[i]);
  } else {
    throw string("Not implemented PassMethod: " + passMethod);
  }

  return elem;

}

void PassMethodFactory::generatePassMethods(std::string& code) {

  string mType;
  string ftype;
  AbstractGPU::getInstance()->getDeviceFunctionQualifier(ftype);
  AbstractGPU::getInstance()->getGlobalQualifier(mType);

  callCode.clear();

  for(size_t i=0;i<NB_PASSMETHOD_TYPE;i++) {
    if (passMethodInfos[i].used) {

      // Pass method code
      code.append(ftype);
      code.append("void " + passMethodInfos[i].name + "(AT_FLOAT* r6," + mType + " ELEMENT* elem,AT_FLOAT turn) {\n");
      passMethodInfos[i].generate(code, passMethodInfos + i, integrator);
      code.append("}\n");

      // Call code (switch/case)
      string classDefine = passMethodInfos[i].name;
      transform(classDefine.begin(), classDefine.end(), classDefine.begin(),
                [](unsigned char c) { return std::toupper(c); });
      callCode.append("      case " + classDefine + ":\n");
      callCode.append("        " + passMethodInfos[i].name + "(r6,elemPtr,fTurn);\n");
      callCode.append("        break;\n");
    }
  }

}

void PassMethodFactory::generatePassMethodsCalls(std::string& code) {
  code.append(callCode);
}

// Recursively calculate the local transverse magnetic field
string PassMethodFactory::polyLoop =
"  AT_FLOAT ReSum = B[max_order];\n"
"  AT_FLOAT ImSum = A[max_order];\n"
"  for(int i = max_order - 1; i >= 0; i--) {\n"
"    AT_FLOAT ReSumTemp = ReSum * r[0] - ImSum * r[2] + B[i];\n"
"    ImSum = ImSum * r[0] + ReSum * r[2] + A[i];\n"
"    ReSum = ReSumTemp;\n"
"  }\n";

void PassMethodFactory::generateUtilsFunctions(std::string& code) {

  // Generate utils for all pass methods
  for(size_t i=0;i<NB_PASSMETHOD_TYPE;i++)
    passMethodInfos[i].ugenerate(code,passMethodInfos+i);

}