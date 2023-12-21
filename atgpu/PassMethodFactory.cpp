#include "PassMethodFactory.h"
#include "DriftPass.h"
#include "StrMPoleSymplectic4RadPass.h"
#include "BndMPoleSymplectic4RadPass.h"
#include "RFCavityPass.h"
#include "AbstractGPU.h"
#include <string.h>
#include <iostream>
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

void PassMethodInfo::Merge(PassMethodInfo *pi) {

  used |= pi->used;
  doR1 |= pi->doR1;
  doR2 |= pi->doR2;
  doT1 |= pi->doT1;
  doT2 |= pi->doT2;
  doEAperture |= pi->doEAperture;
  doRAperture |= pi->doRAperture;
  doQuadEnter |= pi->doQuadEnter;
  doQuadExit |= pi->doQuadExit;

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
    if( elem->getType() != i )
      // An element has been changed from one pass method to an other
      passMethodInfos[elem->getType()].Merge(passMethodInfos+i);
  } else {
    throw string("Not implemented PassMethod: " + passMethod);
  }

  return elem;

}

void PassMethodFactory::generatePassMethods(std::string& code) {

  callCode.clear();

  for(size_t i=0;i<NB_PASSMETHOD_TYPE;i++) {
    if (passMethodInfos[i].used) {

      // Pass method code
      string ftype;
      IdentityPass::getGPUFunctionQualifier(ftype);
      code.append(ftype + "void " + passMethodInfos[i].name + "(AT_FLOAT* r6,ELEMENT* elem,AT_FLOAT turn) {\n");
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
"  int i;\n"
"  AT_FLOAT ReSum = B[max_order];\n"
"  AT_FLOAT ImSum = A[max_order];\n"
"  AT_FLOAT ReSumTemp;\n"
"  for(i = max_order - 1; i >= 0; i--) {\n"
"    ReSumTemp = ReSum * r[0] - ImSum * r[2] + B[i];\n"
"    ImSum = ImSum * r[0] + ReSum * r[2] + A[i];\n"
"    ReSum = ReSumTemp;\n"
"  }\n";

void PassMethodFactory::generateUtilsFunctions(std::string& code) {

  // Generate utils for all pass methods
  for(size_t i=0;i<NB_PASSMETHOD_TYPE;i++)
    passMethodInfos[i].ugenerate(code,passMethodInfos+i);

}