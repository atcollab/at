#ifndef AT_GPU_IDENTITYPASS_H
#define AT_GPU_IDENTITYPASS_H
#include "AbstractElement.h"

class IdentityPass: public AbstractElement {

public:

  // Construct an identity pass
  IdentityPass() noexcept;

  // Retrieve parameters from upper layer (Python, Matlab)
  virtual void getParameters(AbstractInterface *param, PASSMETHOD_INFO *info);

  // GPU code generation
  static void generateGPUKernel(std::string& code, PASSMETHOD_INFO *info) noexcept;
  static void generateEnter(std::string& code, PASSMETHOD_INFO *info) noexcept;
  static void generateExit(std::string& code, PASSMETHOD_INFO *info) noexcept;
  static void generateApertures(std::string& code, PASSMETHOD_INFO *info) noexcept;

  static void generateEAperture(std::string& code) noexcept;
  static void generateRAperture(std::string& code) noexcept;
  static void generateR(std::string& code,const std::string& pname) noexcept;
  static void generateT(std::string& code,const std::string& pname) noexcept;
  static void generateCall(std::string& code) noexcept;

protected:

  ELEMENT elemData;

};


#endif //AT_GPU_IDENTITYPASS_H
