#ifndef ATGPU_CORRECTORPASS_H
#define ATGPU_CORRECTORPASS_H
#include "IdentityPass.h"

class CorrectorPass: public IdentityPass {

public:
    // Construct a drift pass
    CorrectorPass() noexcept;

    // Retrieve parameters from upper layer (Python, Matlab)
    void getParameters(AbstractInterface *param, PassMethodInfo *info) override;
    uint64_t getMemorySize() override;
    void fillGPUMemory(uint8_t *startAdd,ELEMENT *elemMem,uint64_t *offset) override;

    // Generic code generation
    static void generateCode(std::string& code, PassMethodInfo *info,SymplecticIntegrator &integrator) noexcept;
    static void generateUtilsFunction(std::string& code, PassMethodInfo *info) noexcept;

protected:

    AT_FLOAT *KickAngle; // KickAngle

};

#endif //ATGPU_CORRECTORPASS_H
