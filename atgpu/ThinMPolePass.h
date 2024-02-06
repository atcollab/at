#ifndef AT_GPU_THINMPOLEPASS_H
#define AT_GPU_THINMPOLEPASS_H

#include "StrMPoleSymplectic4Pass.h"

class ThinMPolePass: public StrMPoleSymplectic4Pass {

public:
    // Construct a thin multipole pass (single kick)
    explicit ThinMPolePass() noexcept;
    ~ThinMPolePass() noexcept override;

    // Retrieve parameters from upper layer (Python, Matlab)
    void getParameters(AbstractInterface *param, PassMethodInfo *info) override;
    void fillGPUMemory(uint8_t *startAdd,ELEMENT *elemMem,uint64_t *offset) override;

    // Generic code generation
    static void generateCode(std::string& code, PassMethodInfo *info,SymplecticIntegrator &integrator) noexcept;
    static void generateUtilsFunction(std::string& code, PassMethodInfo *info) noexcept;

protected:
    AT_FLOAT *BendingAngle;  // BendingAngle
    int BendingAngleSize;    // 1-> BAx 2-> BAx,BAy

};

#endif //AT_GPU_THINMPOLEPASS_H
