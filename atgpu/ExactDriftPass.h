#ifndef AT_GPU_EXACTDRIFTPASS_H
#define AT_GPU_EXACTDRIFTPASS_H
#include "DriftPass.h"

class ExactDriftPass: public DriftPass {

public:
    // Construct an exact drift pass
    ExactDriftPass() noexcept;

    // Retrieve parameters from upper layer (Python, Matlab)
    void getParameters(AbstractInterface *param, PassMethodInfo *info) override;

    // Generic code generation
    static void generateCode(std::string& code, PassMethodInfo *info,SymplecticIntegrator &integrator) noexcept;
    static void generateUtilsFunction(std::string& code, PassMethodInfo *info) noexcept;

};


#endif //AT_GPU_EXACTDRIFTPASS_H
