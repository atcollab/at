#ifndef AT_GPU_EXACTMULTIPOLERADPASS_H
#define AT_GPU_EXACTMULTIPOLERADPASS_H

#include "ExactMultipolePass.h"

class ExactMultipoleRadPass: public ExactMultipolePass {

public:
    // Construct an exact multipole pass with radiation
    explicit ExactMultipoleRadPass() noexcept;
    ~ExactMultipoleRadPass() noexcept override;

    // Retrieve parameters from upper layer (Python, Matlab)
    void getParameters(AbstractInterface *param, PassMethodInfo *info) override;

    // Generic code generation
    static void generateCode(std::string& code, PassMethodInfo *info,SymplecticIntegrator &integrator) noexcept;
    static void generateUtilsFunction(std::string& code, PassMethodInfo *info) noexcept;

};

#endif //AT_GPU_EXACTMULTIPOLERADPASS_H
