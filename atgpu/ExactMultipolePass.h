#ifndef AT_GPU_EXACTMULTIPOLEPASS_H
#define AT_GPU_EXACTMULTIPOLEPASS_H
#include "StrMPoleSymplectic4Pass.h"
#include "SymplecticIntegrator.h"

class ExactMultipolePass: public StrMPoleSymplectic4Pass {

public:
    // Construct an exact multipole pass
    explicit ExactMultipolePass() noexcept;
    ~ExactMultipolePass() noexcept override;

    // Retrieve parameters from upper layer (Python, Matlab)
    void getParameters(AbstractInterface *param, PassMethodInfo *info) override;

    // Generic code generation
    static void generateCode(std::string& code, PassMethodInfo *info,SymplecticIntegrator &integrator) noexcept;
    static void generateUtilsFunction(std::string& code, PassMethodInfo *info) noexcept;

    static void generateQuadFringeEnter(std::string& code, PassMethodInfo *info) noexcept;
    static void generateQuadFringeExit(std::string& code, PassMethodInfo *info) noexcept;

protected:

    AT_FLOAT *PolynomA;  // PolynomA
    AT_FLOAT *PolynomB;  // PolynomB
    AT_FLOAT *KickAngle; // KickAngle

private:

    bool isQuadrupole();
    bool isSextupole();
    bool isOctupole();

};


#endif //AT_GPU_EXACTMULTIPOLEPASS_H
