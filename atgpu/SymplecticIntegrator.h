#ifndef AT_GPU_SYMPLECTICINTEGRATOR_H
#define AT_GPU_SYMPLECTICINTEGRATOR_H
#include <string>
#include <vector>

class SymplecticIntegrator {

public:

  // Construct a symplectic integrator of the given type
  SymplecticIntegrator(int type);
  ~SymplecticIntegrator();

  // Set type
  void setType(int type);
  int getType() { return type; }

  // Generate integrator code
  void generateCode(std::string& code);

  // Kick/Drift method for all subtype
  void resetMethods();
  void addDriftMethod(const std::string& driftMethod);
  void addKickMethod(const std::string& kickMethod);

  // Get last kick weight
  double getLastKickWeight();

private:

  // Coefficients
  int type;
  int nbCoefficients;
  double *c; // Drift
  double *d; // Kick

  void allocate(int nb);
  void generateLoopCode(std::string& code,size_t subType);
  bool replace(std::string& str, const std::string& from, const std::string& to);

  std::vector<std::string> driftMethods;
  std::vector<std::string> kickMethods;

};

#endif //AT_GPU_SYMPLECTICINTEGRATOR_H
