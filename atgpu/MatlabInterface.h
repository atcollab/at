#ifndef AT_GPU_MATLABINTERFACE_H
#define AT_GPU_MATLABINTERFACE_H
#include "AbstractInterface.h"
#include <mex.h>

class MatlabInterface: public AbstractInterface {

public:

    std::string getString(const std::string& name) override;
    int getInt(const std::string& name) override;
    double getDouble(const std::string& name) override;
    double *getNativeDoubleArray(const std::string& name,std::vector<int64_t>& shape) override;
    float *getNativeFloatArray(const std::string& name,std::vector<int64_t>& shape) override;

    void setObject(mxArray *obj);

private:

    mxArray *getField(const mxArray *pm, const std::string& name);
    mxArray *elem = nullptr;

};

#endif //AT_GPU_MATLABINTERFACE_H
