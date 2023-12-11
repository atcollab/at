//
// AT GPU test file
//

#include "CudaGPU.h"
#include "iostream"
#include "Lattice.h"
#include "AbstractInterface.h"

using namespace std;

class CppElement {
public:
  CppElement() {
  }
  void addField(const string& name,const string& value) {
    fieldNames.push_back(name);
    fieldValue.push_back(value);
  }
  string getField(const string& name) {
    bool found = false;
    size_t i = 0;
    while(!found && i<fieldNames.size()) {
      found = fieldNames[i] == name;
      if(!found) i++;
    }
    if(!found)
      throw string(name + " not found");
    return fieldValue[i];
  }

  vector<string> fieldNames;
  vector<string> fieldValue;
};

class DebugInterface: public AbstractInterface {

public:
  DebugInterface() {
  }

  std::string getString(const std::string& name) {
    return elem->getField(name);
  }

  virtual int getInt(const std::string& name) {
    string val = elem->getField(name);
    return stoi(val);
  }

  AT_FLOAT getDouble(const std::string& name) {
    string val = elem->getField(name);
    return stod(val);
  }

  AT_FLOAT *getNativeDoubleArray(const std::string& name,std::vector<int64_t>& shape) {

    shape.clear();
    string val = elem->getField(name);
    vector<string> tokens;
    split(tokens,val,',');
    size_t elIdx=0;

    if(tokens.size()>1) {
      string shapeStr = tokens[0];
      vector<string> stokens;
      split(stokens,shapeStr,'x');
      for(const auto & stoken : stokens)
        shape.push_back(stoi(stoken));
      elIdx++;
    } else {
      shape.push_back(1);
    }

    AT_FLOAT *d = new AT_FLOAT[tokens.size()];
    for(int i=0;elIdx<tokens.size();elIdx++,i++)
      d[i] = stod(tokens[elIdx]);
    return d;

  }

  void setObject(CppElement *obj) {
    elem = obj;
  }

private:
  CppElement *elem;

  void split(vector<string> &tokens, const string &text, char sep) {

    size_t start = 0, end = 0;
    tokens.clear();

    while ((end = text.find(sep, start)) != string::npos) {
      tokens.push_back(text.substr(start, end - start));
      start = end + 1;
    }

    tokens.push_back(text.substr(start));

  }

};



int main(int argc,char **arv) {

  SymplecticIntegrator integrator(4);
  DebugInterface *dI = new DebugInterface();
  AbstractInterface::setHandler(dI);

  vector<CppElement *> elements;

  CppElement *I = new CppElement();
  I->addField("PassMethod","IdentityPass");
  I->addField("T1","6,0,0,0,0,0,0");
  I->addField("R1","6x6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0");
  I->addField("EApertures","0,0");
  elements.push_back(I);

  CppElement *D1 = new CppElement();
  D1->addField("PassMethod","DriftPass");
  D1->addField("Length","5");
  elements.push_back(D1);

  CppElement *D2 = new CppElement();
  D2->addField("PassMethod","DriftPass");
  D2->addField("Length","5");
  elements.push_back(D2);

  CppElement *Q1 = new CppElement();
  Q1->addField("PassMethod","StrMPoleSymplectic4Pass");
  Q1->addField("Length","5");
  Q1->addField("PolynomA","2,0,0");
  Q1->addField("PolynomB","2,0,1.2");
  Q1->addField("MaxOrder","1");
  Q1->addField("FringeQuadEntrance","1");
  Q1->addField("FringeQuadExit","1");
  elements.push_back(Q1);

  string code;

  try {

    Lattice *l = new Lattice(integrator);
    for(auto & element : elements) {
      dI->setObject(element);
      l->addElement();
    }
    l->generateGPUKernel(code);
    cout << code << endl;
    AbstractGPU::getInstance()->run(code);

  } catch (string& errStr) {
    string err =  "at_gpupass() failed: " + errStr;
    cout << "Error: " << err << endl;
  }

  return 0;

}