#ifndef AT_GPU_REPRLOADER_H
#define AT_GPU_REPRLOADER_H
#include <string>
#include <vector>
#include "AbstractInterface.h"
#include "inttypes.h"

// CppObject (debugging purpose)
class CppObject {
public:
  CppObject();
  void addField(const std::string& name,const std::string& value);
  void removeField(const std::string &name);
  std::string& getField(const std::string& name);
private:
  std::vector<std::string> fieldNames;
  std::vector<std::string> fieldValues;
};

// Abstract interface implementation for CppObject (debugging purpose)
class CppInterface: public AbstractInterface {

public:

  CppInterface();
  std::string getString(const std::string& name) override;
  int getInt(const std::string& name) override;
  double getDouble(const std::string& name) override;
  double *getNativeDoubleArray(const std::string& name,std::vector<int64_t>& shape) override;
  float *getNativeFloatArray(const std::string& name,std::vector<int64_t>& shape) override;
  void setObject(CppObject *obj);

private:

  CppObject *elem;

};

// Load a lattice in repr format (debugging purpose)
class REPRLoader {

public:

  // Construct a .repr file loader
  REPRLoader(const std::string& fileName);
  // Parse a REPR file
  void parseREPR(std::vector<CppObject> &elems);

  // Global item
  CppObject globals;

private:

  size_t currentPos;
  char nextNextChar;
  char nextChar;
  char currentChar;
  bool backSlashed;
  std::string fileName;
  FILE *f;
  std::string fileBuffer;

  void parseExtraParams(CppObject& obj);
  void parseParam(const std::string& name,CppObject& obj);
  void parseIdentity(CppObject& obj);
  void parseDrift(CppObject& obj);
  void parseDipole(CppObject& obj);
  void parseQuadrupole(CppObject& obj);
  void parseSextupole(CppObject& obj);
  void parseMultipole(CppObject& obj);
  void parseRFCavity(CppObject& obj);
  void jumpSpace();
  void readWord(std::string& word);
  void jumpSep(const std::string& sep);
  void jumpSep(char sep);
  void parseParamValue(std::string& value,std::vector<int64_t> *shapePtr= nullptr);
  int getPosMarker() const;
  void parseArray(std::vector<std::string> &ret,std::vector<int64_t>& shape,int level);
  void parseArrayType(std::string& typeStr);
  void getArrayValue(std::string& str,std::vector<int64_t>& shape,std::vector<std::string>& values);
  std::string getErrorLocation(int pos);
  char current();
  bool endOf(char c);
  bool eof();
  bool isSpecialChar(char c);
  char read();
  void readChar();
  void toNextChar();
  std::string getCoord(int pos);

};


#endif //AT_GPU_REPRLOADER_H
