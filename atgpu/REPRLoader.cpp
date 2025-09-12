// Load a lattice in REPR format (debugging purpose)

#include "REPRLoader.h"
#include <errno.h>
#include <string.h>

using namespace std;

// ---------------------------------------------------------------------------------------------------------------

CppObject::CppObject() {
}

void CppObject::addField(const string &name, const string &value) {
  fieldNames.push_back(name);
  fieldValues.push_back(value);
}

string& CppObject::getField(const string &name) {
  // Search from this end
  // overloaded param may appear several times
  bool found = false;
  int i = (int)fieldNames.size() - 1;
  while (!found && i >= 0) {
    found = fieldNames[i] == name;
    if (!found) i--;
  }
  if (!found)
    throw string(name + " not found");
  return fieldValues[i];
}

void CppObject::removeField(const string &name) {
    bool found = false;
    int i = (int)fieldNames.size() - 1;
    while (!found && i >= 0) {
        found = fieldNames[i] == name;
        if (!found) i--;
    }
    if (found) {
        fieldNames.erase(fieldNames.begin()+i);
        fieldValues.erase(fieldValues.begin()+i);
    } else {
        throw string(name + " not found");
    }
}

// ---------------------------------------------------------------------------------------------------------------

CppInterface::CppInterface() {
}

std::string CppInterface::getString(const std::string &name) {
  return elem->getField(name);
}

int CppInterface::getInt(const std::string &name) {
  string val = elem->getField(name);
  return stoi(val);
}

double CppInterface::getDouble(const std::string &name) {
  string val = elem->getField(name);
  return stod(val);
}

double *CppInterface::getNativeDoubleArray(const std::string &name, std::vector<int64_t> &shape) {

  shape.clear();
  string val = elem->getField(name);
  vector<string> tokens;
  split(tokens, val, ',');
  size_t elIdx = 0;

  if (tokens.size() > 1) {
    getShapeFromStr(shape,tokens[0]);
    elIdx++;
  } else {
    shape.push_back(1);
  }

  double *d = new double[tokens.size()];
  for (int i = 0; elIdx < tokens.size(); elIdx++, i++)
    d[i] = stod(tokens[elIdx]);
  return d;

}

float *CppInterface::getNativeFloatArray(const std::string &name, std::vector<int64_t> &shape) {

  shape.clear();
  string val = elem->getField(name);
  vector<string> tokens;
  split(tokens, val, ',');
  size_t elIdx = 0;

  if (tokens.size() > 1) {
    getShapeFromStr(shape,tokens[0]);
    elIdx++;
  } else {
    shape.push_back(1);
  }

  float *d = new float[tokens.size()];
  for (int i = 0; elIdx < tokens.size(); elIdx++, i++)
    d[i] = stof(tokens[elIdx]);
  return d;

}

void CppInterface::setObject(CppObject *obj) {
  elem = obj;
}


// ---------------------------------------------------------------------------------------------------------------

REPRLoader::REPRLoader(const std::string &fileName) {

  this->fileName = fileName;
  f = fopen(fileName.c_str(),"r");
  if( f== nullptr )
    throw string(fileName + ": " + strerror(errno));

  // Determine file size
  fseek(f, 0, SEEK_END);
  size_t fSize = ftell(f);
  if( fSize==0 ) {
    fclose(f);
    throw string(fileName + ": empty file");
  }

  // Read buffer
  char *buff = new char[fSize+1];
  rewind(f);
  size_t nbRead = fread(buff, 1, fSize, f);
  fclose(f);
  if( nbRead != fSize ) {
    throw string(fileName + ": " + strerror(errno));
  }
  buff[fSize]=0;
  fileBuffer.assign(buff);
  delete[] buff;

  // Init parser variables
  currentChar = 0;
  nextChar = 0;
  nextNextChar = 0;
  currentPos = 0;
  backSlashed = false;

  // Done 3 times to initialise nextchars and currentChar
  readChar();
  readChar();
  readChar();

}

void REPRLoader::parseExtraParams(CppObject& obj) {

  string pName;
  string pValue;
  if (current() == ',') {
    jumpSep(',');
    while (!endOf(')')) {
      readWord(pName);
      jumpSep('=');
      parseParamValue(pValue);
      obj.addField(pName, pValue);
      if (current() == ',') jumpSep(',');
    }
  }

}

void REPRLoader::parseParam(const std::string& name,CppObject& obj) {
  string value;
  parseParamValue(value);
  obj.addField(name,value);
}

void REPRLoader::parseIdentity(CppObject& obj) {

  jumpSep('(');
  parseParam("Name",obj);
  obj.addField("PassMethod","IdentityPass");
  parseExtraParams(obj);
  jumpSep(')');

}

void REPRLoader::parseDrift(CppObject& obj) {

  jumpSep('(');
  parseParam("Name",obj);
  jumpSep(',');
  parseParam("Length",obj);
  obj.addField("PassMethod","DriftPass");
  parseExtraParams(obj);
  jumpSep(')');

}

void REPRLoader::parseDipole(CppObject& obj) {

  string K;

  jumpSep('(');
  parseParam("Name",obj);
  jumpSep(',');
  parseParam("Length",obj);
  jumpSep(',');
  parseParam("BendingAngle",obj);
  jumpSep(',');
  readWord(K);
  obj.addField("PolynomA","(2),0.0,0.0");
  obj.addField("PolynomB","(2),0.0," + K);
  obj.addField("MaxOrder","1");
  obj.addField("PassMethod","BndMPoleSymplectic4Pass");
  parseExtraParams(obj);
  jumpSep(')');

}

void REPRLoader::parseQuadrupole(CppObject& obj) {

  string K;

  jumpSep('(');
  parseParam("Name",obj);
  jumpSep(',');
  parseParam("Length",obj);
  jumpSep(',');
  readWord(K);
  obj.addField("PolynomA","(2),0.0,0.0");
  obj.addField("PolynomB","(2),0.0," + K);
  obj.addField("MaxOrder","1");
  obj.addField("PassMethod","StrMPoleSymplectic4Pass");
  parseExtraParams(obj);
  jumpSep(')');

}

void REPRLoader::parseSextupole(CppObject& obj) {

  string K;

  jumpSep('(');
  parseParam("Name",obj);
  jumpSep(',');
  parseParam("Length",obj);
  jumpSep(',');
  readWord(K);
  obj.addField("PolynomA","(3),0.0,0.0,0.0");
  obj.addField("PolynomB","(3),0.0,0.0," + K);
  obj.addField("MaxOrder","2");
  obj.addField("PassMethod","StrMPoleSymplectic4Pass");
  parseExtraParams(obj);
  jumpSep(')');

}

void REPRLoader::parseMultipole(CppObject& obj) {

  jumpSep('(');
  parseParam("Name",obj);
  jumpSep(',');
  parseParam("Length",obj);
  jumpSep(',');

  vector<int64_t> shapeA;
  string polyA;
  parseParamValue(polyA,&shapeA);
  jumpSep(',');
  vector<int64_t> shapeB;
  string polyB;
  parseParamValue(polyB,&shapeB);

  if(shapeA.size()!=1 || shapeA.size() != shapeB.size() || shapeA[0]!=shapeB[0])
    throw string("Unexpected PolynomA " + AbstractInterface::getShapeStr(shapeA) +
                  "or PolynomB shape " + AbstractInterface::getShapeStr(shapeB) );

  obj.addField("PolynomA",polyA);
  obj.addField("PolynomB",polyB);
  obj.addField("MaxOrder", to_string(shapeA[0]-1));
  obj.addField("PassMethod","StrMPoleSymplectic4Pass");
  parseExtraParams(obj);
  jumpSep(')');

}

void REPRLoader::parseRFCavity(CppObject& obj) {

  jumpSep('(');
  parseParam("Name",obj);
  jumpSep(',');
  parseParam("Length",obj);
  jumpSep(',');
  parseParam("Voltage",obj);
  jumpSep(',');
  parseParam("Frequency",obj);
  jumpSep(',');
  parseParam("HarmonicNumber",obj);
  jumpSep(',');
  parseParam("Energy",obj);
  obj.addField("PassMethod","RFCavityPass");
  parseExtraParams(obj);
  jumpSep(')');

}

void REPRLoader::parseREPR(std::vector<CppObject> &elems) {

  elems.clear();

  // Parse global param
  jumpSep('{');
  while (!endOf('}')) {
    string pName;
    string pValue;
    readWord(pName);
    jumpSep(':');
    parseParamValue(pValue);
    globals.addField(pName,pValue);
    if(current() == ',') jumpSep(',');
  }
  jumpSep('}');

  // Parse elements
  while(!eof()) {

    string elemType;
    size_t pos = getPosMarker();
    readWord(elemType);
    if( !elemType.empty() ) {

        CppObject obj;

        if (elemType == "Marker") {
            parseIdentity(obj);
        } else if (elemType == "Monitor") {
            parseIdentity(obj);
        } else if (elemType == "Drift") {
            parseDrift(obj);
        } else if (elemType == "Dipole") {
            parseDipole(obj);
        } else if (elemType == "Quadrupole") {
            parseQuadrupole(obj);
        } else if (elemType == "Sextupole") {
            parseSextupole(obj);
        } else if (elemType == "Multipole") {
            parseMultipole(obj);
        } else if (elemType == "RFCavity") {
            parseRFCavity(obj);
        } else {
            throw (fileName + ": Unsupported element " + elemType + " at " + getCoord((int) pos));
        }
        elems.push_back(obj);

    }

  }

}

void REPRLoader::getArrayValue(std::string& str,std::vector<int64_t>& shape,std::vector<std::string>& values) {

  str.clear();
  str.append(AbstractInterface::getShapeStr(shape));
  str.append(",");
  for(size_t i=0;i<values.size();i++) {
    str.append(values[i]);
    if(i<values.size()-1) str.append(",");
  }

}

void REPRLoader::parseArrayType(std::string& typeStr) {
  if(current()==',') {
    // Jump dtype
    jumpSep(',');
    jumpSep("dtype");
    jumpSep('=');
    readWord(typeStr);
  }
}

void REPRLoader::parseParamValue(std::string& value,std::vector<int64_t> *shapePtr) {

  readWord(value);

  int pos = getPosMarker();

  // Array type
  if(value=="array" && current()=='(') {

    jumpSep('(');
    vector<int64_t> shape;
    vector<string> values;
    try {
      parseArray(values, shape, 0);
    } catch (string& errStr) {
      throw string(fileName + ": " + errStr + " at " + getErrorLocation(pos));
    }
    if(shapePtr) *shapePtr = shape;
    getArrayValue(value,shape,values);
    string typeStr;
    parseArrayType(typeStr);
    jumpSep(')');

  // Particle type
  } else if (value=="Particle" && current()=='(') {
    jumpSep('(');
    readWord(value);
    jumpSep(')');
  }

}

// Return next significant char
char REPRLoader::current() {
  jumpSpace();
  return currentChar;
}

// Get a char from the input string
char REPRLoader::read() {
  if (currentPos < fileBuffer.size())
    return fileBuffer[currentPos++];
  else
    return 0;
}

// Go to next char
void REPRLoader::toNextChar() {

  char c = read();
  currentChar = nextChar;
  nextChar = nextNextChar;
  nextNextChar = c;

}

// Read a char and handle backslash
void REPRLoader::readChar() {

  backSlashed = false;

  toNextChar();

  /* Escape sequence */
  if (currentChar == '\\') {
    switch (nextChar) {
      case '"':
        toNextChar();        // return '"'
        backSlashed = true;  // For string handling
        break;
      case 'n':
        toNextChar();
        currentChar = '\n';  // return '\n'
        break;
      case 'r':
        toNextChar();
        currentChar = '\r';  // return '\r'
        break;
      case 't':
        toNextChar();
        currentChar = '\t'; // return '\t'
        break;
      case '\\':
        toNextChar(); // return '\'
        break;
    }
  }

}

// Jump spaces
void REPRLoader::jumpSpace() {
  if(currentChar>32)
    return;
  while (currentChar <= 32 && currentChar > 0) readChar();
}

// Special character
bool REPRLoader::isSpecialChar(char c) {
  return c == '=' || c == ':' || c == ',' || c == '[' || c == ']' || c == '(' || c == ')' || c == '{' || c == '}';
}

// Read a word
void REPRLoader::readWord(std::string& word) {

  word.clear();

  /* Jump space and comments */
  jumpSpace();

  /* Treat special character */
  if (isSpecialChar(currentChar)) {
    word.push_back(currentChar);
    readChar();
    return;
  }

  /* Treat char */
  if (currentChar == '\'' && nextNextChar == '\'') {
    word.push_back(currentChar);
    readChar();
    word.push_back(currentChar);
    readChar();
    word.push_back(currentChar);
    readChar();
    return;
  }

  /* Treat string */
  if ((currentChar == '\'' || currentChar == '"') && !backSlashed) {
    readChar();
    while (((currentChar != '\'' && currentChar != '"') || backSlashed) && currentChar != 0) {
      word.push_back(currentChar);
      readChar();
    }
    if (currentChar == 0) {
      throw string("Unterminated string");
    }
    readChar();
    return;
  }

  /* Treat other word */
  while (currentChar > 32 && !isSpecialChar(currentChar)) {
    word.push_back(currentChar);
    readChar();
  }

}

int REPRLoader::getPosMarker() const {
  return (int)currentPos-((currentChar==0)?0:1)-((nextNextChar==0)?0:1)-((nextNextChar==0)?0:1);
}

string REPRLoader::getErrorLocation(int pos) {
  return "at " + getCoord(pos) + " in " + fileName;
}

// Jump a separator
void REPRLoader::jumpSep(const string& sep) {
  jumpSpace();
  size_t pos = getPosMarker();
  string w;
  readWord(w);
  if (w!=sep)
    throw string("'" + sep + "' expected but got '" + w + "' " + getErrorLocation((int)pos));
}

void REPRLoader::jumpSep(const char sep) {
  jumpSpace();
  size_t pos = getPosMarker();
  string w;
  readWord(w);
  if (!(w.length()==1 && sep==w[0])) {
    string sepAsStr;
    sepAsStr.push_back(sep);
    throw string("'" + sepAsStr + "' expected but got '" + w + "' " + getErrorLocation((int)pos));
  }
}

// Recursively parse array
void REPRLoader::parseArray(vector<string> &ret,vector<int64_t>& shape,int level) {

  int64_t count = 0;
  jumpSep('[');
  while (!endOf(']')) {
    if(current()=='[') {
      parseArray(ret, shape, level+1);
    } else {
      string w;
      readWord(w);
      if(w=="array" && current()=='(') {
        jumpSep('(');
        vector<int64_t> inShape;
        vector<string> inValues;
        parseArray(inValues,inShape,0);
        string typeStr;
        parseArrayType(typeStr);
        jumpSep(')');
        string inValue;
        getArrayValue(inValue,inShape,inValues);
        ret.push_back(inValue);
      } else {
        ret.push_back(w);
      }
    }
    if (current() == ',') jumpSep(',');
    count++;
  }
  jumpSep(']');

  // Update shape
  if(level>=shape.size()) {
    shape.resize(level+1,-1);
    shape[level] = count;
  } else {
    if( shape[level]==-1 ) {
      shape[level] = count;
    } else if (shape[level] != count) {
      throw string("Malformed array");
    }
  }

}

std::string REPRLoader::getCoord(int pos) {

  // Find line and offset of a position
  size_t p = 0;
  size_t line = 1;
  size_t column = 1;

  if(pos>=(int)fileBuffer.length()-1)
    pos = (int)fileBuffer.length()-1;

  while(p<(size_t)pos) {
    if(fileBuffer[p]=='\n') {
      line++;
      column = 1;
    } else {
      column++;
    }
    p++;
  }

  return to_string(line)+":"+to_string(column);

}

bool REPRLoader::endOf(char c) {
  jumpSpace();
  return (currentChar==c) || (currentChar==0);

}

bool REPRLoader::eof() {
  jumpSpace();
  return (currentChar==0);
}
