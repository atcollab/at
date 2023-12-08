#include "CudaGPU.h"
#include "iostream"

using namespace std;

AbstractGPU *AbstractGPU::gpuHandler = nullptr;

AbstractGPU *AbstractGPU::getInstance() {
  if( gpuHandler== nullptr )
    gpuHandler = new CudaGPU();
  return gpuHandler;
}

void AbstractGPU::outputCode(std::string& code) {

  vector<string> lines;
  split(lines,code,'\n');
  for(size_t i=0;i<lines.size();i++) {
    char tmp[256];
    sprintf(tmp,"%04d: ",(int)(i+1));
    cout << tmp << lines[i] << endl;
  }

}


void AbstractGPU::split(vector<string> &tokens, const string &text, char sep) {

  size_t start = 0, end = 0;
  tokens.clear();

  while ((end = text.find(sep, start)) != string::npos) {
    tokens.push_back(text.substr(start, end - start));
    start = end + 1;
  }

  tokens.push_back(text.substr(start));

}
