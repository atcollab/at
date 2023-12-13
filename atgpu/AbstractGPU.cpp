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

// Add math function
void AbstractGPU::addUtilsFunctions(std::string &code) {

  string ftype;
  getDeviceFunctionQualifier(ftype);
  if(!ftype.empty()) ftype.append(" ");

  // Math constants
  code.append(
          "#define INF   __longlong_as_double(0x7ff0000000000000ULL)\n"
          "#define NAN   __longlong_as_double(0xfff8000000000000ULL)\n"
  );

  // 6D transfrom
  code.append(
          ftype +
          "void translate6(AT_FLOAT* r,AT_FLOAT *t) {\n"
          "  r[0] += t[0];  r[1] += t[1];  r[2] += t[2];\n"
          "  r[3] += t[3];  r[4] += t[4];  r[5] += t[5];\n"
          "}\n"
  );
  code.append(
          ftype +
          "void transform66(AT_FLOAT* r,AT_FLOAT *M) {\n"
          "  int i,j;\n"
          "  AT_FLOAT sum[6];\n"
          "  for(i=0;i<6;i++)\n"
          "  {\n"
          "    sum[i]=0;\n"
          "    for(j=0;j<6;j++)\n"
          "      sum[i]+=M[i+j*6]*r[j];\n"
          "  }\n"
          "  for(i=0;i<6;i++)\n"
          "    r[i]=sum[i];\n"
          "}\n"
  );

  //Drift (small angle)
  code.append(
          ftype +
          "void fastdrift(AT_FLOAT* r,AT_FLOAT NormL,AT_FLOAT p_norm) {\n"
          "  r[0] += NormL * r[1];\n"
          "  r[2] += NormL * r[3];\n"
          "  r[5] += p_norm * NormL * (r[1] * r[1] + r[3] * r[3]) * 0.5;\n"
          "}\n"
  );

  // Bending fringe field correction (edge angle focusing)
  code.append(
          ftype +
    "void edge_fringe(AT_FLOAT* r,AT_FLOAT p_norm,AT_FLOAT edge_angle,AT_FLOAT irho,AT_FLOAT f_corrx,AT_FLOAT f_corry,uint32_t method) {\n"
    "  AT_FLOAT fy;\n"
    "  switch(method) {\n"
    "  case 0:\n"
    "    fy = f_corry;\n"
    "    break;\n"
    "  case 1:\n"
    "    fy = irho * tan(edge_angle - f_corry * p_norm);\n"
    "    break;\n"
    "  case 2:\n"
    "    fy = irho * tan(edge_angle - f_corry * p_norm) * p_norm;\n"
    "    break;\n"
    "  case 3:\n"
    "    fy = irho * tan(edge_angle - f_corry + r[1] * p_norm);\n"
    "    break;\n"
    "  }\n"
    "  r[1] += r[0] * f_corrx;\n"
    "  r[3] -= r[2] * fy;\n"
    "}\n"
  );


  // Kick (No rad)

  // Quad
  code.append(
          ftype +
          "void strthinkick1(AT_FLOAT* r,const AT_FLOAT* A,const AT_FLOAT* B,AT_FLOAT L,int max_order) {\n"
          "  r[1] -= B[1] * L * r[0] + B[0];\n"
          "  r[3] += B[1] * L * r[2] + A[0];\n"
          "}\n"
  );

  // Sextu
  code.append(
          ftype +
          "void strthinkick2(AT_FLOAT* r,const AT_FLOAT* A,const AT_FLOAT* B,AT_FLOAT L,int max_order) {\n"
          "  r[1] -= B[2] * L * (r[0]*r[0]-r[2]*r[2]) + B[0];\n"
          "  r[3] += B[2] * L * (2.0* r[0] * r[2]) + A[0];\n"
          "}\n"
  );

  // Octu
  code.append(
          ftype +
          "void strthinkick3(AT_FLOAT* r,const AT_FLOAT* A,const AT_FLOAT* B,AT_FLOAT L,int max_order) {\n"
          "  AT_FLOAT x2 = r[0]*r[0];\n"
          "  AT_FLOAT y2 = r[2]*r[2];\n"
          "  r[1] -= B[3] * L * r[0] * (x2 - 3.0*y2) + B[0];\n"
          "  r[3] += B[3] * L * r[2] * (3.0*x2 - y2) + A[0];\n"
          "}\n"
  );

  // Recursively calculate the local transverse magnetic field
  string polyLoop =
          "  int i;\n"
          "  AT_FLOAT ReSum = B[max_order];\n"
          "  AT_FLOAT ImSum = A[max_order];\n"
          "  AT_FLOAT ReSumTemp;\n"
          "  for(i = max_order - 1; i >= 0; i--) {\n"
          "    ReSumTemp = ReSum * r[0] - ImSum * r[2] + B[i];\n"
          "    ImSum = ImSum * r[0] + ReSum * r[2] + A[i];\n"
          "    ReSum = ReSumTemp;\n"
          "  }\n";

  code.append(
          ftype +
          "void strthinkick(AT_FLOAT* r,const AT_FLOAT* A,const AT_FLOAT* B,AT_FLOAT L,int max_order) {\n"
          + polyLoop +
          "  r[1] -= L * ReSum;\n"
          "  r[3] += L * ImSum;\n"
          "}\n"
  );

  code.append(
          ftype +
          "void bndthinkick(AT_FLOAT* r,const AT_FLOAT* A,const AT_FLOAT* B,AT_FLOAT L,int max_order,AT_FLOAT irho) {\n"
          + polyLoop +
          "  r[1] -= L * (ReSum - (r[4] - r[0] * irho) * irho);\n"
          "  r[3] += L * ImSum;\n"
          "  r[5] += L * irho * r[0];\n"
          "}\n"
  );

  /*
  code.append(
          ftype +
          "void bndthinkick0(AT_FLOAT* r,const AT_FLOAT* A,const AT_FLOAT* B,AT_FLOAT L,AT_FLOAT irho) {\n"
          "  r[1] -= L * (B[0] - (r[4] - r[0] * irho) * irho);\n"
          "  r[3] += L * A[0];\n"
          "  r[5] += L * irho * r[0];\n"
          "}\n"
  );
  */

  //Lee-Whiting's thin lens limit formula as given in p. 390 of "Beam Dynamics..." by E. Forest
  code.append(
          ftype +
          "void quad_fringe(AT_FLOAT* r, AT_FLOAT b2, AT_FLOAT sign, AT_FLOAT p_norm) {\n"
          "  AT_FLOAT u = p_norm * b2 / 12.0;\n"
          "  AT_FLOAT x2 = r[0] * r[0];\n"
          "  AT_FLOAT z2 = r[2] * r[2];\n"
          "  AT_FLOAT xz = r[0] * r[2];\n"
          "  AT_FLOAT gx = u * (x2 + 3 * z2) * r[0];\n"
          "  AT_FLOAT gz = u * (z2 + 3 * x2) * r[2];\n"
          "  AT_FLOAT r1tmp = 0;\n"
          "  AT_FLOAT r3tmp = 0;\n"
          "  r[0] += sign*gx;\n"
          "  r1tmp = 3 * u * (2 * xz * r[3] - (x2 + z2) * r[1]);\n"
          "  r[2] -= sign*gz;\n"
          "  r3tmp = 3 * u * (2 * xz * r[1] - (x2 + z2) * r[3]);\n"
          "  r[5] -= sign * (gz * r[3] - gx * r[1]) * p_norm;\n"
          "  r[1] += sign*r1tmp;\n"
          "  r[3] -= sign*r3tmp;\n"
          "}\n"
  );

}
