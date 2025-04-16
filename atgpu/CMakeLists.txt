# AT GPU debugging project

set(_CUDA OFF)
set(_OPENCL ON)

cmake_minimum_required(VERSION 3.13)
project(AT_GPU)

set(SOURCE_FILES
        main.cpp
        REPRLoader.cpp
        AbstractGPU.cpp
        AbstractInterface.cpp
        Lattice.cpp
        SymplecticIntegrator.cpp
        PassMethodFactory.cpp
        IdentityPass.cpp
        DriftPass.cpp
        ExactDriftPass.cpp
        ThinMPolePass.cpp
        StrMPoleSymplectic4Pass.cpp
        ExactMultipolePass.cpp
        StrMPoleSymplectic4RadPass.cpp
        ExactMultipoleRadPass.cpp
        BndMPoleSymplectic4Pass.cpp
        BndMPoleSymplectic4RadPass.cpp
        CavityPass.cpp
        RFCavityPass.cpp
)

if(_CUDA)
  set(GPU_CXX_FLAGS "-DCUDA")
  if (MSVC)
    set(LIB ${LIB} cuda.lib nvrtc.lib)
  else()
    set(LIB ${LIB} -lnvrtc -lcuda)
  endif (MSVC)
  set(SOURCE_FILES ${SOURCE_FILES} CudaGPU.cpp)
endif()

if(_OPENCL)
  set(GPU_CXX_FLAGS "-DOPENCL")
  if (MSVC)
    set(LIB OpenCL.lib)
  else()
    set(LIB -lOpenCL)
  endif (MSVC)
  set(SOURCE_FILES ${SOURCE_FILES} OpenCLGPU.cpp)
endif()

message(GPU="${GPU_CXX_FLAGS}")

set(CMAKE_CXX_STANDARD 11)

# Add executable target with source files listed in SOURCE_FILES variable
if(MSVC)
  set(CMAKE_CXX_FLAGS "/EHsc ${GPU_CXX_FLAGS}")
  set(CUDA_PATH "C:\\Program Files\\NVIDIA GPU Computing Toolkit\\CUDA\\v12.3")
  include_directories($ENV{CUDA_PATH}/include)
  link_directories($ENV{CUDA_PATH}/lib/x64)
else()
  set(CMAKE_CXX_FLAGS "-g ${GPU_CXX_FLAGS}")
  set(SOURCE_FILES ${SOURCE_FILES} PyInterface.cpp)
  set(PYTHON_INCL /usr/include/python3.8)
  set(PYTHON_LIB /usr/lib/python3.8)
  set(LIB ${LIB} -lpython3.8)
  include_directories($ENV{CUDA_PATH}/include ${PYTHON_INCL})
  link_directories($ENV{CUDA_PATH}/lib64 ${PYTHON_LIB})
endif (MSVC)

link_libraries(${LIB})
add_executable(atgpu ${SOURCE_FILES})
