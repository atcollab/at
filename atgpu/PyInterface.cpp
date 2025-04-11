#include "PyInterface.h"
#include <string>
#include <Python.h>
#include "AbstractGPU.h"
#include "Lattice.h"
#include "iostream"

using namespace std;

// --------------------------------------------------------------------------------------------------------------------

// Particle type
PyObject *particle_type;

// Methods declaration
static PyObject *at_gpuinfo(PyObject *self, PyObject *args, PyObject *kwargs);
static PyObject *at_gpupass(PyObject *self, PyObject *args, PyObject *kwargs);

// Method table
static PyMethodDef AtGPUMethods[] = {
        {"gpuinfo",  (PyCFunction)at_gpuinfo, METH_VARARGS | METH_KEYWORDS,
                PyDoc_STR("gpuinfo()\n\n"
                          "Return lists of GPU present on the system.\n"
                          "[name,compute capability,compute unit number,platform]\n"
                )},
        {"gpupass",  (PyCFunction)at_gpupass, METH_VARARGS | METH_KEYWORDS,
                PyDoc_STR("gpupass(line: Sequence[Element], r_in, n_turns: int, refpts: Uint32_refs = [], "
                          "reuse: Optional[bool] = False, omp_num_threads: Optional[int] = 0)\n\n"
                          "Track input particles r_in along line for nturns turns.\n"
                          "Record 6D phase space at elements corresponding to refpts for each turn.\n\n"
                          "Parameters:\n"
                          "    line:    list of elements\n"
                          "    rin:     6 x n_particles Fortran-ordered numpy array.\n"
                          "      On return, rin contains the final coordinates of the particles\n"
                          "    n_turns: number of turns to be tracked\n"
                          "    refpts:  numpy array of indices of elements where output is desired\n"
                          "       0 means entrance of the first element\n"
                          "       len(line) means end of the last element\n"
                          "    energy:  nominal energy [eV]\n"
                          "    particle (Optional[Particle]):  circulating particle\n"
                          "    reuse:   if True, use previously cached description of the lattice.\n"
                          "    losses:  if True, process losses\n"
                          "    gpu_pool:  List of GPU id to use\n"
                          "    tracking_starts: numpy array of indices of elements where tracking should start.\n"
                          "       len(tracking_start) must divide the number of particle and it gives the stride size.\n"
                          "       The i-th particle of rin starts at elem tracking_start[i/stride].\n"
                          "       The behavior is similar to lattice.rotate(tracking_starts[i/stride]).\n"
                          "       The stride size should be multiple of 64 for best performance.\n"
                          "    integrator: Type of integrator to use.\n"
                          "       1: Euler 1st order, 1 drift/1 kick per step.\n"
                          "       2: Mclachlan 2nd order, 2 drift/2 kicks per step.\n"
                          "       3: Ruth 3rd order, 3 drifts/3 kicks per step.\n"
                          "       4: Forest/Ruth 4th order, 4 drifts/3 kicks per step (Default).\n"
                          "       5: Optimal 4th order from R. Mclachlan, 4 drifts/4 kicks per step.\n"
                          "       6: Yoshida 6th order, 8 drifts/7 kicks per step.\n\n"
                          "Returns:\n"
                          "    rout:    6 x n_particles x n_refpts x n_turns Fortran-ordered numpy array\n"
                          "         of particle coordinates\n\n"
                          ":meta private:"
                )},
        {NULL, NULL, 0, NULL} // Sentinel
};

PyMODINIT_FUNC PyInit_gpu(void) {

  static struct PyModuleDef moduledef = {
          PyModuleDef_HEAD_INIT,
          "gpu",         /* m_name */
          PyDoc_STR("GPU handler for Accelerator Toolbox"),      /* m_doc */
          -1,           /* m_size */
          AtGPUMethods, /* m_methods */
          NULL,         /* m_slots */
          NULL,         /* m_traverse */
          NULL,         /* m_clear */
          NULL,         /* m_free */
  };

  PyObject *m = PyModule_Create(&moduledef);

  // Import numpy array api
  if( PyArray_API == nullptr ) {
    if (_import_array() < 0) {
      PyErr_Print();
      return nullptr;
    }
  }

  // Init particle type
  PyObject *module = PyImport_ImportModule("at.lattice");
  if (module) {
    particle_type = PyObject_GetAttrString(module, "Particle");
    Py_DECREF(module);
    if(particle_type == nullptr) {
      cerr << "gpu module initialisation error, cannot get at.lattice.Particle type" << endl;
      return nullptr;
    }
  } else {
    PyErr_Print();
    return nullptr;
  }

  // Init Interface handler
  AbstractInterface::setHandler(new PyInterface());

  return m;
}

static PyObject *at_gpuinfo(PyObject *self, PyObject *args, PyObject *kwargs) {

  try {
    vector<GPU_INFO> gpuInfos;
    gpuInfos = AbstractGPU::getInstance()->getDeviceList();
    PyObject *infos = PyList_New(gpuInfos.size());
    for(int i=0;i<gpuInfos.size();i++) {
      PyObject *info = PyList_New(4);
      PyList_SetItem(info, 0, PyUnicode_FromString(gpuInfos[i].name.c_str()));
      PyList_SetItem(info, 1, PyUnicode_FromString(gpuInfos[i].version.c_str()));
      PyList_SetItem(info, 2, PyLong_FromLong(gpuInfos[i].mpNumber));
      PyList_SetItem(info, 3, PyUnicode_FromString(gpuInfos[i].platform.c_str()));
      PyList_SetItem(infos, i, info);
    }
    return infos;
  } catch (string& errorStr) {
    cout << "at_gpuinfo(): " << errorStr << endl;
    return PyList_New(0);
  }

}

static PyObject *at_gpupass(PyObject *self, PyObject *args, PyObject *kwargs) {

  // Default symplectic integrator (4th order)
  static SymplecticIntegrator integrator(4);
  // Lattice object
  static Lattice *gpuLattice = nullptr;

  static const char *kwlist[] = {"line","rin","nturns","refpts","turn",
                                 "energy", "particle", "keep_counter",
                                 "reuse","losses",
                                 "bunch_spos", "bunch_currents", "gpu_pool",
                                 "tracking_starts","integrator","verbose",
                                 nullptr};

  NPY_TYPES floatType;
  string floatTypeStr;
  if(sizeof(AT_FLOAT)==8) {
    floatType = NPY_DOUBLE;
    floatTypeStr = "NPY_DOUBLE";
  } else if(sizeof(AT_FLOAT)==4) {
    floatType = NPY_FLOAT;
    floatTypeStr = "NPY_FLOAT";
  } else
    return PyErr_Format(PyExc_ValueError, "PyAT GPU is compatible only with NPY_FLOAT or NPY_DOUBLE");

  PyObject *lattice;
  PyObject *particle;
  PyObject *energy;
  PyArrayObject *rin;
  PyArrayObject *refs = nullptr;
  PyArrayObject *trackstarts = nullptr;
  PyArrayObject *bcurrents;
  PyArrayObject *bspos;
  PyObject *gpupool = nullptr;
  int num_turns;
  int keep_lattice=0;
  int keep_counter=0;
  int verbose=0;
  int counter=0;
  int losses=0;
  int integratorType=4;
  double t0,t1;

  // Get input args
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!i|O!$iO!O!pppO!O!O!O!ip", const_cast<char **>(kwlist),
                                   &PyList_Type, &lattice,
                                   &PyArray_Type, &rin,
                                   &num_turns,
                                   &PyArray_Type, &refs,
                                   &counter,
                                   &PyFloat_Type, &energy,
                                   particle_type, &particle,
                                   &keep_counter,
                                   &keep_lattice,
                                   &losses,
                                   &PyArray_Type, &bspos,
                                   &PyArray_Type, &bcurrents,
                                   &PyList_Type, &gpupool,
                                   &PyArray_Type, &trackstarts,
                                   &integratorType,
                                   &verbose
                                   )) {
    return nullptr;
  }

  // Input particles
  if (PyArray_DIM(rin,0) != 6) {
    return PyErr_Format(PyExc_ValueError, "rin is not 6D");
  }
  if (PyArray_TYPE(rin) != floatType) {
    return PyErr_Format(PyExc_ValueError, ("rin is not a "+floatTypeStr+" array").c_str());
  }
  if ((PyArray_FLAGS(rin) & NPY_ARRAY_FARRAY_RO) != NPY_ARRAY_FARRAY_RO) {
    return PyErr_Format(PyExc_ValueError, "rin is not Fortran-aligned");
  }

  uint32_t num_particles = (PyArray_SIZE(rin)/6);
  AT_FLOAT *drin = (AT_FLOAT *)PyArray_DATA(rin);

  // Reference points
  uint32_t *ref_pts;
  uint32_t num_refs;

  if (refs) {
    if (PyArray_TYPE(refs) != NPY_UINT32) {
      return PyErr_Format(PyExc_ValueError, "refpts is not a numpy uint32 array");
    }
    ref_pts = (uint32_t *)PyArray_DATA(refs);
    num_refs = PyArray_SIZE(refs);
  } else {
    ref_pts = nullptr;
    num_refs = 0;
  }

  // Starting elements
  uint32_t *track_starts = nullptr;
  uint32_t num_starts = 0;

  if (trackstarts) {
    if (PyArray_TYPE(refs) != NPY_UINT32) {
      return PyErr_Format(PyExc_ValueError, "tracking_starts is not a numpy uint32 array");
    }
    track_starts = (uint32_t *)PyArray_DATA(trackstarts);
    num_starts = PyArray_SIZE(trackstarts);
    if( num_particles % num_starts != 0 ) {
      return PyErr_Format(PyExc_ValueError, "len(tracking_starts) must divide number of particle");
    }
  }

  // GPU
  int gpuId = 0;
  if( gpupool ) {
    size_t nGPU = PyList_Size(gpupool);
    if(nGPU==0) {
      string err = "at_gpupass() gpu_pool is empty";
      return PyErr_Format(PyExc_RuntimeError, err.c_str());
    } else if(nGPU>1) {
      string err = "at_gpupass() Multi GPU support not implemented";
      return PyErr_Format(PyExc_RuntimeError, err.c_str());
    }
    gpuId = (int)PyLong_AsLong(PyList_GET_ITEM(gpupool, 0));
  }

  // Integrator
  if( integratorType!=integrator.getType() ) {
    if( verbose && keep_lattice )
      cout << "Warning, lattice is recreated when integrator type is changed" << endl;
    delete gpuLattice;
    gpuLattice = nullptr;
    integrator.setType(integratorType);
  }

  // Set up lattice and run tracking
  if( !keep_lattice ) {
    delete gpuLattice;
    gpuLattice = nullptr;
  }

  if( !gpuLattice ) {

    // Create the GPU lattice
    try {

      PyInterface *pyI = (PyInterface *) AbstractInterface::getInstance();
      size_t nElements = PyList_Size(lattice);
      gpuLattice = new Lattice(nElements,integrator, 0.0, gpuId);
      for (size_t i = 0; i < nElements; i++) {
        PyObject *elem = PyList_GET_ITEM(lattice, i);
        pyI->setObject(elem);
        gpuLattice->addElement();
      }
      gpuLattice->generateGPUKernel();

      if( verbose )
        cout << "Lattice created on " << gpuLattice->getGPUContext()->name() << " [" << nElements << " elements]" << endl;

    } catch (string& errStr) {
      delete gpuLattice;
      gpuLattice = nullptr;
      string err =  "at_gpupass() build lattice failed: " + errStr;
      return PyErr_Format(PyExc_RuntimeError,err.c_str());
    }

  }

  // Load lattice on the GPU
  try {
    uint64_t size = gpuLattice->fillGPUMemory();
    if( verbose )
      cout << "Lattice successfully loaded on " << gpuLattice->getGPUContext()->name() << " [" << size << " bytes]" << endl;
  } catch (string& errStr) {
    string err =  "at_gpupass() fill GPU memory failed: " + errStr;
    return PyErr_Format(PyExc_RuntimeError,err.c_str());
  }

  // Turn counter
  if( !keep_counter )
    gpuLattice->setTurnCounter(counter);

  try {

    if( verbose )
      cout << "Tracking " << num_particles << " particles for " << num_turns << " turns on " << gpuLattice->getGPUContext()->name() << " #" << gpuId << endl;

    npy_intp outdims[4] = {6,(npy_intp)(num_particles),num_refs,num_turns};
    PyObject *rout = PyArray_EMPTY(4, outdims, floatType, 1);
    if( rout==nullptr )
      return PyErr_Format(PyExc_RuntimeError,"Not enough memory while trying to allocate particle output coordinates");
    AT_FLOAT *drout = (AT_FLOAT *)PyArray_DATA((PyArrayObject *)rout);

    if(losses) {

      npy_intp pdims[] = {(npy_intp)(num_particles)};
      npy_intp lxdims[] = {6,(npy_intp)(num_particles)};
      PyObject *xnturn = PyArray_EMPTY(1, pdims, NPY_UINT32, 1);
      PyObject *xnelem = PyArray_EMPTY(1, pdims, NPY_UINT32, 1);
      PyObject *xlost = PyArray_EMPTY(1, pdims, NPY_BOOL, 1);
      PyObject *xlostcoord = PyArray_EMPTY(2, lxdims, floatType, 1);
      uint32_t *xnturnPtr = (uint32_t *)PyArray_DATA((PyArrayObject *)xnturn);
      uint32_t *xnelemPtr = (uint32_t *)PyArray_DATA((PyArrayObject *)xnelem);
      bool *xlostPtr = (bool *)PyArray_DATA((PyArrayObject *)xlost);
      AT_FLOAT *xlostcoordPtr = (AT_FLOAT *)PyArray_DATA((PyArrayObject *)xlostcoord);

      gpuLattice->run(num_turns,num_particles,drin,drout,num_refs,ref_pts,num_starts,track_starts,
                      xnturnPtr,xnelemPtr,xlostcoordPtr,true);

      // Format result for AT
      for(uint32_t i=0;i<num_particles;i++) {
        xlostPtr[i] = (xnturnPtr[i] != num_turns);
        if(!xlostPtr[i]) xnturnPtr[i] = 0;
      }

      PyObject *tout = PyTuple_New(2);
      PyObject *dict = PyDict_New();
      PyDict_SetItemString(dict,(char *)"islost",(PyObject *)xlost);
      PyDict_SetItemString(dict,(char *)"turn",(PyObject *)xnturn);
      PyDict_SetItemString(dict,(char *)"elem",(PyObject *)xnelem);
      PyDict_SetItemString(dict,(char *)"coord",(PyObject *)xlostcoord);
      PyTuple_SetItem(tout, 0, rout);
      PyTuple_SetItem(tout, 1, dict);
      Py_DECREF(xlost);
      Py_DECREF(xnturn);
      Py_DECREF(xnelem);
      Py_DECREF(xlostcoord);
      return tout;

    } else {

      gpuLattice->run(num_turns,num_particles,drin,drout,num_refs,ref_pts,num_starts,track_starts,
                      nullptr,nullptr,nullptr,true);
      return rout;

    }

  } catch (string& errStr) {
    string err =  "at_gpupass() run failed: " + errStr;
    return PyErr_Format(PyExc_RuntimeError,err.c_str());
  }

}

// --------------------------------------------------------------------------------------------------------------------

void PyInterface::setObject(PyObject *obj) {
  self = obj;
}

int PyInterface::getInt(const std::string& name) {

  PyObject *attr = PyObject_GetAttrString(self, name.c_str());
  if (!attr) {
    if (PyErr_Occurred()) PyErr_Clear(); // Reset python error
    throw (name + " attribute not found");
  }
  Py_DECREF(attr);
  return (int)PyLong_AsLong(attr);

}

std::string PyInterface::getString(const std::string& name) {

  PyObject *attr = PyObject_GetAttrString(self, name.c_str());
  if (!attr) {
    if (PyErr_Occurred()) PyErr_Clear(); // Reset python error
    throw (name + " attribute not found");
  }
  Py_DECREF(attr);
  return PyUnicode_AsUTF8(attr);

}


double PyInterface::getDouble(const std::string& name) {

  PyObject *attr = PyObject_GetAttrString(self, name.c_str());
  if (!attr) {
    if (PyErr_Occurred()) PyErr_Clear(); // Reset python error
    throw string(name + " attribute not found");
  }


  Py_DECREF(attr);
  return (AT_FLOAT)PyFloat_AsDouble(attr);

}

double *PyInterface::getNativeDoubleArray(const std::string& name,std::vector<int64_t>& shape) {

  PyArrayObject *array = (PyArrayObject *) PyObject_GetAttrString(self, name.c_str());
  if (array == nullptr) {
    if (PyErr_Occurred()) PyErr_Clear(); // Reset python error
    throw string(name + " array attribute not found");
  }

  if ((PyArray_FLAGS(array) & NPY_ARRAY_FARRAY_RO) != NPY_ARRAY_FARRAY_RO) {
    Py_DECREF(array);
    throw string(name + " array attribute is not Fortran-aligned");
  }

  if (PyArray_TYPE(array) != NPY_DOUBLE) {
    Py_DECREF(array);
    throw string(name + " array attribute is not a double array");
  }

  size_t nDim = PyArray_NDIM(array);
//int64_t *dims = PyArray_SHAPE(array);
  npy_intp *dims = PyArray_SHAPE(array);

  shape.resize(nDim);
  uint32_t nbItem = 1;
  for(int i=0;i<nDim;i++) {
    shape[i] = dims[i];
    nbItem *= dims[i];
  }

  // Map python memory
  double *pyPtr = (double *) PyArray_DATA(array);
  Py_DECREF(array);
  return (double *)pyPtr;

}

float *PyInterface::getNativeFloatArray(const std::string& name,std::vector<int64_t>& shape) {

  PyArrayObject *array = (PyArrayObject *) PyObject_GetAttrString(self, name.c_str());
  if (array == nullptr) {
    if (PyErr_Occurred()) PyErr_Clear(); // Reset python error
    throw string(name + " array attribute not found");
  }

  if ((PyArray_FLAGS(array) & NPY_ARRAY_FARRAY_RO) != NPY_ARRAY_FARRAY_RO) {
    Py_DECREF(array);
    throw string(name + " array attribute is not Fortran-aligned");
  }

  if (PyArray_TYPE(array) != NPY_FLOAT) {
    Py_DECREF(array);
    throw string(name + " array attribute is not a float array");
  }

  size_t nDim = PyArray_NDIM(array);
//int64_t *dims = PyArray_SHAPE(array);
  npy_intp *dims = PyArray_SHAPE(array);

  shape.resize(nDim);
  uint32_t nbItem = 1;
  for(int i=0;i<nDim;i++) {
    shape[i] = dims[i];
    nbItem *= dims[i];
  }

  // Map python memory
  float *pyPtr = (float *) PyArray_DATA(array);
  Py_DECREF(array);
  return (float *)pyPtr;

}

