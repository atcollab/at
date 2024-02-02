#include <mex.h>
#include <MatlabInterface.h>
#include <AbstractGPU.h>
#include <Lattice.h>

using namespace std;

#define LATTICE prhs[0]
#define RIN prhs[1]
#define NEWLATTICE prhs[2]
#define NTURNS prhs[3]
#define REFPTS prhs[4]
#define TURN prhs[5]
#define KEEPCOUNTER prhs[6]
#define GPUPOOL prhs[7]
#define INTEGRATOR prhs[8]

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  if (nlhs > 2)
    mexErrMsgIdAndTxt("Atpass:WrongParameter","GPU does not support more than 2 output arguments [Rout,lossinfo]");

  // Init Interface handler
  if( !AbstractInterface::isValidHandler() ) {
    AbstractInterface::setHandler(new MatlabInterface());
  }

  // Default symplectic integrator (4th order)
  static SymplecticIntegrator integrator(4);
  // Lattice object
  static Lattice *gpuLattice = nullptr;

  int num_turns=(int)mxGetScalar(NTURNS);
  int keep_lattice=(mxGetScalar(NEWLATTICE) == 0) ? 0 : 1;
  int keep_counter=(int)mxGetScalar(KEEPCOUNTER);
  int counter=(int)mxGetScalar(TURN);
  int losses=(nlhs == 2);
  int integratorType=(int)mxGetScalar(INTEGRATOR);

  // Input particles
  if (!mxIsDouble(RIN) || mxIsComplex(RIN)) {
    mexErrMsgIdAndTxt("Atpass:WrongType","RIN must be real");
  }
  if (mxGetM(RIN) != 6) {
    mexErrMsgIdAndTxt("Atpass:WrongSize","RIN must have 6 rows");
  }

  uint32_t num_particles =  mxGetN(RIN);
  AT_FLOAT *drin = (AT_FLOAT *)mxGetDoubles(RIN);

  // Reference points
  uint32_t *ref_pts;
  uint32_t num_refs = (uint32_t)mxGetNumberOfElements(REFPTS);

  if( num_refs==0 ) {
    // One ref at the end of the turn
    num_refs = 1;
    ref_pts = (uint32_t *) mxCalloc(num_refs, sizeof(uint32_t));
    ref_pts[0] = mxGetNumberOfElements(LATTICE);
  } else {
    // Convert indices to uint32_t
    mxDouble *dblrefpts = mxGetDoubles(REFPTS);
    ref_pts = (uint32_t *) mxCalloc(num_refs, sizeof(uint32_t));
    for (int i = 0; i < num_refs; i++)
      ref_pts[i] = ((int) dblrefpts[i]) - 1;
  }

  // Starting elements (not supported in Matlab)
  uint32_t *track_starts = nullptr;
  uint32_t num_starts = 0;

  // GPU (pool not supported)
  int gpuId = (int)mxGetScalar(GPUPOOL)-1;

  // Integrator
  if( integratorType!=integrator.getType() ) {
    if( keep_lattice )
      mexWarnMsgIdAndTxt("AT:PassWarning", "Lattice is recreated when integrator type is changed");
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

      MatlabInterface *mI = (MatlabInterface *) AbstractInterface::getInstance();
      size_t nElements = mxGetNumberOfElements(LATTICE);
      gpuLattice = new Lattice(nElements,integrator, 0.0, gpuId);
      for (size_t i = 0; i < nElements; i++) {
        mxArray *elem = mxGetCell(LATTICE,i);
        mI->setObject(elem);
        gpuLattice->addElement();
      }
      gpuLattice->generateGPUKernel();

    } catch (string& errStr) {
      delete gpuLattice;
      gpuLattice = nullptr;
      mxFree(ref_pts);
      string err =  "at_gpupass() build lattice failed: " + errStr;
      mexErrMsgIdAndTxt("Atpass:RuntimeError",err.c_str());
    }

  }

  // Load lattice on the GPU
  try {
    gpuLattice->fillGPUMemory();
  } catch (string& errStr) {
    mxFree(ref_pts);
    string err =  "at_gpupass() fill GPU memory failed: " + errStr;
    mexErrMsgIdAndTxt("Atpass:RuntimeError",err.c_str());
  }

  // Turn counter
  if( !keep_counter )
    gpuLattice->setTurnCounter(counter);

  // Buffer for lost info
  uint32_t *xnturnPtr = nullptr;
  uint32_t *xnelemPtr = nullptr;
  bool *xlostPtr = nullptr;
  mxArray *mxLostCoord = nullptr;

  try {

    //mexPrintf("Tracking %d particles on %s #%d\n",num_particles,gpuLattice->getGPUContext()->name().c_str(),gpuId);

    uint32_t outsize=num_particles*num_refs*num_turns;
    plhs[0] = mxCreateDoubleMatrix(6,outsize,mxREAL);
    if( plhs[0]==nullptr ) {
      mxFree(ref_pts);
      mexErrMsgIdAndTxt("Atpass:RuntimeError","Not enough memory while trying to allocate particle output coordinates");
    }
    AT_FLOAT *drout = (AT_FLOAT *)mxGetDoubles(plhs[0]);

    if(losses) {

      // Coordinates history not supported
      const char *lossinfo[] = {"lost", "turn", "element", "coordinates_at_loss"/*,"coordinates"*/};
      mwSize LostCoordDims[2] = {6,num_particles};
      mxLostCoord=mxCreateNumericArray(2,LostCoordDims,mxDOUBLE_CLASS,mxREAL);   // Coordinates when particle is lost

      xnturnPtr = new uint32_t[num_particles];
      xnelemPtr = new uint32_t[num_particles];
      xlostPtr = new bool[num_particles];
      AT_FLOAT *xlostcoordPtr = (AT_FLOAT *)mxGetDoubles(mxLostCoord);

      gpuLattice->run(num_turns,num_particles,drin,drout,num_refs,ref_pts,num_starts,track_starts,xnturnPtr,xnelemPtr,xlostcoordPtr,false);

      // Format result for AT
      for(uint32_t i=0;i<num_particles;i++) {
        xlostPtr[i] = (xnturnPtr[i] != num_turns);
        if(!xlostPtr[i]) xnturnPtr[i] = 0;
      }

      mxArray *mxLost=mxCreateLogicalMatrix(1,num_particles);          // lost particle flag
      mxArray *mxNturn=mxCreateDoubleMatrix(1,num_particles,mxREAL);   // Turn number when lost
      mxArray *mxNelem=mxCreateDoubleMatrix(1,num_particles,mxREAL);   // Element number when lost

      // Fill returned array
      mxLogical *xlost=mxGetLogicals(mxLost);
      double *xnturn=mxGetDoubles(mxNturn);
      double *xnelem=mxGetDoubles(mxNelem);
      for(uint32_t i=0;i<num_particles;i++) {
        xlost[i] = xlostPtr[i];
        xnturn[i] = xnturnPtr[i];
        xnelem[i] = xnelemPtr[i];
      }

      // Argout
      mxArray *mxLoss=mxCreateStructMatrix(1,1,4,lossinfo);
      mxSetField(mxLoss, 0, lossinfo[0], mxLost);
      mxSetField(mxLoss, 0, lossinfo[1], mxNturn);
      mxSetField(mxLoss, 0, lossinfo[2], mxNelem);
      mxSetField(mxLoss, 0, lossinfo[3], mxLostCoord);
      plhs[1]=mxLoss;

    } else {

      gpuLattice->run(num_turns,num_particles,drin,drout,num_refs,ref_pts,num_starts,track_starts,nullptr,nullptr,nullptr,false);

    }

  } catch (string& errStr) {
    mxFree(ref_pts);
    mxDestroyArray(plhs[0]);
    if(mxLostCoord) mxDestroyArray(mxLostCoord);
    delete[] xnturnPtr;
    delete[] xnelemPtr;
    delete[] xlostPtr;
    string err =  "at_gpupass() run failed: " + errStr;
    mexErrMsgIdAndTxt("Atpass:RuntimeError",err.c_str());
    return; // Avoid warning (delete non allocated memory)
  }

  delete[] xnturnPtr;
  delete[] xnelemPtr;
  delete[] xlostPtr;
  mxFree(ref_pts);

}
