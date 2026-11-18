/*
 * Intra-beam scattering pass method.
 *
 * Each pass, the eigen emittances, momentum spread and bunch length of each
 * bunch are computed from the 6D sigma matrix of its particles. Every
 * UpdateTurns turns, the IBS emittance growth rates are computed from these
 * beam parameters and the optics sampled around the ring, with the
 * Bjorken-Mtingwa formalism and the Coulomb logarithm of MAD-X. Each particle
 * then receives random momentum kicks weighted by the longitudinal line
 * density (R. Bruce et al., PRSTAB 13, 091001 (2010)). The momentum kick is
 * compensated by the local dispersion so that it does not change the
 * betatron coordinates.
 *
 * IBS is not defined for a bunch with a singular sigma matrix (for instance
 * a plane with zero spread): such a bunch receives no kick and its growth
 * rates are set to NaN.
 */

#include "atelem.c"
#include "atlalib.c"
#include "atrandom.c"
#include "atibslib.c"
#include <float.h>
#ifdef MPI
#include <mpi.h>
#include <mpi4py/mpi4py.h>
#endif

#define NMOMENTS 28            /* Moments per bunch: count, 6 means, 21 products */

struct elem
{
    long update_turns;
    long nslice;
    int npoints;
    double *optics;         /* (9, npoints) */
    double *local_beta;     /* (2,) */
    double *local_disp;     /* (4,) */
    double npart_per_amp;   /* Particles per ampere of bunch current */
    double *rates;          /* (nbunch, 3) emittance growth rates [1/s] */
};

static long ibs_slice(double ct, double ctmin, double ctmax, long nslice)
{
    long k = (long)(nslice*(ct - ctmin)/(ctmax - ctmin));
    return (k < nslice) ? k : nslice-1;
}

void IBSPass(double *r_in, int num_particles, struct elem *Elem, struct parameters *Param)
{
    int nbunch = Param->nbunch;
    long nslice = Elem->nslice;
    double *disp = Elem->local_disp;
    double *rates = Elem->rates;
    double mass = (Param->rest_energy > 0.0) ? Param->rest_energy : __E0*1.0e9;
    double gamma = Param->energy/mass;
    struct ibs_ring ring = {Param->energy, mass, Param->charge,
        (Param->rest_energy > 0.0) ? sqrt(1.0 - 1.0/(gamma*gamma)) : 1.0};
    double dt = Param->T0;   /* Time between passes */
    int update = (Param->nturn % Elem->update_turns) == 0;
    int i, j, ib;

    double *buffer = atCalloc(nbunch*(NMOMENTS + 2 + nslice + 3), sizeof(double));
    double *moments = buffer;
    double *ctrange = moments + nbunch*NMOMENTS;   /* -min(ct), then max(ct) */
    double *hist = ctrange + 2*nbunch;
    double *kick = hist + nbunch*nslice;

    /* First and second moments of the 6D distribution */
    for (ib = 0; ib < 2*nbunch; ib++) ctrange[ib] = -DBL_MAX;
    for (i = 0; i < num_particles; i++) {
        double *r6 = r_in + 6*i;
        ib = i%nbunch;
        if (!atIsNaN(r6[0])) {
            double *m = moments + NMOMENTS*ib;
            int a, c, n = 7;
            m[0] += 1.0;
            for (a = 0; a < 6; a++) {
                m[1+a] += r6[a];
                for (c = a; c < 6; c++) m[n++] += r6[a]*r6[c];
            }
            if (-r6[5] > ctrange[ib]) ctrange[ib] = -r6[5];
            if (r6[5] > ctrange[nbunch+ib]) ctrange[nbunch+ib] = r6[5];
        }
    }
    #ifdef MPI
    MPI_Allreduce(MPI_IN_PLACE, moments, NMOMENTS*nbunch, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, ctrange, 2*nbunch, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    #endif

    /* Longitudinal line density */
    for (i = 0; i < num_particles; i++) {
        double *r6 = r_in + 6*i;
        ib = i%nbunch;
        double ctmin = -ctrange[ib], ctmax = ctrange[nbunch+ib];
        if (!atIsNaN(r6[0]) && ctmax > ctmin)
            hist[nslice*ib + ibs_slice(r6[5], ctmin, ctmax, nslice)] += 1.0;
    }
    #ifdef MPI
    MPI_Allreduce(MPI_IN_PLACE, hist, nbunch*nslice, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    #endif

    /* Beam parameters, growth rates, kick amplitudes and line density weights */
    for (ib = 0; ib < nbunch; ib++) {
        double *m = moments + NMOMENTS*ib;
        double *h = hist + nslice*ib;
        double *k = kick + 3*ib;
        double n = m[0];
        if (n < 2.0) continue;
        double sigma[36], emit[3];
        int a, c, nm = 7;
        for (a = 0; a < 6; a++) {
            for (c = a; c < 6; c++) {
                sigma[6*a+c] = sigma[6*c+a] = m[nm++]/n - (m[1+a]/n)*(m[1+c]/n);
            }
        }
        if (!eigen_emittances(sigma, emit)) {
            for (j = 0; j < 3; j++) rates[ib + nbunch*j] = atGetNaN();
            continue;
        }
        struct ibs_beam b;
        b.npart = Param->bunch_currents[ib]*Elem->npart_per_amp;
        b.emitx = emit[0];
        b.emity = emit[1];
        b.sigma_e = sqrt(sigma[6*4+4]);
        b.bunch_length = sqrt(sigma[6*5+5]);
        if (update || atIsNaN(rates[ib])) {
            double r[3];
            ibs_growth_rates(Elem->npoints, Elem->optics, &ring, &b, r);
            for (j = 0; j < 3; j++) rates[ib + nbunch*j] = r[j];
        }
        double h2 = 0.0;
        for (j = 0; j < nslice; j++) h2 += h[j]*h[j];
        if (h2 <= 0.0) continue;
        /* Weights sqrt(line density) normalised to a mean density of 1 */
        for (j = 0; j < nslice; j++) h[j] = sqrt(h[j]*n/h2);
        /* <dp^2> per pass = 2 d(emittance) / beta in the transverse planes */
        k[0] = sqrt(2.0*fmax(rates[ib], 0.0)*dt*b.emitx/Elem->local_beta[0]);
        k[1] = sqrt(2.0*fmax(rates[ib + nbunch], 0.0)*dt*b.emity/Elem->local_beta[1]);
        k[2] = sqrt(2.0*fmax(rates[ib + 2*nbunch], 0.0)*dt)*b.sigma_e;
    }

    /* Random kicks */
    #pragma omp parallel for if (num_particles > OMP_PARTICLE_THRESHOLD*10) default(shared) private(i, ib)
    for (i = 0; i < num_particles; i++) {
        double *r6 = r_in + 6*i;
        ib = i%nbunch;
        double *k = kick + 3*ib;
        if (!atIsNaN(r6[0]) && (k[0] > 0.0 || k[1] > 0.0 || k[2] > 0.0)) {
            long s = ibs_slice(r6[5], -ctrange[ib], ctrange[nbunch+ib], nslice);
            double w = hist[nslice*ib + s];
            double dpx = k[0]*w*atrandn_r(Param->thread_rng, 0.0, 1.0);
            double dpy = k[1]*w*atrandn_r(Param->thread_rng, 0.0, 1.0);
            double dd = k[2]*w*atrandn_r(Param->thread_rng, 0.0, 1.0);
            r6[0] += disp[0]*dd;
            r6[1] += disp[1]*dd + dpx;
            r6[2] += disp[2]*dd;
            r6[3] += disp[3]*dd + dpy;
            r6[4] += dd;
        }
    }
    atFree(buffer);
}

#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData, struct elem *Elem,
                                      double *r_in, int num_particles, struct parameters *Param)
{
    if (!Elem) {
        int msz, nsz;
        long update_turns = atGetLong(ElemData, "UpdateTurns"); check_error();
        long nslice = atGetLong(ElemData, "NSlice"); check_error();
        double npart_per_amp = atGetDouble(ElemData, "_npart_per_amp"); check_error();
        double *optics = atGetDoubleArraySz(ElemData, "_optics", &msz, &nsz); check_error();
        double *local_beta = atGetDoubleArray(ElemData, "_local_beta"); check_error();
        double *local_disp = atGetDoubleArray(ElemData, "_local_disp"); check_error();
        double *rates = atGetDoubleArray(ElemData, "_rates"); check_error();
        int dims[] = {Param->nbunch, 3};
        atCheckArrayDims(ElemData, "_rates", 2, dims); check_error();
        if (msz != IBS_NOPTICS) atError("_optics must have 9 rows."); check_error();
        Elem = (struct elem *)atMalloc(sizeof(struct elem));
        Elem->update_turns = (update_turns > 0) ? update_turns : 1;
        Elem->nslice = (nslice > 0) ? nslice : 1;
        Elem->npoints = nsz;
        Elem->optics = optics;
        Elem->local_beta = local_beta;
        Elem->local_disp = local_disp;
        Elem->npart_per_amp = npart_per_amp;
        Elem->rates = rates;
    }
    IBSPass(r_in, num_particles, Elem, Param);
    return Elem;
}

MODULE_DEF(IBSPass)        /* Dummy module initialisation */

#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#if defined(MATLAB_MEX_FILE)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs >= 2) {
        /* Without the ring parameters, the bunch current is zero: no kick */
        if (mxGetM(prhs[1]) != 6) mexErrMsgIdAndTxt("AT:WrongArg","Second argument must be a 6 x N matrix");
        plhs[0] = mxDuplicateArray(prhs[1]);
    }
    else if (nrhs == 0) {
        /* list of required fields */
        plhs[0] = mxCreateCellMatrix(7,1);
        mxSetCell(plhs[0],0,mxCreateString("UpdateTurns"));
        mxSetCell(plhs[0],1,mxCreateString("NSlice"));
        mxSetCell(plhs[0],2,mxCreateString("_npart_per_amp"));
        mxSetCell(plhs[0],3,mxCreateString("_optics"));
        mxSetCell(plhs[0],4,mxCreateString("_local_beta"));
        mxSetCell(plhs[0],5,mxCreateString("_local_disp"));
        mxSetCell(plhs[0],6,mxCreateString("_rates"));
    }
    else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
    }
}
#endif /*defined(MATLAB_MEX_FILE)*/
