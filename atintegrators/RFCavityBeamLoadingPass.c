/*
 *  RFCavityBeamLoadingPass.c
 *  Accelerator Toolbox
 *  07/03/2016
 *  Nicola Carmignani
 *
 *  modifications for new trackfunction in 07/04/2017 by N.C.
 *  modifications for matlab 2020 atGetComplex in 05/10/2020 by N.C.
 *
 */

#include "atelem.c"
#include "atlalib.c"
#include <math.h>
#include <complex.h>
#define TWOPI  6.28318530717959
#define C0  	2.99792458e8

struct elem
{
    double Length;
    double Voltage;
    double Energy;
    double Frequency;
    double HarmNumber;
    double _Complex *VoltBeam;
    double *ct_ave_p;
    double f_res;
    double Q;
    double R;
    double TotalCharge;
    /* optional fields */
    double TimeLag;
};


void RFCavityBeamLoadingPass(double *r_in, double le, double nv,
        double freq, double h, double lag, int nturn, double T0,
        int num_particles, double energy,double _Complex *VoltBeam,
        double *ct_ave_p, double f_res,
        double Q, double R, double TotalCharge)
        /* le - physical length
         * nv - peak voltage (V) normalized to the design enegy (eV)
         * r is a 6-by-N matrix of initial conditions reshaped into
         * 1-d array of 6*N elements
         */
{
    int c, c6;
    double halflength , p_norm, NormL, k, psi;
    double ct_ave_new=0.0;
    double _Complex volt_beam_old = *VoltBeam;
    double _Complex volt_beam_new, volt_beam_kick;
    double deltat;
    int debug=0;
    
    
    k = TWOPI * f_res * R / Q;
    psi = atan( 2 * Q * (f_res - freq) / freq );
    
    for(c = 0;c<num_particles;c++)
    {
        c6 = c*6;
        if(!atIsNaN(r_in[c6]))
            ct_ave_new += r_in[c6+5];
    }
    ct_ave_new = ct_ave_new/num_particles;
    
    deltat = T0 + ct_ave_new/C0 - *ct_ave_p/C0;
    
    volt_beam_new = volt_beam_old * cexp(deltat * TWOPI * f_res * _Complex_I - deltat * TWOPI * f_res/2/Q) - k * TotalCharge;
    volt_beam_kick = volt_beam_new + k * TotalCharge / 2;
    
    if(debug)
    {
        printf("volt_beam_old= %f + %f i \n",creal(volt_beam_old), cimag(volt_beam_old));
        printf("volt_beam_new= %f + %f i \n",creal(volt_beam_new), cimag(volt_beam_new));
    }
    
    *VoltBeam = creal(volt_beam_new) + cimag(volt_beam_new) * _Complex_I;
    
#if defined(MATLAB_MEX_FILE)
    if(debug)
    {
        mxArray *Vbeam_8k, *Vgen_8k, *Vtot_8k, *ct_8k;
        mxArray *Vbeam_16k, *Vgen_16k, *Vtot_16k, *ct_16k;
        int i;
        if (nturn==8000)
        {
            Vbeam_8k=mxCreateDoubleMatrix(850,1,mxREAL);
            double *Vbeam_8k_d=mxGetPr(Vbeam_8k);
            Vgen_8k=mxCreateDoubleMatrix(850,1,mxREAL);
            double *Vgen_8k_d=mxGetPr(Vgen_8k);
            Vtot_8k=mxCreateDoubleMatrix(850,1,mxREAL);
            double *Vtot_8k_d=mxGetPr(Vtot_8k);
            ct_8k=mxCreateDoubleMatrix(850,1,mxREAL);
            double *ct_8k_d=mxGetPr(ct_8k);
            for (i=0; i<850; i++)
            {
                Vbeam_8k_d[i]=cabs(volt_beam_kick)*
                        cos(TWOPI*freq*(0.001*(i-400) - ct_ave_new)/C0 + carg(volt_beam_kick));
                Vgen_8k_d[i]=energy*nv*
                        sin(psi+TWOPI*freq*((0.001*(i-400) -lag)/C0 /*- (h/freq-T0)*nturn */));
                Vtot_8k_d[i]=Vbeam_8k_d[i]+Vgen_8k_d[i];
                ct_8k_d[i]=0.001*(i-400);
            }
            int hnum=(int)(freq*T0);
            char name1[80];
            char name2[80];
            char name3[80];
            char name4[80];
            sprintf(name1,"V_beam_8k_h%d",hnum);
            sprintf(name2,"V_gen_8k_h%d",hnum);
            sprintf(name3,"V_tot_8k_h%d",hnum);
            sprintf(name4,"ct_8k_h%d",hnum);
            mexPutVariable("base",name1,Vbeam_8k);
            mexPutVariable("base",name2,Vgen_8k);
            mexPutVariable("base",name3,Vtot_8k);
            mexPutVariable("base",name4,ct_8k);
        }
        if (nturn==16000)
        {
            Vbeam_16k=mxCreateDoubleMatrix(850,1,mxREAL);
            double *Vbeam_16k_d=mxGetPr(Vbeam_16k);
            Vgen_16k=mxCreateDoubleMatrix(850,1,mxREAL);
            double *Vgen_16k_d=mxGetPr(Vgen_16k);
            Vtot_16k=mxCreateDoubleMatrix(850,1,mxREAL);
            double *Vtot_16k_d=mxGetPr(Vtot_16k);
            ct_16k=mxCreateDoubleMatrix(850,1,mxREAL);
            double *ct_16k_d=mxGetPr(ct_16k);
            for (i=0; i<850; i++)
            {
                Vbeam_16k_d[i]=cabs(volt_beam_kick)*cos(TWOPI*freq*(0.001*(i-400) - ct_ave_new)/C0 + carg(volt_beam_kick));
                Vgen_16k_d[i]=energy*nv*sin(psi+TWOPI*freq*((0.001*(i-400) -lag)/C0 /*- (h/freq-T0)*nturn */));
                Vtot_16k_d[i]=Vbeam_16k_d[i]+Vgen_16k_d[i];
                ct_16k_d[i]=0.001*(i-400);
            }
            int hnum=(int)(freq*T0);
            char name1[80];
            char name2[80];
            char name3[80];
            char name4[80];
            sprintf(name1,"V_beam_16k_h%d",hnum);
            sprintf(name2,"V_gen_16k_h%d",hnum);
            sprintf(name3,"V_tot_16k_h%d",hnum);
            sprintf(name4,"ct_16k_h%d",hnum);
            mexPutVariable("base",name1,Vbeam_16k);
            mexPutVariable("base",name2,Vgen_16k);
            mexPutVariable("base",name3,Vtot_16k);
            mexPutVariable("base",name4,ct_16k);
        }
    }
#endif /* MATLAB_MEX_FILE */
    if(le == 0)
    {
        for(c = 0;c<num_particles;c++)
        {	c6 = c*6;
            if(!atIsNaN(r_in[c6]))
            {
                r_in[c6+4] += nv*sin(psi+TWOPI*freq*((r_in[c6+5]-lag)/C0 /*- (h/freq-T0)*nturn */))
                + (cabs(volt_beam_kick)/energy)*cos(TWOPI*freq*(r_in[c6+5] - ct_ave_new)/C0 + carg(volt_beam_kick));
            }
        }
    }
    else
    {	halflength = le/2;
        for(c = 0;c<num_particles;c++)
        {	c6 = c*6;
            if(!atIsNaN(r_in[c6]))
            {   p_norm = 1/(1+r_in[c6+4]);
                NormL  = halflength*p_norm;
                /* Propagate through a drift equal to half cavity length */
                r_in[c6+0]+= NormL*r_in[c6+1];
                r_in[c6+2]+= NormL*r_in[c6+3];
                r_in[c6+5]+= NormL*p_norm*(r_in[c6+1]*r_in[c6+1]+r_in[c6+3]*r_in[c6+3])/2;
                /* Longitudinal momentum kick */
                /*r_in[c6+4] += -nv*sin(TWOPI*freq*((r_in[c6+5]-lag)/C0 - (h/freq-T0)*nturn ));*/
                r_in[c6+4] += nv*sin(psi+TWOPI*freq*((r_in[c6+5]-lag)/C0 /*- (h/freq-T0)*nturn */))
                + (cabs(volt_beam_kick)/energy)*cos(TWOPI*freq*(r_in[c6+5] - ct_ave_new)/C0 + carg(volt_beam_kick));
                p_norm = 1/(1+r_in[c6+4]);
                NormL  = halflength*p_norm;
                /* Propagate through a drift equal to half cavity length */
                r_in[c6+0]+= NormL*r_in[c6+1];
                r_in[c6+2]+= NormL*r_in[c6+3];
                r_in[c6+5]+= NormL*p_norm*(r_in[c6+1]*r_in[c6+1]+r_in[c6+3]*r_in[c6+3])/2;
            }
        }
    }
    *ct_ave_p=ct_ave_new;
}

#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
        double *r_in, int num_particles, struct parameters *Param)
{
    int nturn=Param->nturn;
    double T0=Param->T0;
    if (!Elem) {
        double Length, Voltage, Energy, Frequency, HarmNumber, TimeLag, 
                f_res, Q, R, TotalCurrent;
        double *ct_ave_p;
        double _Complex *Volt_beam;
        Length=atGetDouble(ElemData,"Length"); check_error();
        Voltage=atGetDouble(ElemData,"Voltage"); check_error();
        Energy=atGetDouble(ElemData,"Energy"); check_error();
        Frequency=atGetDouble(ElemData,"Frequency"); check_error();
        HarmNumber=atGetDouble(ElemData,"HarmNumber"); check_error();
        f_res=atGetDouble(ElemData,"f_res"); check_error();
        Q=atGetDouble(ElemData,"Q"); check_error();
        R=atGetDouble(ElemData,"R"); check_error();
        TotalCurrent=atGetDouble(ElemData,"TotalCurrent"); check_error();
        Volt_beam=atGetComplexArray(ElemData,"VoltageBeam"); check_error();
        ct_ave_p=atGetDoubleArray(ElemData,"ct_ave"); check_error();
        
        TimeLag=atGetOptionalDouble(ElemData,"TimeLag",0); check_error();
        
        /*printf("volt_beam= %f + %f i \n",creal(* Volt_beam),
         * cimag(*Volt_beam));*/
        
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->Length=Length;
        Elem->Voltage=Voltage;
        Elem->Energy=Energy;
        Elem->Frequency=Frequency;
        Elem->HarmNumber=HarmNumber;
        Elem->f_res=f_res;
        Elem->Q=Q;
        Elem->R=R;
        Elem->TotalCharge=TotalCurrent*T0;
        Elem->ct_ave_p=ct_ave_p;
        Elem->VoltBeam=Volt_beam;
        Elem->TimeLag=TimeLag;
    }
    
    RFCavityBeamLoadingPass(r_in, Elem->Length, Elem->Voltage/Elem->Energy,
            Elem->Frequency, Elem->HarmNumber, Elem->TimeLag, nturn, T0,
            num_particles, Elem->Energy, Elem->VoltBeam,
            Elem->ct_ave_p, Elem->f_res, Elem->Q,
            Elem->R, Elem->TotalCharge);
    return Elem;
}

MODULE_DEF(RFCavityBeamLoadingPass)        /* Dummy module initialisation */

#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#if defined(MATLAB_MEX_FILE)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if(nrhs == 2)
    {
        double Length, Voltage, Energy, Frequency, HarmNumber, TimeLag, 
                f_res, Q, R, TotalCurrent, TotalCharge;
        double *ct_ave_p, *r_in;
        double _Complex *Volt_beam;
        int nturn=1;
        double T0=3e-6;
        const mxArray *ElemData = prhs[0];
        int num_particles = mxGetN(prhs[1]);
        Length=atGetDouble(ElemData,"Length"); check_error();
        Voltage=atGetDouble(ElemData,"Voltage"); check_error();
        Energy=atGetDouble(ElemData,"Energy"); check_error();
        Frequency=atGetDouble(ElemData,"Frequency"); check_error();
        HarmNumber=atGetDouble(ElemData,"HarmNumber"); check_error();
        f_res=atGetDouble(ElemData,"f_res"); check_error();
        Q=atGetDouble(ElemData,"Q"); check_error();
        R=atGetDouble(ElemData,"R"); check_error();
        TotalCurrent=atGetDouble(ElemData,"TotalCurrent"); check_error();
        Volt_beam=atGetComplexArray(ElemData,"VoltageBeam"); check_error();
        ct_ave_p=atGetDoubleArray(ElemData,"ct_ave"); check_error();
        
        TimeLag=atGetOptionalDouble(ElemData,"TimeLag",0); check_error();
        TotalCharge=TotalCurrent*T0;
        
        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetPr(plhs[0]);
        
        RFCavityBeamLoadingPass(r_in, Length, Voltage/Energy, Frequency, 
                HarmNumber, TimeLag, nturn, T0, num_particles, Energy, 
                Volt_beam, ct_ave_p, f_res, Q, R, TotalCharge);
        
    }
    else if (nrhs == 0)
    {
        /* return list of required fields */
        plhs[0] = mxCreateCellMatrix(11,1);
        mxSetCell(plhs[0],0,mxCreateString("Length"));
        mxSetCell(plhs[0],1,mxCreateString("Voltage"));
        mxSetCell(plhs[0],2,mxCreateString("Energy"));
        mxSetCell(plhs[0],3,mxCreateString("Frequency"));
        mxSetCell(plhs[0],4,mxCreateString("HarmNumber"));
        mxSetCell(plhs[0],5,mxCreateString("f_res"));
        mxSetCell(plhs[0],6,mxCreateString("Q"));
        mxSetCell(plhs[0],7,mxCreateString("R"));
        mxSetCell(plhs[0],8,mxCreateString("TotalCurrent"));
        mxSetCell(plhs[0],9,mxCreateString("VoltageBeam"));
        mxSetCell(plhs[0],10,mxCreateString("ct_ave"));
        
        if(nlhs>1) /* optional fields */
        {
            plhs[1] = mxCreateCellMatrix(1,1);
            mxSetCell(plhs[1],0,mxCreateString("TimeLag"));
        }
    }
    else
    {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
    }
}
#endif /* MATLAB_MEX_FILE */
