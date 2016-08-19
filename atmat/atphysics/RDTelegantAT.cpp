/*
 *  This code is based on the elegant function computeDrivingTerms
 *  it has been modified to be used from the matlab accelerator toolbox 
 *  using the function computeRDT()
 *  The formulas have not been changed
 *  
 */

#include "mex.h"
#include <complex>
#include <cmath>
#include <cstring>

typedef struct {
  double betax, betay;
  double rbetax, rbetay;           /* rbetax = sqrt(betax) */
  double betax2, betay2;           /* betax2 = betax^2 */
  double phix, phiy;
  std::complex <double> px[5], py[5];    /* px[j]=(exp(i*phix))^j j>0 */
  double b2L, b3L, s;
} ELEMDATA;

typedef struct {
  /* First-order geometric terms (abs, real, imag) */
  double h21000[3], h30000[3], h10110[3], h10020[3], h10200[3];
  /* First order chromatic terms */
  double h11001[3], h00111[3], h20001[3], h00201[3], h10002[3];
  /* First order coupling terms */
  double h10010[3], h10100[3];
  /* Second-order geometric terms */
  double h22000[3], h11110[3], h00220[3], h31000[3], h40000[3];
  double h20110[3], h11200[3], h20020[3], h20200[3], h00310[3], h00400[3];
  /* tune shifts with amplitude */
  double dnux_dJx, dnux_dJy, dnuy_dJy;

  /* resonance driving terms, obtained from hamiltonian driving terms */
  /* First-order geometric terms (abs, real, imag) */
  double f21000[3], f30000[3], f10110[3], f10020[3], f10200[3];
  /* First order coupling terms */
  double f10010[3], f10100[3];
  /* Second-order geometric terms */
  double f31000[3], f40000[3];
  double f20110[3], f11200[3], f20020[3], f20200[3], f00310[3], f00400[3];
  
  
} DRIVING_TERMS;

void computeDrivingTerms(DRIVING_TERMS *d, int NumElem,
							double *s, double *betax, double *betay, double *phix, double *phiy, double *etax,
							double *Lista2L, double *Listb2L, double *Listb3L, double *Listb4L, double *tune,
							char Geometric1, char Geometric2, char Chromatic1, char Coupling1, char TuneShifts, long nPeriods);

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]);

/* square function */
double sqr(double input)
{
	double square=0;
	square = input * input;
	return square;
}
/* sign function */
double SIGN(double number)
{
	if (number>0)
		return 1;
	else
		if (number<0)
			return -1;
		else
			return 0;
}

/* void computeDrivingTerms(DRIVING_TERMS *d, ELEMENT_LIST *elem, TWISS *twiss0, double *tune, long nPeriods)*/
void computeDrivingTerms(DRIVING_TERMS *d, int NumElem,
							double *s, double *betax, double *betay, double *phix, double *phiy, double *etax,
							double *Lista2L, double *Listb2L, double *Listb3L, double *Listb4L, double *tune, 
							char Chromatic1, char Coupling1, char Geometric1, char Geometric2, char TuneShifts,
							long nPeriods)
/* Based on J. Bengtsson, SLS Note 9/97, March 7, 1997, with corrections per W. Guo (NSLS) */
/* Revised to follow C. X. Wang AOP-TN-2009-020 for second-order terms */
{
  std::complex <double> h11001, h00111, h20001, h00201, h10002, h10010, h10100;
  std::complex <double> h21000, h30000, h10110, h10020, h10200;
  std::complex <double> h22000, h11110, h00220, h31000, h40000;
  std::complex <double> h20110, h11200, h20020, h20200, h00310, h00400;
  std::complex <double> t1, t2, t3, t4;
  std::complex <double> ii;
  std::complex <double> periodicFactor[9][9];
#define PF(i,j) (periodicFactor[4+i][4+j])
  double betax1, betay1, phix1, phiy1, etax1, termSign;
  double coef, b2L, a2L, b3L, b4L, nux, nuy;
  double PIx2 = 6.283185307179586;
  double PI = 3.141592653589793;
//  ELEMENT_LIST *eptr1;
  double two=2, three=3, four=4;
  ELEMDATA *ed = NULL;
  long nE=0, iE, jE, i, j;
  double sqrt8, sqrt2, tilt;

  sqrt8 = sqrt((double)8);
  sqrt2 = sqrt((double)2);
  ii = std::complex<double>(0,1);

//printf("leading_order_driving_terms_only=%d\n",leading_order_driving_terms_only);

  /* accumulate real and imaginary parts */
  h11001 = h00111 = h20001 = h00201 = h10002 = std::complex<double>(0,0);
  h21000 = h30000 = h10110 = h10020 = h10200 = std::complex<double>(0,0);
  h22000 = h11110 = h00220 = h31000 = h40000 = std::complex<double>(0,0);
  h20110 = h11200 = h20020 = h20200 = h00310 = h00400 = std::complex<double>(0,0);
  h10100 = h10010 = std::complex<double>(0,0);

  d->dnux_dJx = d->dnux_dJy = d->dnuy_dJy = 0;

  if (nPeriods!=1) {
    double a1, a2;
    for (i=0; i<9; i++) {
      for (j=0; j<9; j++) {
        a1 = PIx2*(tune[0]*(i-4)+tune[1]*(j-4));
        a2 = a1/nPeriods;
        periodicFactor[i][j] = (exp(ii*a1)-1.0)/(exp(ii*a2)-1.0);
      }
    }
  } else {
    for (i=0; i<9; i++)
      for (j=0; j<9; j++)
        periodicFactor[i][j] = 1;
  }

ed = (ELEMDATA*)malloc(NumElem*sizeof(ELEMDATA));

/*  eptr1 = elem;*/
/*  while (eptr1) {*/
  for (nE=0; nE<NumElem; nE++)
	{
    a2L = b2L = b3L = b4L = 0;
    a2L = Lista2L[nE];
    
    b2L = Listb2L[nE];
    
    b3L = Listb3L[nE];
    
    b4L = Listb4L[nE];
        
    //printf("a2L=%f,b2L=%f,b3L=%f,b4L=%f\n",a2L,b2L,b3L,b4L);
    
//    if (a2L || b2L || b3L || b4L) 
    	{
        betax1 = betax[nE];
        etax1  = etax[nE];
        phix1  = phix[nE];
        betay1 = betay[nE];
        phiy1  = phiy[nE];
    	
    
/*      ed = (ELEMDATA*)SDDS_Realloc(ed, sizeof(*ed)*(nE+1));*/
      ed[nE].s = s[nE];
      ed[nE].b2L = Listb2L[nE];
      ed[nE].b3L = Listb3L[nE]; //printf("el %d. b3L= %f \n",nE,ed[nE].b3L);
      ed[nE].betax = betax1;
      ed[nE].betax2 = sqr(betax1);
      ed[nE].rbetax = sqrt(betax1);
      ed[nE].phix = phix1;
      ed[nE].px[1] = exp(ii*phix1);
      ed[nE].px[2] = ed[nE].px[1]*ed[nE].px[1];
      ed[nE].px[3] = ed[nE].px[1]*ed[nE].px[2];
      ed[nE].px[4] = ed[nE].px[1]*ed[nE].px[3];
      ed[nE].betay = betay1;
      ed[nE].betay2 = sqr(betay1);
      ed[nE].rbetay = sqrt(betay1);
      ed[nE].phiy = phiy1;
      ed[nE].py[1] = exp(ii*phiy1);
      ed[nE].py[2] = ed[nE].py[1]*ed[nE].py[1];
      ed[nE].py[3] = ed[nE].py[1]*ed[nE].py[2];
      ed[nE].py[4] = ed[nE].py[1]*ed[nE].py[3];

      if (fabs(a2L)>1e-6) {
		if(Coupling1)
		{	/* linear coupling terms */
			h10010 += (a2L/4)*ed[nE].rbetax*ed[nE].rbetay*ed[nE].px[1]/ed[nE].py[1]*PF(1, -1);
			h10100 += (a2L/4)*ed[nE].rbetax*ed[nE].rbetay*ed[nE].px[1]*ed[nE].py[1]*PF(1, 1);
		}
      }
      if (fabs(b2L)>1e-6 || fabs(b3L)>1e-6) 
      {
		if(Chromatic1)
		{
			/* first-order chromatic terms */
			/* h11001 and h00111 */
			h11001 += (b3L*betax1*etax1/2-b2L*betax1/4)*nPeriods;
			h00111 += (b2L*betay1/4-b3L*betay1*etax1/2)*nPeriods;
			/* h20001, h00201 */
			h20001 += (b3L*betax1*etax1/2-b2L*betax1/4)/2*ed[nE].px[2]*PF(2,0);
			h00201 += (b2L*betay1/4-b3L*betay1*etax1/2)/2*ed[nE].py[2]*PF(0,2);

			/* h10002 */
			h10002 += (b3L*ed[nE].rbetax*pow(etax1,2)-b2L*ed[nE].rbetax*etax1)/2*ed[nE].px[1]*PF(1,0);
			//	h10002 += (b3L*ed[nE].rbetax*ipow(etax1,2)-b2L*ed[nE].rbetax*etax1)/2*ed[nE].px[1]*PF(1,0);
		}
      }
      if (fabs(ed[nE].b3L)>1e-6) 
      {
      	if(Geometric1)
        {
        	/* first-order geometric terms from sextupoles */
        	/* h21000 */
			h21000 += b3L*ed[nE].rbetax*betax1/8*ed[nE].px[1]*PF(1,0);
        	/* h30000 */
			h30000 += b3L*ed[nE].rbetax*betax1/24*ed[nE].px[3]*PF(3,0);
            /* h10110 */
			h10110 += -b3L*ed[nE].rbetax*betay1/4*ed[nE].px[1]*PF(1,0);
	        /* h10020 and h10200 */
			h10020 += -b3L*ed[nE].rbetax*betay1/8*ed[nE].px[1]*conj(ed[nE].py[2])*PF(1,-2);
			h10200 += -b3L*ed[nE].rbetax*betay1/8*ed[nE].px[1]*ed[nE].py[2]*PF(1,2);
		}
      }
      if (fabs(b4L)>1e-6) 
      {
      	if (TuneShifts)
      	{
        	/* second-order terms from leading order effects of octupoles */
			/* Ignoring a large number of terms that are not also driven by sextupoles */
        	d->dnux_dJx += 3*b4L*ed[nE].betax2/(8*PI)*nPeriods;
        	d->dnux_dJy -= 3*b4L*betax1*betay1/(4*PI)*nPeriods;
        	d->dnuy_dJy += 3*b4L*ed[nE].betay2/(8*PI)*nPeriods;
		}
	if (Geometric2)
	{
		h22000 += 3*b4L*ed[nE].betax2/32*nPeriods;
		h11110 += -3*b4L*betax1*betay1/8*nPeriods;
		h00220 += 3*b4L*ed[nE].betay2/32*nPeriods;
		h31000 += b4L*ed[nE].betax2/16*ed[nE].px[2]*PF(2,0);
		h40000 += b4L*ed[nE].betax2/64*ed[nE].px[4]*PF(4,0);
		h20110 += -3*b4L*betax1*betay1/16*ed[nE].px[2]*PF(2,0);
		h11200 += -3*b4L*betax1*betay1/16*ed[nE].py[2]*PF(0,2);
		h20020 += -3*b4L*betax1*betay1/32*ed[nE].px[2]*conj(ed[nE].py[2])*PF(2,-2);
		h20200 += -3*b4L*betax1*betay1/32*ed[nE].px[2]*ed[nE].py[2]*PF(2,2);
		h00310 += b4L*ed[nE].betay2/16*ed[nE].py[2]*PF(0,2);
		h00400 += b4L*ed[nE].betay2/64*ed[nE].py[4]*PF(0,4);
	}
      }
  }
  }
	
  /* Done with the leading-order quad and sext terms */
  d->h11001[0] = std::abs<double>(h11001);
  d->h11001[1] = h11001.real();
  d->h11001[2] = h11001.imag();
  d->h00111[0] = std::abs<double>(h00111);
  d->h00111[1] = h00111.real();
  d->h00111[2] = h00111.imag();
  d->h20001[0] = std::abs<double>(h20001);
  d->h20001[1] = h20001.real();
  d->h20001[2] = h20001.imag();
  d->h00201[0] = std::abs<double>(h00201);
  d->h00201[1] = h00201.real();
  d->h00201[2] = h00201.imag();
  d->h10002[0] = std::abs<double>(h10002);
  d->h10002[1] = h10002.real();
  d->h10002[2] = h10002.imag();

  d->h10100[0] = std::abs<double>(h10100);
  d->h10100[1] = h10100.real();
  d->h10100[2] = h10100.imag();
  d->h10010[0] = std::abs<double>(h10010);
  d->h10010[1] = h10010.real();
  d->h10010[2] = h10010.imag();

  d->h21000[0] = std::abs<double>(h21000);
  d->h21000[1] = h21000.real();
  d->h21000[2] = h21000.imag();
  d->h30000[0] = std::abs<double>(h30000);
  d->h30000[1] = h30000.real();
  d->h30000[2] = h30000.imag();
  d->h10110[0] = std::abs<double>(h10110);
  d->h10110[1] = h10110.real();
  d->h10110[2] = h10110.imag();
  d->h10020[0] = std::abs<double>(h10020);
  d->h10020[1] = h10020.real();
  d->h10020[2] = h10020.imag();
  d->h10200[0] = std::abs<double>(h10200);
  d->h10200[1] = h10200.real();
  d->h10200[2] = h10200.imag();
//	printf("d->dnux_dJx=%f, d->dnux_dJy=%f, d->dnuy_dJy=%f\n",d->dnux_dJx,d->dnux_dJy,d->dnuy_dJy);
 
//  if (!leading_order_driving_terms_only) 
  if (Geometric2 || TuneShifts) 
  {
    /* compute sextupole contributions to second-order terms */
    if (nPeriods!=1) 
    {
      printf("ERROR: computating of higher-order driving terms not available when n_periods!=1");
//      return 0;
	}
	//printf("entered in !leading_order_driving_terms_only\n");
    nux = tune[0];
    nuy = tune[1];

    for (iE=0; iE<NumElem; iE++) {
      if (fabs(ed[iE].b3L)>1e-6) {
	for (jE=0; jE<NumElem; jE++) {
	  if (fabs(ed[jE].b3L)>1e-6) {
	  if(TuneShifts)
	  {
	    d->dnux_dJx += ed[iE].b3L*ed[jE].b3L/(-16*PI)*pow(ed[iE].betax*ed[jE].betax, 1.5)*
	      (3*cos(fabs(ed[iE].phix-ed[jE].phix)-PI*nux)/sin(PI*nux) + cos(fabs(3*(ed[iE].phix-ed[jE].phix))-3*PI*nux)/sin(3*PI*nux));
	    d->dnux_dJy += ed[iE].b3L*ed[jE].b3L/(8*PI)*sqrt(ed[iE].betax*ed[jE].betax)*ed[iE].betay*
	      (2*ed[jE].betax*cos(fabs(ed[iE].phix-ed[jE].phix)-PI*nux)/sin(PI*nux) 
	       - ed[jE].betay*cos(fabs(ed[iE].phix-ed[jE].phix)+2*fabs(ed[iE].phiy-ed[jE].phiy)-PI*(nux+2*nuy))/sin(PI*(nux+2*nuy))
	       + ed[jE].betay*cos(fabs(ed[iE].phix-ed[jE].phix)-2*fabs(ed[iE].phiy-ed[jE].phiy)-PI*(nux-2*nuy))/sin(PI*(nux-2*nuy)));
	    d->dnuy_dJy += ed[iE].b3L*ed[jE].b3L/(-16*PI)*sqrt(ed[iE].betax*ed[jE].betax)*ed[iE].betay*ed[jE].betay*
	      (4*cos(fabs(ed[iE].phix-ed[jE].phix)-PI*nux)/sin(PI*nux) 
	       + cos(fabs(ed[iE].phix-ed[jE].phix)+2*fabs(ed[iE].phiy-ed[jE].phiy)-PI*(nux+2*nuy))/sin(PI*(nux+2*nuy)) 
	       + cos(fabs(ed[iE].phix-ed[jE].phix)-2*fabs(ed[iE].phiy-ed[jE].phiy)-PI*(nux-2*nuy))/sin(PI*(nux-2*nuy)));
	   }
	   if (Geometric2)
	   {
	    termSign = SIGN(ed[iE].s - ed[jE].s);
            if (termSign) {
              /* geometric terms */
              h22000 += (1./64)*termSign*ii*ed[iE].b3L*ed[jE].b3L*
                ed[iE].rbetax*ed[jE].rbetax*ed[iE].betax*ed[jE].betax*
                (ed[iE].px[3]*conj(ed[jE].px[3]) + three*ed[iE].px[1]*conj(ed[jE].px[1]));
              h31000 += (1./32)*termSign*ii*ed[iE].b3L*ed[jE].b3L*
                ed[iE].rbetax*ed[jE].rbetax*ed[iE].betax*ed[jE].betax*
                ed[iE].px[3]*conj(ed[jE].px[1]);
	      t1 = conj(ed[iE].px[1])*ed[jE].px[1];
	      t2 = ed[iE].px[1]*conj(ed[jE].px[1]);
              h11110 += (1./16)*termSign*ii*ed[iE].b3L*ed[jE].b3L*
                ed[iE].rbetax*ed[jE].rbetax*ed[iE].betay*
                (ed[jE].betax*(t1 - conj(t1) )
                 + 
		 ed[jE].betay*ed[iE].py[2]*conj(ed[jE].py[2])*(conj(t1) + t1)
		 );
	      t1 = exp(-ii*(ed[iE].phix-ed[jE].phix));
	      t2 = conj(t1);
              h11200 += (1./32)*termSign*ii*ed[iE].b3L*ed[jE].b3L*
                ed[iE].rbetax*ed[jE].rbetax*ed[iE].betay*exp(ii*(2*ed[iE].phiy))*
                (ed[jE].betax*(t1 - t2)
                 + 
		 two*ed[jE].betay*(t2 + t1)
		 );
              h40000 += (1./64)*termSign*ii*ed[iE].b3L*ed[jE].b3L*
                ed[iE].rbetax*ed[jE].rbetax*ed[iE].betax*ed[jE].betax*
                ed[iE].px[3]*ed[jE].px[1];
              h20020 += (1./64)*termSign*ii*ed[iE].b3L*ed[jE].b3L*
                ed[iE].rbetax*ed[jE].rbetax*ed[iE].betay*
                (ed[jE].betax*conj(ed[iE].px[1]*ed[iE].py[2])*ed[jE].px[3]
                 -(ed[jE].betax+four*ed[jE].betay)*ed[iE].px[1]*ed[jE].px[1]*conj(ed[iE].py[2])
		 );
              h20110 += (1./32)*termSign*ii*ed[iE].b3L*ed[jE].b3L*
                ed[iE].rbetax*ed[jE].rbetax*ed[iE].betay*
                (ed[jE].betax*(
			       conj(ed[iE].px[1])*ed[jE].px[3]
			       - 
			       ed[iE].px[1]*ed[jE].px[1]
			       ) 
                 + two*ed[jE].betay*ed[iE].px[1]*ed[jE].px[1]*ed[iE].py[2]*conj(ed[jE].py[2])
		 );
              h20200 += (1./64)*termSign*ii*ed[iE].b3L*ed[jE].b3L*
                ed[iE].rbetax*ed[jE].rbetax*ed[iE].betay*
                (ed[jE].betax*
		 conj(ed[iE].px[1])*ed[jE].px[3]*ed[iE].py[2]
                 -(ed[jE].betax-four*ed[jE].betay)*
		 ed[iE].px[1]*ed[jE].px[1]*ed[iE].py[2]
		 );
              h00220 += (1./64)*termSign*ii*ed[iE].b3L*ed[jE].b3L*
                ed[iE].rbetax*ed[jE].rbetax*ed[iE].betay*ed[jE].betay*
		(ed[iE].px[1]*ed[iE].py[2]*conj(ed[jE].px[1]*ed[jE].py[2])
                 + four*ed[iE].px[1]*conj(ed[jE].px[1])
		 - conj(ed[iE].px[1]*ed[jE].py[2])*ed[jE].px[1]*ed[iE].py[2]
		 );
              h00310 += (1./32)*termSign*ii*ed[iE].b3L*ed[jE].b3L*
                ed[iE].rbetax*ed[jE].rbetax*ed[iE].betay*ed[jE].betay*ed[iE].py[2]*
                (ed[iE].px[1]*conj(ed[jE].px[1])
		 -conj(ed[iE].px[1])*ed[jE].px[1]
		 );
              h00400 += (1./64)*termSign*ii*ed[iE].b3L*ed[jE].b3L*
                ed[iE].rbetax*ed[jE].rbetax*ed[iE].betay*ed[jE].betay*
		ed[iE].px[1]*conj(ed[jE].px[1])*ed[iE].py[2]*ed[jE].py[2];
	    }
	    }
    }
    }
    }
    }
  }
//	printf("d->dnux_dJx=%f, d->dnux_dJy=%f, d->dnuy_dJy=%f\n",d->dnux_dJx,d->dnux_dJy,d->dnuy_dJy);
  d->h22000[0] = std::abs<double>(h22000);
  d->h22000[1] = h22000.real();
  d->h22000[2] = h22000.imag();
  d->h11110[0] = std::abs<double>(h11110);
  d->h11110[1] = h11110.real();
  d->h11110[2] = h11110.imag();
  d->h00220[0] = std::abs<double>(h00220);
  d->h00220[1] = h00220.real();
  d->h00220[2] = h00220.imag();
  d->h31000[0] = std::abs<double>(h31000);
  d->h31000[1] = h31000.real();
  d->h31000[2] = h31000.imag();
  d->h40000[0] = std::abs<double>(h40000);
  d->h40000[1] = h40000.real();
  d->h40000[2] = h40000.imag();
  d->h20110[0] = std::abs<double>(h20110);
  d->h20110[1] = h20110.real();
  d->h20110[2] = h20110.imag();
  d->h11200[0] = std::abs<double>(h11200);
  d->h11200[1] = h11200.real();
  d->h11200[2] = h11200.imag();
  d->h20020[0] = std::abs<double>(h20020);
  d->h20020[1] = h20020.real();
  d->h20020[2] = h20020.imag();
  d->h20200[0] = std::abs<double>(h20200);
  d->h20200[1] = h20200.real();
  d->h20200[2] = h20200.imag();
  d->h00310[0] = std::abs<double>(h00310);
  d->h00310[1] = h00310.real();
  d->h00310[2] = h00310.imag();
  d->h00400[0] = std::abs<double>(h00400);
  d->h00400[1] = h00400.real();
  d->h00400[2] = h00400.imag();

  if (ed)
    free(ed);
}





void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    double *outMatrixRe;		/* output matrix Real part of hamiltonian driving terms */
    double *outMatrixIm;		/* output matrix Imaginary part of hamiltonian driving terms */
    double *outMatrixTSwA;		/* output matrix Tune Shifts with Amplitude */
	DRIVING_TERMS *d;			/* output struct, to be converted in arrays */
	double *s;
	double *betax;
    double *betay;
    double *etax;
    double *phix;
    double *phiy;
    double *Lista2L;
    double *Listb2L;
    double *Listb3L;
    double *Listb4L;
    int leading_order_driving_terms_only;
	double Tunex;
	double Tuney;
	double tune[2];
	int nPeriods = 1;
	int NumElem, ncols;
	char Geometric1, Geometric2, Chromatic1, Coupling1, TuneShifts;	/*boolean values, taken as input*/
	
	d = (DRIVING_TERMS*)malloc(sizeof(DRIVING_TERMS));
    /* check for proper number of arguments */
    if(nrhs!=18) {
        mexErrMsgIdAndTxt("MyToolbox:RDTelegantAT:nrhs","18 inputs required.");
    }
    if(nlhs!=3) {
        mexErrMsgIdAndTxt("MyToolbox:RDTelegantAT:nlhs","3 outputs required.");
    }
	
//	printf("entered\n");
	
	/* the number of rows of all the inputs must be 1 and the number of columns NumElem, i.e. the 13th input element */
    ncols = mxGetN(prhs[0]);
    
    /* get the value of the 4 scalar inputs: 11,12,13,14 (Tunex, Tuney, NumElem,leading_order_driving_terms_only)  */
    Tunex = mxGetScalar(prhs[10]);
    Tuney = mxGetScalar(prhs[11]);
    NumElem = mxGetScalar(prhs[12]);
    //leading_order_driving_terms_only=mxGetScalar(prhs[13]);
    Geometric1 = mxGetScalar(prhs[13]);
    Geometric2 = mxGetScalar(prhs[14]);
    Chromatic1 = mxGetScalar(prhs[15]);
    Coupling1 = mxGetScalar(prhs[16]);
    TuneShifts = mxGetScalar(prhs[17]);
    
	tune[0]=Tunex;
	tune[1]=Tuney;
		
    /* create a pointer to the real data in the input matrix  */
    s = mxGetPr(prhs[0]);
    betax = mxGetPr(prhs[1]);
    betay = mxGetPr(prhs[2]);
    etax = mxGetPr(prhs[3]);
    phix = mxGetPr(prhs[4]);
    phiy = mxGetPr(prhs[5]);
    Lista2L = mxGetPr(prhs[6]);
    Listb2L = mxGetPr(prhs[7]);
    Listb3L = mxGetPr(prhs[8]);
    Listb4L = mxGetPr(prhs[9]);
	
	//printf("s0=%f, s1=%f, s2=%f \n",s[0],s[1],s[2]);	
	//printf("NumElem=%d\n",NumElem);
    /* get dimensions of the input matrix */
    ncols = mxGetN(prhs[0]);

    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(1,(mwSize)23,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1,(mwSize)23,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1,(mwSize)3,mxREAL);
	plhs[3] = mxCreateDoubleMatrix(1,(mwSize)15,mxREAL);
    plhs[4] = mxCreateDoubleMatrix(1,(mwSize)15,mxREAL);
    
    /* call the computational routine */
    computeDrivingTerms(d, NumElem,
							s, betax, betay, phix, phiy, etax,
							Lista2L, Listb2L, Listb3L, Listb4L, tune, 
							Geometric1, Geometric2, Chromatic1, Coupling1, TuneShifts, nPeriods);

	outMatrixRe = mxGetPr(plhs[0]);
	outMatrixIm = mxGetPr(plhs[1]);
	outMatrixTSwA = mxGetPr(plhs[2]);
	/* Compute the resonance driving terms from the hamiltonian driving terms */
	
	    
	/* put the output d in the vectors outMatrixRe, outMatrixIm and outMatrixTSwA that will be sent to matlab */
	/* First-order geometric terms (real, imag) */
	outMatrixRe[0] = d->h21000[1];
	outMatrixIm[0] = d->h21000[2];
	outMatrixRe[1] = d->h30000[1];
	outMatrixIm[1] = d->h30000[2];
	outMatrixRe[2] = d->h10110[1];
	outMatrixIm[2] = d->h10110[2];
	outMatrixRe[3] = d->h10020[1];
	outMatrixIm[3] = d->h10020[2];
	outMatrixRe[4] = d->h10200[1];
	outMatrixIm[4] = d->h10200[2];
	/* First order chromatic terms */
	outMatrixRe[5] = d->h11001[1];
	outMatrixIm[5] = d->h11001[2];
	outMatrixRe[6] = d->h00111[1];
	outMatrixIm[6] = d->h00111[2];
	outMatrixRe[7] = d->h20001[1];
	outMatrixIm[7] = d->h20001[2];
	outMatrixRe[8] = d->h00201[1];
	outMatrixIm[8] = d->h00201[2];
	outMatrixRe[9] = d->h10002[1];
	outMatrixIm[9] = d->h10002[2];
	/* First order coupling terms */	
	outMatrixRe[10] = d->h10010[1];
	outMatrixIm[10] = d->h10010[2];
	outMatrixRe[11] = d->h10100[1];
	outMatrixIm[11] = d->h10100[2];
	/* Second-order geometric terms */
	outMatrixRe[12] = d->h22000[1];
	outMatrixIm[12] = d->h22000[2];
	outMatrixRe[13] = d->h11110[1];
	outMatrixIm[13] = d->h11110[2];
	outMatrixRe[14] = d->h00220[1];
	outMatrixIm[14] = d->h00220[2];
	outMatrixRe[15] = d->h31000[1];
	outMatrixIm[15] = d->h31000[2];
	outMatrixRe[16] = d->h40000[1];
	outMatrixIm[16] = d->h40000[2];
	outMatrixRe[17] = d->h20110[1];
	outMatrixIm[17] = d->h20110[2];
	outMatrixRe[18] = d->h11200[1];
	outMatrixIm[18] = d->h11200[2];	
	outMatrixRe[19] = d->h20020[1];
	outMatrixIm[19] = d->h20020[2];	
	outMatrixRe[20] = d->h20200[1];
	outMatrixIm[20] = d->h20200[2];
	outMatrixRe[21] = d->h00310[1];
	outMatrixIm[21] = d->h00310[2];
	outMatrixRe[22] = d->h00400[1];
	outMatrixIm[22] = d->h00400[2];
	/* tune shifts with amplitude */
	outMatrixTSwA[0] = d->dnux_dJx;
	outMatrixTSwA[1] = d->dnux_dJy;
	outMatrixTSwA[2] = d->dnuy_dJy;




/*	outMatrix[0] = d->h21000[0];
	outMatrix[1] = d->h21000[1];
	outMatrix[2] = d->h21000[2];
	outMatrix[3] = d->h30000[0];
	outMatrix[4] = d->h30000[1];
	outMatrix[5] = d->h30000[2];
	outMatrix[6] = d->h10110[0];
	outMatrix[7] = d->h10110[1];
	outMatrix[8] = d->h10110[2];
	outMatrix[9] = d->h10020[0];
	outMatrix[10] = d->h10020[1];
	outMatrix[11] = d->h10020[2];
	outMatrix[12] = d->h10200[0];
	outMatrix[13] = d->h10200[1];
	outMatrix[14] = d->h10200[2];*/

/*	outMatrix[15] = d->h11001[0];
	outMatrix[16] = d->h11001[1];
	outMatrix[17] = d->h11001[2];
	outMatrix[18] = d->h00111[0];
	outMatrix[19] = d->h00111[1];
	outMatrix[20] = d->h00111[2];
	outMatrix[21] = d->h20001[0];
	outMatrix[22] = d->h20001[1];
	outMatrix[23] = d->h20001[2];
	outMatrix[24] = d->h00201[0];
	outMatrix[25] = d->h00201[1];
	outMatrix[26] = d->h00201[2];
	outMatrix[27] = d->h10002[0];
	outMatrix[28] = d->h10002[1];
	outMatrix[29] = d->h10002[2];*/

/*	outMatrix[30] = d->h10010[0];
	outMatrix[31] = d->h10010[1];
	outMatrix[32] = d->h10010[2];
	outMatrix[33] = d->h10100[0];
	outMatrix[34] = d->h10100[1];
	outMatrix[35] = d->h10100[2];*/
	
/*	outMatrix[36] = d->h22000[0];
	outMatrix[37] = d->h22000[1];
	outMatrix[38] = d->h22000[2];
	outMatrix[39] = d->h11110[0];
	outMatrix[40] = d->h11110[1];
	outMatrix[41] = d->h11110[2];
	outMatrix[42] = d->h00220[0];
	outMatrix[43] = d->h00220[1];
	outMatrix[44] = d->h00220[2];
	outMatrix[45] = d->h31000[0];
	outMatrix[46] = d->h31000[1];
	outMatrix[47] = d->h31000[2];
	outMatrix[48] = d->h40000[0];
	outMatrix[49] = d->h40000[1];
	outMatrix[50] = d->h40000[2];
	outMatrix[51] = d->h20110[0];
	outMatrix[52] = d->h20110[1];
	outMatrix[53] = d->h20110[2];
	outMatrix[54] = d->h11200[0];
	outMatrix[55] = d->h11200[1];
	outMatrix[56] = d->h11200[2];
	outMatrix[57] = d->h20020[0];
	outMatrix[58] = d->h20020[1];
	outMatrix[59] = d->h20020[2];
	outMatrix[60] = d->h20200[0];
	outMatrix[61] = d->h20200[1];
	outMatrix[62] = d->h20200[2];
	outMatrix[63] = d->h00310[0];
	outMatrix[64] = d->h00310[1];
	outMatrix[65] = d->h00310[2];
    outMatrix[66] = d->h00400[0];
	outMatrix[67] = d->h00400[1];
	outMatrix[68] = d->h00400[2];*/
	
	
	/* tune shifts with amplitude */
/*	outMatrix[69] = d->dnux_dJx;
	outMatrix[70] = d->dnux_dJy;
	outMatrix[71] = d->dnuy_dJy;*/

}                                    
