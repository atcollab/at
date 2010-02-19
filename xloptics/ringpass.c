#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef WIN32
#ifndef __MWERKS__
#include <float.h>
#define isfinite _finite
#define isnan _isnan
#endif
#elif defined __SUNPRO_C
#include <ieeefp.h>
#define isfinite finite
#endif

#include "attrack.h"

#define NTRACKS 10

static const char revision[] = "$Rev$";

static const double zerodef[] = { 0.0, 0.0, 0.0, };
static const double A0[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,};
static double E0; /* energy in eV */
static double f0; /* revolution frequency for RF cav */

/*  Tracking engine RINGPASS for Accelerator Toolbox 1.2    */

typedef struct elistelement elistelement;
typedef const struct elistelement *elist;
typedef void (*ftrack)(double *r_in, int np, const double* args, elist elem);
typedef long (*fstore)(elistelement *elem, const double *source, int availargs, const char *xtype, int nargs);

static long storeargs(elistelement *elem, const double *source, int availargs, const char *xtype, int nargs);
static long bendargs(elistelement *elem, const double *source, int availargs, const char *xtype, int nargs);
static long thinmpargs(elistelement *elem, const double *source, int availargs, const char *xtype, int nargs);
static long thickmpargs(elistelement *elem, const double *source, int availargs, const char *xtype, int nargs);
static void cleanargs(elistelement *elem);

static void tdrift(double *r_in, int np, const double *args, elist elem);
static void tqpole(double *r_in, int np, const double *args, elist elem);
static void tdpole(double *r_in, int np, const double *args, elist elem);
static void tthinmag(double *r_in, int np, const double *args, elist elem);
static void tthinmpoleN(double *r_in, int np, const double *args, elist elem);
static void tthinmpoleS(double *r_in, int np, const double *args, elist elem);
static void tthincavity(double *r_in, int np, const double *args, elist elem);
static void tthincorrector(double *r_in, int np, const double *args, elist elem);
static void tsymplecticbend(double *r_in, int np, const double *args, elist elem);
static void tsymplecticE2bend(double *r_in, int np, const double *args, elist elem);
static void tsymplecticE2Radbend(double *r_in, int np, const double *args, elist elem);
static void tsymplecticRadbend(double *r_in, int np, const double *args, elist elem);
static void tsymplecticstr(double *r_in, int np, const double *args, elist elem);
static void tsymplecticQuadFringe(double *r_in, int np, const double *args, elist elem);
static void twigglinear(double *r_in, int np, const double *args, elist elem);
static void taperture(double *r_in, int np, const double *args, elist elem);

struct elistelement {
/* double args[12];*/
   double *args;
   double *T1;
   double *T2;
   double *R1;
   double *R2;
   ftrack track;
   const char *passmethod;
};

typedef struct tracksel {
   ftrack track;
   fstore store;
   const char *passmethod;
   int v5availargs;
   int nargs;
} tracksel;

static tracksel sel[] = {
/* func,	  -,	      passname,			v5 nargs,	nb. args */
{NULL,		NULL,		"IdentityPass",		4,	0,},
{tdrift,	storeargs,	"DriftPass",		4,	1,},
{tqpole,	storeargs,	"QuadLinearPass",	4,	2,},
{tdpole,	storeargs,	"BendFringeTiltPass",	5,	5,},
{tdpole,	storeargs,	"BendLinearPass",	4,	5,},
{tthinmag,	thinmpargs,	"ThinMPolePass",	4,	8,},
{tthinmpoleN,	storeargs,	"ThinNMPolePass",	4,	2,},
{tthinmpoleS,	storeargs,	"ThinSMPolePass",	4,	2,},
{tthincorrector,storeargs,	"ThinCorrectorPass",    4,	2,},
{tsymplecticbend, bendargs,	"BndMPoleSymplectic4Pass", 5,	 12,},
{tsymplecticE2bend, bendargs,	"BndMPoleSymplectic4E2Pass", 5,	 12,},
{tsymplecticRadbend, bendargs,	"BndMPoleSymplectic4RadPass", 5,	 12,},
{tsymplecticE2Radbend, bendargs,	"BndMPoleSymplectic4E2RadPass", 5,	 12,},
{tsymplecticstr, thickmpargs,	"StrMPoleSymplectic4Pass", 4,	 12,},
{tsymplecticQuadFringe, thickmpargs,	"QuadMPoleFringePass", 4,	 12,},
{tthincavity,	storeargs,		"ThinCavityPass",	4,	4,},
{twigglinear,	storeargs,	"WiggLinearPass",	4,	2,},
{taperture,	storeargs,	"AperturePass",		4,	4,},
};

#define NB_METHODS (sizeof(sel)/sizeof(sel[0]))

static double threshold[6*NTRACKS];
static elistelement *tbl[NTRACKS];
static int n_elements[NTRACKS], mx_elements[NTRACKS];

static tracksel *trackmethod(const char *methstr)
{
   int i;
   tracksel *method = sel;
   for (i=0; i<NB_METHODS; i++) {
      if (!strcmp(methstr, method->passmethod)) return method;
      method++;
   }
   return NULL;
}

static void tdrift(double *r_in, int np, const double *args, elist elem)
{
   DriftPass(r_in, args[0], np);
}
static void tqpole(double *r_in, int np, const double *args, elist elem)
{
   QuadLinearPass(r_in,
   args[0],			/* length */
   args[1],			/* strength */
   elem->T1, elem->T2,		/* displacement */
   elem->R1, elem->R2,		/* rotation */
   np);				/* number of particles */
}
static void tdpole(double *r_in, int np, const double *args, elist elem)
{
   BendLinearPass(r_in,
   args[1]*args[2],		/* length */
   -args[3]/args[2]/args[2],	/* Ky */
   args[1],			/* angle */
   0.0,				/* deltaB/B */
   args[0],			/* entrance pole angle */
   args[4],			/* exit angle */
   0.0,				/* Fint1 */
   0.0,				/* Fint2 */
   0.0,				/* Gap */
   elem->T1, elem->T2,		/* displacement */
   elem->R1, elem->R2,		/* rotation */
   np);				/* number of particles */
}
static void tthinmag(double *r_in, int np, const double *args, elist elem)
{
   int order = (int) args[0];
   const double *B = args+1;
   ThinMPolePass(r_in,
   A0,					/* skew field */
   B,					/* normal field */
   order, 0.0, 0.0,		/* order, bending */
   elem->T1, elem->T2,	/* displacement */
   elem->R1, elem->R2,	/* rotation */
   np);					/* number of particles */
}
static void tthinmpoleN(double *r_in, int np, const double *args, elist elem)
{
   int order = (int) args[0]-1;
   double *A = (double *)calloc(order+1, sizeof(double));
   double *B = (double *)calloc(order+1, sizeof(double));
   B[order] = args[1];
   ThinMPolePass(r_in,
   A,					/* skew field */
   B,					/* normal field */
   order, 0.0, 0.0,		/* order, bending */
   elem->T1, elem->T2,	/* displacement */
   elem->R1, elem->R2,	/* rotation */
   np);					/* number of particles */
   free(A);
   free(B);
}
static void tthinmpoleS(double *r_in, int np, const double *args, elist elem)
{
   int order = (int) args[0]-1;
   double *A = (double *)calloc(order+1, sizeof(double));
   double *B = (double *)calloc(order+1, sizeof(double));
   A[order]=args[1];
   ThinMPolePass(r_in,
   A,					/* skew field */
   B,					/* normal field */
   order, 0.0, 0.0,		/* order, bending */
   elem->T1, elem->T2,	/* displacement */
   elem->R1, elem->R2,	/* rotation */
   np);					/* number of particles */
   free(A);
   free(B);
}
static void tthincavity(double *r_in, int np, const double *args, elist elem)/*(l,v,h)->(l,v/E,h f_0)*/
{
   CavityPass(r_in, args[0], args[1]/(1000*E0), args[2]*f0, 0, np);
}
static void tthincorrector(double *r_in, int np, const double *args, elist elem)
{
   CorrectorPass(r_in, args[0], args[1], 0.0, np);
}
static void tsymplecticbend(double *r_in, int np, const double *args, elist elem)
{
   int order = (int) args[4];
   const double *B = args+5;
   BndMPoleSymplectic4Pass(r_in,
   args[0],			/* length */
   args[1],			/* 1/rho */
   A0,				/* multipoles */
   B,				/* multipoles */
   order,			/* max order */
   10,				/* number of slices */
   args[2],			/* entrance pole angle */
   args[3],			/* exit angle */
   0.0,				/* Fint1 */
   0.0,				/* Fint2 */
   0.0,				/* Gap */
   elem->T1, elem->T2,		/* displacement */
   elem->R1, elem->R2,		/* rotation */
   np);				/* number of particles */
}
static void tsymplecticE2bend(double *r_in, int np, const double *args, elist elem)
{
   int order = (int) args[4];
   const double *B = args+5;
   BndMPoleSymplectic4E2Pass(r_in,
   args[0],			/* length */
   args[1],			/* 1/rho */
   A0,				/* multipoles */
   B,				/* multipoles */
   order,			/* max order */
   10,				/* number of slices */
   args[2],			/* entrance pole angle */
   args[3],			/* exit angle */
   0.0,				/* Fint1 */
   0.0,				/* Fint2 */
   0.0,				/* Gap */
   0.0,				/* h1 */
   0.0,				/* h2 */
   elem->T1, elem->T2,		/* displacement */
   elem->R1, elem->R2,		/* rotation */
   np);				/* number of particles */
}

static void tsymplecticRadbend(double *r_in, int np, const double *args, elist elem)
{
   int order = (int) args[4];
   const double *B = args+5;
   BndMPoleSymplectic4RadPass(r_in,
   args[0],			/* length */
   args[1],			/* 1/rho */
   A0,				/* multipoles */
   B,				/* multipoles */
   order,			/* max order */
   10,				/* number of slices */
   args[2],			/* entrance pole angle */
   args[3],			/* exit angle */
   0.0,				/* Fint1 */
   0.0,				/* Fint2 */
   0.0,				/* Gap */
   elem->T1, elem->T2,		/* displacement */
   elem->R1, elem->R2,		/* rotation */
   1e9*E0,             /* energy  */
   np);				/* number of particles */
}

static void tsymplecticE2Radbend(double *r_in, int np, const double *args, elist elem)
{
   int order = (int) args[4];
   const double *B = args+5;
   BndMPoleSymplectic4E2RadPass(r_in,
   args[0],			/* length */
   args[1],			/* 1/rho */
   A0,				/* multipoles */
   B,				/* multipoles */
   order,			/* max order */
   10,				/* number of slices */
   args[2],			/* entrance pole angle */
   args[3],			/* exit angle */
   0.0,				/* Fint1 */
   0.0,				/* Fint2 */
   0.0,				/* Gap */
   0.0,				/* h1 */
   0.0,				/* h2 */
   elem->T1, elem->T2,		/* displacement */
   elem->R1, elem->R2,		/* rotation */
   1e9*E0,             /* energy  */
   np);				/* number of particles */
}

static void tsymplecticstr(double *r_in, int np, const double *args, elist elem)
{
   int order = (int) args[1];
   const double *B = args+2;
   StrMPoleSymplectic4Pass(r_in,
   args[0],			/* length */
   A0,				/* multipoles */
   B,				/* multipoles */
   order,			/* max order */
   10,				/* number of slices */
   elem->T1, elem->T2,		/* displacement */
   elem->R1, elem->R2,		/* rotation */
   np);				/* number of particles */
}

static void tsymplecticQuadFringe(double *r_in, int np, const double *args, elist elem)
{
   int order = (int) args[1];
   const double *B = args+2;
   QuadMPoleFringePass(r_in,
   args[0],			/* length */
   A0,				/* multipoles */
   B,				/* multipoles */
   order,			/* max order */
   10,				/* number of slices */
   elem->T1, elem->T2,		/* displacement */
   elem->R1, elem->R2,		/* rotation */
   np);				/* number of particles */
}

static void twigglinear(double *r_in, int np, const double *args, elist elem)
{
   WiggLinearPass(r_in,
   args[0],			/* length */
   args[1],			/* 1/rho */
   0.0,			        /* kxkz */
   elem->T1, elem->T2,		/* displacement */
   elem->R1, elem->R2,		/* rotation */
   np);				/* number of particles */
}

static void taperture(double *r_in, int np, const double *args, elist elem)
{
   AperturePass(r_in, args, np);
}

static void setshift(elistelement *elem, double dx, double dz)
{
   double *t1, *t2;
   if ((dx != 0.0) || (dz != 0.0)) {
      t1 = (double *) calloc(6, sizeof(double));
      t2 = (double *) calloc(6, sizeof(double));
      t1[0] = -dx;
      t1[2] = -dz;
      t2[0] = dx;
      t2[2] = dz;
   }
   else {
      t1 = NULL;
      t2 = NULL;
   }
   elem->T1 = t1;
   elem->T2 = t2;
}
static void settilt(elistelement *elem, double tilt)
{
   double *r1, *r2;
   if (tilt != 0.0) {
      double C = cos(tilt);
      double S = sin(tilt);
      r1=calloc(36, sizeof(double));
      r2=calloc(36, sizeof(double));
      r1[0] = C; r2[0] = C;
      r1[7] = C; r2[7] = C;
      r1[14] = C; r2[14] = C;
      r1[21] = C; r2[21] = C;
      r1[28] = 1.0; r2[28] = 1.0;
      r1[35] = 1.0; r2[35] = 1.0;
      r1[2] = -S; r2[2] = S;
      r1[9] = -S; r2[9] = S;
      r1[12] = S; r2[12] = -S;
      r1[19] = S; r2[19] = -S;
   }
   else {
      r1 = NULL;
      r2 = NULL;
   }
   elem->R1 = r1;
   elem->R2 = r2;
}
static void zeroargs(elistelement *elem, const double *args, int nargs)
{
   elem->args = (nargs <= 0) ? NULL : (double*)calloc(nargs, sizeof(double));
   setshift(elem, 0.001*args[0], 0.001*args[1]);
   settilt(elem, 0.001*args[2]);
}
static long fillargs(const double *source, int availargs, double *dest, int nargs)
{
   if (availargs > nargs) availargs = nargs;
   memcpy(dest, source, availargs*sizeof(double));
   return 0;
}

static long storeargs(elistelement *elem, const double *source, int availargs, const char *xtype, int nargs)
{
   return fillargs(source, availargs, elem->args, nargs);
}

static int max_order(const double *args, int nargs)
{
	int ordermax = nargs-1;
	while ((args[ordermax] == 0) && (ordermax-- > 0)) ;
	if (ordermax < 0) ordermax = 0;
	return ordermax;
}
static long bendargs(elistelement *elem, const double *source, int availargs, const char *xtype, int nargs)
{
   long ret = 0;
   double *dest = elem->args;
   double *B = dest+5;
   register double ro = source[2];
   dest[0] = ro * source[1];	/* length */
   dest[1] = 1.0 / ro;			/* radius */
   dest[2] = source[0];			/* entrance face */
   dest[3] = source[4];			/* exit face */
   B[1] =  -source[3]/ro/ro;	/* 4 poles */
   B[2] = source[5]/ro/ro/ro;	/* 6 poles */
   if (availargs > 7)  {
      B[4] = source[6];			/* 10 poles */
      B[6] = source[7];			/* 14 poles */
   }
   dest[4] = max_order(B, nargs-5);
   return ret;
}
static long thinmpargs(elistelement *elem, const double *source, int availargs, const char *xtype, int nargs)
{
   long ret;
   double *dest = elem->args;
   double *B = dest+1;
   if ((xtype == NULL) || (strcmp(xtype, "SEXT") == 0)) {	/* SEXT */
      B[2] = *source++;
      ret = 0;
   }
   else {													/* THINMAG */
      ret = fillargs(source, availargs, B+1, nargs-2);	/* remaining args */
   }
   dest[0] = max_order(B, nargs-1);
   return ret;
}
static long thickmpargs(elistelement *elem, const double *source, int availargs, const char *xtype, int nargs)
{
   long ret;
   double *dest = elem->args;
   double *B = dest+2;
   dest[0] = source[0];			/* length */
   ret = fillargs(source+1, availargs-1, B+1, 3);	/* remaining args */
   if (availargs > 5)  {
      B[5] = source[4];			/* 12 poles */
      B[9] = source[5];			/* 20 poles */
   }
   dest[1] = max_order(B, nargs-2);
   return ret;
}

static void cleanargs(elistelement *elem)
{
   free(elem->args);
   free(elem->T1);
   free(elem->T2);
   free(elem->R1);
   free(elem->R2);
}

static long addelem(long tptr, const char *xtype, const tracksel* method, const double *args, int availargs, const double *defs)
{
/* tptr = 0;*/
   long ret = -1;
   if (method && (tptr >= 0) && (tptr < NTRACKS)) {
	elistelement *next = tbl[tptr] + n_elements[tptr];
	next->track = method->track;
	next->passmethod = method->passmethod;
	zeroargs(next, defs, method->nargs);
	if (method->store) ret = method->store(next, args, availargs, xtype, method->nargs);
	n_elements[tptr]++;
  }
   return ret;
}

long STDCALL ataddtypeelem(long tptr, long nargs, const double *args, const char* passmethod, const char *xtype)
{
   return addelem(tptr,  xtype, trackmethod(passmethod), args+3, nargs-3, args);
}

long STDCALL dbgatarg(long tptr, const char* passmethod, const double *args)
{
   const tracksel* method = trackmethod(passmethod);
   int availargs = method ? method->v5availargs : 0;
   return addelem(tptr, NULL, method, args, method->v5availargs, zerodef);
}

double STDCALL atgetarg(long tptr, long nel, long np)
{
/* tptr = 0;*/
   if ((tptr >= 0) && (tptr < NTRACKS))
      return (nel < n_elements[tptr]) ? tbl[tptr][nel].args[np] : 888;
   else
      return 888;
}

static long ptrack(elist firstel, int nel, const double *thresh, double *r_in, long np, long nturns)
{
	long turn;
	int el,j;
	for (turn=0; turn<nturns; turn++) {
		elist next=firstel;
		for (el=0; el<nel; el++) {
			if (next->track) next->track(r_in, np, next->args, next);
			for (j=0; j<6; j++) {
				if (isnan(r_in[j])) return -1;
				if (!isfinite(r_in[j])) return -2;
				if (fabs(r_in[j]) >= thresh[j]) return turn+1;
			}
			next++;
		}
	}
	return 0;
}

long STDCALL attrack(long tptr, double *r_in, long np, long nturns)
{
	long ret = 0;
	if ((tptr >= 0) && (tptr < NTRACKS)) {
		ret=ptrack(tbl[tptr], n_elements[tptr], threshold+6*tptr, r_in, np, nturns);
	}
	return ret;
}

long STDCALL atstable(long tptr, const double *xl_in, long nturns)
{
	long ret = 0;
	if ((tptr >= 0) && (tptr < NTRACKS)) {
		double r_in[6];
		memcpy(r_in, xl_in, 6*sizeof(double));
		ret=ptrack(tbl[tptr], n_elements[tptr], threshold+6*tptr, r_in, 1, nturns);
	}
	return ret;
}


long STDCALL atinitialize2(long tptr, long nelems, const double *thresh, const double *energy, const double *freq)
{
   static long cpt = 0;
   if ((tptr >= 0) && (tptr < NTRACKS)) {
	  int i, n;
      double *threshptr = threshold + 6*tptr;
	  elistelement *next=tbl[tptr];
      for (n=0; n<n_elements[tptr]; n++) cleanargs(next++);
      tbl[tptr]=realloc(tbl[tptr], nelems*sizeof(elistelement));
      mx_elements[tptr] = nelems;
      for (i=0; i<6; i++) *threshptr++ = *thresh++;
      n_elements[tptr] = 0;
   }
	E0= *energy;
	f0= *freq;
   return cpt++;
}


/* Kept for backward compatibilty */

long STDCALL atinitialize(long tptr, long nelems, const double *thresh)
{
   return atinitialize2(tptr, nelems, thresh, zerodef, zerodef);
}

long STDCALL xatsize(long nelems, double thresh)
{
   double t6[6];
   int i;
   for (i=0; i<6; i++) t6[i] = thresh;
   return atinitialize2(0L, nelems, t6, zerodef, zerodef);
}

long STDCALL atsize(long tptr, long nelems, double thresh)
{
   double t6[6];
   int i;
   for (i=0; i<6; i++) t6[i] = thresh;
   return atinitialize2(tptr, nelems, t6, zerodef, zerodef);
}

long STDCALL xataddelem(long code, const double *args)
{
   const tracksel* method = ((code >= 0) && (code < NB_METHODS)) ? sel+code : NULL;
   int availargs = method ? method->v5availargs : 0;
   return addelem(0L, NULL, method, args, availargs, zerodef);
}

long STDCALL ataddcodeelem(long tptr, long code, const double *args)
{
   const tracksel* method = ((code >= 0) && (code < NB_METHODS)) ? sel+code : NULL;
   int availargs = method ? method->v5availargs : 0;
   return addelem(tptr, NULL, method, args,  availargs, zerodef);
}

long STDCALL ataddpasselem(long tptr, const char* passmethod, const double *args)
{
   const tracksel* method = trackmethod(passmethod);
   int availargs = method ? method->v5availargs : 0;
   return addelem(tptr, NULL, method, args, availargs, zerodef);
}

long STDCALL xattrack(double *r_in, long np, long nturns)
{
   return attrack(0L, r_in, np, nturns);
}

long STDCALL atversion(void)
{
   long revnumber = 0L;
   int ok = sscanf(revision, "$Rev: %ld", &revnumber);
   return revnumber;
}
