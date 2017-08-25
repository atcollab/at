#include <stdlib.h>
#include <stdio.h>
#define _USE_MATH_DEFINES	/* For Visual Studio */
#include <math.h>
#include <float.h>

/*--------*/
/* define */
/*--------*/
#define BOOL int
#define FALSE 0
#define TRUE  1
#define Myyerror printf
#define MAX(x,y) ((x)<(y)?(y):(x))


/*#ifdef FUNCINLINEEXTERN_MUSTBEINLIB*/
#define INLINE_EXTERN static inline
/*#else 
#define INLINE_EXTERN inline extern
#endif *//*FUNCINLINEEXTERN_MUSTBEINLIB*/

 
/*---------------------*/
/*routine d'allocation */
/*---------------------*/
#define SYSCHECKMALLOCSIZE(variable, type, size) \
  if ((variable=(type*)malloc(sizeof(type)*(size)))==NULL)\
  { printf("error : malloc failed!\n"); exit(1);}

#define SYSCHECKMALLOC(variable, type) \
  if ((variable=(type*)malloc(sizeof(type)))==NULL)\
  { printf("error : malloc failed!\n"); exit(1);}

#define DIM2(prow,row,col,type,l) \
  {\
  register type *pdata;\
  int I;\
  SYSCHECKMALLOCSIZE(pdata, type,(row) * (col));\
  SYSCHECKMALLOCSIZE(prow, type *,(row));\
  for (I = 0; I < (row); I++)\
     {\
     prow[I] = pdata;\
     pdata += col;\
     }\
  }

#define SYSFREE(variable) free(variable)
#define HFREE2(variable) {SYSFREE(*variable); SYSFREE(variable);}

/*---------*/
/*structure*/
/*--------*/
struct complexe {
                double reel,imag;
                };
typedef struct complexe t_complexe;

/* v0.96 M. GASTINEAU 18/12/98 : ajout */
/*liste des frequences pour NAF */
struct list_fenetre_naf
{
 double                   dFreqMin; /* frequence minimale */
 double                   dFreqMax; /* frequence maximale */
 int                      iNbTerme; /* nombre de termes a rechercher */
 struct list_fenetre_naf *suivant; /*fenetre suivante */
};
/* v0.96 M. GASTINEAU 18/12/98 : fin ajout */

typedef struct list_fenetre_naf t_list_fenetre_naf; /* v0.96 M. GASTINEAU 18/12/98 : ajout */   

/*v0.96 M. GASTINEAU 04/09/98 : ajout pour la gestion de naf */
/* pour le role de ces champs, cf. modnaff.c */
struct stnaf 
{
 /*champ utilise par modnaff.c */
 FILE *NFPRT;
 double EPSM; 
 int NTERM,KTABS,NFS;
 int ICPLX,IW,ISEC;
 int NERROR,IPRT;

 double *TFS;
 t_complexe *ZAMP;
 t_complexe **ZALP; /*tableau a deux dimensions*/
 t_complexe *ZTABS;

 double DTOUR,UNIANG,FREFON;
 double XH,T0;
 
 /*autre champ utilise en tant que flag */
 double m_dneps; /*equivaut a DNEPS */
 int m_iNbLineToIgnore; /*nombre de lignes a ignorer en debut du fichier des solutions */
 BOOL m_bFSTAB; /*=TRUE => sauve le tableau ZTABS.*/
 /* v0.96 M. GASTINEAU 06/01/99 : ajout */
 t_list_fenetre_naf *m_pListFen; /*liste des fenetres */
 /* v0.96 M. GASTINEAU 06/01/99 : fin ajout */
};

typedef struct stnaf t_naf;
/*v0.96 M. GASTINEAU 04/09/98 : fin ajout */


/*-----------------*/
/*variable globale */
/*-----------------*/
extern t_naf g_NAFVariable;
extern double pi;


/*-----------------*/
/* public functions*/
/*-----------------*/
void naf_initnaf_notab();
void naf_cleannaf_notab();
void naf_initnaf();
void naf_cleannaf();
BOOL naf_mftnaf(int NBTERM, double EPS);
void naf_prtabs(int KTABS, t_complexe *ZTABS, int IPAS);
