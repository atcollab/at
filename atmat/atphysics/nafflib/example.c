/*
 *  Hello World for the CodeWarrior
 *  � 1997-1998 Metrowerks Corp.
 *
 *  Questions and comments to:
 *       <mailto:support@metrowerks.com>
 *       <http://www.metrowerks.com/>
 */

#include "modnaff.h"
#include "complexe.h"

int main(void)
{
    int i, iCpt;
    const int ndata=9996; /* multiple of 6 */
    
 g_NAFVariable.DTOUR=2*M_PI; /* size of a "cadran" */
 g_NAFVariable.XH=1;         /* step */
 g_NAFVariable.T0=0;         /* time t0 */
 g_NAFVariable.NTERM=10;     /* max term to find */
 g_NAFVariable.KTABS=ndata;  /* number of data : must be a multiple of 6 */
 g_NAFVariable.m_pListFen=NULL; /*no window*/
 g_NAFVariable.TFS=NULL;    /* will contain frequency */
 g_NAFVariable.ZAMP=NULL;   /* will contain amplitude */
 g_NAFVariable.ZTABS=NULL;  /* will contain data to analyze */

 /*internal use in naf */
 g_NAFVariable.NERROR=0;
 g_NAFVariable.ICPLX=1;
 g_NAFVariable.IPRT=-1; /*0*/
 g_NAFVariable.NFPRT=stdout; /*NULL;*/
 g_NAFVariable.NFS=0;
 g_NAFVariable.IW=1;
 g_NAFVariable.ISEC=1;
 g_NAFVariable.EPSM=0; 
 g_NAFVariable.UNIANG=0;
 g_NAFVariable.FREFON=0;
 g_NAFVariable.ZALP=NULL; 
 g_NAFVariable.m_iNbLineToIgnore=1; /*unused*/
 g_NAFVariable.m_dneps=1.E100;
 g_NAFVariable.m_bFSTAB=FALSE; /*unused*/
 /*end of interl use in naf */
 
 
    naf_initnaf();
    
    /*remplit les donnees initiales*/
    for(i=0;i<ndata;i++)
    {
     g_NAFVariable.ZTABS[i].reel=2.E0+0.1*cos(M_PI*i)+0.00125*cos(M_PI/3*i);
     g_NAFVariable.ZTABS[i].imag=2.E0+0.1*sin(M_PI*i)+0.00125*sin(M_PI/3*i);
     fprintf(stdout,"%2d = % .15f % .15f\n",i,g_NAFVariable.ZTABS[i].reel
     ,g_NAFVariable.ZTABS[i].imag);
    }
    
    /*analyse en frequence*/
    /* recherche de 5 termes */
    printf("cte=%g\n",fabs(g_NAFVariable.FREFON)/g_NAFVariable.m_dneps);
    naf_mftnaf(5,fabs(g_NAFVariable.FREFON)/g_NAFVariable.m_dneps);

   /* affichage des resultats */

   printf("NFS=%d\n",g_NAFVariable.NFS);
   for(iCpt=1;iCpt<=g_NAFVariable.NFS; iCpt++)
   {
    printf("AMPL=% .15E+i*% .15E abs(AMPL)=% .15E arg(AMPL)=% .15E FREQ=% .15E\n",
           g_NAFVariable.ZAMP[iCpt].reel,g_NAFVariable.ZAMP[iCpt].imag, 
           i_compl_module(g_NAFVariable.ZAMP[iCpt]), 
           i_compl_angle(g_NAFVariable.ZAMP[iCpt]),
           g_NAFVariable.TFS[iCpt]);
   }
    /*liberation de la memoire*/
	naf_cleannaf();
	return 0;
}

