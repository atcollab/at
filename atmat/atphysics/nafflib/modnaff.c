/* MODNAFF.C version 0.96                                                 */
/* Astronomie et systemes dynamiques                                      */
/* J. LASKAR Bureau des Longitudes, Paris (version originale fortran)     */
/* M. GASTINEAU, portage en C, 03/09/98                                   */
/* v0.96 M. GASTINEAU 09/09/98 : modification dans la fonction naf_ztder  */
/*                               dans les boucles de I en I-1.            */
/* v0.96 M. GASTINEAU 09/11/98 : remplacement des variables PI2 par PICARRE*/
/*                               car PI2 est un define sous LINUX.        */
/* v0.96 M. GASTINEAU 01/12/98 : utilisation des fonctions  complexes     */
/*                               inlines.                                 */
/* v0.96 M. GASTINEAU 07/12/98 : modification dans naf_frefin car bug lors*/
/*                               de l'optimisation sur Mac.               */
/* v0.96 M. GASTINEAU 18/12/98 : modification de naf_iniwin               */
/* v0.96 M. GASTINEAU 06/01/99 : ajout de la liberation de la liste des   */
/*                        fenetres dans naf_cleannaf et naf_cleannaf_notab*/
/* v0.96 M. GASTINEAU 12/01/99 : correction dans le cas ou ICPLX=0 dans   */
/*                        naf_modfre et naf_gramsc.                       */
/*v0.96 M. GASTINEAU 14/01/99 : ajout du support des fenetres dans        */
/*                        naf_mfttab et naf_fftmax.                       */
/*v0.97 M. GASTINEAU 26/05/99 : correction bug (0.97/99/05/26/A) dans     */
/*                        naf_mftnaf                                      */
/* v0.97 M. GASTINEAU 27/05/99 : correction bug (0.97/99/05/27/A) dans le */
/*                  cas ou ICPLX=0 et FS=0 dans naf_modfre et naf_gramsc. */

#define FUNCINLINEEXTERN_MUSTBEINLIB
#include "modnaff.h"
#include "complexe.h"

/*variable globale de naf */
t_naf g_NAFVariable;

/*legere difference entre NAF_USE_OPTIMIZE=0 et NAF_USE_OPTIMIZE=1*/
/* car division differente dans naf_gramsc */
#define NAF_USE_OPTIMIZE 1

/*
!-------------------------------------------------------------------------
!   NAFF.0.84 NUMERICAL ANALYSIS OF FUNDAMENTAL FREQUENCIES
!            (27 JANVIER 1996)
!  
!   (C) JACQUES LASKAR 
!       ASTRONOMIE ET SYSTEMES DYNAMIQUES
!       BUREAU DES LONGITUDES
!       75006 PARIS
!       EMAIL : LASKAR@BDL.FR
!
!    MAIN REFERENCES : 
!    LASKAR, J.: THE CHAOTIC MOTION OF THE SOLAR SYSTEM. A NUMERICAL
!    ESTIMATE OF THE SIZE OF THE CHAOTIC ZONES,ICARUS,88,(1990),266--291
!
!   (NAFF.MAC.1.6 - 27/5/98)
!
!**************************************************************************
!  THIS PROGRAMM  CANNOT BE COPIED, 
!  DISTRIBUTED NOR MODIFIED WITHOUT THE AGREEMENT OF THE AUTHOR.
!**************************************************************************
! 
!                        PROGRAMME D'ANALYSE DE FOURIER AUTOMATIQUE
!   -- NAFF --           CALCULE UN TERME A LA FOIS ET LE RETRANCHE
!                        A TABS
!                        NOUVELLE VERSION 26/9/87
!                        MODIFIEE POUR NFS LE 18/10/87
!                        MODIFIEE LE 10/4/88 POUR TRAITER
!                        A LA FOIS UN TABLEAU REEL OU COMPLEXE
!                        CALCUL BATCH SUR CIRCE
! 8/10/90 MODIF POUR REMETTRE SUR VP (FORTRAN STANDARD)
! 13/12/90   MODIF POUR RENDRE PLUS MODULAIRE
! 18/12/90   MODIF IL N'Y A PLUS DE DIMENSIONS A REGLER DANS
!            LES SOUS PROGRAMMES. TABS(2,*) AU LIEU DE TABS(2,0:KTABS)
!            NFR EST ENCORE EN PARAMETER, MAIS EST RAREMENT A CHANGER.
!
!            IL FAUT REMPLACER LE * DES TABLEAUX
!            DIMENSIONES T(*) OU T(N,*) PAR T(1) OU T(N,1)
! 19/12/90
!  16/4/91   COMPILATION SEPAREE DE NAFF EN UNE BIBLIOTHEQUE
!  30/4/91   TRAITEMENT DU CAS OU ON RETROUVE LA MEME FREQUENCE
!  27/1/96   PLUSIEURS FENETRES
!  27/2/97   MODIF POUR DIMENSIONNER TOUT EXPLICITEMENT AVEC 
!            INCLUDE NAFF.INC POUR LES PARAMETRES
!  23/45/1997    CORRECTION DE MAXIQUA.F (F. JOUTEL)
!  22/4/1998 MODIF POUR QUADRUPLE PRECISION (J. LASKAR)
!            CHANGEMENT POUR UN EPSILON MACHINE DONNE PAR NAFF.INC
!  27/5/98   AMELIORATION DE LA RECHERCHE DU MAXIMUM PAR
!            LA METHODE DES SECANTES
!            ROUTINE DE CORRECTION DE LA FREQUENCE PRINCIPALE PAR
!            DEVELOPPEMENT ASYMTOTIQUE (F. JOUTEL)
!***********************************************************************
!          MODULE  NAFF
!----------------------------------------------------------------------
!  PROCEDURE D'UTILISATION:
!  - PHASE D'INITIALISATION
!      1)  INITIALISER LES VARIABLES :
!              DTOUR,KTABS,XH,NTERM,IW,T0,ICPLX,ISEC,NFPRT,IPRT
!      2)  CALL INITNAF
!             CALCULE UNIANG,FREFON, EPSM, PI 
!             CREE LES TABLEAUX DYNAMIQUES ((ZTABS)), ((TFS)), ((ZAMP)),
!	     ((ZALP)) et((TWIN))
!             APPELLE INIWIN
!  - PHASE D'UTILISATION
!      1)  EMPLIR LE TABLEAU ZTABS
!      2)  CALL MFTNAF(NBTERM,EPS)
!          ANALYSE ((ZTABS)) ET ESSAIE DE RECHERCHER (NBTERM) FREQUENCES, 
!          (EPS) ETANT L'ERREUR TOLEREE
!          --! MODIFIE ((ZTABS))
!  - PHASE D'EXPLOITATION OU DE MODIFICATION
!     * CALL PRTABS(KTABS,ZTABS,IPAS)
!         IMPRIME LES ELEMENTS DU TABLEAU COMPLEXE (ZTABS(0:KTABS)) TOUS 
!         LES (IPAS), ((IPRT)) DOIT ETRE MIS A 1 AVANT LA PROCEDURE POUR
!         IMPRIMER QUELQUE CHOSE 
!     * CALL SMOY(ZM)
!         CALCULE LA MOYENNNE DES ELEMENTS DE ((ZTABS)) ET MODIFIE CE  
!         TABLEAU EN SOUSTRAYANT DE CHACUN DES ELEMENTS CETTE VALEUR 
!         MOYENNE (ZM)
!         --! MODIFIE ((ZTABS))
!     * CALL TESSOL (EPS,TFSR,ZAMPR)
!         TESTE LA COHERENCE DES SOLUTIONS OBTENUES PAR NAFF:
!          -  RECONSTRUIT LE SIGNAL A PARTIR DES ((NFS)) TERMES TROUVES 
!             QUI SONT DANS ((TFS)) ET ((ZAMP))
!          -  ANALYSE CE NOUVEAU SIGNAL, (EPS) ETANT L'ERREUR TOLEREE,
!             ET RANGE LES NOUVEAUX TERMES DANS (TFSR) ET (ZAMPR)
!          --! MODIFIE ((ZTABS))
!     * CALL INIFRE
!         REMET A ZERO ((TFS)), ((ZAMP)),((ZALP)) et ((NFS))
!         UTILE QUAND ON BOUCLE SUR PLUSIEURS CAS
!     * CALL CORRECTION (FREQ)
!        UTILISE LES TERMES TROUVES POUR CORRIGER LA FREQUENCE
!        PRINCIPALE, LA NOUVELLE VALEUR CORRIGEE EST (FREQ)
!  - PHASE DE SORTIE (POUR RECUPERER DE LA PLACE OU POUR REINITIALISER)
!      1) CALL CLEANNAF
!         DESALLOUE LES TABLEAUX DYNAMIQUES ((ZTABS)), ((TFS)), ((ZAMP)),
!	  ((ZALP)) et ((TWIN))
!-------------------------------------------------------------------------
!  VARIABLES INDISPENSABLES
!  NTERM                : LE NOMBRE MAXIMAL DE FREQUENCES A RECHERCHER
!  KTABS                : LE NOMBRE D'INTERVALLES ENTRE LES DONNEES
!                         (IL Y A DONC KTABS+1 DONNEES)
!  XH                   : PAS ENTRE DEUX DONNEES 
!                         DANS UNE UNITE DE TEMPS QUELCONQUE
!                         LA DUREE DE L'INTERVALLE EST ALORS XH!KTABS
!  T0                   : DATE DE DEPART DES DONNEES
!  DTOUR                : LONGUEUR D'UN TOUR DE CADRAN
!  ICPLX                :  0 SI LA FONCTION EST REELLE, 1 SINON
!  IW                   : PARAMETRE DE LA FENETRE D'INTEGRATION
!                            IW=0  PAS DE FENETRE
!                            IW=1  FENETRE DE HANNING  etc
!  ISEC                 : DRAPEAU DE LA METHODE DES SECANTES
!                           ISEC=1 oui,  ISEC=0 non
!  NFPRT                : NUMERO DE LA CONNECTION PUR L'IMPRESSION
!  IPRT                 : INDIQUE LE TYPE D'IMPRESSION, EN GENERAL :
!                         -1 PAS D'IMPRESSION
!                         O  MESSAGES D'ERREURS EVENTUELLES
!                         1  RESULTATS
!                         2  (DEBUG)
!  TABLEAU OU RANGER LE SIGNAL
!  ZTABS(0:KTABS)       : TABLEAU DE KTABS+1 COMPLEXES
!
!  VARIABLES CALCULEES PAR NAFF
!  UNIANG              : UNITE D'ANGLE = 1RD=UNIANG UNITE D'ANGLE
!                          (SI L'UNITE EST EN SECONDES D'ARC
!                              UNIANG=206264.80624709D0 )
!                        CALCULE EN FONCTION DE DTOUR
!  FREFON              : FREQUENCE FONDAMENTALE
!                           FREFON=2*PI/(KTABS*XH) OU EN SECONDES
!                           FREFON=360.D0*3600.D0/(KTABS*XH)  
!  NFS                 : NOMBRE DE FREQUENCES ENTIEREMENT DETERMINEES
!                        AINSI QUE LEUR AMPLITUDE
!  TFS(NBTERM)         : TABLEAU DES FREQUENCES
!  ZAMP(NBTERM)        : AMPLITUDE COMPLEXE DES TERMES
!  ZALP(NBTERM,NBTERM) : TABLEAU DE CHANGEMENT DE BASE
!
!-----------------------------------------------------------------------*/
static double AF,BF;
static double *TWIN=NULL;

/*!-----------------------------------------------------------------------	  
! VARIABLES A INITIALISER PAR L'UTILISATEUR AVANT DE LANCER INITNAF
          PUBLIC :: DTOUR,KTABS,XH,NTERM,IW,NFPRT,IPRT;
! VARIABLES PUBLIQUES INITIALISEES PAR INITNAF
          PUBLIC :: EPSM,PI,UNIANG,FREFON;
! TABLEAUX PUBLICS CREES PAR INITNAF ET DETRUITES PAR CLEANNAFF
          PUBLIC ::  TFS,ZAMP,ZTABS,ZALP;
!-----------------------------------------------------------------------	  
! VARIABLES  A INITIALISER PAR L'UTILISATEUR AVANT DE LANCER NAFF
          PUBLIC :: T0,ICPLX,ISEC
! VARIABLES INITIALISES ET UTILISEES PAR NAFF
          PUBLIC :: NERROR,NFS
!-----------------------------------------------------------------------
! ROUTINES PUBLIQUES
          PUBLIC :: INITNAF,CLEANNAF,MFTNAF,PRTABS,SMOY,TESSOL
          PUBLIC :: INIFRE,CORRECTION*/
/*v0.96 M. GASTINEAU 06/10/98 : ajout */
void naf_initnaf_notab();
void naf_cleannaf_notab();
/*v0.96 M. GASTINEAU 06/10/98 : fin ajout */
void naf_initnaf();
void naf_inifre();
void naf_cleannaf();
BOOL naf_mftnaf(int NBTERM, double EPS);
void naf_prtabs(int KTABS, t_complexe *ZTABS, int IPAS);
void naf_smoy(t_complexe *ZM);
BOOL naf_tessol(double EPS, double *TFSR, t_complexe *ZAMPR);
void naf_correction(double *FREQ);
void naf_four1(double *DATA /*tableau commencant a l'indice 1 */,
               int NN, int ISIGN);
void naf_puiss2(int NT, int *N2);
/*v0.96 M. GASTINEAU 18/12/98 : modification du prototype */
/*void naf_iniwin();*//*remplacee par: */
void naf_iniwin(double *p_pardTWIN);
/*v0.96 M. GASTINEAU 18/12/98 : fin modification */
void delete_list_fenetre_naf(t_list_fenetre_naf *p_pListFenNaf);
t_list_fenetre_naf *concat_list_fenetre_naf(t_list_fenetre_naf *p_pListFenHead,
                                            t_list_fenetre_naf *p_pListFenEnd);
t_list_fenetre_naf *cree_list_fenetre_naf(const double p_dFreqMin, 
                                          const double p_dFreqMax, 
                                          const int p_iNbTerm);

/*!-----------------------------------------------------------------------	  
! VARIABLES PRIVEES
          PRIVATE :: TWIN,AF,BF                   
! ROUTINES PRIVEES
          PRIVATE :: INIWIN,FRETES,ZTPOW,FFTMAX,FOUR1,PUISS2,MAXX        
          PRIVATE :: MODFRE, GRAMSC,PROSCA,SECANTES,MAXIQUA
          PRIVATE :: FUNC,FUNCP,PRODER,ZTDER,FREFIN,PROFRE
          PRIVATE :: ZTPOW2,ZARDYD,PROSCAA,ZTPOW2A,MODTAB*/
static void naf_fretes(double FR, int *IFLAG, double TOL, int * NUMFR);
static void naf_ztpow(int N, int N1, t_complexe *ZT, t_complexe ZA, t_complexe ZAST);
/* v0.96 M. GASTINEAU 12/01/99 : optimisation */
#if NAF_USE_OPTIMIZE==0
static void naf_maxx(int N, double *T, int *INDX);
static void naf_fftmax(double *FR);
#else /*remplacee par:*/
static int naf_maxx(int N, double *T);
static double naf_fftmax(int p_iFrMin, int p_iFrMax, double FREFO2, int KTABS2); 
#endif /*NAF_USE_OPTIMIZE*/
/* v0.96 M. GASTINEAU 12/01/99 : fin optimisation */
static void naf_modtab(int N, double *T);
static void naf_modfre(int NUMFR, double *A, double *B);
static BOOL naf_gramsc(double FS, double A, double B);
static void naf_prosca(double F1, double F2, t_complexe* ZP);
static void naf_proder(double FS, double *DER, double *A, double *B,double *RM);
static double naf_funcp(double X);
static void naf_ztder(int N, int N1, t_complexe *ZTF, t_complexe *ZTA, double *TW, t_complexe ZA, t_complexe ZAST, double T0, double XH);
static void naf_secantes(double X, double PASS, double EPS, double *XM, int IPRT, FILE *NFPRT);
static double naf_func(double X);
static void naf_maxiqua(double X, double PASS, double EPS, double *XM, double *YM, int IPRT, FILE *NFPRT);
static void naf_frefin(double *FR, double *A, double *B, double *RM, const double RPAS0, const double RPREC);
static void naf_ztpow2(int N, int N1, t_complexe *ZTF, t_complexe *ZTA, double *TW, t_complexe ZA, t_complexe ZAST);
static BOOL naf_profre(double FS, double *A, double *B, double *RMD);
static BOOL naf_proscaa(double F1, double F2, t_complexe *ZP);
static BOOL naf_zardyd(t_complexe *ZT, int N, double H, t_complexe *ZOM);
static void naf_ztpow2a(int N, int N1, t_complexe *ZTF, double *TW, t_complexe ZA, t_complexe ZAST);
          
/*!------------------------------------------------------------------------          
          CONTAINS*/
 
 
/*!*/

/*!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!           ROUTINES DE NAFF
!-----------------------------------------------------------------------
      SUBROUTINE INITNAF
      IMPLICIT NONE
!-----------------------------------------------------------------------
!        INITNAFF
!  - effectue les initialisations necessaires :  PI,UNIANG,FREFON,EPSM
!  - alloue les tableaux dynamiques TFS,ZAMP,ZALP,ZTABS,TWIN
!  - appelle INIWIN
!-----------------------------------------------------------------------*/
void naf_initnaf()
{  
/*!----------------- PREMIERES INITIALISATIONS*/
      g_NAFVariable.EPSM = DBL_EPSILON;
      /*PI = ATAN2(1.D0,0.D0)*2*/
      g_NAFVariable.UNIANG = g_NAFVariable.DTOUR/(2*M_PI) ;
      g_NAFVariable.FREFON = g_NAFVariable.DTOUR/(g_NAFVariable.KTABS*g_NAFVariable.XH)	;      
      SYSCHECKMALLOCSIZE(g_NAFVariable.TFS, double, g_NAFVariable.NTERM+1);/*allocate(TFS(1:NTERM),stat = NERROR)*/
      SYSCHECKMALLOCSIZE(g_NAFVariable.ZAMP, t_complexe, g_NAFVariable.NTERM+1); /*allocate(ZAMP(1:NTERM),stat = NERROR)*/
      DIM2(g_NAFVariable.ZALP, (g_NAFVariable.NTERM+1), (g_NAFVariable.NTERM+1), t_complexe,"ZALP"); /*allocate(ZALP(1:NTERM,1:NTERM),stat = NERROR)*/
      SYSCHECKMALLOCSIZE(g_NAFVariable.ZTABS, t_complexe, g_NAFVariable.KTABS+1);/* allocate(ZTABS(0:KTABS),stat = NERROR)*/
      SYSCHECKMALLOCSIZE(TWIN, double, g_NAFVariable.KTABS+1); /*allocate(TWIN(0:KTABS),stat = NERROR)*/
      /*v0.96 M. GASTINEAU 18/12/98 : modification du prototype */
      /*naf_iniwin();  */ naf_iniwin(TWIN);     
}/*      end SUBROUTINE INITNAF*/
/*!-----------------------------------------------------------------------
      subroutine CLEANNAF*/
/* v0.96 M. GASTINEAU 06/01/99 : ajout de la liberation de la liste des fenetres */
void naf_cleannaf()
{
/*!-----------------------------------------------------------------------
!     desalloue les tableaux dynamiques
!-----------------------------------------------------------------------*/
 
      SYSFREE(g_NAFVariable.TFS);
      SYSFREE(g_NAFVariable.ZAMP);
      HFREE2(g_NAFVariable.ZALP);
      SYSFREE(g_NAFVariable.ZTABS);
      SYSFREE(TWIN);
      /* v0.96 M. GASTINEAU 06/01/99 : ajout */
      delete_list_fenetre_naf(g_NAFVariable.m_pListFen);
      g_NAFVariable.m_pListFen =NULL;
      /* v0.96 M. GASTINEAU 06/01/99 : fin ajout */
}/*      end  subroutine CLEANNAF  */
      
/*-----------------------------------------------------------------------*/
/* fonction identique a naf_initnaf mais :                               */
/*n'alloue pas de tableau pour ZTABS, TFS et ZAMP.                       */
/*-----------------------------------------------------------------------*/
/*v0.96 M. GASTINEAU 06/10/98 : ajout */
void naf_initnaf_notab()
{  
/*!----------------- PREMIERES INITIALISATIONS*/
      g_NAFVariable.EPSM = DBL_EPSILON;
      g_NAFVariable.UNIANG = g_NAFVariable.DTOUR/(2*M_PI) ;
      g_NAFVariable.FREFON = g_NAFVariable.DTOUR/(g_NAFVariable.KTABS*g_NAFVariable.XH)	;      
      DIM2(g_NAFVariable.ZALP, (g_NAFVariable.NTERM+1), (g_NAFVariable.NTERM+1), t_complexe,"ZALP"); /*allocate(ZALP(1:NTERM,1:NTERM),stat = NERROR)*/
      SYSCHECKMALLOCSIZE(TWIN, double, g_NAFVariable.KTABS+1); /*allocate(TWIN(0:KTABS),stat = NERROR)*/
      /*v0.96 M. GASTINEAU 18/12/98 : modification du prototype */
      /*naf_iniwin();*/ naf_iniwin(TWIN);    
}/*      end SUBROUTINE naf_initnaf_notab*/

/*-----------------------------------------------------------------------*/
/* fonction identique a naf_initnaf mais :                               */
/*ne libere pas  ZTABS, TFS et ZAMP.                                     */
/*-----------------------------------------------------------------------*/
/*v0.96 M. GASTINEAU 06/10/98 : ajout */
/* v0.96 M. GASTINEAU 06/01/99 : ajout de la liberation de la liste des fenetres */
void naf_cleannaf_notab()
{
/*!-----------------------------------------------------------------------
!     desalloue les tableaux dynamiques
!-----------------------------------------------------------------------*/
 
      HFREE2(g_NAFVariable.ZALP);
      SYSFREE(TWIN);
      /* v0.96 M. GASTINEAU 06/01/99 : ajout */
      delete_list_fenetre_naf(g_NAFVariable.m_pListFen);
      g_NAFVariable.m_pListFen =NULL;
      /* v0.96 M. GASTINEAU 06/01/99 : fin ajout */
}/*      end  subroutine naf_initnaf_notab  */

         
/*!-----------------------------------------------------------------------
      SUBROUTINE MFTNAF(NBTERM,EPS)*/
/*v0.97 M. GASTINEAU 26/05/99 : correction bug 0.97/99/05/26/A  pour pas<0 (XH<0) */
BOOL naf_mftnaf(int NBTERM, double EPS)
{
/*!-----------------------------------------------------------------------
!     MFTNAF
!                CALCULE UNE APPROXIMATION QUASI PERIODIQUE
!                DE LA FONCTION TABULEE DANS g_NAFVariable.ZTABS(0:KTABS)
!
!     NBTERM               : NOMBRE DE TERMES RECHERCHES (<= NTERM)
!
!
!     EPS                 : PRECISION AVEC LAQUELLE ON RECHERCHE 
!                           LES FREQUENCES
!
!-----------------------------------------------------------------------

      IMPLICIT NONE
! (EPS)
      integer :: NBTERM
      REAL (8) :: EPS     
!
      integer :: I,IFLAG,NUMFR
      REAL (8) :: TOL,STAREP,FR,A,B,RM     
!-----------------------------------------------------------------------*/
      int I, IFLAG, NUMFR;
      double TOL,STAREP,FR,A,B,RM;
/*v0.96 M. GASTINEAU 14/01/99 : ajout - support des fenetres*/
#if NAF_USE_OPTIMIZE>0
      int iFrMin;
      int iFrMax;
      double FREFO2;
      int KTABS2;
#endif /*NAF_USE_OPTIMIZE*/
/*v0.96 M. GASTINEAU 14/01/99 : fin ajout*/
      
      if (NBTERM >g_NAFVariable.NTERM)
      {
       Myyerror("Nbre de termes cherches trop grand");
       return FALSE;
      }
      TOL = 1.E-4;
      STAREP=fabs(g_NAFVariable.FREFON)/3;
      naf_inifre();
/*v0.96 M. GASTINEAU 14/01/99 : ajout support des fenetres*/
#if NAF_USE_OPTIMIZE>0
     naf_puiss2(g_NAFVariable.KTABS+1,&KTABS2);
     FREFO2=(g_NAFVariable.FREFON*g_NAFVariable.KTABS)/KTABS2;
     do
     {
      if (g_NAFVariable.m_pListFen==NULL)
      {/*pas de fenetre => on recherche dans toutes les frequences */
       iFrMin=KTABS2;
       iFrMax=KTABS2;
      }
      else
      {
       const int iMaxValue=KTABS2/2-1;
       NBTERM=g_NAFVariable.m_pListFen->iNbTerme;
       iFrMin=g_NAFVariable.m_pListFen->dFreqMin/FREFO2;
       iFrMax=g_NAFVariable.m_pListFen->dFreqMax/FREFO2;
       /*v0.97 M. GASTINEAU 26/05/99 : ajout - correction bug (cas pas <0) */
       if (iFrMin>iFrMax) 
       {/*swap*/
        double temp = iFrMin;
        iFrMin = iFrMax;
        iFrMax = temp;
       }
       /*v0.97 M. GASTINEAU 26/05/99 : fin ajout  */
       if (iFrMin<-iMaxValue) 
       { iFrMin=-iMaxValue; }
       else if (iFrMin>iMaxValue) 
       { iFrMin=+iMaxValue; }
       if (iFrMax<-iMaxValue) 
       { iFrMax=-iMaxValue; }
       else if (iFrMax>iMaxValue) 
       { iFrMax=+iMaxValue; }
      }
#endif /*NAF_USE_OPTIMIZE*/
/*v0.96 M. GASTINEAU 14/01/99 : fin ajout */
      for(I=1;I<=NBTERM; I++)
      {
/*v0.96 M. GASTINEAU 14/01/99 : modification support des fenetres*/
#if NAF_USE_OPTIMIZE==0
         naf_fftmax(&FR);
#else /*remplacee par:*/
         FR=naf_fftmax(iFrMin,iFrMax,FREFO2,KTABS2);
#endif /*NAF_USE_OPTIMIZE*/
/*v0.96 M. GASTINEAU 14/01/99 : fin modification*/
         naf_frefin(&FR,&A,&B,&RM,STAREP,EPS);
         naf_fretes(FR,&IFLAG,TOL,&NUMFR);
         if (IFLAG == 0) break; /*GOTO 999*/
         if (IFLAG == 1)
         {
            if(naf_gramsc(FR,A,B)==FALSE)
            {
             return FALSE;
            }
            g_NAFVariable.NFS++; /*g_NAFVariable.NFS=g_NAFVariable.NFS+1;*/
         }
         if (IFLAG==-1)
         {
          naf_modfre(NUMFR,&A,&B);
         }
      }
/*v0.96 M. GASTINEAU 14/01/99 : ajout - support des fenetres*/
#if NAF_USE_OPTIMIZE>0
      /*passage a la fenetre suivante */
      if (g_NAFVariable.m_pListFen!=NULL)
      {
       t_list_fenetre_naf *pListTemp=g_NAFVariable.m_pListFen;
       g_NAFVariable.m_pListFen=pListTemp->suivant;
       SYSFREE(pListTemp);
      }
     } while (g_NAFVariable.m_pListFen!=NULL);
#endif /*NAF_USE_OPTIMIZE*/
/*v0.96 M. GASTINEAU 14/01/99 : fin ajout*/

 /*999   CONTINUE*/
 return TRUE;
} /*      END SUBROUTINE MFTNAF */
/*!
      SUBROUTINE PRTABS(KTABS,g_NAFVariable.ZTABS,IPAS)*/
void naf_prtabs(int KTABS, t_complexe *ZTABS, int IPAS)
{
/*
!-----------------------------------------------------------------------
!     IMPRESSION DE g_NAFVariable.ZTABS
!-----------------------------------------------------------------------
  
      IMPLICIT NONE
! (KTABS,g_NAFVariable.ZTABS,IPAS)
      integer :: KTABS,IPAS
      complex (8) :: g_NAFVariable.ZTABS(0:KTABS)
!
      integer :: I         
!*/
      int I;
      
      if  (g_NAFVariable.IPRT==1)
      {
         for (I=0;I<=KTABS;I+=IPAS)
         {
          fprintf(g_NAFVariable.NFPRT,"%6d  %-+20.15E %-+20.15E\n",I, ZTABS[I].reel, ZTABS[I].imag);
          /* WRITE(g_NAFVariable.NFPRT,1000) I,DREAL(ZTABS(I)),DIMAG(ZTABS(I))*/
         }
      }
/*1000  FORMAT (1X,I6,2X,2D20.6)*/
}/*      END SUBROUTINE PRTABS */
/*!
!
      SUBROUTINE SMOY(ZM)*/
void naf_smoy(t_complexe *ZM)
{
/*
!-----------------------------------------------------------------------
!    CALCUL DES MOYENNNES DE g_NAFVariable.ZTABS ET SOUSTRAIT DE CHACUNE
!    LA VALEUR MOYENNE
!-----------------------------------------------------------------------
 
      IMPLICIT NONE 
! (ZM)
      complex (8) :: ZM 
!
      integer :: I          
!-----------------------------------------------------------------------
*/
      int I;
/* v0.96 M. GASTINEAU 01/12/98 : utilisation des fonctions  complexes inlines et optimisation */
#if NAF_USE_OPTIMIZE==0
       *ZM=cmplx(0.E0,0.E0);
      for (I=0; I<=g_NAFVariable.KTABS; I++)
      {
         *ZM = addcomplexe(*ZM,g_NAFVariable.ZTABS[I]);
      }
      *ZM=muldoublcomplexe(1E0/((double)(g_NAFVariable.KTABS+1)), *ZM); /*ZM=ZM/(g_NAFVariable.KTABS+1)*/
      for (I=0; I<=g_NAFVariable.KTABS; I++)
      {
         g_NAFVariable.ZTABS[I]=subcomplexe(g_NAFVariable.ZTABS[I],*ZM);
      }
#else /*remplacees par: */
      t_complexe *pzarTabs;
      double dNbKTabs1=g_NAFVariable.KTABS+1;
      i_compl_cmplx(ZM,0.E0,0.E0);
      for (I=0, pzarTabs = g_NAFVariable.ZTABS; 
           I<=g_NAFVariable.KTABS;
           I++, pzarTabs++)
      {
         i_compl_padd(ZM,pzarTabs);
      }
      i_compl_pdivdoubl(ZM,&dNbKTabs1);/*ZM=ZM/(g_NAFVariable.KTABS+1)*/
      for (I=0, pzarTabs = g_NAFVariable.ZTABS;
           I<=g_NAFVariable.KTABS; 
           I++, pzarTabs++)
      {
         i_compl_psub(pzarTabs,ZM);
      }
#endif /*NAF_USE_OPTIMIZE*/
/* v0.96 M. GASTINEAU 01/12/98 : fin modification */
}/*      END SUBROUTINE SMOY*/
/*!
!
      SUBROUTINE FRETES (FR,IFLAG,TOL,NUMFR)*/
void naf_fretes(double FR, int *IFLAG, double TOL, int * NUMFR)
{
/*!**********************************************************************
!     TEST DE LA NOUVELLE FREQUENCE TROUVEE PAR RAPPORT AUX ANCIENNES.
!     LA DISTANCE ENTRE DEUX FREQ  DOIT ETRE DE g_NAFVariable.FREFON
!
!     RENVOIE IFLAG =  1 SI LE TEST REUSSIT (ON PEUT CONTINUER)
!             IFLAG =  0 SI LE TEST ECHOUE (IL VAUT MIEUX S'ARRETER)
!             IFLAG = -1 SI TEST < ECART, MAIS TEST/ECART < TOL
!                        (ON RETROUVE PRATIQUEMENT LA MEME FREQUENCE 
!                         D'INDICE NFR)
!     TOL (ENTREE) TOLERANCE ( 1.D-7) EST UN BON EXEMPLE
!     NFR (SORTIE) INDICE DE LA FREQUENCE RETROUVEE
!
!***********************************************************************
  
      IMPLICIT NONE
! (FR,IFLAG,TOL,NUMFR)
      integer :: IFLAG,NUMFR
      REAL (8) :: FR,TOL
!
      integer :: I
      REAL (8) :: ECART,TEST
! */
      int I;
      double ECART, TEST;                
      *IFLAG = 1;
      ECART = fabs(g_NAFVariable.FREFON) ; 
      for (I = 1; I<=g_NAFVariable.NFS; I++)
      {
        TEST = fabs(g_NAFVariable.TFS[I] - FR);
        if (TEST<ECART)
        {
            if ((TEST/ECART)<TOL)
            {
               *IFLAG = -1;
               *NUMFR   = I;
               if (g_NAFVariable.IPRT>=1)
               {
                fprintf(g_NAFVariable.NFPRT, "TEST/ECART = %g   ON CONTINUE\n", TEST/ECART);
               }
               break; /*GOTO 999*/
            }
            else
            {
               *IFLAG = 0 ;
               if (g_NAFVariable.IPRT>=0)
               {
                fprintf(g_NAFVariable.NFPRT,"TEST = %g ECART = %g \n", TEST, ECART);
                fprintf(g_NAFVariable.NFPRT,"FREQUENCE   FR = %g  TROP PROCHE  DE  %g\n", FR, g_NAFVariable.TFS[I]);
               }
            }
         }
      }
/*999   CONTINUE*/
}/*      END SUBROUTINE FRETES*/
/*!
      SUBROUTINE ZTPOW (N,N1,ZT,ZA,ZAST)*/
void naf_ztpow(int N, int N1, t_complexe *ZT, t_complexe ZA, t_complexe ZAST)
{
/*!-----------------------------------------------------------------------
!     ZTPOW   CALCULE  ZT(I) = ZAST*ZA**I EN VECTORIEL
!             ZT(N)
!-----------------------------------------------------------------------
      IMPLICIT NONE
! (N,N1,ZT,ZA,ZAST)
      integer :: N,N1
      complex (8) :: ZT(0:N),ZA,ZAST
!           
      integer :: I,INC,NX,IT,NT
      complex (8) ZT1,ZINC
! */
      int  I,INC,NX,IT,NT;
      t_complexe  ZT1,ZINC;   
/* v0.96 M. GASTINEAU 01/12/98 : utilisation des fonctions  complexes inlines et optimisation */
#if NAF_USE_OPTIMIZE>0
      t_complexe *pzarZT;
#endif     /*NAF_USE_OPTIMIZE*/ 

      if (N<N1-1)
      {
         printf("DANS ZTPOW, N = %d\n", N);
         return;
      }
/*!----------! */
/* v0.96 M. GASTINEAU 01/12/98 : utilisation des fonctions  complexes inlines et optimisation */
#if NAF_USE_OPTIMIZE==0
      ZT[0] = mulcomplexe(ZAST,ZA);
      for(I = 1; I<N1; I++)
      {
         ZT[I] = mulcomplexe(ZT[I-1],ZA);
      }
      ZT1=divcomplexe(ZT[N1-1], ZAST); /*ZT1 = ZT(N1-1)/ZAST*/
      ZINC=cmplx(1E0, 0E0); /*ZINC= 1;*/
      INC =0;
      NT = (N+1)/N1;  
      for (IT = 2; IT<= NT; IT++)
      {
         ZINC = mulcomplexe(ZINC,ZT1);
         INC += N1; /*INC  = INC + N1*/
         for (I = 0; I<N1; I++)
         {
            ZT[INC+I]=mulcomplexe(ZT[I], ZINC); /*ZT(INC +I) = ZT(I)*ZINC*/
         }
      }
      ZINC = mulcomplexe(ZINC, ZT1); /*ZINC = ZINC*ZT1*/
      INC += N1; /*INC  = INC + N1;*/
      NX = N+1-NT*N1;
      for (I = 0; I<NX; I++)
      {
         ZT[INC+I]=mulcomplexe(ZT[I], ZINC); /*ZT(INC +I) = ZT(I)*ZINC*/
      }
#else /*remplacee par: */
      ZT[0] = i_compl_mul(ZAST,ZA);
      for(I = 1, pzarZT=ZT+1 ; I<N1; I++,pzarZT++)
      {
         *pzarZT = i_compl_mul(*(pzarZT-1),ZA);
      }
      ZT1=i_compl_div(ZT[N1-1], ZAST); /*ZT1 = ZT(N1-1)/ZAST*/
      i_compl_cmplx(&ZINC,1E0, 0E0); /*ZINC= 1;*/
      INC =0;
      NT = (N+1)/N1;  
      for (IT = 2; IT<= NT; IT++)
      {
         i_compl_pmul(&ZINC,&ZT1);
         INC += N1; /*INC  = INC + N1*/
         for (I = 0; I<N1; I++)
         {
            ZT[INC+I]=i_compl_mul(ZT[I], ZINC); /*ZT(INC +I) = ZT(I)*ZINC*/
         }
      }
      i_compl_pmul(&ZINC, &ZT1); /*ZINC = ZINC*ZT1*/
      INC += N1; /*INC  = INC + N1;*/
      NX = N+1-NT*N1;
      for (I = 0; I<NX; I++)
      {
         ZT[INC+I]=i_compl_mul(ZT[I], ZINC); /*ZT(INC +I) = ZT(I)*ZINC*/
      }
#endif /*NAF_USE_OPTIMIZE*/ 
/* v0.96 M. GASTINEAU 01/12/98 : fin modification */
}/*      END SUBROUTINE ZTPOW*/
/*!
      SUBROUTINE TESSOL (EPS,TFSR,ZAMPR)*/
BOOL naf_tessol(double EPS, double *TFSR, t_complexe *ZAMPR)
{      
/*!-----------------------------------------------------------------------
!     TESSOL
!                SOUS PROGRAMME POUR VERIFIER LA PRECISION DES
!                SOLUTIONS OBTENUES PAR ANALYSE DE FOURIER
!                ON ANALYSE A NOUVEAU LA SOLUTION ET ON COMPARE 
!                LES TERMES
!
!     J. LASKAR 25/5/89
!-----------------------------------------------------------------------
      
      IMPLICIT NONE
! (EPS,TFSR,ZAMPR)
      REAL (8) :: EPS, TFSR(g_NAFVariable.NTERM)
      complex (8) :: ZAMPR(g_NAFVariable.NTERM)
!          
      integer :: IT,IFR,JFR, NVTERM
      REAL (8) :: OFFSET
      complex (8) :: ZI,ZA,ZOM,ZEX
      REAL (8), dimension(:),allocatable :: TFST
      complex (8),allocatable,dimension(:) :: ZAMPT
      complex (8),allocatable,dimension(:,:) :: ZALPT
      complex (8) :: ZINC
      complex (8), dimension (:), allocatable :: ZT
 */
      int IT,IFR,JFR, NVTERM;
      double OFFSET;
      t_complexe ZI,ZA,ZOM,ZEX;
      double *TFST=NULL;
      t_complexe *ZAMPT=NULL;
      t_complexe **ZALPT=NULL;
      t_complexe ZINC;
      t_complexe *ZT=NULL;
/* v0.96 M. GASTINEAU 01/12/98 : utilisation des fonctions  complexes inlines et optimisation */
#if NAF_USE_OPTIMIZE>0
      t_complexe *pzarTab;
#endif /**/      
      SYSCHECKMALLOCSIZE(ZAMPT, t_complexe, g_NAFVariable.NTERM+1); /* allocate(ZAMPT(1:g_NAFVariable.NTERM),stat = NERROR)*/
      DIM2(ZALPT, (g_NAFVariable.NTERM+1), (g_NAFVariable.NTERM+1), t_complexe, "ZALPT"); /*allocate(ZALPT(1:g_NAFVariable.NTERM,1:g_NAFVariable.NTERM),stat = NERROR)*/
      SYSCHECKMALLOCSIZE(TFST, double, g_NAFVariable.NTERM+1);/* allocate(TFST(1:g_NAFVariable.NTERM),stat = NERROR)*/
      SYSCHECKMALLOCSIZE(ZT, t_complexe, g_NAFVariable.KTABS+1);/* allocate (ZT (0:g_NAFVariable.KTABS),stat = NERROR)*/

/*!*************************************************************/
      OFFSET=g_NAFVariable.KTABS*g_NAFVariable.XH;
/* v0.96 M. GASTINEAU 01/12/98 : utilisation des fonctions  complexes inlines et optimisation */
#if NAF_USE_OPTIMIZE==0
/*!----------! INITIALISATION DE g_NAFVariable.ZTABS*/
      for ( IT=0;IT<=g_NAFVariable.KTABS; IT++)
       {
          g_NAFVariable.ZTABS[IT]=cmplx(0.E0,0.E0);
       }
/*!----------! CALCUL DE LA NOUVELLE SOLUTION*/
      ZI=cmplx(0.E0,1.E0);
      for ( IFR = 1; IFR<=g_NAFVariable.NFS; IFR++)
      {
         ZOM = muldoublcomplexe(g_NAFVariable.TFS[IFR]/g_NAFVariable.UNIANG,ZI); /* ZOM=g_NAFVariable.TFS(IFR)/g_NAFVariable.UNIANG*ZI */
         ZA=g_NAFVariable.ZAMP[IFR];
/*!-----------! ATTENTION ICI (CAS REEL) ON AURA AUSSI LE TERME CONJUGUE
!-----------! QU'ON NE CALCULE PAS. LE TERME TOTAL EST
!-----------!       2*RE(g_NAFVariable.ZAMP(I)*EXP(ZI*g_NAFVariable.TFS(I)*T) )
!-----------! ON RAJOUTE LA CONTRIBUTION TOTALE DU TERME g_NAFVariable.TFS(I) DANS TAB2*/
         if (g_NAFVariable.ICPLX==1)
         {
            ZEX = mulcomplexe(ZA,expcomplexe(muldoublcomplexe(g_NAFVariable.T0+OFFSET-g_NAFVariable.XH,ZOM))); /*ZEX = ZA*EXP ((g_NAFVariable.T0+OFFSET-g_NAFVariable.XH)*ZOM)*/
            ZINC= expcomplexe(muldoublcomplexe(g_NAFVariable.XH,ZOM)); /*ZINC= EXP (g_NAFVariable.XH*ZOM)*/
            naf_ztpow(g_NAFVariable.KTABS,64,ZT,ZINC,ZEX); /*naf_ztpow(g_NAFVariable.KTABS+1,64,ZT,ZINC,ZEX);*/
            for (IT=0; IT<=g_NAFVariable.KTABS; IT++)
            {
               g_NAFVariable.ZTABS[IT]=addcomplexe(g_NAFVariable.ZTABS[IT], ZT[IT]); /*g_NAFVariable.ZTABS(IT)=g_NAFVariable.ZTABS(IT)+ZT(IT)*/
            }
         }
         else
         {
            ZEX = mulcomplexe(ZA,expcomplexe(muldoublcomplexe(g_NAFVariable.T0+OFFSET-g_NAFVariable.XH,ZOM))); /*ZEX = ZA*EXP ((g_NAFVariable.T0+OFFSET-g_NAFVariable.XH)*ZOM)*/
            ZINC= expcomplexe(muldoublcomplexe(g_NAFVariable.XH,ZOM)); /*ZINC= EXP (g_NAFVariable.XH*ZOM)*/
            naf_ztpow(g_NAFVariable.KTABS,64,ZT,ZINC,ZEX); /*naf_ztpow(g_NAFVariable.KTABS+1,64,ZT,ZINC,ZEX);*/
            for (IT=0; IT<=g_NAFVariable.KTABS; IT++)
            {
               g_NAFVariable.ZTABS[IT].reel += ZT[IT].reel; /*g_NAFVariable.ZTABS(IT)=g_NAFVariable.ZTABS(IT)+DREAL(ZT(IT))*/
            }
         }
      }
#else /*remplacee par:*/
/*!----------! INITIALISATION DE g_NAFVariable.ZTABS*/
      for ( IT=0, pzarTab=g_NAFVariable.ZTABS;
            IT<=g_NAFVariable.KTABS;
            IT++, pzarTab++)
       {
          i_compl_cmplx(pzarTab, 0.E0, 0.E0);
       }
/*!----------! CALCUL DE LA NOUVELLE SOLUTION*/
      i_compl_cmplx(&ZI,0.E0,1.E0);
      for ( IFR = 1; IFR<=g_NAFVariable.NFS; IFR++)
      {
         ZOM = i_compl_muldoubl(g_NAFVariable.TFS[IFR]/g_NAFVariable.UNIANG,ZI); /* ZOM=g_NAFVariable.TFS(IFR)/g_NAFVariable.UNIANG*ZI */
         ZA=g_NAFVariable.ZAMP[IFR];
/*!-----------! ATTENTION ICI (CAS REEL) ON AURA AUSSI LE TERME CONJUGUE
!-----------! QU'ON NE CALCULE PAS. LE TERME TOTAL EST
!-----------!       2*RE(g_NAFVariable.ZAMP(I)*EXP(ZI*g_NAFVariable.TFS(I)*T) )
!-----------! ON RAJOUTE LA CONTRIBUTION TOTALE DU TERME g_NAFVariable.TFS(I) DANS TAB2*/
         if (g_NAFVariable.ICPLX==1)
         {
            ZEX = i_compl_mul(ZA,i_compl_exp(i_compl_muldoubl(g_NAFVariable.T0+OFFSET-g_NAFVariable.XH,ZOM))); /*ZEX = ZA*EXP ((g_NAFVariable.T0+OFFSET-g_NAFVariable.XH)*ZOM)*/
            ZINC= i_compl_exp(i_compl_muldoubl(g_NAFVariable.XH,ZOM)); /*ZINC= EXP (g_NAFVariable.XH*ZOM)*/
            naf_ztpow(g_NAFVariable.KTABS,64,ZT,ZINC,ZEX); /*naf_ztpow(g_NAFVariable.KTABS+1,64,ZT,ZINC,ZEX);*/
            for (IT=0, pzarTab=g_NAFVariable.ZTABS;
                 IT<=g_NAFVariable.KTABS;
                 IT++, pzarTab++)
            {
               i_compl_padd(pzarTab,ZT+IT); /*g_NAFVariable.ZTABS(IT)=g_NAFVariable.ZTABS(IT)+ZT(IT)*/
            }
         }
         else
         {
            ZEX = i_compl_mul(ZA,i_compl_exp(i_compl_muldoubl(g_NAFVariable.T0+OFFSET-g_NAFVariable.XH,ZOM))); /*ZEX = ZA*EXP ((g_NAFVariable.T0+OFFSET-g_NAFVariable.XH)*ZOM)*/
            ZINC= i_compl_exp(i_compl_muldoubl(g_NAFVariable.XH,ZOM)); /*ZINC= EXP (g_NAFVariable.XH*ZOM)*/
            naf_ztpow(g_NAFVariable.KTABS,64,ZT,ZINC,ZEX); /*naf_ztpow(g_NAFVariable.KTABS+1,64,ZT,ZINC,ZEX);*/
            for (IT=0; IT<=g_NAFVariable.KTABS; IT++)
            {
               g_NAFVariable.ZTABS[IT].reel += ZT[IT].reel; /*g_NAFVariable.ZTABS(IT)=g_NAFVariable.ZTABS(IT)+DREAL(ZT(IT))*/
            }
         }
      }
#endif /*NAF_USE_OPTIMIZE==0*/
/* v0.96 M. GASTINEAU 01/12/98 : utilisation des fonctions  complexes inlines et optimisation */
/*!----------------------*/
      SYSFREE(ZT); /*deallocate(ZT)*/
/*!-----------! SAUVEGARDE DANS DES TABLEAUX ANNEXES DE LA SOLUTION*/
       for( IFR = 1; IFR<=g_NAFVariable.NFS; IFR++)
       {
         TFST[IFR] = g_NAFVariable.TFS[IFR];
         ZAMPT[IFR] = g_NAFVariable.ZAMP[IFR];
          for(JFR = 1;JFR <=g_NAFVariable.NFS; JFR++)
          {
            ZALPT[IFR][JFR] = g_NAFVariable.ZALP[IFR][JFR];
          }
      }
/*!-----------! NOUVEAU CALCUL DE LA SOLUTION*/
      NVTERM=g_NAFVariable.NFS;
/*!***DEBUG      CALL PRTABS(g_NAFVariable.KTABS,g_NAFVariable.ZTABS,g_NAFVariable.KTABS/10)*/
      if (naf_mftnaf (NVTERM,EPS)==FALSE)
      {
       HFREE2(ZALPT);
       SYSFREE(ZAMPT);
       SYSFREE(TFST);
       return FALSE;
      }
/*!-----------! RESULTATS*/
      for( IFR = 1; IFR<=g_NAFVariable.NFS; IFR++)
      {
         TFSR[IFR] = g_NAFVariable.TFS[IFR];
         ZAMPR[IFR] = g_NAFVariable.ZAMP[IFR];
      }
/*!-----------! RESTITUTION  DE LA SOLUTION*/
      for( IFR = 1; IFR<=g_NAFVariable.NFS; IFR++)
      {
         g_NAFVariable.TFS[IFR] = TFST[IFR];
         g_NAFVariable.ZAMP[IFR] = ZAMPT[IFR];
         for(JFR = 1;JFR <=g_NAFVariable.NFS; JFR++)
         {
             g_NAFVariable.ZALP[IFR][JFR] = ZALPT[IFR][JFR];
         }
      }
      HFREE2(ZALPT);
      SYSFREE(ZAMPT);
      SYSFREE(TFST);
      return TRUE;
}/*      END SUBROUTINE TESSOL*/
/*!
!
      SUBROUTINE INIFRE*/
void  naf_inifre()
{
/*
!************************************************************************
!     REMET A ZERO g_NAFVariable.TFS, g_NAFVariable.ZAMP,g_NAFVariable.ZALP et g_NAFVariable.NFS
!     UTILE QUAND ON BOUCLE SUR PLUSIEURS CAS
!
!***********************************************************************
       
      IMPLICIT NONE
*/
      int I,J;
/* v0.96 M. GASTINEAU 12/01/99 : optimisation */
#if NAF_USE_OPTIMIZE==0
      t_complexe cZERO;
      cZERO=cmplx(0.E0,0.E0);
      for(I = 1;I<= g_NAFVariable.NTERM; I++)
      {
         g_NAFVariable.TFS[I] = 0.E0;
         g_NAFVariable.ZAMP[I] = cZERO;
         for(J = 1; J<= g_NAFVariable.NTERM; J++)
         {
            g_NAFVariable.ZALP[I][J] = cZERO;
         }
      }
#else /* remplace par: */
      const int iNterm=g_NAFVariable.NTERM;
      const double dZero=0.E0;
      double *pdTFS=g_NAFVariable.TFS+1;
      double *pdZAMP=(double*)(g_NAFVariable.ZAMP+1);
      double *pdZALP;
      for(I = 1;I<= iNterm; I++)
      {
         pdZALP=(double*)(g_NAFVariable.ZALP[I]);
         *pdTFS++ = dZero;
         *pdZAMP++ = dZero;
         *pdZAMP++ = dZero;
         for(J = 1; J<= iNterm; J++)
         {
            *pdZALP++ = dZero;
            *pdZALP++ = dZero;
         }
      }

#endif /*NAF_USE_OPTIMIZE*/
/* v0.96 M. GASTINEAU 12/01/99 : fin optimisation */
      g_NAFVariable.NFS = 0;
} /*  END SUBROUTINE INIFRE*/
/*!
!
      SUBROUTINE FFTMAX(FR)*/
/*
!-----------------------------------------------------------------------
!     FFTMAX          MODIF NICE 90 (J. LASKAR)
!                     MODIF 13/12/90
!                     VERSION STANDARD. (FOUR1)
!                     ANALYSE DE FOURIER RAPIDE
!                     RECHERCHE D'UN MAXIMUM DE SPECTRE POUR UNE FREQ. F
!***********************************************************************
!     ON UTILISERA CELA ULTERIEUREMENT (13/12/90)
!     NFR11, NFR22 NUMERO MIN ET MAX DES FREQUENCES RECHERCHEES
!                   NFR1 ET NFR2 INCLUSES
!                -NRTAB<NFR1<NFR2<NRTAB
!
!     RTAB(NRTAB,2) EN 1 AMPLITUDE DE SFREQUENCES POSITIVES
!                      2 AMPLITUDE DES FREQUENCES NEGATIVES
!     DISTANCE ENTRE DEUX LIGNES  g_NAFVariable.FREFON
!-----------------------------------------------------------------------
      
      IMPLICIT NONE
! (FR)
      REAL (8) :: FR      
!
      integer :: KTABS2,ISG,IPAS,I,IDV,INDX,IFR
      REAL (8) :: FREFO2
      REAL (8), dimension(:),allocatable :: TAB
      REAL (8), dimension(:),allocatable :: RTAB
! */
/* v0.96 M. GASTINEAU 14/01/99 : support des fenetres */
#if NAF_USE_OPTIMIZE==0
void naf_fftmax(double *FR)
{
      int   KTABS2,ISG,IPAS,I,IDV,INDX,IFR;
      double  FREFO2;
      double *pdTAB=NULL; /*=TAB*/
      double *RTAB=NULL;

      naf_puiss2(g_NAFVariable.KTABS+1,&KTABS2);
/*! */     
      SYSCHECKMALLOCSIZE(pdTAB,double,2*KTABS2); pdTAB--;/*  allocate(TAB(2*KTABS2),stat = NERROR)*/
      SYSCHECKMALLOCSIZE(RTAB,double,KTABS2);/*allocate(RTAB(0:KTABS2-1),stat = NERROR)*/
/*!*/      
      FREFO2=(g_NAFVariable.FREFON*g_NAFVariable.KTABS)/KTABS2;
/*!****************** */  
      if (g_NAFVariable.IPRT==1)
      {
       fprintf(g_NAFVariable.NFPRT,"KTABS2= %d  FREFO2= %g\n",KTABS2, FREFO2);
      }
/*!****************! CALCUL DES FREQUENCES */
      ISG=-1;
      IDV=KTABS2;
      IPAS=1;
      for(I=0;I<=KTABS2-1; I++)
      {
         pdTAB[2*I+1]=g_NAFVariable.ZTABS[I].reel*TWIN[I];/*pdTAB(2*I+1)=DREAL(g_NAFVariable.ZTABS(I))*TWIN(I)*/
         pdTAB[2*I+2]=g_NAFVariable.ZTABS[I].imag*TWIN[I]; /*pdTAB(2*I+2)=DIMAG(g_NAFVariable.ZTABS(I))*TWIN(I)*/          
      }
      naf_four1(pdTAB,KTABS2,ISG);
      for(I=0;I<=KTABS2-1; I++)
      {
         RTAB[I]=sqrt(pdTAB[2*I+1]*pdTAB[2*I+1]+pdTAB[2*I+2]*pdTAB[2*I+2])/IDV;
      }
      pdTAB++; SYSFREE(pdTAB);
/*!**********************        CALL MODTAB(KTABS2,RTAB)*/
      naf_maxx(KTABS2-1, RTAB, &INDX);/*naf_maxx(KTABS2,RTAB,&INDX);*/
      /* car naf_maxx travaille sur des tableaux de 1 a N */
      
      IFR=((INDX+1)<=(KTABS2/2))?INDX:INDX-KTABS2; /*IF (INDX.LE.KTABS2/2) THEN
          IFR = INDX-1
      ELSE
          IFR = INDX-1-KTABS2
      ENDIF*/
      
      *FR = IFR*FREFO2*IPAS;
      if (g_NAFVariable.IPRT==1)
      {
       fprintf(g_NAFVariable.NFPRT,"IFR=%d FR=%g RTAB=%g INDX=%d KTABS2=%d\n",IFR,*FR, RTAB[INDX],INDX,KTABS2); 
      }
      SYSFREE(RTAB);
}
#else /*remplacee par:*/
double naf_fftmax(int p_iFrMin, int p_iFrMax, double FREFO2, int KTABS2) 
{
/* la fonction retourne la frequence FR dertminee */
/* On suppose p_iFrMin < p_iFrMax */
      double FR;
      int   ISG,IPAS,I,INDX,IFR;
      /*double  FREFO2;*/
      double *pdTAB=NULL; /*=TAB*/
      double *RTAB=NULL;
      int iKTABS2m1, iKTABS2;
      double dDIV;
      double *pdTABTemp1,*pdTABTemp2;
      
      
      /*naf_puiss2(g_NAFVariable.KTABS+1,&KTABS2);*//*v0.96 M. GASTINEAU 14/01/99 */
      iKTABS2 = KTABS2;
      iKTABS2m1 = iKTABS2-1;
/*! */     
      SYSCHECKMALLOCSIZE(pdTAB,double,2*iKTABS2); /*  allocate(TAB(2*KTABS2),stat = NERROR)*/
      SYSCHECKMALLOCSIZE(RTAB,double,iKTABS2);/*allocate(RTAB(0:KTABS2-1),stat = NERROR)*/
/*!*/      
     /* FREFO2=(g_NAFVariable.FREFON*g_NAFVariable.KTABS)/iKTABS2;*//*v0.96 M. GASTINEAU 14/01/99 */
/*!****************** */  
      if (g_NAFVariable.IPRT==1)
      {
       fprintf(g_NAFVariable.NFPRT,"KTABS2= %d  FREFO2= %g\n",iKTABS2, FREFO2);
      }
/*!****************! CALCUL DES FREQUENCES */
      ISG=-1;
      dDIV=iKTABS2;
      IPAS=1;
      for(I=0, pdTABTemp1=pdTAB, pdTABTemp2=(double*)(g_NAFVariable.ZTABS);
          I<=iKTABS2m1;
          I++)
      {
         *pdTABTemp1++ = (*pdTABTemp2++) * TWIN[I];/*pdTAB(2*I+1)=DREAL(g_NAFVariable.ZTABS(I))*TWIN(I)*/
         *pdTABTemp1++ = (*pdTABTemp2++) * TWIN[I]; /*pdTAB(2*I+2)=DIMAG(g_NAFVariable.ZTABS(I))*TWIN(I)*/          
      }
      naf_four1(pdTAB-1,iKTABS2,ISG);
      for(I=0, pdTABTemp1=pdTAB, pdTABTemp2=pdTABTemp1+1;
          I<=iKTABS2m1;
          I++, pdTABTemp1+=2, pdTABTemp2+=2)
      {
         RTAB[I]=sqrt((*pdTABTemp1)*(*pdTABTemp1)+(*pdTABTemp2)*(*pdTABTemp2))/dDIV;
      }
      SYSFREE(pdTAB);
/*!**********************        CALL MODTAB(KTABS2,RTAB)*/
      /*v0.96 M. GASTINEAU 14/01/99 : modification pour le support des fenetres */
      /*INDX=naf_maxx(iKTABS2m1, RTAB);*//*naf_maxx(KTABS2,RTAB,&INDX);*/
      /*IFR=((INDX+1)<=(iKTABS2/2))?INDX:INDX-iKTABS2;*/ /*IF (INDX.LE.KTABS2/2) THEN
          IFR = INDX-1
      ELSE
          IFR = INDX-1-KTABS2
      ENDIF*/
      /* *FR = IFR*FREFO2*IPAS; 
      if (g_NAFVariable.IPRT==1)
      {
       fprintf(g_NAFVariable.NFPRT,"IFR=%d FR=%g RTAB=%g INDX=%d KTABS2=%d\n",IFR,*FR, RTAB[INDX],INDX,iKTABS2); 
      }
      SYSFREE(RTAB); */
      /*remplacee par:*/
      if (p_iFrMin==iKTABS2)
      {
       INDX=naf_maxx(iKTABS2m1, RTAB);
       IFR=((INDX+1)<=(iKTABS2/2))?INDX:INDX-iKTABS2;
      }
      else  if (p_iFrMin>=0)
      {/* p_iFrMin et p_iFrMax positif */
       IFR=INDX=p_iFrMin+naf_maxx(p_iFrMax-p_iFrMin, RTAB+p_iFrMin);
      }
      else  if (p_iFrMax<0)
      {/* p_iFrMin et p_iFrMax negatif */
       INDX=naf_maxx(p_iFrMax-p_iFrMin, RTAB+iKTABS2+p_iFrMin);
       IFR=INDX+p_iFrMin;
      }
      else 
      {/* p_iFrMin et p_iFrMax de signe differents */
       int INDXNeg;
       /* frequence negative */
       INDXNeg=p_iFrMin+naf_maxx(-1-p_iFrMin, RTAB+iKTABS2+p_iFrMin);
       /* frequence positive */
       INDX=naf_maxx(p_iFrMax, RTAB);
       if (RTAB[INDX]<RTAB[INDXNeg+iKTABS2])
       {
        IFR=INDXNeg;
       }
       else
       {
        IFR=INDX;
       }
      }
      FR = IFR*FREFO2*IPAS; 
      if (g_NAFVariable.IPRT==1)
      {
       fprintf(g_NAFVariable.NFPRT,"IFRMIN=%d IFRMAX=%d IFR=%d FR=%g RTAB=%g INDX=%d KTABS2=%d\n",p_iFrMin, p_iFrMax,
               IFR,FR, RTAB[INDX],INDX,iKTABS2); 
      }
      SYSFREE(RTAB);
      return FR;
      /* v0.96 M. GASTINEAU 14/01/99 : fin de modification */
}/*      END SUBROUTINE FFTMAX*/
#endif /* NAF_USE_OPTIMIZE */
/*!
!
      SUBROUTINE FOUR1(DATA,NN,ISIGN)*/
/*v0.96 M. GASTINEAU 14/12/98 : modification du calcul de la puissance */
void naf_four1(double *DATA, int NN, int ISIGN)
{
/*
!************** FOURIER RAPIDE DE NUMERICAL RECIPES
      implicit none
! (DATA,NN,ISIGN)
       integer :: NN,ISIGN
       REAL (8) :: DATA(2*NN)
! 
      integer :: N,I,J,M,MMAX,ISTEP
      REAL (8) :: THETA,WPR,WPI,WR,WI,TEMPR,TEMPI,WTEMP
*/
      int  N,I,J,M,MMAX,ISTEP;
      double THETA,WPR,WPI,WR,WI,TEMPR,TEMPI,WTEMP;
      N=2*NN ;
      J=1;
      for(I=1;I<=N; I+=2)
      {
/*! en fait N est pair donc I<= N-1    */  
        if(J>I)
        {
          TEMPR=DATA[J];
          TEMPI=DATA[J+1];
          DATA[J]=DATA[I];
          DATA[J+1]=DATA[I+1];
          DATA[I]=TEMPR;
          DATA[I+1]=TEMPI;
        }
        M=N/2;
        while ((M>=2) && (J>M))
        {
          J-=M; /*J=J-M;*/
          M >>=1; /*M=M/2;*/
        };
        J +=M; /*J=J+M*/
    }
      MMAX=2;
      while (N>MMAX)
      {
        ISTEP=2*MMAX;
        THETA=6.28318530717959E0/(ISIGN*MMAX);
        /*v0.96 M. GASTINEAU 14/12/98 : modification */
        /*WPR=-2.E0*pow(sin(0.5E0*THETA),2);*//*remplacee par:*/
        WTEMP=sin(0.5E0*THETA);
        WPR=-2.E0*WTEMP*WTEMP;
        /*v0.96 M. GASTINEAU 14/12/98 : fin modification */
        WPI=sin(THETA);
        WR=1.E0;
        WI=0.E0;
        for ( M=1; M<=MMAX; M+=2)
        {
          for (I=M; I<=N; I+=ISTEP)
          {
            J=I+MMAX;
            TEMPR=WR*DATA[J]-WI*DATA[J+1];
            TEMPI=WR*DATA[J+1]+WI*DATA[J];
/*!**            TEMPR=SNGL[WR]*DATA[J]-SNGL[WI]*DATA[J+1]
!**            TEMPI=SNGL[WR]*DATA[J+1]+SNGL[WI]*DATA[J]*/
            DATA[J]=DATA[I]-TEMPR;
            DATA[J+1]=DATA[I+1]-TEMPI;
            DATA[I]=DATA[I]+TEMPR;
            DATA[I+1]=DATA[I+1]+TEMPI;
          }
          WTEMP=WR;
          WR=WR*WPR-WI*WPI+WR;
          WI=WI*WPR+WTEMP*WPI+WI;
        }
        MMAX=ISTEP;
      }
}/*      END SUBROUTINE FOUR1*/


/*      SUBROUTINE PUISS2(NT,N2)*/
void naf_puiss2(int NT, int *N2)
{
/*!*******************************************
!    CALCULE LA PLUS GRANDE PUISSANCE DE 2
!    CONTENUE DANS NT
!*******************************************
      implicit none
!  (NT,N2)    
      INTEGER ::  N2,NT
!    
      integer :: n
! */
      int N;           
      N=NT;
      if (N==0)
      {
         *N2=0;
      return;
      }
      *N2=1;
      while (N>=2)
      {
         *N2 <<=1;/*N2=N2*2*/
         N >>=1; /*N=N/2*/
         
      }
}/*      END SUBROUTINE PUISS2 */

/*      SUBROUTINE MAXX(N,T,INDX)*/
/*!*******************************************
!    CALCULE LE MAX ET L'INDICE DU MAX
!*******************************************
       implicit none
!  (N,T,INDX)
       integer :: N,INDX
       REAL (8) :: T(0:N)
!
      integer :: J
      REAL (8) VMAX             
!*/
/* v0.96 M. GASTINEAU 12/01/99 : optimisation  passage par retour et non par la pile */
#if NAF_USE_OPTIMIZE==0
void naf_maxx(int N, double *T, int *INDX)
{
      int J;
      double VMAX;
      VMAX=T[0];
      *INDX=0;
      for(J=1; J<=N; J++)
      {
         if (T[J]>VMAX)
         {
            VMAX=T[J];
            *INDX=J;
         }
      }
}/*      END SUBROUTINE MAXX*/
#else /*remplacee par:*/
int naf_maxx(int N, double *T)
{
      int J;
      double VMAX;
      int INDX=0;
      VMAX=T[0];
      for(J=1; J<=N; J++)
      {
         if (T[J]>VMAX)
         {
            VMAX=T[J];
            INDX=J;
         }
      }
      return INDX;
}
#endif /*NAF_USE_OPTIMIZE*/
/* v0.96 M. GASTINEAU 12/01/99 : fin optimisation */
/*      END SUBROUTINE MAXX*/

/*      SUBROUTINE MODTAB(N,T)*/
void naf_modtab(int N, double *T)
{
/*!****************************************
!    SUPPRESSION DE CERTAINS TERMES
!    A MODIFIER ULTERIEUREMENT (13/12/90)
!****************************************
      implicit none
!  (N,T)
      integer :: N
      REAL (8) :: T(0:N) 
!          
      integer NOK,I
! */
      int NOK, I;     
      NOK = 40;
      for(I = 0; I<=NOK; I++)
      {   
       T[I] = 0;
       T[N-I] = 0;
      }
}    /*  END SUBROUTINE MODTAB*/


/*      SUBROUTINE MODFRE(NUMFR,A,B)*/
/*v0.96 M. GASTINEAU 12/01/99 : modification dans le cas ou ICPLX=0 */
/*v0.97 M. GASTINEAU 27/05/99 : modification dans le cas ou ICPLX=0 et FS=0 */
void naf_modfre(int NUMFR, double *A, double *B)
{
/*!-----------------------------------------------------------------------
!     PERMET DE MODIFIER UNE AMPLITUDE DEJA CALCULEE QUAND
!     ON RETROUVE LE MEME TERME, A TOL PRES DANS LE DEVELOPPEMENT
! 
!     NUMFR     INDICE DE LA FREQUENCE EXISTANT DEJA
!     A,B       PARTIES REELLES ET IMAGINAIRES DE L'AMPLITUDE DE LA MODIF 
!               A APPORTER A L'AMPLITUDE DE LA FREQUENCE NUMFR
!
!     AUTRES PARAMETRES TRANSMIS PAR LE COMMON/CGRAM/g_NAFVariable.ZAMP,g_NAFVariable.ZALP,g_NAFVariable.TFS,g_NAFVariable.NFS
!
!     30/4/91
!-----------------------------------------------------------------------
      
!----------! g_NAFVariable.TFS EST LE TABLEAU FINAL DES FREQUENCES
!----------! g_NAFVariable.ZAMP   TABLEAU DES AMPLITUDES COMPLEXES
!----------! g_NAFVariable.ZALP   MATRICE DE PASSAGE DE LA BASE ORTHONORMALISEE
!----------! g_NAFVariable.NFS    NOMBRE DE FREQ.DEJA DETERMINEES
      IMPLICIT NONE
! (NUMFR,A,B)
      integer :: NUMFR
      REAL (8) :: A,B
!            
      integer :: IT
      complex (8) :: ZI,ZOM,ZA,ZEX,ZINC
      complex (8), dimension (:),  allocatable :: ZT
! */
      int IT;
      t_complexe ZI,ZOM,ZA,ZEX,ZINC;
      t_complexe *ZT=NULL;
           
/* v0.96 M. GASTINEAU 01/12/98 : utilisation des fonctions  complexes inlines et optimisation */
#if NAF_USE_OPTIMIZE==0
      SYSCHECKMALLOCSIZE(ZT,t_complexe, g_NAFVariable.KTABS+1); /*allocate(ZT(0:g_NAFVariable.KTABS),stat = NERROR)*/
      ZI = cmplx(0.E0,1.E0);
      ZOM=muldoublcomplexe(g_NAFVariable.TFS[NUMFR]/g_NAFVariable.UNIANG,ZI); /*ZOM=g_NAFVariable.TFS[NUMFR]/g_NAFVariable.UNIANG*ZI*/
      ZA=cmplx(*A,*B);
      fprintf(g_NAFVariable.NFPRT,"CORRECTION DE  IFR = %d AMPLITUDE  = %g",NUMFR, module(ZA));
/*!-----------! L' AMPLITUDES DU TERMES EST CORRIGEES
!-----------! ATTENTION ICI (CAS REEL) ON AURA AUSSI LE TERME CONJUGUE
!-----------! QU'ON NE CALCULE PAS. LE TERME TOTAL EST
!-----------!       2*RE(g_NAFVariable.ZAMP(I)*EXP(ZI*g_NAFVariable.TFS(I)*T) )*/
      g_NAFVariable.ZAMP[NUMFR]=addcomplexe(g_NAFVariable.ZAMP[NUMFR],ZA);
      if (g_NAFVariable.IPRT==1)
      { 
       fprintf(g_NAFVariable.NFPRT," %+-20.15E %+-20.15E %+-20.15E %+-20.15E %+-20.15E\n",g_NAFVariable.TFS[NUMFR],
                     module(g_NAFVariable.ZAMP[NUMFR]),g_NAFVariable.ZAMP[NUMFR].reel,
                     g_NAFVariable.ZAMP[NUMFR].imag,atan2(g_NAFVariable.ZAMP[NUMFR].imag,g_NAFVariable.ZAMP[NUMFR].reel));
      }
/*!-----------! ON RETIRE LA CONTRIBUTION DU TERME g_NAFVariable.TFS(NUMFR) DANS TABS*/
      if (g_NAFVariable.ICPLX==1)
      {
         ZEX=mulcomplexe(ZA,expcomplexe(muldoublcomplexe(g_NAFVariable.T0-g_NAFVariable.XH,ZOM))); /*ZEX = ZA*EXP(ZOM*(g_NAFVariable.T0-g_NAFVariable.XH))*/
         ZINC=expcomplexe(muldoublcomplexe(g_NAFVariable.XH,ZOM)); /*ZINC=EXP(ZOM*g_NAFVariable.XH)*/
         naf_ztpow(g_NAFVariable.KTABS,64,ZT,ZINC,ZEX); /*CALL  ZTPOW (KTABS+1,64,ZT,ZINC,ZEX)*/
         for(IT=0;IT<=g_NAFVariable.KTABS;IT++)
         {
            g_NAFVariable.ZTABS[IT]=subcomplexe(g_NAFVariable.ZTABS[IT], ZT[IT]); 
         }
      }
      else
      {
         ZEX=mulcomplexe(ZA,expcomplexe(muldoublcomplexe(g_NAFVariable.T0-g_NAFVariable.XH,ZOM))); /*ZEX = ZA*EXP(ZOM*(g_NAFVariable.T0-g_NAFVariable.XH))*/
         ZINC=expcomplexe(muldoublcomplexe(g_NAFVariable.XH,ZOM)); /*ZINC=EXP(ZOM*XH)*/
         naf_ztpow(g_NAFVariable.KTABS,64,ZT,ZINC,ZEX); /*CALL  ZTPOW (KTABS+1,64,ZT,ZINC,ZEX)*/
         /*v0.97 M. GASTINEAU 27/05/99 : ajout - correction bug si FS==0 */
         if (g_NAFVariable.TFS[NUMFR]==0.E0)
         {
          for(IT=0;IT<=g_NAFVariable.KTABS;IT++)
          {
           g_NAFVariable.ZTABS[IT].reel -= ZT[IT].reel; 
           /*g_NAFVariable.ZTABS(IT+1)=g_NAFVariable.ZTABS(IT+1)- DREAL(ZT(IT+1)) */
          }
         }
         else
         /*v0.97 M. GASTINEAU 27/05/99 : fin ajout */
         for(IT=0;IT<=g_NAFVariable.KTABS;IT++)
         {
            /*v0.96 M. GASTINEAU 12/01/99 : modification */
            /*g_NAFVariable.ZTABS[IT].reel -=ZT[IT].reel;*/
            /*remplacee par:*/
            g_NAFVariable.ZTABS[IT].reel -= 2*ZT[IT].reel;
            /*v0.96 M. GASTINEAU 12/01/99 : fin modification */
            /* g_NAFVariable.ZTABS(IT)=g_NAFVariable.ZTABS(IT)- DREAL(ZT(IT)) */
         }
      }
      SYSFREE(ZT); /*deallocate(ZT)*/
#else /*remplacee par:*/
      t_complexe *pzarTabs, *pzarZT;
      const int ikTabs=g_NAFVariable.KTABS; /*v0.96 M. GASTINEAU 12/01/99 : optimisation*/
      
      SYSCHECKMALLOCSIZE(ZT,t_complexe, g_NAFVariable.KTABS+1); /*allocate(ZT(0:g_NAFVariable.KTABS),stat = NERROR)*/
      i_compl_cmplx(&ZI,0.E0,1.E0);
      ZOM=i_compl_muldoubl(g_NAFVariable.TFS[NUMFR]/g_NAFVariable.UNIANG,ZI); /*ZOM=g_NAFVariable.TFS[NUMFR]/g_NAFVariable.UNIANG*ZI*/
      i_compl_cmplx(&ZA,*A,*B);
      fprintf(g_NAFVariable.NFPRT,"CORRECTION DE  IFR = %d AMPLITUDE  = %g",NUMFR, i_compl_module(ZA));
/*!-----------! L' AMPLITUDES DU TERMES EST CORRIGEES
!-----------! ATTENTION ICI (CAS REEL) ON AURA AUSSI LE TERME CONJUGUE
!-----------! QU'ON NE CALCULE PAS. LE TERME TOTAL EST
!-----------!       2*RE(g_NAFVariable.ZAMP(I)*EXP(ZI*g_NAFVariable.TFS(I)*T) )*/
      i_compl_padd(g_NAFVariable.ZAMP+NUMFR,&ZA);
      if (g_NAFVariable.IPRT==1)
      { 
       fprintf(g_NAFVariable.NFPRT," %+-20.15E %+-20.15E %+-20.15E %+-20.15E %+-20.15E\n",g_NAFVariable.TFS[NUMFR],
                     i_compl_module(g_NAFVariable.ZAMP[NUMFR]),g_NAFVariable.ZAMP[NUMFR].reel,
                     g_NAFVariable.ZAMP[NUMFR].imag,atan2(g_NAFVariable.ZAMP[NUMFR].imag,g_NAFVariable.ZAMP[NUMFR].reel));
      }
/*!-----------! ON RETIRE LA CONTRIBUTION DU TERME g_NAFVariable.TFS(NUMFR) DANS TABS*/
      if (g_NAFVariable.ICPLX==1)
      {
         ZEX=i_compl_mul(ZA,i_compl_exp(i_compl_muldoubl(g_NAFVariable.T0-g_NAFVariable.XH,ZOM))); /*ZEX = ZA*EXP(ZOM*(g_NAFVariable.T0-g_NAFVariable.XH))*/
         ZINC=i_compl_exp(i_compl_muldoubl(g_NAFVariable.XH,ZOM)); /*ZINC=EXP(ZOM*g_NAFVariable.XH)*/
         naf_ztpow(g_NAFVariable.KTABS,64,ZT,ZINC,ZEX); /*CALL  ZTPOW (KTABS+1,64,ZT,ZINC,ZEX)*/
         /*v0.96 M. GASTINEAU 12/01/99 : optimisation*/
         /*for(IT=0, pzarTabs=g_NAFVariable.ZTABS, pzarZT = ZT;
             IT<=g_NAFVariable.KTABS;
             IT++, pzarTabs++, pzarZT++)*//*remplacee par:*/
         for(IT=0, pzarTabs=g_NAFVariable.ZTABS, pzarZT = ZT;
             IT<=ikTabs;
             IT++, pzarTabs++, pzarZT++)/*v0.96 M. GASTINEAU 12/01/99 : fin optimisation*/
         {
            i_compl_psub(pzarTabs,pzarZT);
         }
      }
      else
      {
         ZEX=i_compl_mul(ZA,i_compl_exp(i_compl_muldoubl(g_NAFVariable.T0-g_NAFVariable.XH,ZOM))); /*ZEX = ZA*EXP(ZOM*(g_NAFVariable.T0-g_NAFVariable.XH))*/
         ZINC=i_compl_exp(i_compl_muldoubl(g_NAFVariable.XH,ZOM)); /*ZINC=EXP(ZOM*XH)*/
         naf_ztpow(g_NAFVariable.KTABS,64,ZT,ZINC,ZEX); /*CALL  ZTPOW (KTABS+1,64,ZT,ZINC,ZEX)*/
         /*v0.96 M. GASTINEAU 12/01/99 : optimisation*/
         /*for(IT=0;IT<=g_NAFVariable.KTABS;IT++)*/
         /*v0.97 M. GASTINEAU 27/05/99 : ajout - correction bug si FS==0 */
         if (g_NAFVariable.TFS[NUMFR]==0.E0)
         {
          for(IT=0;IT<=ikTabs;IT++)
          {
           g_NAFVariable.ZTABS[IT].reel -= ZT[IT].reel; 
           /*g_NAFVariable.ZTABS(IT+1)=g_NAFVariable.ZTABS(IT+1)- DREAL(ZT(IT+1)) */
          }
         }
         else
         /*v0.97 M. GASTINEAU 27/05/99 : fin ajout */
         for(IT=0;IT<=ikTabs;IT++)/*v0.96 M. GASTINEAU 12/01/99 : fin optimisation*/
         {
           /*v0.96 M. GASTINEAU 12/01/99 : modification */
           /* g_NAFVariable.ZTABS[IT].reel -=ZT[IT].reel;*//* g_NAFVariable.ZTABS(IT)=g_NAFVariable.ZTABS(IT)- DREAL(ZT(IT)) */
           /*remplacee par:*/
           g_NAFVariable.ZTABS[IT].reel -=2*ZT[IT].reel;/* g_NAFVariable.ZTABS(IT)=g_NAFVariable.ZTABS(IT)- DREAL(ZT(IT)) */
           /*v0.96 M. GASTINEAU 12/01/99 : fin modification */        
         }
      }
      SYSFREE(ZT); /*deallocate(ZT)*/
#endif /*NAF_USE_OPTIMIZE==0*/
/* v0.96 M. GASTINEAU 01/12/98 : fin optimisation */
}/*      END SUBROUTINE MODFRE*/


/*      SUBROUTINE GRAMSC(FS,A,B)*/
/*v0.96 M. GASTINEAU 12/01/99 : modification dans le cas ou ICPLX=0 */
/*v0.97 M. GASTINEAU 27/05/99 : modification dans le cas ou ICPLX=0 et FS=0 */
BOOL naf_gramsc(double FS, double A, double B)
{
/*
!-----------------------------------------------------------------------
!     GRAMSC   ORTHONORMALISATION  PAR GRAM-SCHIMDT DE LA BASE DE FONCTI
!               CONSTITUEE DES TERMES PERIODIQUES TROUVES PAR ANALYSE
!               SPECTRALE .
!     FS        FREQUENCE EN "/AN
!     A,B       PARTIES REELLES ET IMAGINAIRES DE L'AMPLITUDE
!
!     AUTRES PARAMETRES TRANSMIS PAR LE COMMON/CGRAM/g_NAFVariable.ZAMP,g_NAFVariable.ZALP,g_NAFVariable.TFS,g_NAFVariable.NFS
!
!     MODIFIE LE 26 9 87 POUR LES FONCTIONS REELLES   J. LASKAR
!-----------------------------------------------------------------------
      
!----------! g_NAFVariable.TFS EST LE TABLEAU FINAL DES FREQUENCES
!----------! g_NAFVariable.ZAMP   TABLEAU DES AMPLITUDES COMPLEXES
!----------! g_NAFVariable.ZALP   MATRICE DE PASSAGE DE LA BASE ORTHONORMALISEE
!----------! ZTEE   TABLEAU DE TRAVAIL
!----------! g_NAFVariable.NFS    NOMBRE DE FREQ.DEJA DETERMINEES
      IMPLICIT NONE
! (g_NAFVariable.KTABS,FS,A,B)
      REAL (8) :: FS,A,B     
!
      integer :: I, J, K, NF, IT
      REAL (8) DIV
      complex (8), dimension(:),allocatable :: ZTEE
      complex (8) :: ZDIV,ZMUL,ZI,ZEX,ZINC,ZA,ZOM
      complex (8), dimension (:),allocatable :: ZT
! */
      int   I, J, K, NF, IT;
      double DIV;
      t_complexe *ZTEE=NULL;
      t_complexe ZDIV,ZMUL,ZI,ZEX,ZINC,ZA,ZOM;
      t_complexe *ZT=NULL;  
/* v0.96 M. GASTINEAU 01/12/98 : utilisation des fonctions  complexes inlines et optimisation */
#if NAF_USE_OPTIMIZE>0
      t_complexe *pzarTabs, *pzarZT;
      const int ikTabs=g_NAFVariable.KTABS; /*v0.96 M. GASTINEAU 12/01/99 : optimisation*/
#endif /**/
      SYSCHECKMALLOCSIZE(ZT, t_complexe, g_NAFVariable.KTABS+1); /*allocate(ZT(0:g_NAFVariable.KTABS),stat = NERROR)*/
      SYSCHECKMALLOCSIZE(ZTEE, t_complexe, g_NAFVariable.NTERM+1); /*allocate(ZTEE(1:g_NAFVariable.NTERM),stat = NERROR)*/

/*!----------! CALCUL DE ZTEE(I)=<EN,EI>*/
      for(I =1;I<=g_NAFVariable.NFS;I++)
      {
        if(naf_proscaa(FS,g_NAFVariable.TFS[I],ZTEE+I)==FALSE)
        {
         SYSFREE(ZT);
         SYSFREE(ZTEE);
         return FALSE;
        }
      }
/* v0.96 M. GASTINEAU 01/12/98 : utilisation des fonctions  complexes inlines et optimisation */
#if NAF_USE_OPTIMIZE==0
/*!----------! NF EST LE NUMERO DE LA NOUVELLE FREQUENCE*/
      NF=g_NAFVariable.NFS+1;
      ZTEE[NF]=cmplx(1.E0,0.E0);
/*!----------! CALCUL DE FN = EN - SOM(<EN,FI>FI) QUI EST ORTHOGONAL AUX F*/
      g_NAFVariable.TFS[NF]=FS;
      for( K=1;K<=g_NAFVariable.NFS;K++)
      {
       for(I=K;I<=g_NAFVariable.NFS;I++)
       {
        for(J=1;J<=I;J++)
        {
          /*g_NAFVariable.ZALP(NF,K)=g_NAFVariable.ZALP(NF,K)-DCONJG(g_NAFVariable.ZALP(I,J))*g_NAFVariable.ZALP(I,K)*ZTEE(J)*/
          g_NAFVariable.ZALP[NF][K]=subcomplexe(g_NAFVariable.ZALP[NF][K],mulcomplexe(mulcomplexe(conjcomplexe(g_NAFVariable.ZALP[I][J]),g_NAFVariable.ZALP[I][K]),ZTEE[J]));
        }
       }
      }
      g_NAFVariable.ZALP[NF][NF]=cmplx(1.E0,0.E0);
/*!----------! ON REND LA NORME DE FN = 1*/
      DIV=1.E0;
      ZDIV=cmplx(0.E0,0.E0);
      for( I=1; I<=NF; I++)
      {
         /*ZDIV=ZDIV+DCONJG(g_NAFVariable.ZALP(NF,I))*ZTEE(I)*/
         ZDIV=addcomplexe(ZDIV, mulcomplexe(conjcomplexe(g_NAFVariable.ZALP[NF][I]),ZTEE[I]));
      }
      DIV=sqrt(module(ZDIV));
      if (g_NAFVariable.IPRT==1)
      {
       /*v0.96 M. GASTINEAU 19/11/98 : modification*/
       /*fprintf(g_NAFVariable.NFPRT,"ZDIV= %g+i%g DIV=%g\n",ZDIV,DIV);*//*remplacee par:*/
       fprintf(g_NAFVariable.NFPRT,"ZDIV= %g+i%g DIV=%g\n",ZDIV.reel,ZDIV.imag,DIV);
       /*v0.96 M. GASTINEAU 19/11/98 : fin modification*/
      }
      for(I=1; I<=NF; I++)
      {
         g_NAFVariable.ZALP[NF][I]=muldoublcomplexe(1.E0/DIV,g_NAFVariable.ZALP[NF][I]); /*g_NAFVariable.ZALP(NF,I) = g_NAFVariable.ZALP(NF,I)/DIV*/
      }
/*!-----------! F1,F2....., FN EST UNE BASE ORTHONORMEE
!-----------! ON RETIRE MAINTENANT A F  <F,FN>FN  (<F,FN>=<F,EN>)*/
      ZMUL=muldoublcomplexe(1.E0/DIV, cmplx(A,B));
      ZI=cmplx(0.E0,1.E0);
      for(I=1; I<=NF; I++)
      {
         ZOM=muldoublcomplexe(g_NAFVariable.TFS[I]/g_NAFVariable.UNIANG,ZI); /*ZOM=g_NAFVariable.TFS(I)/g_NAFVariable.UNIANG*ZI*/
         ZA=mulcomplexe(g_NAFVariable.ZALP[NF][I],ZMUL);
/*!-----------! LES AMPLITUDES DES TERMES SONT CORRIGEES
!-----------! ATTENTION ICI (CAS REEL) ON AURA AUSSI LE TERME CONJUGUE
!-----------! QU'ON NE CALCULE PAS. LE TERME TOTAL EST
!-----------!       2*RE(g_NAFVariable.ZAMP(I)*EXP(ZI*g_NAFVariable.TFS(I)*T) )*/
         g_NAFVariable.ZAMP[I]=addcomplexe(g_NAFVariable.ZAMP[I],ZA);
         if (g_NAFVariable.IPRT==1)
         {
           fprintf(g_NAFVariable.NFPRT," %g %g %g %g %g\n", g_NAFVariable.TFS[I],module(g_NAFVariable.ZAMP[I]),g_NAFVariable.ZAMP[I].reel,
                   g_NAFVariable.ZAMP[I].imag,atan2(g_NAFVariable.ZAMP[I].imag,g_NAFVariable.ZAMP[I].reel));
         }
/*!-----------! ON RETIRE LA CONTRIBUTION TOTALE DU TERME g_NAFVariable.TFS(I) DANS TABS*/
       if (g_NAFVariable.ICPLX==1)
       {
         ZEX=mulcomplexe(ZA, expcomplexe(muldoublcomplexe(g_NAFVariable.T0-g_NAFVariable.XH,ZOM))); /*ZEX = ZA*EXP(ZOM*(g_NAFVariable.T0-g_NAFVariable.XH))*/
         ZINC=expcomplexe(muldoublcomplexe(g_NAFVariable.XH,ZOM)); /*ZINC=EXP(ZOM*g_NAFVariable.XH)*/
         naf_ztpow(g_NAFVariable.KTABS,64,ZT,ZINC,ZEX); /*naf_ztpow(g_NAFVariable.KTABS+1,64,ZT,ZINC,ZEX);*/
         for(IT=0;IT<=g_NAFVariable.KTABS;IT++)
         {
            g_NAFVariable.ZTABS[IT]=subcomplexe(g_NAFVariable.ZTABS[IT],ZT[IT]); 
         }
       }
       else
       {
         ZEX=mulcomplexe(ZA, expcomplexe(muldoublcomplexe(g_NAFVariable.T0-g_NAFVariable.XH,ZOM))); /*ZEX = ZA*EXP(ZOM*(g_NAFVariable.T0-g_NAFVariable.XH))*/
         ZINC=expcomplexe(muldoublcomplexe(g_NAFVariable.XH,ZOM)); /*ZINC=EXP(ZOM*g_NAFVariable.XH)*/
         naf_ztpow(g_NAFVariable.KTABS,64,ZT,ZINC,ZEX);/*naf_ztpow(g_NAFVariable.KTABS+1,64,ZT,ZINC,ZEX);*/
         /*v0.97 M. GASTINEAU 27/05/99 : ajout - correction bug si FS==0 */
         if (FS==0.E0)
         {
         for(IT=0;IT<=g_NAFVariable.KTABS;IT++)
         {
           g_NAFVariable.ZTABS[IT].reel -= ZT[IT].reel; 
           /*g_NAFVariable.ZTABS(IT+1)=g_NAFVariable.ZTABS(IT+1)- DREAL(ZT(IT+1)) */
         }
         }
         else
         /*v0.97 M. GASTINEAU 27/05/99 : fin ajout */
         for(IT=0;IT<=g_NAFVariable.KTABS;IT++)
         {
           /*v0.96 M. GASTINEAU 12/01/99 : modification*/
           /*g_NAFVariable.ZTABS[IT].reel -=ZT[IT].reel;*/ /*g_NAFVariable.ZTABS(IT+1)=g_NAFVariable.ZTABS(IT+1)- DREAL(ZT(IT+1)) */
           /*remplacee par: */
           g_NAFVariable.ZTABS[IT].reel -= 2*ZT[IT].reel; /*g_NAFVariable.ZTABS(IT+1)=g_NAFVariable.ZTABS(IT+1)- DREAL(ZT(IT+1)) */
           /*v0.96 M. GASTINEAU 12/01/99 : fin modification*/
         }
       }
      }
#else /*remplacee par:*/
/*!----------! NF EST LE NUMERO DE LA NOUVELLE FREQUENCE*/
      NF=g_NAFVariable.NFS+1;
      i_compl_cmplx(ZTEE+NF,1.E0,0.E0);
/*!----------! CALCUL DE FN = EN - SOM(<EN,FI>FI) QUI EST ORTHOGONAL AUX F*/
      g_NAFVariable.TFS[NF]=FS;
      for( K=1;K<=g_NAFVariable.NFS;K++)
      {
       for(I=K;I<=g_NAFVariable.NFS;I++)
       {
        for(J=1;J<=I;J++)
        {
          /*g_NAFVariable.ZALP(NF,K)=g_NAFVariable.ZALP(NF,K)-DCONJG(g_NAFVariable.ZALP(I,J))*g_NAFVariable.ZALP(I,K)*ZTEE(J)*/
          t_complexe zSubExp;
          zSubExp=i_compl_mul(i_compl_mul(i_compl_conj(&(g_NAFVariable.ZALP[I][J])),g_NAFVariable.ZALP[I][K]),ZTEE[J]);
          i_compl_psub(&(g_NAFVariable.ZALP[NF][K]),&zSubExp);
        }
       }
      }
      i_compl_cmplx(&(g_NAFVariable.ZALP[NF][NF]),1.E0,0.E0);
/*!----------! ON REND LA NORME DE FN = 1*/
      DIV=1.E0;
      i_compl_cmplx(&ZDIV,0.E0,0.E0);
      for( I=1; I<=NF; I++)
      {
         /*ZDIV=ZDIV+DCONJG(g_NAFVariable.ZALP(NF,I))*ZTEE(I)*/
         t_complexe zSubExp;
         zSubExp = i_compl_mul(i_compl_conj(&(g_NAFVariable.ZALP[NF][I])),ZTEE[I]);
         i_compl_padd(&ZDIV,&zSubExp);
      }
      DIV=sqrt(i_compl_module(ZDIV));
      if (g_NAFVariable.IPRT==1)
      {
       /*v0.96 M. GASTINEAU 19/11/98 : modification*/
       /*fprintf(g_NAFVariable.NFPRT,"ZDIV= %g+i%g DIV=%g\n",ZDIV,DIV);*//*remplacee par:*/
       fprintf(g_NAFVariable.NFPRT,"ZDIV= %g+i%g DIV=%g\n",ZDIV.reel,ZDIV.imag,DIV);
       /*v0.96 M. GASTINEAU 19/11/98 : fin modification*/
      }
      for(I=1; I<=NF; I++)
      {
         i_compl_pdivdoubl(&(g_NAFVariable.ZALP[NF][I]),&DIV); /*g_NAFVariable.ZALP(NF,I) = g_NAFVariable.ZALP(NF,I)/DIV*/
      }
/*!-----------! F1,F2....., FN EST UNE BASE ORTHONORMEE
!-----------! ON RETIRE MAINTENANT A F  <F,FN>FN  (<F,FN>=<F,EN>)*/
      /*v0.96 M. GASTINEAU 01/12/98 : modification */
      i_compl_cmplx(&ZMUL,A/DIV,B/DIV);
      i_compl_cmplx(&ZI,0.E0,1.E0);
      for(I=1; I<=NF; I++)
      {
         ZOM=i_compl_muldoubl(g_NAFVariable.TFS[I]/g_NAFVariable.UNIANG,ZI); /*ZOM=g_NAFVariable.TFS(I)/g_NAFVariable.UNIANG*ZI*/
         ZA=i_compl_mul(g_NAFVariable.ZALP[NF][I],ZMUL);
/*!-----------! LES AMPLITUDES DES TERMES SONT CORRIGEES
!-----------! ATTENTION ICI (CAS REEL) ON AURA AUSSI LE TERME CONJUGUE
!-----------! QU'ON NE CALCULE PAS. LE TERME TOTAL EST
!-----------!       2*RE(g_NAFVariable.ZAMP(I)*EXP(ZI*g_NAFVariable.TFS(I)*T) )*/
         i_compl_padd(g_NAFVariable.ZAMP+I,&ZA);
         if (g_NAFVariable.IPRT==1)
         {
           fprintf(g_NAFVariable.NFPRT," %g %g %g %g %g\n", g_NAFVariable.TFS[I],i_compl_module(g_NAFVariable.ZAMP[I]),g_NAFVariable.ZAMP[I].reel,
                   g_NAFVariable.ZAMP[I].imag,atan2(g_NAFVariable.ZAMP[I].imag,g_NAFVariable.ZAMP[I].reel));
         }
/*!-----------! ON RETIRE LA CONTRIBUTION TOTALE DU TERME g_NAFVariable.TFS(I) DANS TABS*/
       if (g_NAFVariable.ICPLX==1)
       {
         ZEX=i_compl_mul(ZA, i_compl_exp(i_compl_muldoubl(g_NAFVariable.T0-g_NAFVariable.XH,ZOM))); /*ZEX = ZA*EXP(ZOM*(g_NAFVariable.T0-g_NAFVariable.XH))*/
         ZINC=i_compl_exp(i_compl_muldoubl(g_NAFVariable.XH,ZOM)); /*ZINC=EXP(ZOM*g_NAFVariable.XH)*/
         naf_ztpow(g_NAFVariable.KTABS,64,ZT,ZINC,ZEX); /*naf_ztpow(g_NAFVariable.KTABS+1,64,ZT,ZINC,ZEX);*/
         /*v0.96 M. GASTINEAU 12/01/99 : optimisation*/
         /*for(IT=0, pzarTabs=g_NAFVariable.ZTABS, pzarZT=ZT;
             IT<=g_NAFVariable.KTABS;
             IT++,pzarTabs++,pzarZT++)*//*remplacee par:*/
         for(IT=0, pzarTabs=g_NAFVariable.ZTABS, pzarZT=ZT;
             IT<=ikTabs;
             IT++,pzarTabs++,pzarZT++) /*v0.96 M. GASTINEAU 12/01/99 : optimisation*/
         {
            i_compl_psub(pzarTabs, pzarZT);
         }
       }
       else
       {
         ZEX=i_compl_mul(ZA, i_compl_exp(i_compl_muldoubl(g_NAFVariable.T0-g_NAFVariable.XH,ZOM))); /*ZEX = ZA*EXP(ZOM*(g_NAFVariable.T0-g_NAFVariable.XH))*/
         ZINC=i_compl_exp(i_compl_muldoubl(g_NAFVariable.XH,ZOM)); /*ZINC=EXP(ZOM*g_NAFVariable.XH)*/
         naf_ztpow(g_NAFVariable.KTABS,64,ZT,ZINC,ZEX);/*naf_ztpow(g_NAFVariable.KTABS+1,64,ZT,ZINC,ZEX);*/
         /*v0.96 M. GASTINEAU 12/01/99 : optimisation*/
         /*for(IT=0;IT<=g_NAFVariable.KTABS;IT++)*//*remplacee par:*/
         /*v0.97 M. GASTINEAU 27/05/99 : ajout - correction bug si FS==0 */
         if (FS==0.E0)
         {
         for(IT=0;IT<=ikTabs;IT++)
         {
           g_NAFVariable.ZTABS[IT].reel -= ZT[IT].reel; 
           /*g_NAFVariable.ZTABS(IT+1)=g_NAFVariable.ZTABS(IT+1)- DREAL(ZT(IT+1)) */
         }
         }
         else
         /*v0.97 M. GASTINEAU 27/05/99 : fin ajout */
         for(IT=0;IT<=ikTabs;IT++)
          /*v0.96 M. GASTINEAU 12/01/99 : fin optimisation*/
         {
           /*v0.96 M. GASTINEAU 12/01/99 : modification*/
           /* g_NAFVariable.ZTABS[IT].reel -=ZT[IT].reel; *//*g_NAFVariable.ZTABS(IT+1)=g_NAFVariable.ZTABS(IT+1)- DREAL(ZT(IT+1)) */
           /*remplacee par:*/
           g_NAFVariable.ZTABS[IT].reel -= 2*ZT[IT].reel; /*g_NAFVariable.ZTABS(IT+1)=g_NAFVariable.ZTABS(IT+1)- DREAL(ZT(IT+1)) */
           /*v0.96 M. GASTINEAU 12/01/99 : fin modification*/
         }
       }
      }
#endif /*NAF_USE_OPTIMIZE==0*/
/* v0.96 M. GASTINEAU 01/12/98 : fin optimisation */
      SYSFREE(ZT);
      SYSFREE(ZTEE);
      return TRUE;
}/*      END SUBROUTINE GRAMSC*/


/*      SUBROUTINE PROSCA(F1,F2,ZP)*/
void naf_prosca(double F1, double F2, t_complexe* ZP)
{
/*!-----------------------------------------------------------------------
!     PROSCA   CALCULE LE PRODUIT SCALAIRE  DE EXP(I*F1*T) PAR EXP (-I*F
!               SUR L'INTERVALLE [0:g_NAFVariable.KTABS]  T=g_NAFVariable.T0+g_NAFVariable.XH*IT
!              CALCUL ANALYTIQUE
!     ZP        PRODUIT SCALAIRE COMPLEXE
!     F1,F2     FREQUENCES EN "/AN
!               REVU LE 26/9/87 J. LASKAR
!-----------------------------------------------------------------------
      
      IMPLICIT NONE
! (g_NAFVariable.KTABS,F1,F2,ZP)
      REAL (8) :: F1,F2
      complex (8) :: ZP
!
      complex (8) ZI,ZF,ZF1,ZF2
      REAL (8) PICARRE, T1,T2, XT,DIV,FAC,T,FR1,FR2
* */
      t_complexe ZI,ZF,ZF1,ZF2;
      double  PICARRE, T1,T2, XT,DIV,FACTEUR,T,FR1,FR2;               
/* v0.96 M. GASTINEAU 01/12/98 : utilisation des fonctions  complexes inlines et optimisation */
#if NAF_USE_OPTIMIZE==0
      ZI=cmplx(0.E0,1.E0);
#else /*remplacee par:*/
      i_compl_cmplx(&ZI,0.E0,1.E0);
#endif /*NAF_USE_OPTIMIZE==0*/
/* v0.96 M. GASTINEAU 01/12/98 : fin modification*/
/*!----------! FREQUENCES EN UNITE D'ANGLE PAR UNITE DE TEMPS*/
      FR1=F1/g_NAFVariable.UNIANG;
      FR2=F2/g_NAFVariable.UNIANG;
      T1 =g_NAFVariable.T0;
      T2 =g_NAFVariable.T0+g_NAFVariable.KTABS*g_NAFVariable.XH;
/* v0.96 M. GASTINEAU 01/12/98 : utilisation des fonctions  complexes inlines et optimisation */
#if NAF_USE_OPTIMIZE==0
      /*v0.96 M. GASTINEAU 01/12/98 : correction bug */
      /*ZF=muldoublcomplexe((FR1-FR2)*(T2-T1), ZF);*/ /*ZF=ZI*(FR1-FR2)*(T2-T1)*/
      /*remplacee par: */ZF=muldoublcomplexe((FR1-FR2)*(T2-T1), ZI);
      ZF1=muldoublcomplexe((FR1-FR2)*T1, ZI); /*ZF1=ZI*(FR1-FR2)*T1*/
      ZF2=muldoublcomplexe((FR1-FR2)*T2, ZI); /*ZF2=ZI*(FR1-FR2)*T2*/
      *ZP=divcomplexe(subcomplexe(expcomplexe(ZF2),expcomplexe(ZF1)),ZF); /*ZP=(EXP(ZF2)-EXP(ZF1))/ZF*/
#else /*remplacee par:*/
      ZF=i_compl_muldoubl((FR1-FR2)*(T2-T1), ZI); /*ZF=ZI*(FR1-FR2)*(T2-T1)*/
      ZF1=i_compl_muldoubl((FR1-FR2)*T1, ZI); /*ZF1=ZI*(FR1-FR2)*T1*/
      ZF2=i_compl_muldoubl((FR1-FR2)*T2, ZI); /*ZF2=ZI*(FR1-FR2)*T2*/
      *ZP=i_compl_div(i_compl_sub(i_compl_exp(ZF2),i_compl_exp(ZF1)),ZF); /*ZP=(EXP(ZF2)-EXP(ZF1))/ZF*/
#endif /*NAF_USE_OPTIMIZE==0*/
/* v0.96 M. GASTINEAU 01/12/98 : fin optimisation */
/*!
!      PI=2.D0*ATAN2(1.D0,0.D0)*/
      PICARRE=M_PI*M_PI;
      T=(T2-T1)/2.E0;
      XT=(FR1-FR2)*T;
      DIV=XT*XT-PICARRE;
      if(fabs(DIV)<1.E-10)
      {
         /*v0.96 M. GASTINEAU 19/11/98 : modification*/
         /*fprintf(stdout, "g_NAFVariable.T0= %g, g_NAFVariable.XH=%g, g_NAFVariable.KTABS==%g\n",
                   g_NAFVariable.T0,g_NAFVariable.XH,g_NAFVariable.KTABS);*/
         /*remplacee par:*/
         fprintf(stdout, "g_NAFVariable.T0= %g, g_NAFVariable.XH=%g, g_NAFVariable.KTABS=%d\n",
                 g_NAFVariable.T0, g_NAFVariable.XH, g_NAFVariable.KTABS);
         /*v0.96 M. GASTINEAU 19/11/98 : fin modification*/
         fprintf(stdout, "F1=%g ,F2=%g\n",F1,F2);
         fprintf(stdout, "T1= %g , T2= %g, T= %g\n",T1, T2, T);
         fprintf(stdout, "FR1= %g ,FR2= %g ,XT= %g \n",FR1,FR2,XT);
         fprintf(stdout, "DIV= %g\n",DIV);
      } else
      {
         FACTEUR=-PICARRE/DIV;
/* v0.96 M. GASTINEAU 01/12/98 : utilisation des fonctions  complexes inlines et optimisation */
#if NAF_USE_OPTIMIZE==0
         *ZP=muldoublcomplexe(FACTEUR,*ZP);
#else /*remplacee par:*/
         i_compl_pmuldoubl(ZP,&FACTEUR);
#endif /*NAF_USE_OPTIMIZE==0*/
/* v0.96 M. GASTINEAU 01/12/98 : fin optimisation */
      }
}/*      END  SUBROUTINE PROSCA*/
      
/*      FUNCTION FUNC(X)  */
double naf_func(double X)
{
/*      IMPLICIT NONE
! (X)
      REAL (8) :: X,FUNC
!
      REAL (8) ::RMD
!            */
      double RMD;
      naf_profre(X,&AF,&BF,&RMD);
      return RMD;
}/*      END FUNCTION FUNC*/
     
/*      FUNCTION FUNCP(X)*/
double naf_funcp(double X)
{
/*      IMPLICIT NONE
!  (X)
      REAL (8) :: X,FUNCP
!          
      REAL (8) :: DER, RMF
! */
      double DER, RMF;     
      naf_proder(X,&DER,&AF,&BF,&RMF);
      return DER;
}/*      END FUNCTION FUNCP*/

/*      SUBROUTINE PRODER(FS,DER,A,B,RM)*/
void naf_proder(double FS, double *DER, double *A, double *B,double *RM)
{
/*!***********************************************************************
!  CE SOUS PROG. CALCULE LA DERIVEE DU CARRE DU MODULE DU PRODUIT 
!  SCALAIRE DE EXP(-I*FS)*T PAR TAB(0:7
!  UTILISE LA METHODE D'INTEGRATION DE HARDY
!  LE RESULTAT EST DANS DER 
!  A +IB LE PRODUIT SCALAIRE , RM SON MODULE
!  FS EST DONNEE EN SECONDES PAR AN
!               REVU LE 26/9/87 J. LASKAR
!                       26/5/98 (CHGT DE SIGNE) (*m/4)
!***********************************************************************
      
      IMPLICIT NONE
! (FS,DER,A,B,RM)
      REAL (8) :: FS,DER,A,B,RM      
! 
      integer,parameter :: NVECT=64
      REAL (8) :: OM,ANG0,ANGI,H
      complex (8) :: ZI,ZAC,ZINC,ZEX,ZB,ZT,ZA
      integer :: LTF
      complex (8), dimension (:), allocatable :: ZTF
! */
#define NVECT 64
      double  OM,ANG0,ANGI,H;    
      t_complexe ZI,ZAC,ZINC,ZEX,ZB,/*ZT,*/ZA;
      int LTF;
      t_complexe *ZTF=NULL;/*tableau de 1 a KTABS+1 */
      SYSCHECKMALLOCSIZE(ZTF,t_complexe, g_NAFVariable.KTABS+1); /*allocate(ZTF(0:g_NAFVariable.KTABS),stat = NERROR)*/
/* v0.96 M. GASTINEAU 01/12/98 : utilisation des fonctions  complexes inlines et optimisation */
#if NAF_USE_OPTIMIZE==0
/*!
!------------------ CONVERSION DE FS EN RD/AN*/
      OM=FS/g_NAFVariable.UNIANG;
      ZI=cmplx(0.E0,1.E0);
/*!------------------ LONGUEUR DU TABLEAU A INTEGRER
!------------------ DANS LE CAS REEL MEME CHOSE */
      LTF=g_NAFVariable.KTABS;
      ANG0=OM*g_NAFVariable.T0;
      ANGI=OM*g_NAFVariable.XH;
      ZAC = expcomplexe (muldoublcomplexe(-ANG0,ZI));
      ZINC=  expcomplexe (muldoublcomplexe(-ANGI,ZI));
      ZEX = divcomplexe(ZAC,ZINC); /*ZEX = ZAC/ZINC*/
      naf_ztpow2(g_NAFVariable.KTABS,NVECT,ZTF,g_NAFVariable.ZTABS,TWIN,ZINC,ZEX); /*CALL  ZTPOW2(g_NAFVariable.KTABS+1,NVECT,ZTF,g_NAFVariable.ZTABS,TWIN,ZINC,ZEX)*/ 
/*!------------------ TAILLE DU PAS*/
      H=1.E0/((double)LTF);
      naf_zardyd(ZTF,LTF,H,&ZA);
      *A = ZA.reel;
      *B = ZA.imag;
      *RM=module(ZA);
      naf_ztder(g_NAFVariable.KTABS,NVECT,ZTF,g_NAFVariable.ZTABS,TWIN,ZINC,ZEX,g_NAFVariable.T0,g_NAFVariable.XH); /*CALL ZTDER(g_NAFVariable.KTABS+1,NVECT,ZTF,g_NAFVariable.ZTABS,TWIN,ZINC,ZEX,g_NAFVariable.T0,g_NAFVariable.XH)*/
      naf_zardyd(ZTF,LTF,H,&ZB);
      *DER=(mulcomplexe(conjcomplexe(ZA),ZB)).imag*2.E0;
#else /*remplacee par: */
/*!
!------------------ CONVERSION DE FS EN RD/AN*/
      OM=FS/g_NAFVariable.UNIANG;
      i_compl_cmplx(&ZI,0.E0,1.E0);
/*!------------------ LONGUEUR DU TABLEAU A INTEGRER
!------------------ DANS LE CAS REEL MEME CHOSE */
      LTF=g_NAFVariable.KTABS;
      ANG0=OM*g_NAFVariable.T0;
      ANGI=OM*g_NAFVariable.XH;
      ZAC = i_compl_exp (i_compl_muldoubl(-ANG0,ZI));
      ZINC=  i_compl_exp (i_compl_muldoubl(-ANGI,ZI));
      ZEX = i_compl_div(ZAC,ZINC); /*ZEX = ZAC/ZINC*/
      naf_ztpow2(g_NAFVariable.KTABS,NVECT,ZTF,g_NAFVariable.ZTABS,TWIN,ZINC,ZEX); /*CALL  ZTPOW2(g_NAFVariable.KTABS+1,NVECT,ZTF,g_NAFVariable.ZTABS,TWIN,ZINC,ZEX)*/ 
/*!------------------ TAILLE DU PAS*/
      H=1.E0/((double)LTF);
      naf_zardyd(ZTF,LTF,H,&ZA);
      *A = ZA.reel;
      *B = ZA.imag;
      *RM=i_compl_module(ZA);
      naf_ztder(g_NAFVariable.KTABS,NVECT,ZTF,g_NAFVariable.ZTABS,TWIN,ZINC,ZEX,g_NAFVariable.T0,g_NAFVariable.XH); /*CALL ZTDER(g_NAFVariable.KTABS+1,NVECT,ZTF,g_NAFVariable.ZTABS,TWIN,ZINC,ZEX,g_NAFVariable.T0,g_NAFVariable.XH)*/
      naf_zardyd(ZTF,LTF,H,&ZB);
      *DER=(i_compl_mul(i_compl_conj(&ZA),ZB)).imag*2.E0;
#endif /*NAF_USE_OPTIMIZE==0*/
/* v0.96 M. GASTINEAU 01/12/98 : fin optimisation */
      SYSFREE(ZTF);
#undef NVECT
}/*      END SUBROUTINE PRODER*/

/*      SUBROUTINE ZTDER (N,N1,ZTF,ZTA,TW,ZA,ZAST,T0,XH)*/
/*v0.96 M. GASTINEAU 09/09/98 : modification dans les boucles de I en I-1 */
void naf_ztder(int N, int N1, t_complexe *ZTF, t_complexe *ZTA, double *TW, t_complexe ZA, t_complexe ZAST, double T0, double XH)
{
/*!-----------------------------------------------------------------------
!ZTPOW   CALCULE  ZTF(I) = ZTA(I)* TW(I)*ZAST*ZA**I *(T0+(I-1)*XH EN VECTORIEL
!          
!-----------------------------------------------------------------------
      IMPLICIT NONE
!  (N,N1,ZT,ZTF,ZTA,TW,ZA,ZAST,T0,XH)
      integer :: N,N1
      complex (8) :: ZTF(0:N),ZTA(0:N),ZA,ZAST
      REAL (8) :: TW(0:N),T0,XH
!
      integer :: I,INC,NT,IT,NX
      complex (8) :: ZT1,ZINC
      complex (8),dimension (:), allocatable :: ZT 
!      */
      int I,INC,NT,IT,NX;
      t_complexe ZT1, ZINC;
      t_complexe *ZT=NULL;
      
      SYSCHECKMALLOCSIZE(ZT,t_complexe,N1); /*allocate(ZT(0:N1-1),stat = NERROR)*/
      if (N<N1-1)
      {
         fprintf(stdout,"DANS ZTDER, N = %d\n", N);
         return;
      }
/*!----------! */
/* v0.96 M. GASTINEAU 01/12/98 : utilisation des fonctions  complexes inlines et optimisation */
#if NAF_USE_OPTIMIZE==0
      ZT[0] = mulcomplexe(ZAST,ZA);
      for( I = 1; I<N1; I++)
      {
         ZT[I] = mulcomplexe(ZT[I-1],ZA);
      }
      for(I = 0; I<N1; I++)
      {
         /*v0.96M. GASTINEAU 09/09/98 : modification */
         /*ZTF[I] = muldoublcomplexe((T0+(I-1)*XH)*TW[I], mulcomplexe(ZTA[I],ZT[I]));*/ /*ZTF(I) = ZTA(I)*TW(I)*ZT(I)*(T0+(I-1)*XH)*/
         ZTF[I] = muldoublcomplexe((T0+(I)*XH)*TW[I], mulcomplexe(ZTA[I],ZT[I])); /*ZTF(I) = ZTA(I)*TW(I)*ZT(I)*(T0+(I-1)*XH)*/
      }
      ZT1 = divcomplexe(ZT[N1-1],ZAST);
      ZINC= cmplx(1,0);
      INC =0;
      NT = (N+1)/N1;  
      for( IT = 2; IT<= NT; IT++)
      {
         ZINC = mulcomplexe(ZINC,ZT1);
         INC += N1; /*INC  = INC + N1*/
         for(I = 0; I<N1; I++)
         {
           /*v0.96M. GASTINEAU 09/09/98 : modification */
           /*ZTF[INC+I]=muldoublcomplexe((T0+(INC+I-1)*XH)*TW[INC+I],mulcomplexe(ZTA[INC+I],mulcomplexe(ZT[I],ZINC)));*/
           ZTF[INC+I]=muldoublcomplexe((T0+(INC+I)*XH)*TW[INC+I],mulcomplexe(ZTA[INC+I],mulcomplexe(ZT[I],ZINC)));
           /*ZTF(INC+I)=ZTA(INC+I)*TW(INC+I)*ZT(I)*ZINC*(T0+(INC+I-1)*XH)*/
         }
      }
      ZINC = mulcomplexe(ZINC,ZT1);
      INC += N1; /*INC  = INC + N1*/
      NX = (N+1)-NT*N1;
      for (I = 0; I<NX; I++)
      {
        /*v0.96M. GASTINEAU 09/09/98 : modification */
        /*ZTF[INC+I]=muldoublcomplexe(TW[INC+I]*(T0+(INC+I-1)*XH), mulcomplexe(ZTA[INC+I], mulcomplexe(ZT[I],ZINC)));*/
        ZTF[INC+I]=muldoublcomplexe(TW[INC+I]*(T0+(INC+I)*XH), mulcomplexe(ZTA[INC+I], mulcomplexe(ZT[I],ZINC)));
        /*ZTF(INC+I)=ZTA(INC+I)*TW(INC+I)*ZT(I)*ZINC*(T0+(INC+I-1)*XH)*/
      }
#else /*remplacee par:*/
      ZT[0] = i_compl_mul(ZAST,ZA);
      for( I = 1; I<N1; I++)
      {
         ZT[I] = i_compl_mul(ZT[I-1],ZA);
      }
      for(I = 0; I<N1; I++)
      {
         ZTF[I] = i_compl_muldoubl((T0+(I)*XH)*TW[I], i_compl_mul(ZTA[I],ZT[I])); 
         /*ZTF(I) = ZTA(I)*TW(I)*ZT(I)*(T0+(I-1)*XH)*/
      }
      ZT1 = i_compl_div(ZT[N1-1],ZAST);
      i_compl_cmplx(&ZINC,1,0);
      INC =0;
      NT = (N+1)/N1;  
      for( IT = 2; IT<= NT; IT++)
      {
         i_compl_pmul(&ZINC,&ZT1);
         INC += N1; /*INC  = INC + N1*/
         for(I = 0; I<N1; I++)
         {
           ZTF[INC+I]=i_compl_muldoubl((T0+(INC+I)*XH)*TW[INC+I],i_compl_mul(ZTA[INC+I],i_compl_mul(ZT[I],ZINC)));
           /*ZTF(INC+I)=ZTA(INC+I)*TW(INC+I)*ZT(I)*ZINC*(T0+(INC+I-1)*XH)*/
         }
      }
      i_compl_pmul(&ZINC,&ZT1);
      INC += N1; /*INC  = INC + N1*/
      NX = (N+1)-NT*N1;
      for (I = 0; I<NX; I++)
      {
        ZTF[INC+I]=i_compl_muldoubl(TW[INC+I]*(T0+(INC+I)*XH), i_compl_mul(ZTA[INC+I], i_compl_mul(ZT[I],ZINC)));
        /*ZTF(INC+I)=ZTA(INC+I)*TW(INC+I)*ZT(I)*ZINC*(T0+(INC+I-1)*XH)*/
      }
#endif /*NAF_USE_OPTIMIZE==0*/
/* v0.96 M. GASTINEAU 01/12/98 : fin optimisation */
      SYSFREE(ZT);
}/*      END SUBROUTINE ZTDER*/
      
/*      SUBROUTINE SECANTES(X,PASS,EPS,XM,IPRT,NFPRT)*/
void naf_secantes(double X, double PASS, double EPS, double *XM, int IPRT, FILE *NFPRT)
{
/*!      SUBROUTINE SECANTES(X,PASS,EPS,XM,FONC,IPRT,NFPRT)
!***********************************************************
!      CALCUL PAR  LA METHODE DES SECANTES D'UN ZERO
!      D'UNE FONCTION FUNC FOURNIE PAR L'UTILISATEUR
!      
!     X, PASS : LE ZERO EST SUPPOSE ETRE DANS X-PASS,X+PASS
!     EPS     : PRECISION AVEC LAQUELLE ON VA CALCULER
!               LE MAXIMUM
!     XM      : ABSCISSE DU  ZERO
!    g_NAFVariable.IPRT     : -1 PAS D'IMPRESSION
!               0: IMPRESSION D'UNE EVENTUELLE ERREUR FATALE
!               1: IMPRESSION SUR LE FICHIER g_NAFVariable.NFPRT
!               2: IMPRESSION DES RESULTATS INTERMEDIAIRES
!     
!    FONC(X) :FONCTION DONNEE PAR L'UTILISATEUR
!    PARAMETRES
!       NMAX  : NOMBRE D'ESSAIS
!       IENC  : 1 : test d'encadrement 0 sinon
!   VARIABLE NERROR
!       0: A PRIORI TOUT S'EST BIEN PASSE
!       1: LE ZERO N'EST PAS TROUVE A EPSI PRES
!          -->   SEULEMENT SI IENC=1
!       2: PENTRE TROP FAIBLE (DIVISION PAR PRESQUE ZERO)
!       3: ECHEC TOTAL 
!  F. Joutel 26/5/98
!  Modif 26/8/98 Enleve l'appel de FONC car dans un module
! il ne semble pas etre possible de passer une fonction
! interne en parametre, remplace FONC par FUNCP
!***********************************************************
       
       IMPLICIT NONE
! (X,PASS,EPS,XM,FUNCP,IPRT,NFPRT)
       real (8) :: X, PASS,EPS,XM
       integer :: IPRT,NFPRT
!       interface
!         function fonc(x)
! 	   real (8) :: x,fonc
! 	 end function fonc
!       end interface
!              
       integer, parameter :: NMAX=30, IENC=1
       integer :: I
       REAL (8) :: EPSI,A,B,FA,FB,DELTA,COR
!*/
#define SECANTES_NMAX (30)
#define SECANTES_IENC (1)
       int I;
       double EPSI,A,B,FA,FB,DELTA,COR;
       
       g_NAFVariable.NERROR=0;
       EPSI=MAX(g_NAFVariable.EPSM,EPS);
       if (fabs(PASS)>EPSI)
       {
         if (IPRT>=1)
         {
           fprintf(NFPRT," AMELIORATION PAR LES SECANTES\n");
         }
       }
       I=0;
       A=X-PASS;
       B=X;
       FB=naf_funcp(A);
/*!      FB=FONC(A)*/
/*10    CONTINUE*/
SECANTES_10 :
       if (fabs(B-A)>EPSI)
       { 
          FA=FB;
          FB=naf_funcp(B);
/*! 	  FB=FONC(B)*/
          DELTA= FB-FA;
          if (IPRT>=2)
          {
            fprintf(NFPRT,"SEC: A=%g, B=%g, abs(B-A)=%g \n", A, B,fabs(B-A));
            fprintf(NFPRT,"SEC: F(A)=%g, F(B)=%g\n", FA, FB);
          }
          if (fabs(DELTA)<=g_NAFVariable.EPSM)
          {
            g_NAFVariable.NERROR=2;
            if (IPRT>=1)
            {
              fprintf(NFPRT,"ECHEC DE LA METHODE DES SECANTES\n");
              fprintf(NFPRT," DIVISION PAR PRESQUEZERO\n");
              fprintf(NFPRT," ON CONTINUE AVEC LA VALEUR TROUVEE\n");
            }
            *XM=B;
            return;
          }
          COR = FB*(A-B)/DELTA;
          A = B;
          B = B+COR;
          I++; /*I=I+1 */
          if (I>SECANTES_NMAX)
          {
             g_NAFVariable.NERROR=3;
             if (IPRT>=0)
             {
               fprintf(NFPRT,"ECHEC DE LA METHODE DES SECANTES\n");
               fprintf(NFPRT," BEAUCOUP TROP D ITERATIONS\n");
               fprintf(NFPRT," ON CONTINUE AVEC LA VALEUR INITIALE\n");
             }
             *XM=X;
             return;
          }
          goto SECANTES_10;  /*GOTO 10*/
       }
/*!----- !fin de la boucle*/
       *XM=B;
       if (SECANTES_IENC==1)
       {
            if (naf_funcp(B-EPSI)*naf_funcp(B+EPSI)>0.E0)
            {
/*!           if (FONC(B-EPSI)*FONC(B+EPSI)>0.D0) {*/
             g_NAFVariable.NERROR=1;
           }
       }
       if (IPRT==1)
       {
          /*v0.96 M. GASTINEAU 19/11/98 : modification*/
          /*fprintf(NFPRT,"POSITION SECANTES %g TROUVEE A %g",XM, EPSI);*/
          /*remplacee par:*/
          fprintf(NFPRT,"POSITION SECANTES %g TROUVEE A %g",*XM, EPSI);
          /*v0.96 M. GASTINEAU 19/11/98 : fin modification*/
          if (SECANTES_IENC==1)
          {
             if (g_NAFVariable.NERROR==1)
             {
                fprintf(NFPRT," SANS GARANTIE\n");
             }
             if (g_NAFVariable.NERROR==0)
             {
                fprintf(NFPRT," AVEC GARANTIE\n");
             }
           }
          /*v0.96 M. GASTINEAU 19/11/98 : modification*/
          /* fprintf(NFPRT,"\n %19.9E %12.8E\n", PASS,XM);*/
          /*remplacee par:*/
           fprintf(NFPRT,"\n %19.9E %12.8E\n", PASS,*XM);
          /*v0.96 M. GASTINEAU 19/11/98 : fin modification*/
       }
#undef SECANTES_NMAX
#undef SECANTES_IENC
}/*       END SUBROUTINE SECANTES*/

/*      SUBROUTINE MAXIQUA(X,PASS,EPS,XM,YM,IPRT,NFPRT)*/
void naf_maxiqua(double X, double PASS, double EPS, double *XM, double *YM, int IPRT, FILE *NFPRT)
{
/*
!      SUBROUTINE MAXIQUA(X,PASS,EPS,XM,YM,FONC,IPRT,NFPRT)
!***********************************************************
!      CALCUL PAR INTERPOLATION QUADRATIQUE DU MAX
!      D'UNE FONCTION FONC FOURNIE PAR L'UTILISATEUR
!      
!     X, PASS : LE MAX EST SUPPOSE ETRE SI POSSIBLE 
!               DANS X-PASS,X+PASS
!     EPS     : PRECISION AVEC LAQUELLE ON VA CALCULER
!               LE MAXIMUM
!     XM      : ABSCISSE DU MAX
!     YM      : YM=FUNC(XM)
!     IPRT     : -1 PAS D'IMPRESSION 
!               0: IMPRESSION D'UNE EVENTUELLE ERREUR FATALE
!               1: IMPRESSION SUR LE FICHIER NFPRT
!               2: IMPRESSION DES RESULTATS INTERMEDIAIRES
! 
!    FONC(X)  :FONCTION DONNEE PAR L'UTILISATEUR
! 
!***********************************************************
!  PARAMETRES
!   NMAX     : NOMBRE MAXIMAL D'ITERATIONS
!  VARIABLE (COMMON ASD)
!   g_NAFVariable.NERROR:
!        0 : OK
!        1 : COURBURE BIEN FAIBLE !
!        2 : LA FONCTION EST BIEN PLATE !
!        3 : OSCILLATION  DE L'ENCADREMENT
!        4 : PAS D'ENCADREMENT TROUVE
!        5 : ERREUR FATALE
!***********************************************************
!  ASD VERSION 0.9  J.LASKAR & F.JOUTEL 16/5/97
!  CHANGEMENT POUR UN EPSILON MACHINE CALCULE PAR CALCEPSM
!  JXL 22/4/98
!  Modif 26/8/98 Enleve l'appel de FONC car dans un module
! il ne semble pas etre possible de passer une fonction
! interne en parametre, remplace FONC par FUNCP
!***********************************************************
      
      IMPLICIT NONE
!  (X,PASS,EPS,XM,YM,FONC,IPRT,NFPRT)
      integer :: IPRT,NFPRT
      REAL (8) :: X,PASS,EPS,XM,YM
!      interface
!        function fonc (x)
! 	  real (8) :: x,fonc
! 	end function fonc
!      end interface
!          
      integer, parameter :: NMAX=200, NFATAL=100 
      REAL (8), parameter :: RABS=1.D0
      integer :: NITER,NCALCUL
      REAL (8)  :: PAS,EPSLOC,X1,X2,X3,Y1,Y2,Y3,A,ERR,DX,TEMP
      REAL (8) :: D1,D2   
*/
#define MAXIQUA_NMAX (200)
#define MAXIQUA_NFATAL (100)
#define MAXIQUA_RABS (1.E0) 
      int NITER, NCALCUL;
      double PAS,EPSLOC,X1,X2,X3,Y1,Y2,Y3,A,ERR,DX,TEMP;
      double D1=0E0,D2=0E0;
     
/*! 
!* ENTRE 0 ET RABS LES ERREURS SONT ENVISAGEES DE FACON ABSOLUE
!* APRES RABS, LES ERREURS SONT ENVISAGEES DE FACON RELATIVE.
!*/
      EPSLOC=MAX(EPS,sqrt(g_NAFVariable.EPSM));
/*!      
!*  AU SOMMET D'UNE PARABOLE Y=Y0-A*(X-X0)**2, UN ECART DE X AUTOUR DE X0
!*  INFERIEUR A EPSLOC DONNE UN ECART DE Y INFERIEUR A A*EPSLOC**2*/
      if (IPRT>=1)
      {
          fprintf(NFPRT," ROUTINE DE RECHERCHE DU MAXIMUM :\n");
      }
      PAS = PASS;
      NITER =0;
      NCALCUL=0;
      g_NAFVariable.NERROR=0;
      X2 = X;
      Y2 = naf_func (X2);
      X1 = X2 - PAS;
      X3 = X2 + PAS;
      Y1 = naf_func (X1);
      Y3 = naf_func (X3);
/*!* DECALAGE POUR ENCADRER LE MAXIMUM*/
MAXIQUA_10:
      if ((NITER>MAXIQUA_NMAX)||(NCALCUL>MAXIQUA_NFATAL))
      {
         if (NCALCUL>MAXIQUA_NFATAL)
         {
            g_NAFVariable.NERROR =5;
            if (IPRT>=0)
            {
               fprintf(NFPRT," ERREUR FATALE\n");
            }
         } 
         else
         { 
           if (PAS>=PASS)
           {
              g_NAFVariable.NERROR = 4;
              if (IPRT>=0)
              {
                fprintf(NFPRT,"  PAS D''ENCADREMENT TROUVE\n");
              }
           } 
           else 
           {
              g_NAFVariable.NERROR =3;
              if (IPRT>=0)
              {
                fprintf(NFPRT," OSCILLATION DE L''ENCADREMENT ?\n");
              }
           }
         }
         ERR=PAS;
         goto MAXIQUA_EXIT;
      }
      if ((Y1>Y2)||(Y3>Y2))
      {
            NITER ++; /*NITER = NITER +1*/
            if (Y1>Y3)
            {
                  X2 = X1;
                  Y2 = Y1;
            }
            else
            {
                  X2 = X3;
                  Y2 = Y3;
            }
            X1 = X2 - PAS;
            X3 = X2 + PAS;
            Y1 = naf_func (X1);
            Y3 = naf_func (X3);
            goto MAXIQUA_10;
      }
      D1 = Y2-Y1;
      D2 = Y3-Y2;
      A = fabs(0.5E0*(D2-D1)/(PAS*PAS));
      if (IPRT>=2)
      {
        if (sqrt(A)!=0)
        { 
         fprintf(NFPRT,"%12.8E %12.8E %25.16E %25.16E %25.16E %25.16E\n", PAS,X2,Y1,Y2,Y3,EPSLOC/sqrt(A));
        }
        else
        { 
         fprintf(NFPRT,"%12.8E %12.8E %25.16E %25.16E %25.16E INF\n", PAS,X2,Y1,Y2,Y3);
        }        
      }
/*!* TEST POUR UN CALCUL SIGNifICATif DES VALEURS DE LA FONCTION*/
      if ((fabs (D1-D2))<2.E0*g_NAFVariable.EPSM*MAX(MAXIQUA_RABS,fabs(Y2)))
      {
        g_NAFVariable.NERROR= 2;
        if (IPRT>=1)
        {
          fprintf(NFPRT," PLATEAU DE LA FONCTION\n");
        }
        ERR=PAS;
        goto MAXIQUA_EXIT;
      }

/*!* TEST SUR LA FORME DE LA COURBE */
      if (A<g_NAFVariable.EPSM*MAX(MAXIQUA_RABS,fabs(Y2)))
      {
         g_NAFVariable.NERROR= 1;
         if (IPRT>=1)
         {
          fprintf(NFPRT,"COURBE TROP PLATE\n");
         }
         ERR=PAS;
        goto MAXIQUA_EXIT;
      }
/*!* TEST POUR UN CALCUL SIGNifICATif AU SOMMET D'UNE PARABOLE
!* ON PEUT IGNORER LE TEST MAIS ON N'A PLUS DE GARANTIE SUR
!* LES VALEURS DE LA FONCTION PLUS TARD */
      if (PAS<=EPSLOC*sqrt(MAX(MAXIQUA_RABS,fabs(Y2))/A))
      {
          if (IPRT>=1)
          {
             fprintf(NFPRT,"SORTIE AU SOMMET DE LA PARABOLE\n"); 
          }
          ERR=EPSLOC*sqrt(MAX(MAXIQUA_RABS,Y2)/A);
          goto MAXIQUA_EXIT;
        } 
      DX= 0.5E0*PAS*(Y1-Y3)/(D2-D1);
/*!* TEST POUR UN NOUVEAU PAS SIGNifICATif*/
      if (fabs(DX)<g_NAFVariable.EPSM*MAX(MAXIQUA_RABS,fabs(X2)))
      {
          if (IPRT>=1)
          {
           fprintf(NFPRT,"CORRECTION MINUSCULE\n");
          } 
          ERR=PAS;
          goto MAXIQUA_EXIT;
      }
      PAS=fabs(DX); 
      if (DX>0.E0)
      {
         X1 = X2;
         Y1 = Y2;
         X2 +=PAS; /*X2 = X2+PAS*/
         Y2 = naf_func(X2);
         X3 = X2 + PAS;
         Y3 = naf_func(X3);
      }
      else
      {
         X3 = X2;
         Y3 = Y2;
         X2 -=PAS; /*X2 = X2-PAS*/
         Y2 = naf_func(X2);
         X1 = X2 - PAS;
         Y1 = naf_func(X1);
      }
/*!* A LA SORTIE DE CE CALCUL LE PAS DOIT ETRE DIVISE AU MOINS PAR 2 
!! MAIS ON N'EST PAS ASSURE QUE LES VALEURS FINALES Y1,Y2,Y3 
!* VERifIENT Y1<=Y2>=Y3, AUSSI, ON REPART AU DEBUT*/
      NCALCUL++; /*NCALCUL=NCALCUL+1*/
      goto MAXIQUA_10;
MAXIQUA_EXIT :
      *XM = X2;
      *YM = Y2;
      if (IPRT>=1)
      {
         fprintf(NFPRT,"%12.8E %12.8E %12.8E\n", PAS,*XM,*YM);
         fprintf(NFPRT," POSITION TROUVEE A %g PRES\n",ERR); 
         TEMP =  fabs (Y3-Y2)+ fabs(Y2-Y1);
         if (TEMP<2.E0*g_NAFVariable.EPSM*MAX(MAXIQUA_RABS,fabs(Y2)))
         {
            fprintf(NFPRT," TROUVEE SOUS UN PLATEAU DE LA FONCTION\n");
         } 
      }
#undef MAXIQUA_NMAX 
#undef MAXIQUA_NFATAL 
#undef MAXIQUA_RABS    
}/*      END SUBROUTINE MAXIQUA*/

/*      SUBROUTINE FREFIN (FR,A,B,RM, RPAS0,RPREC)*/
/*v0.96 M. GASTINEAU 07/12/98 : modification car bug lors de l'optimisation sur Mac */
void naf_frefin(double *FR, double *A, double *B, double *RM, const double RPAS0, const double RPREC)
{
/*!-----------------------------------------------------------------------
!     FREFIN            RECHERCHE FINE DE LA FREQUENCE
!                       FR     FREQUENCE DE BASE EN "/AN  ENTREE ET SORT
!                       A      PARTIE REELLE DE L'AMPLITUDE     (SORTIE)
!                       B      PARTIE IMAGINAIRE DE L'AMPLITUDE (SORTIE)
!                       RM     MODULE DE L'AMPLITUDE            (SORTIE)
!                       RPAS0  PAS INITIAL EN "/AN
!                       RPREC  PRECISION DEMANDEE EN "/AN
!               REVU LE 26/9/87 J. LASKAR
!               MODif LE 15/11/90 MAX PAR INTERP. QUAD.
! Modif 26/8/98 Enleve les appels aux fonctions func et funcp internes
! au module
!-----------------------------------------------------------------------
      
      IMPLICIT NONE
! (FR,A,B,RM, RPAS0,RPREC)
      REAL (8) FR,A,B,RM,RPAS0,RPREC
!
      REAL (8) X,PASS,EPS,XM,YM
! */
      double X,PASS,EPS,XM,YM;
   /*v0.96 M. GASTINEAU 07/12/98 : modification car bug quand optimisation sur Mac */
   /*   X    = *FR;
      PASS = RPAS0;
      EPS  = RPREC;*/ /*remplacee par:*/
      PASS = RPAS0;
      EPS  = RPREC;
      X    = *FR;
   /*v0.96 M. GASTINEAU 07/12/98 : fin modification */
      naf_maxiqua(X,PASS,EPS,&XM,&YM,g_NAFVariable.IPRT,g_NAFVariable.NFPRT);
/*!     CALL MAXIQUA(X,PASS,EPS,XM,YM,FUNC,g_NAFVariable.IPRT,g_NAFVariable.NFPRT)
!************! AMELIORATION DE LA RECHERCHE PAR LE METHODE DES SECANTES
!************!  APPLIQUEE A LA DERIVEE DU MODULE.*/
      if (g_NAFVariable.ISEC==1)
      {
         X=XM;
         naf_secantes(X,PASS,EPS,&XM,g_NAFVariable.IPRT,g_NAFVariable.NFPRT);
/*! 	 CALL SECANTES(X,PASS,EPS,XM,FUNCP,g_NAFVariable.IPRT,g_NAFVariable.NFPRT)*/
         YM=naf_func(XM);
      }   
      *FR=XM;
      *RM=YM;
      *A=AF;
      *B=BF; 
      if (g_NAFVariable.IPRT==1)
      {
       /*v0.96 M. GASTINEAU 19/11/98 : modification*/
       /*fprintf(g_NAFVariable.NFPRT,"\n%10.6E %19.9E %19.9E %19.9E\n", FR,RM,AF,BF);*/
       /*remplacee par:*/
       fprintf(g_NAFVariable.NFPRT,"\n%10.6E %19.9E %19.9E %19.9E\n", *FR,*RM,AF,BF);
       /*v0.96 M. GASTINEAU 19/11/98 : fin modification*/
      }
/*1000  FORMAT (1X,F10.6,3D19.9)*/
}/*      END SUBROUTINE FREFIN*/



/*      SUBROUTINE PROFRE(FS,A,B,RMD)*/
BOOL naf_profre(double FS, double *A, double *B, double *RMD)
{
/*!***********************************************************************
!  CE SOUS PROG. CALCULE LE PRODUIT SCALAIRE DE EXP(-I*FS)*T PAR TAB(0:7
!  UTILISE LA METHODE D'INTEGRATION DE HARDY
!  LES RESULTATS SONT DANS A,B,RMD PARTIE REELLE, IMAGINAIRE, MODULE
!  FS EST DONNEE EN SECONDES PAR AN
!               REVU LE 26/9/87 J. LASKAR
!***********************************************************************
      
      IMPLICIT NONE
! (g_NAFVariable.KTABS,FS,A,B,RMD)
      REAL (8) :: FS,A,B,RMD
!
      integer :: LTF
      REAL (8) :: OM,ANG0,ANGI,H
      complex (8) ZI,ZAC,ZINC,ZEX,ZA
      complex (8), dimension (:), allocatable :: ZTF      
! */
      int LTF;
      double OM,ANG0,ANGI,H;
      t_complexe ZI,ZAC,ZINC,ZEX,ZA;
      t_complexe *ZTF=NULL;

      SYSCHECKMALLOCSIZE(ZTF, t_complexe, g_NAFVariable.KTABS+1); /*allocate(ZTF(0:g_NAFVariable.KTABS),stat = g_NAFVariable.NERROR)*/
/*!
!------------------ CONVERSION DE FS EN RD/AN*/
      OM=FS/g_NAFVariable.UNIANG ;
/* v0.96 M. GASTINEAU 01/12/98 : utilisation des fonctions  complexes inlines et optimisation */
#if NAF_USE_OPTIMIZE==0
      ZI=cmplx(0.E0,1.E0);
#else /*remplacee par:*/
      i_compl_cmplx(&ZI,0.E0,1.E0);
#endif /*NAF_USE_OPTIMIZE==0*/
/* v0.96 M. GASTINEAU 01/12/98 : fin optimisation */
/*!------------------ LONGUEUR DU TABLEAU A INTEGRER
!------------------ DANS LE CAS REEL MEME CHOSE */
      LTF=g_NAFVariable.KTABS;
      ANG0=OM*g_NAFVariable.T0;
      ANGI=OM*g_NAFVariable.XH;
/* v0.96 M. GASTINEAU 01/12/98 : utilisation des fonctions  complexes inlines et optimisation */
#if NAF_USE_OPTIMIZE==0
      ZAC = expcomplexe (muldoublcomplexe(-ANG0,ZI));
      ZINC= expcomplexe (muldoublcomplexe(-ANGI,ZI));
      ZEX = divcomplexe(ZAC,ZINC);
#else /*remplacee par:*/
      ZAC = i_compl_exp(i_compl_muldoubl(-ANG0,ZI));
      ZINC= i_compl_exp(i_compl_muldoubl(-ANGI,ZI));
      ZEX = i_compl_div(ZAC,ZINC);
#endif /*NAF_USE_OPTIMIZE==0*/
/* v0.96 M. GASTINEAU 01/12/98 : fin optimisation */
      naf_ztpow2(g_NAFVariable.KTABS,64,ZTF,g_NAFVariable.ZTABS,TWIN,ZINC,ZEX); /*CALL  ZTPOW2(g_NAFVariable.KTABS+1,64,ZTF,g_NAFVariable.ZTABS,TWIN,ZINC,ZEX)*/
/*!------------------ TAILLE DU PAS*/
      H=1.E0/((double)LTF);
      if (naf_zardyd(ZTF,LTF,H,&ZA)==FALSE)
      {
       SYSFREE(ZTF);
       return FALSE;
      }
/* v0.96 M. GASTINEAU 01/12/98 : utilisation des fonctions  complexes inlines et optimisation */
#if NAF_USE_OPTIMIZE==0
      *RMD=module(ZA);
#else /*remplacee par:*/
      *RMD=i_compl_module(ZA);
#endif /*NAF_USE_OPTIMIZE==0*/
/* v0.96 M. GASTINEAU 01/12/98 : fin optimisation */
      *A = ZA.reel;
      *B = ZA.imag;
      SYSFREE(ZTF);
      return TRUE;
}/*      END SUBROUTINE PROFRE*/

/*      SUBROUTINE ZTPOW2 (N,N1,ZTF,ZTA,TW,ZA,ZAST)*/
void naf_ztpow2(int N, int N1, t_complexe *ZTF, t_complexe *ZTA, double *TW, t_complexe ZA, t_complexe ZAST)
{
/*!-----------------------------------------------------------------------
!     ZTPOW   CALCULE  ZTF(I) = ZTA(I)* TW(I)*ZAST*ZA**I EN VECTORIEL
!           
!-----------------------------------------------------------------------
      IMPLICIT NONE
! (N,N1,ZTF,ZTA,TW,ZA,ZAST)
      integer :: N,N1
      complex (8) :: ZTF(0:N),ZTA(0:N),ZA,ZAST
      REAL (8) :: TW(0:N)
!
      integer :: INC,NT,IT,I,NX
      complex (8) :: ZT1,ZINC
      complex (8), dimension(:),allocatable :: ZT
!               */
      int INC,NT,IT,I,NX;
      t_complexe ZT1,ZINC;
      t_complexe *ZT=NULL;
/* v0.96 M. GASTINEAU 01/12/98 : utilisation des fonctions  complexes inlines et optimisation */
#if NAF_USE_OPTIMIZE>0
      t_complexe *pzarZT;
#endif /**/      
      SYSCHECKMALLOCSIZE(ZT,t_complexe, N1); /*allocate(ZT(0:N1-1),stat = nerror)*/
      if (N<N1-1)
      {
         fprintf(stdout,"DANS ZTPOW, N = %d\n", N);
         return;
      }
/*!----------! */
/* v0.96 M. GASTINEAU 01/12/98 : utilisation des fonctions  complexes inlines et optimisation */
#if NAF_USE_OPTIMIZE==0
      ZT[0] = mulcomplexe(ZAST,ZA);
      for(I = 1; I<N1; I++)
      {
         ZT[I]=mulcomplexe(ZT[I-1], ZA); /*ZT(I) = ZT(I-1)*ZA*/
      }
      for(I = 0; I<N1; I++)
      {
         ZTF[I]=muldoublcomplexe(TW[I], mulcomplexe(ZTA[I], ZT[I])); /*ZTF(I) = ZTA(I)*TW(I)*ZT(I)*/
      }
      ZT1= divcomplexe(ZT[N1-1], ZAST); /*ZT1 = ZT(N1-1)/ZAST*/
      ZINC= cmplx(1E0,0E0);
      INC =0;
      NT = (N+1)/N1;  
      for(IT = 2;IT<= NT; IT++)
      {
         ZINC=mulcomplexe(ZINC,ZT1); /*ZINC = ZINC*ZT1*/
         INC += N1; /*INC  = INC + N1*/
         for(I = 0; I< N1;I++)
         {
            /*ZTF(INC +I) = ZTA(INC+I)*TW(INC+I)*ZT(I)*ZINC*/
            ZTF[INC +I] = muldoublcomplexe(TW[INC+I], mulcomplexe(mulcomplexe(ZTA[INC+I], ZT[I]),ZINC));
         }
      }
      ZINC = mulcomplexe(ZINC, ZT1); /*ZINC = ZINC*ZT1*/
      INC += N1; /*INC  = INC + N1*/
      NX = N+1-NT*N1;
      for(I = 0; I< NX; I++)
      {
          /*  ZTF(INC +I) = ZTA(INC+I)*TW(INC+I)*ZT(I)*ZINC*/
       ZTF[INC +I] = muldoublcomplexe(TW[INC+I], mulcomplexe(mulcomplexe(ZTA[INC+I], ZT[I]),ZINC));
      }
#else /*remplacee par:*/
      ZT[0] = i_compl_mul(ZAST,ZA);
      for(I = 1, pzarZT=ZT+1; I<N1; I++, pzarZT++)
      {
         *pzarZT=i_compl_mul(*(pzarZT-1), ZA); /*ZT(I) = ZT(I-1)*ZA*/
      }
      for(I = 0; I<N1; I++)
      {
         ZTF[I]=i_compl_muldoubl(TW[I], i_compl_mul(ZTA[I], ZT[I])); /*ZTF(I) = ZTA(I)*TW(I)*ZT(I)*/
      }
      ZT1= i_compl_div(ZT[N1-1], ZAST); /*ZT1 = ZT(N1-1)/ZAST*/
      i_compl_cmplx(&ZINC,1E0,0E0);
      INC =0;
      NT = (N+1)/N1;  
      for(IT = 2;IT<= NT; IT++)
      {
         i_compl_pmul(&ZINC,&ZT1); /*ZINC = ZINC*ZT1*/
         INC += N1; /*INC  = INC + N1*/
         for(I = 0; I< N1;I++)
         {
            /*ZTF(INC +I) = ZTA(INC+I)*TW(INC+I)*ZT(I)*ZINC*/
            ZTF[INC +I] = i_compl_muldoubl(TW[INC+I], i_compl_mul(i_compl_mul(ZTA[INC+I], ZT[I]),ZINC));
         }
      }
      i_compl_pmul(&ZINC,&ZT1); /*ZINC = ZINC*ZT1*/
      INC += N1; /*INC  = INC + N1*/
      NX = N+1-NT*N1;
      for(I = 0; I< NX; I++)
      {
          /*  ZTF(INC +I) = ZTA(INC+I)*TW(INC+I)*ZT(I)*ZINC*/
       ZTF[INC +I] = i_compl_muldoubl(TW[INC+I], i_compl_mul(i_compl_mul(ZTA[INC+I], ZT[I]),ZINC));
      }
#endif /*NAF_USE_OPTIMIZE==0*/
/* v0.96 M. GASTINEAU 01/12/98 : fin optimisation */
      SYSFREE(ZT);
}/*      END SUBROUTINE ZTPOW2*/

/*      SUBROUTINE INIWIN*/
/*v0.96 M. GASTINEAU 18/12/98 : modification du prototype */
/*void naf_iniwin()*//*remplacee par: */
void naf_iniwin(double *p_pardTWIN)
/*v0.96 M. GASTINEAU 18/12/98 : fin modification */
{
/*!******************************************************************
!     INITIALISE LE TABLEAU TWIN POUR UTILISATION D'UNE FENETRE
!     DANS L'ANALYSE DE FOURIER.
!     FENETRE DE HANNING:
!      PHI(T) = (1+COS(PI*T))
!     MODif LE 13/12/90 (J. LASKAR)
!     MODif LE 27/1/96 FOR VARIOUS WINDOWS  (J. LASKAR)
!
!     WORKING AREA TWIN(*) SHOULD BE GIVEN IN 
!
!     IW IS THE WINDOW FLAG
!     IW = 0 : NO WINDOW
!     IW = N > 0   PHI(T) = CN*(1+COS(PI*T))**N
!                  WITH CN = 2^N*(N!)^2/(2N)!
!     IW = -1  EXPONENTIAL WINDOW PHI(T) = 1/CE*EXP(-1/(1-T^2))
!
!     MODif LE 22/4/98 POUR CALCUL EN PLUS DE L'EPSILON MACHINE 
!     g_NAFVariable.EPSM 
!      26/5/98 CORRECTION DE L'ORIGINE DES TABLEAUX (*m/4)
!******************************************************************
!
       
      IMPLICIT NONE
!
      integer :: I,IT           
      REAL (8) :: CE,T1,T2,TM,T,CN,PIST
!
!---------------------------------------------------------------- */
      int I, IT;
      double CE,T1,T2,TM,T,CN,PIST;
      
      CE= 0.22199690808403971891E0;
/*!
!      PI=ATAN2(1.D0,0.E0)*2.E0*/
      T1=g_NAFVariable.T0;
      T2=g_NAFVariable.T0+g_NAFVariable.KTABS*g_NAFVariable.XH;
      TM=(T2-T1)/2;
      PIST=M_PI/TM;
      if (g_NAFVariable.IW==0)
      {
         for(IT=0;IT<=g_NAFVariable.KTABS;IT++)
         {
           /*v0.96 M. GASTINEAU 18/12/98 : modification */
           /* TWIN[IT]=1.E0; *//*remplacee par:*/
           p_pardTWIN[IT]=1.E0;           
           /*v0.96 M. GASTINEAU 18/12/98 : fin modification */
         }
      }
      else if(g_NAFVariable.IW>=0)
       {
         CN = 1.E0;
         for(I = 1;I<=g_NAFVariable.IW;I++)
         {
            CN *= I*(2.E0/((double)(g_NAFVariable.IW+I))); /*CN = CN*2.D0*I*1.D0/(g_NAFVariable.IW+I)*/
         };       
         for(IT=0;IT<=g_NAFVariable.KTABS;IT++)
         {
            T=IT*g_NAFVariable.XH-TM;
           /*v0.96 M. GASTINEAU 18/12/98 : modification */
            /*TWIN[IT]=CN*pow((1.E0+cos(T*PIST)),g_NAFVariable.IW);*/
            /*remplacee par:*/
            p_pardTWIN[IT]=CN*pow((1.E0+cos(T*PIST)),g_NAFVariable.IW);
           /*v0.96 M. GASTINEAU 18/12/98 : fin modification */
         }
      }
      else if(g_NAFVariable.IW==-1)
       {
         /*v0.96 M. GASTINEAU 18/12/98 : modification */
        /* TWIN[0] =0.E0;
         TWIN[g_NAFVariable.KTABS] =0.E0;*/
         /*remplacee par:*/
         p_pardTWIN[0] =0.E0;
         p_pardTWIN[g_NAFVariable.KTABS] =0.E0;
         /*v0.96 M. GASTINEAU 18/12/98 : fin modification */
         for(IT=1; IT<=g_NAFVariable.KTABS-1; IT++)
         {
            T=(IT*g_NAFVariable.XH-TM)/TM;
            /*v0.96 M. GASTINEAU 18/12/98 : modification */
            /* TWIN[IT]= exp(-1.E0/(1.E0-T*T))/CE;*/
            /*remplacee par:*/
            p_pardTWIN[IT]= exp(-1.E0/(1.E0-T*T))/CE;
            /*v0.96 M. GASTINEAU 18/12/98 : fin modification */
         }
      }
      if (g_NAFVariable.IPRT==1)
      {
/*!----------------   IMPRESSION TEMOIN*/
         for(IT=0; IT<=g_NAFVariable.KTABS; (g_NAFVariable.KTABS>20?IT+=g_NAFVariable.KTABS/20:IT++))
         {
            /*v0.96 M. GASTINEAU 18/12/98 : modification */
            /*fprintf(g_NAFVariable.NFPRT,"%20.3E\n", TWIN[IT]*1.E6);*/
            /*remplacee par:*/
            fprintf(g_NAFVariable.NFPRT,"%20.3E\n", p_pardTWIN[IT]*1.E6);
            /*v0.96 M. GASTINEAU 18/12/98 : fin modification */
         }
      }
}/*      END SUBROUTINE INIWIN*/

/*      SUBROUTINE ZARDYD(ZT,N,H,ZOM)*/
BOOL naf_zardyd(t_complexe *ZT, int N, double H, t_complexe *ZOM)
{
/*!***********************************************************************
!     CALCULE L'INTEGRALE D'UNE FONCTION TATULEE PAR LA METHODE DE HARDY
!      T(N) TABLEAU DOUBLE PRECISION DES VALEURS DE LA FONCTION
!      N = 6*K  ENTIER
!      H PAS ENTRE DEUX VALEURS
!      SOM VALEUR DE L'INTEGRALE SUR L'INTERVALLE [X1,XN]
!      LE PROGRAMME EST EN DOUBLE PRECISION
!               REVU LE 26/9/87 J. LASKAR
!      
!************************************************************************
      IMPLICIT NONE
!  (ZT,N,H,ZOM)
      integer :: N
      complex (8) :: ZT(0:N),ZOM
      REAL (8) :: H
!
      integer :: ITEST,K,I 
!               */
      int ITEST, K, I;
/* v0.96 M. GASTINEAU 13/01/99 : optimisation */
#if NAF_USE_OPTIMIZE==0
      int INC;
#else /*NAF_USE_OPTIMIZE*/
      double zomreel, zomimag;
      double zomreeltemp, zomimagtemp;
      double *pdZT=(double*)(ZT);
#endif /*NAF_USE_OPTIMIZE*/
/* v0.96 M. GASTINEAU 13/01/99 : fin optimisation */

      ITEST=N%6; /*ITEST=MOD(N,6);*/
      if (ITEST!=0)
      {
         Myyerror("naf_zradyd - N N'EST PAS UN MULTIPLE DE 6\n");
         return FALSE;
      }
      K=N/6;
/* v0.96 M. GASTINEAU 01/12/98 : utilisation des fonctions  complexes inlines et optimisation */
#if NAF_USE_OPTIMIZE==0
      /*ZOM=41*ZT(1)+216*ZT(2)+27*ZT(3)+272*ZT(4)
     $          +27*ZT(5)+216*ZT(6)+41*ZT(N)*/
      *ZOM = addcomplexe(
             addcomplexe(
                addcomplexe( muldoublcomplexe(41,ZT[0]),  muldoublcomplexe(216,ZT[1]) )
                , addcomplexe( muldoublcomplexe(27,ZT[2]),  muldoublcomplexe(272,ZT[3]) ) )
             , addcomplexe(              
                addcomplexe( muldoublcomplexe(27,ZT[4]),  muldoublcomplexe(216,ZT[5]) )
                , muldoublcomplexe(41,ZT[N]) )
              );
      INC=0;
      for(I=1; I<=K-1; I++)
      {
/*         ZOM=ZOM+82*ZT(6*I+1)+216*ZT(6*I+2)+27*ZT(6*I+3)+272*ZT(6*I+4)
     $           +27*ZT(6*I+5)+216*ZT(6*I+6)*/
       INC +=6;
       *ZOM = addcomplexe(
             addcomplexe(
                addcomplexe( muldoublcomplexe(82,ZT[INC]),  muldoublcomplexe(216,ZT[INC+1]) )
                , addcomplexe( muldoublcomplexe(27,ZT[INC+2]),  muldoublcomplexe(272,ZT[INC+3]) ) )
             , addcomplexe(              
                addcomplexe( muldoublcomplexe(27,ZT[INC+4]),  muldoublcomplexe(216,ZT[INC+5]) )
                , *ZOM )
              );
      }
      *ZOM= muldoublcomplexe(H*6.E0/840.E0, *ZOM); /*ZOM =ZOM*H*6.D0/840.D0*/
#else /*remplacee par:*/
     /* v0.96 M. GASTINEAU 13/01/99 : optimisation */
#if 0
      /*ZOM=41*ZT(1)+216*ZT(2)+27*ZT(3)+272*ZT(4)
     $          +27*ZT(5)+216*ZT(6)+41*ZT(N)*/
      *ZOM = i_compl_add(
             i_compl_add(
                i_compl_add( i_compl_muldoubl(41,ZT[0]),  i_compl_muldoubl(216,ZT[1]) )
                , i_compl_add( i_compl_muldoubl(27,ZT[2]),  i_compl_muldoubl(272,ZT[3]) ) )
             , i_compl_add(              
                i_compl_add( i_compl_muldoubl(27,ZT[4]),  i_compl_muldoubl(216,ZT[5]) )
                , i_compl_muldoubl(41,ZT[N]) )
              );
      INC=0;
      for(I=1; I<=K-1; I++)
      {
/*         ZOM=ZOM+82*ZT(6*I+1)+216*ZT(6*I+2)+27*ZT(6*I+3)+272*ZT(6*I+4)
     $           +27*ZT(6*I+5)+216*ZT(6*I+6)*/
       INC +=6;
       *ZOM = i_compl_add(
             i_compl_add(
                i_compl_add( i_compl_muldoubl(82,ZT[INC]),  i_compl_muldoubl(216,ZT[INC+1]) )
                , i_compl_add( i_compl_muldoubl(27,ZT[INC+2]),  i_compl_muldoubl(272,ZT[INC+3]) ) )
             , i_compl_add(              
                i_compl_add( i_compl_muldoubl(27,ZT[INC+4]),  i_compl_muldoubl(216,ZT[INC+5]) )
                , *ZOM )
              );
      }
      *ZOM= i_compl_muldoubl(H*6.E0/840.E0, *ZOM); /*ZOM =ZOM*H*6.D0/840.D0*/
#endif /*0*/ /*remplacee par:*/
     /*ZOM=41*ZT(1)+216*ZT(2)+27*ZT(3)+272*ZT(4)
     $          +27*ZT(5)+216*ZT(6)+41*ZT(N)*/
      zomreel=41 * ( *pdZT++); zomimag=41 * ( *pdZT++);
      zomreel+=216 * ( *pdZT++); zomimag+=216 * ( *pdZT++);
      zomreel+=27 * ( *pdZT++); zomimag+=27 * ( *pdZT++);
      zomreel+=272 * ( *pdZT++); zomimag+=272 * ( *pdZT++);
      zomreel+=27 * ( *pdZT++); zomimag+=27 * ( *pdZT++);
      zomreel+=216 * ( *pdZT++); zomimag+=216 * ( *pdZT++);
      zomreel+=41 * ( ZT[N].reel); zomimag+=41 * ( ZT[N].imag);

      for(I=1; I<=K-1; I++)
      {
/*         ZOM=ZOM+82*ZT(6*I+1)+216*ZT(6*I+2)+27*ZT(6*I+3)+272*ZT(6*I+4)
     $           +27*ZT(6*I+5)+216*ZT(6*I+6)*/
       zomreeltemp=82 * ( *pdZT++); zomimagtemp=82 * ( *pdZT++);
       zomreeltemp+=216 * ( *pdZT++); zomimagtemp+=216 * ( *pdZT++);
       zomreeltemp+=27 * ( *pdZT++); zomimagtemp+=27 * ( *pdZT++);
       zomreeltemp+=272 * ( *pdZT++); zomimagtemp+=272 * ( *pdZT++);
       zomreeltemp+=27 * ( *pdZT++); zomimagtemp+=27 * ( *pdZT++);
       zomreeltemp+=216 * ( *pdZT++); zomimagtemp+=216 * ( *pdZT++);
       zomreel += zomreeltemp; zomimag += zomimagtemp;
      }
      ZOM->reel = H*6.E0/840.E0*zomreel;
      ZOM->imag = H*6.E0/840.E0*zomimag;
     /* v0.96 M. GASTINEAU 13/01/99 : fin d'optimisation */
#endif /*NAF_USE_OPTIMIZE==0*/
/* v0.96 M. GASTINEAU 01/12/98 : fin optimisation */
  return TRUE;
}/*      END SUBROUTINE ZARDYD*/


/*      SUBROUTINE PROSCAA(F1,F2,ZP)*/
BOOL naf_proscaa(double F1, double F2, t_complexe *ZP)
{
/*!-----------------------------------------------------------------------
!     PROSCAA   CALCULE LE PRODUIT SCALAIRE 
!           < EXP(I*F1*T),  EXP(I*F2*T)>
!               SUR L'INTERVALLE [0:g_NAFVariable.KTABS]  T=g_NAFVariable.T0+g_NAFVariable.XH*IT
!              CALCUL NUMERIQUE
!     ZP        PRODUIT SCALAIRE COMPLEXE
!     F1,F2     FREQUENCES EN "/AN
!               REVU LE 26/9/87 J. LASKAR
!-----------------------------------------------------------------------
      
      IMPLICIT NONE
! (g_NAFVariable.KTABS,F1,F2,ZP)
      REAL (8) :: F1,F2
      complex (8) :: ZP
!
      integer :: LTF
      real *8 :: OM,ANG0,ANGI,H
      complex (8) :: ZI,ZAC,ZINC,ZEX
      complex (8), dimension (:),allocatable :: ZTF
      
!*/
      int LTF;
      double OM,ANG0,ANGI,H;
      t_complexe ZI,ZAC,ZINC,ZEX;
      t_complexe *ZTF=NULL;
      
      SYSCHECKMALLOCSIZE(ZTF, t_complexe, g_NAFVariable.KTABS+1); /*allocate(ZTF(0:g_NAFVariable.KTABS),stat = nerror)*/
/* v0.96 M. GASTINEAU 01/12/98 : utilisation des fonctions  complexes inlines et optimisation */
#if NAF_USE_OPTIMIZE==0
      ZI=cmplx(0.E0,1.E0);
/*!----------! FREQUENCES EN UNITE D'ANGLE PAR UNITE DE TEMPS*/
      OM = (F1-F2)/g_NAFVariable.UNIANG;
      LTF=g_NAFVariable.KTABS;
      ANG0=OM*g_NAFVariable.T0;
      ANGI=OM*g_NAFVariable.XH;
      ZAC= expcomplexe(muldoublcomplexe(-ANG0, ZI)); /*ZAC = EXP (-ZI*ANG0)*/
      ZINC=expcomplexe(muldoublcomplexe(-ANGI, ZI)); /*ZINC= EXP (-ZI*ANGI)*/
      ZEX= divcomplexe(ZAC, ZINC); /*ZEX = ZAC/ZINC*/      
      naf_ztpow2a(g_NAFVariable.KTABS,64,ZTF,TWIN,ZINC,ZEX);/*CALL  ZTPOW2A(g_NAFVariable.KTABS+1,64,ZTF,TWIN,ZINC,ZEX)*/
#else/*remplacee par:*/
      i_compl_cmplx(&ZI,0.E0,1.E0);
/*!----------! FREQUENCES EN UNITE D'ANGLE PAR UNITE DE TEMPS*/
      OM = (F1-F2)/g_NAFVariable.UNIANG;
      LTF=g_NAFVariable.KTABS;
      ANG0=OM*g_NAFVariable.T0;
      ANGI=OM*g_NAFVariable.XH;
      ZAC= i_compl_exp(i_compl_muldoubl(-ANG0, ZI)); /*ZAC = EXP (-ZI*ANG0)*/
      ZINC=i_compl_exp(i_compl_muldoubl(-ANGI, ZI)); /*ZINC= EXP (-ZI*ANGI)*/
      ZEX= i_compl_div(ZAC, ZINC); /*ZEX = ZAC/ZINC*/      
      naf_ztpow2a(g_NAFVariable.KTABS,64,ZTF,TWIN,ZINC,ZEX);/*CALL  ZTPOW2A(g_NAFVariable.KTABS+1,64,ZTF,TWIN,ZINC,ZEX)*/
#endif /*NAF_USE_OPTIMIZE==0*/
/* v0.96 M. GASTINEAU 01/12/98 : fin optimisation */

/*!------------------ TAILLE DU PAS*/
      H=1.E0/LTF;
      /*CALL ZARDYD(ZTF,LTF+1,H,ZP)*/
      if (naf_zardyd(ZTF,LTF,H,ZP)==FALSE)
      { 
       SYSFREE(ZTF);
       return FALSE;
      }
      SYSFREE(ZTF);
      return TRUE; 
}/*      END SUBROUTINE PROSCAA*/

/*      SUBROUTINE ZTPOW2A (N,N1,ZTF,TW,ZA,ZAST)*/
void naf_ztpow2a(int N, int N1, t_complexe *ZTF, double *TW, t_complexe ZA, t_complexe ZAST)
{
/*!-----------------------------------------------------------------------
!     ZTPOW   CALCULE  ZTF(I) = TW(I)*ZAST*ZA**I EN VECTORIEL
!             ZT(N1) ZONE DE TRAVAIL
!-----------------------------------------------------------------------
      IMPLICIT NONE
!  (N,N1,ZT,ZTF,TW,ZA,ZAST)
      integer :: N,N1
      real *8 :: TW(0:N)
      complex (8) :: ZTF(0:N),ZA,ZAST
!
      integer :: I,INC,IT, NX,NT
      complex (8) :: ZT1,ZINC     
      complex (8), dimension (:),allocatable :: ZT
!      */
      int I,INC,IT, NX,NT;
      t_complexe ZT1,ZINC;
      t_complexe *ZT=NULL;
      
      SYSCHECKMALLOCSIZE(ZT, t_complexe, N1); /*allocate(ZT(0:N1-1),stat = nerror)*/
      if (N<N1-1)
      {
         fprintf(stdout,"DANS ZTPOW, N = %d\n",N);
         return;
      }
/*!----------! */
/* v0.96 M. GASTINEAU 01/12/98 : utilisation des fonctions  complexes inlines et optimisation */
#if NAF_USE_OPTIMIZE==0
      ZT[0] = mulcomplexe(ZAST,ZA);
      for(I = 1; I<N1; I++)
      {
         ZT[I] = mulcomplexe(ZT[I-1],ZA);
      }
      for(I = 0;I<N1; I++)
      {
         ZTF[I] = muldoublcomplexe(TW[I],ZT[I]);
      }
      ZT1 = divcomplexe(ZT[N1-1],ZAST); /*ZT1 = ZT(N1-1)/ZAST*/
      ZINC= cmplx(1E0,0E0);
      INC =0;
      NT = (N+1)/N1  ;
      for(IT = 2;IT<= NT; IT++)
      {
         ZINC=mulcomplexe(ZINC, ZT1); /*ZINC = ZINC*ZT1*/
         INC += N1; /*INC  = INC + N1*/
         for(I = 0; I<N1; I++)
         {
            ZTF[INC +I] = muldoublcomplexe(TW[INC+I],mulcomplexe(ZT[I],ZINC));/*ZTF(INC +I) = TW(INC+I)*ZT(I)*ZINC*/
            
         }
      }
      ZINC=mulcomplexe(ZINC, ZT1); /*ZINC = ZINC*ZT1*/
      INC += N1; /*INC  = INC + N1*/
      NX = N+1-NT*N1;
      for(I = 0;I< NX; I++)
      {
            ZTF[INC +I] = muldoublcomplexe(TW[INC+I],mulcomplexe(ZT[I],ZINC));/*ZTF(INC +I) = TW(INC+I)*ZT(I)*ZINC*/
      }
#else /*remplacee par:*/
      ZT[0] = i_compl_mul(ZAST,ZA);
      for(I = 1; I<N1; I++)
      {
         ZT[I] = i_compl_mul(ZT[I-1],ZA);
      }
      for(I = 0;I<N1; I++)
      {
         ZTF[I] = i_compl_muldoubl(TW[I],ZT[I]);
      }
      ZT1 = i_compl_div(ZT[N1-1],ZAST); /*ZT1 = ZT(N1-1)/ZAST*/
      i_compl_cmplx(&ZINC,1E0,0E0);
      INC =0;
      NT = (N+1)/N1  ;
      for(IT = 2;IT<= NT; IT++)
      {
         i_compl_pmul(&ZINC,&ZT1); /*ZINC = ZINC*ZT1*/
         INC += N1; /*INC  = INC + N1*/
         for(I = 0; I<N1; I++)
         {
            ZTF[INC +I] = i_compl_muldoubl(TW[INC+I],i_compl_mul(ZT[I],ZINC));/*ZTF(INC +I) = TW(INC+I)*ZT(I)*ZINC*/
            
         }
      }
      i_compl_pmul(&ZINC,&ZT1); /*ZINC = ZINC*ZT1*/
      INC += N1; /*INC  = INC + N1*/
      NX = N+1-NT*N1;
      for(I = 0;I< NX; I++)
      {
            ZTF[INC +I] = i_compl_muldoubl(TW[INC+I],i_compl_mul(ZT[I],ZINC));/*ZTF(INC +I) = TW(INC+I)*ZT(I)*ZINC*/
      }
#endif /*NAF_USE_OPTIMIZE==0*/
/* v0.96 M. GASTINEAU 01/12/98 : fin optimisation */
      SYSFREE(ZT);
}/*      END SUBROUTINE ZTPOW2A*/

/*      SUBROUTINE CORRECTION(FREQ)*/
void naf_correction(double *FREQ)
{
/*!------------------------------------------------------------------
!   RETOURNE LA PREMIERE FREQUENCE (FREQ) CORRIGE PAR LES SUIVANTES
!
!  F. Joutel 26/5/98
!  Ref : J. Laskar : Introduction to frequency map analysis
!------------------------------------------------------------------
       
      IMPLICIT NONE
! (FREQ)
      real *8 :: FREQ
!
      integer :: I           
      real *8 :: TL,TM2,PICARRE,TM,FAC,A,COR,OMEGA,DELTA
      complex (8) :: ZI,ZALPHA
!      
!      PI = ATAN2(1.D0,1.D0)*4.D0*/
      int I;
      double TL,TM2,PICARRE,TM,FACTEUR,A,COR,OMEGA,DELTA;
      t_complexe ZI,ZALPHA;
      
      TL=g_NAFVariable.KTABS*g_NAFVariable.XH*0.5E0;
      TM2=1.E0/(TL*TL);
      PICARRE=M_PI*M_PI;
/*!  g_NAFVariable.IW >= 1*/
      TM=g_NAFVariable.T0+TL;
      i_compl_cmplx(&ZI,0.E0,1.E0);
      if (g_NAFVariable.IW>=0)
      {
         FACTEUR = -TM2;
         A=-1.E0/3.E0;
         for(I=1;I<=g_NAFVariable.IW;I++)
         {
            FACTEUR *= -PICARRE*TM2*(I*I);  /*FACTEUR=-FACTEUR*PICARRE*TM2*(I**2)*/
            A+=2.E0/(PICARRE*I*I); /*A=A+2.D0/(PICARRE*I**2) */
         }
         FACTEUR /=A; /*FACTEUR = FACTEUR/A*/
         COR=0.E0;
         for(I=2;I<=g_NAFVariable.NFS;I++)
         {
          OMEGA=(g_NAFVariable.TFS[I]-g_NAFVariable.TFS[1])/g_NAFVariable.UNIANG;
/* v0.96 M. GASTINEAU 01/12/98 : utilisation des fonctions  complexes inlines et optimisation */
#if NAF_USE_OPTIMIZE==0
          ZALPHA = mulcomplexe(divcomplexe(g_NAFVariable.ZAMP[I],g_NAFVariable.ZAMP[1]), expcomplexe(muldoublcomplexe(OMEGA*TM,ZI))); /*ZALPHA=g_NAFVariable.ZAMP[I]/g_NAFVariable.ZAMP[1]*EXP(ZI*OMEGA*TM)*/
#else /*remplacee par:*/
          ZALPHA = i_compl_mul(i_compl_div(g_NAFVariable.ZAMP[I],g_NAFVariable.ZAMP[1]), 
                               i_compl_exp(i_compl_muldoubl(OMEGA*TM,ZI))); 
          /*ZALPHA=g_NAFVariable.ZAMP[I]/g_NAFVariable.ZAMP[1]*EXP(ZI*OMEGA*TM)*/
#endif /*NAF_USE_OPTIMIZE==0*/
/* v0.96 M. GASTINEAU 01/12/98 : fin optimisation */
          DELTA=FACTEUR*ZALPHA.reel/(pow(OMEGA,(2*g_NAFVariable.IW+1)))*cos(OMEGA*TL); /*DELTA=FACTEUR*DREAL(ZALPHA)/OMEGA**(2*g_NAFVariable.IW+1)*COS(OMEGA*TL)*/
          COR += DELTA; /*COR=COR+DELTA*/
          if (g_NAFVariable.IPRT>=2)
          {
             fprintf(g_NAFVariable.NFPRT,"OMEGA=%g, DELTA=%g ,COR=%g\n",OMEGA,DELTA,COR);
          }
         }
         COR *= g_NAFVariable.UNIANG; /*COR=COR*g_NAFVariable.UNIANG*/
         *FREQ=g_NAFVariable.TFS[1]+COR;
         if (g_NAFVariable.IPRT>=1)
         {
           fprintf(g_NAFVariable.NFPRT,"CORRECTION DE LA 1ERE FREQUENCE\n");  
           fprintf(g_NAFVariable.NFPRT, "FREQ. DE DEPART=%20.15E, CORRECTION=%20.15E\n",g_NAFVariable.TFS[1],COR);
           fprintf(g_NAFVariable.NFPRT, "FREQUENCE CORRIGEE:%20.15E\n", *FREQ);
         }
      } 
      else 
      {
          if (g_NAFVariable.IPRT>=0) 
          {
             fprintf(g_NAFVariable.NFPRT, "PAS DE CALCUL\n");
          }
      }
}/*      END SUBROUTINE CORRECTION*/
      
/*      end module naff
*/
/*-----------------------------------------------------------*/
/* cree et retourneun element de liste de fenetre pour naf.  */
/* Les champs sont remplis avec les arguments de la fonction.*/
/* L'element retourne doit etre libere avec                  */
/* delete_list_fenetre_naf.                                  */
/*-----------------------------------------------------------*/
/* v0.96 M. GASTINEAU 18/12/98 : ajout */
t_list_fenetre_naf *cree_list_fenetre_naf(const double p_dFreqMin, const double p_dFreqMax, 
                                          const int p_iNbTerm)
{
 t_list_fenetre_naf *pListRes;
#ifdef DBUG
 if (niveau>=4) 
 {
  printf("cree_list_fenetre_naf(%g, %g, %d) - entree \n",
          p_dFreqMin, p_dFreqMax, p_iNbTerm);
 }
#endif /*DBUG*/
 SYSCHECKMALLOC(pListRes, t_list_fenetre_naf);
 pListRes->suivant = NULL;
 pListRes->dFreqMin = p_dFreqMin;
 pListRes->dFreqMax = p_dFreqMax;
 pListRes->iNbTerme = p_iNbTerm;
#ifdef DBUG
 if (niveau>=4) 
 {
  printf("cree_list_fenetre_naf() - sortie et retourne 0x%p \n",pListRes);
 }
#endif /*DBUG*/
 return pListRes;
}

/*-----------------------------------------------------------*/
/*Detruit une liste cree avec cree_list_fenetre_naf.         */
/* En sortie, p_pListFenNaf pointe vers une adresse invalide.*/
/*-----------------------------------------------------------*/
/* v0.96 M. GASTINEAU 18/12/98 : ajout */
void delete_list_fenetre_naf(t_list_fenetre_naf *p_pListFenNaf)
{
#ifdef DBUG
 if (niveau>=4) 
  printf("delete_list_fenetre_naf(0x%p) - entree \n",p_pListFenNaf);
#endif /*DBUG*/
 while(p_pListFenNaf!=NULL)
 {
  t_list_fenetre_naf *next=p_pListFenNaf->suivant;
  SYSFREE(p_pListFenNaf);
  p_pListFenNaf=next;
 }

#ifdef DBUG
 if (niveau>=4) 
  printf("delete_list_fenetre_naf() - sortie \n");
#endif /*DBUG*/ 
}

/*-----------------------------------------------------------*/
/* Concatene 2 listes de fenetre de naf.                     */
/* En sortie, p_pListFenHead et p_pListFenEnd ne doivent pas */
/* etre liberee. Seul la liste retournee doit etre liberee.  */
/* Les deux parametres peuvent etre NULL.                    */
/*-----------------------------------------------------------*/
/* v0.96 M. GASTINEAU 18/12/98 : ajout */
t_list_fenetre_naf *concat_list_fenetre_naf(t_list_fenetre_naf *p_pListFenHead,
                                            t_list_fenetre_naf *p_pListFenEnd)
{
 t_list_fenetre_naf *pListTemp;
#ifdef DBUG
 if (niveau>=4) 
  printf("concat_list_fenetre_naf(0x%p, 0x%p) - entree \n",p_pListFenHead,p_pListFenEnd);
#endif /*DBUG*/
 if (p_pListFenHead==NULL)
 {
  return p_pListFenEnd;
 }
 else if (p_pListFenEnd==NULL)
 {
  return p_pListFenHead;
 }
 else
 { /* p_pListFenHead is not NULL */
  for(pListTemp=p_pListFenHead; 
      pListTemp->suivant!=NULL;
      pListTemp = pListTemp->suivant)
      { /*do nothing */}
  /*concatene */
  pListTemp->suivant = p_pListFenEnd;
  return p_pListFenHead;
 }
}

