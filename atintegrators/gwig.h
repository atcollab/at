/*
 *----------------------------------------------------------------------------
 * Modification Log:
 * -----------------
 * .02  2003-04-29      YK Wu, Duke University, wu@fel.duke.edu
 *                      use scientific notation for constants
 *
 * .01  2003-04-28      YK Wu, Duke University, wu@fel.duke.edu
 *                      header file for gwig.c
 *----------------------------------------------------------------------------
 *  Accelerator Physics Group, Duke FEL Lab, www.fel.duke.edu
*/

#define WHmax 200
#define GWIG_EPS (1e-6)

static const int Elem_Entrance = 1;
static const int Elem_Exit =-1;

static const double q_e    = 1.602176462e-19; /* electron charge, [C] */
static const double m_e    = 9.10938188e-31;  /* electron mass, [kg] */
static const double clight = 2.99792458e8;    /* speed of light [m/s] */
static const double r_e    = 2.817940285e-15; /* electron classic radius,[m]*/
static const double XMC2   = 0.510998902e-03; /* mc^2 in GeV */
static const double PI     = 3.141592653589793238462643383279502884197e0;

struct gwig{
  int Pmethod;      /* Integration Method */
  int PN;           /* Number of integration steps */
  double E0;        /* Energy of ring, [GeV] */
  double PB0;       /* B0 in [Tesla] */
  int Nw;           /* Number of periods */
  double Lw;        /* Wiggler Period [m] */
  int NHharm;       /* No. of horizontal harmonics */
  int NVharm;       /* No. of vertical harmonics */
  double Aw;        /* Wiggler parameter */
  double Zw;        /* Longitudinal variable [m] */

  double HCw[WHmax];
  double VCw[WHmax];
  double HCw_raw[WHmax];
  double VCw_raw[WHmax];
  double Hkx[WHmax];
  double Hky[WHmax];
  double Hkz[WHmax];
  double Htz[WHmax];
  double Vkx[WHmax];
  double Vky[WHmax];
  double Vkz[WHmax];
  double Vtz[WHmax];
};

