#include "AlphaDll.h"

EXP long STDCALL atinitialize2(long tptr, long nelems, const double *thresh,const double *energy, const double *freq);
EXP long STDCALL ataddtypeelem(long tptr, long nargs, const double *args, const char* passmethod, const char *xtype);
EXP long STDCALL attrack(long tptr, double *r_in, long np, long nturns);
EXP long STDCALL atstable(long tptr, const double *r_in,long nturns);
EXP double STDCALL atgetarg(long tptr, long nel, long np);
EXP const char* STDCALL atgetmethod(long tptr, long nel);
EXP long STDCALL atversion(void);

/* Kept for backward compatibilty */
EXP long STDCALL atinitialize(long tptr, long nelems, const double *threshold);
EXP long STDCALL xatsize(long nelems, double threshold);
EXP long STDCALL atsize(long tptr, long nelems, double threshold);
EXP long STDCALL xataddelem(long code, const double *args);
EXP long STDCALL ataddcodeelem(long tptr, long code, const double *args);
EXP long STDCALL ataddpasselem(long tptr, const char *passmethod, const double *args);
EXP long STDCALL xattrack(double *r_in, long np, long nturns);

EXP long STDCALL dbgatarg(long tptr, const char* passmethod, const double *args);
EXP long STDCALL dbgatstring(long code, const char* etype, const double *args);

void AperturePass(double *r_in, const double *limitsptr, int num_particles);
void DriftPass(double *r_in, double le,
               const double *T1, const double *T2,
               const double *R1, const double *R2,
               double *RApertures, double *EApertures,
               int num_particles);
void QuadLinearPass(double *r, double le, double kv, double *T1, double *T2, double *R1, double *R2, int num_particles);
void BendLinearPass(double *r, double le, double grd ,double ba, double bye,
	double entrance_angle, double exit_angle, double fint1, double fint2, double gap,
	double *T1, double *T2,	double *R1, double *R2, int num_particles);
void CavityPass(double *r_in, double le, double nv, double freq, double lag, int num_particles);
void ThinMPolePass(double *r, const double *A, const double *B, int max_order, double bax, double bay,
	double *T1, double *T2, double *R1, double *R2, int num_particles);
void CorrectorPass(double *r_in, double xkick, double ykick, double le, int num_particles);
void BndMPoleSymplectic4Pass(double *r, double le, double irho, double *A, double *B,
                             int max_order, int num_int_steps,
                             double entrance_angle, 	double exit_angle,
                             double fint1, double fint2, double gap,
                             double *T1, double *T2,
                             double *R1, double *R2,
                             double *RApertures, double *EApertures,int num_particles);
void BndMPoleSymplectic4RadPass(double *r, double le, double irho, double *A, double *B,
                                int max_order, int num_int_steps,
                                double entrance_angle, 	double exit_angle,
                                double fint1, double fint2, double gap,
                                double *T1, double *T2,
                                double *R1, double *R2,
                                double *RApertures, double *EApertures,
                                double E0, int num_particles);
void BndMPoleSymplectic4FrgFPass(double *r, double le, double irho, double *A, double *B,
                                 int max_order, int num_int_steps,
                                 double entrance_angle, 	double exit_angle,
                                 double fint1, double fint2, double gap,
                                 double *T1, double *T2,
                                 double *R1, double *R2,
                                 double *RApertures, double *EApertures, int num_particles);
void BndMPoleSymplectic4FrgFRadPass(double *r, double le, double irho, double *A, double *B,
                                    int max_order, int num_int_steps,
                                    double entrance_angle, 	double exit_angle,
                                    double fint1, double fint2, double gap,
                                    double *T1, double *T2,
                                    double *R1, double *R2,
                                    double *RApertures, double *EApertures,
                                    double E0, int num_particles);
void BndMPoleSymplectic4E2Pass(double *r, double le, double irho, const double *A, const double *B,
	int max_order, int num_int_steps,
	double entrance_angle, 	double exit_angle,
	double fint1, double fint2, double gap,
	double h1, double h2,
	double *T1, double *T2,	
	double *R1, double *R2, int num_particles);
void BndMPoleSymplectic4E2RadPass(double *r, double le, double irho, const double *A, const double *B,
	int max_order, int num_int_steps,
	double entrance_angle, 	double exit_angle,
	double fint1, double fint2, double gap,
	double h1, double h2,
	double *T1, double *T2,	
	double *R1, double *R2,
	double E0,
	int num_particles);
void StrMPoleSymplectic4Pass(double *r, double le, double *A, double *B,
                             int max_order, int num_int_steps,
                             double *T1, double *T2,
                             double *R1, double *R2,
                             double *RApertures, double *EApertures, int num_particles);
void QuadMPoleFringePass(double *r, double le, const double *A, const double *B,
                         int max_order, int num_int_steps,
                         const double *T1, const double *T2,
                         const double *R1, const double *R2,
                         double *fringeIntM0,  /* I0m/K1, I1m/K1, I2m/K1, I3m/K1, Lambda2m/K1 */
                         double *fringeIntP0,  /* I0p/K1, I1p/K1, I2p/K1, I3p/K1, Lambda2p/K1 */
                         double *RApertures, double *EApertures, int num_particles);
void WiggLinearPass(double *r, double le, double invrho, double kxkz,
        double *T1, double *T2, double *R1, double *R2, int num_particles);

void ATmultmv(double *r, const double* A);
void ATaddvv(double *r, const double *dr);
void ATdrift6(double* r, double L);
