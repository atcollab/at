/* Astronomie et Systemes dynamiques                            */
/* M. GASTINEAU  Bureau des Longitudes, Paris 15/07/98          */
/* ensemble de fonctions inlines traitant uniquement les complexes */
/* v0.97 M. GASTINEAU 15/02/99: ajout des fonctions trigonometrique*/
/* v0.97 M. GASTINEAU 22/03/99: ajout de i_compl_powreel        */

/* Modification Laurent */
#ifndef __COMPLEXEHEADER_
#define __COMPLEXEHEADER_

/* prototypes des fonctions publiques */
double     i_compl_module(t_complexe c);              /* niveau 1 */
double     i_compl_angle(t_complexe c);               /* niveau 1 */
t_complexe i_compl_add(t_complexe const c1,const t_complexe c2); /* niveau 2 */ 
t_complexe i_compl_mul(const t_complexe c1, const t_complexe c2); /* niveau 2 */ 
void       i_compl_padd(t_complexe* c1,t_complexe* c2); /* niveau 2 */ 
void       i_compl_paddconst(t_complexe *c1,t_complexe c2); /* niveau 2 */ 
void       i_compl_pmul(t_complexe *c1,t_complexe *c2); /* niveau 2 */ 
t_complexe i_compl_muldoubl(const double c1, const t_complexe c2); /* niveau 2 */ 
void       i_compl_pmuldoubl(t_complexe *c2, double *c1); /* niveau 2 */ 
void       i_compl_pdivdoubl(t_complexe *c2, double *c1); /* niveau 2 */ 
void       i_compl_cmplx(t_complexe *c, double a,double b);         /* niveau 2 */ 
void       i_compl_pdiv(t_complexe* const c1,const t_complexe * const c2); /* niveau 2 */ 
t_complexe i_compl_div2d(const double c1, const double  c2_r, const double c2_i); /* niveau 2 */ 
t_complexe i_compl_div4d(register const double c1_r, register const double c1_i,
                                       register const double  c2_r, register const double c2_i); /* niveau 2 */ 
t_complexe i_compl_conj(t_complexe *c1);             /* niveau 2 */
t_complexe i_compl_pow(const t_complexe z,int k);          /* niveau 2 */
/* v0.97 M. GASTINEAU 22/03/99: ajout de i_compl_powreel        */
t_complexe i_compl_powreel(const t_complexe z,double p_dK);          /* niveau 2 */
/* v0.97 M. GASTINEAU 22/03/99: fin d'ajout */
t_complexe i_compl_pow2d(const double p_r, const double p_i,int k);          /* niveau 2 */
t_complexe i_compl_exp(t_complexe c1);               /* niveau 2 */
void       i_compl_psub(t_complexe *c1,t_complexe *c2); /* niveau 2 */ 

t_complexe i_compl_div(const t_complexe c1,const t_complexe c2); /* niveau 2 */ 
t_complexe i_compl_sub(t_complexe c1,t_complexe c2); /* niveau 2 */ 

/* v0.97 M. GASTINEAU 15/02/99: ajout */
t_complexe i_compl_cos(t_complexe c1);               /* niveau 2 */
t_complexe i_compl_sin(t_complexe c1);               /* niveau 2 */
t_complexe i_compl_cosh(t_complexe c1);              /* niveau 2 */
t_complexe i_compl_sinh(t_complexe c1);              /* niveau 2 */
t_complexe i_compl_tan(t_complexe c1);               /* niveau 2 */
t_complexe i_compl_tanh(t_complexe c1);              /* niveau 2 */
t_complexe i_compl_log(t_complexe c1);               /* niveau 2 */
t_complexe i_compl_log10(t_complexe c1);             /* niveau 2 */
/* v0.97 M. GASTINEAU 15/02/99: fin ajout */
#endif
