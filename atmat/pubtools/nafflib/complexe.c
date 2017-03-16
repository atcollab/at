/* Astronomie et Systemes dynamiques                            */
/* M. GASTINEAU  Bureau des Longitudes, Paris 15/07/98          */
/* ensemble de fonctions inlines traitant uniquement les complexes */
/* v0.97 M. GASTINEAU 15/02/99: ajout des fonctions trigonometrique*/
/* v0.97 M. GASTINEAU 22/03/99: ajout de i_compl_powreel        */

#include "modnaff.h"

#ifndef __COMPLEXE_H
#define __COMPLEXE_H
#include "complexe.h"

/*----------------IMPLEMENTATION--------------------------------*/

/*----------------i_compl_cmplx---------------------------------*/
/*Construit un complexe a partir de deux doubles(retourne c=a+i*b)*/
/*--------------------------------------------------------------*/
 void i_compl_cmplx(t_complexe *c, double a,double b)
{
 c->reel=a;
 c->imag=b;
}


/*----------------module----------------------------------------*/
/* retourne le module du nombre complexe c                      */
/*--------------------------------------------------------------*/
 double    i_compl_module(t_complexe c)
{
 return hypot(c.reel,c.imag);
}

 double    i_compl_angle(t_complexe c)
{
 /*return atan2(c.imag/module(c),c.reel/module(c));*/
 return atan2(c.imag,c.reel);
}

/*----------------i_compl_add-----------------------------------*/
/* addition de deux nombres complexes c1 et c2 : c1+c2          */
/*--------------------------------------------------------------*/
 t_complexe i_compl_add(const t_complexe c1,const t_complexe c2)
{
 t_complexe c;
 i_compl_cmplx(&c, c1.reel+c2.reel,c1.imag+c2.imag);
 return c;
}

/*----------------i_compl_padd----------------------------------*/
/* addition de deux nombres complexes c1 et c2 : c1+=c2         */
/*--------------------------------------------------------------*/
 void i_compl_padd(t_complexe *c1,t_complexe *c2)
{
 c1->reel += c2->reel;
 c1->imag += c2->imag;
}

/*----------------i_compl_paddconst-----------------------------*/
/* addition de deux nombres complexes c1 et c2 : c1+=c2         */
/*--------------------------------------------------------------*/
 void i_compl_paddconst(t_complexe *c1,t_complexe c2)
{
 c1->reel += c2.reel;
 c1->imag += c2.imag;
}

/*----------------i_compl_mul----------------------------------*/
/* multiplication de deux nombres complexes c1 et c2 : c1*c2   */
/*-------------------------------------------------------------*/
 t_complexe i_compl_mul(const t_complexe c1,const t_complexe c2)
{
 t_complexe c;
 i_compl_cmplx(&c, c1.reel*c2.reel-c1.imag*c2.imag,c1.reel*c2.imag+c1.imag*c2.reel);
 return c;
}

/*----------------i_compl_pmul---------------------------------*/
/* multiplication de deux nombres complexes c1 et c2 : c1*=c2  */
/*-------------------------------------------------------------*/
 void i_compl_pmul(t_complexe *c1,t_complexe *c2)
{
 double creel=c1->reel;
 c1->reel=creel*c2->reel-c1->imag*c2->imag;
 c1->imag =creel*c2->imag+c1->imag*c2->reel;
}

/*----------------i_compl_muldoubl-----------------------------*/
/* multiplication d'un double c1 d'un nombre complexe c2: c1*c2*/
/*-------------------------------------------------------------*/
 t_complexe i_compl_muldoubl(const double c1, const t_complexe c2)
{
 t_complexe c;
 c.reel = c1*c2.reel;
 c.imag = c1*c2.imag;
 return c;
}

/*----------------i_compl_pmuldoubl----------------------------*/
/* multiplication d'un double c1 d'un nombre complexe c2: c2*=c1*/
/*-------------------------------------------------------------*/
 void i_compl_pmuldoubl(t_complexe *c2, double* const c1)
{
 c2->reel *= *c1;
 c2->imag *= *c1;
}

/*----------------i_compl_pdivdoubl----------------------------*/
/* division d'un double c1 d'un nombre complexe c2: c2/=c1     */
/*-------------------------------------------------------------*/
 void i_compl_pdivdoubl(t_complexe *c2, double *c1)
{
 c2->reel /= *c1;
 c2->imag /= *c1;
}

/*----------------i_compl_conj---------------------------------*/
/*retourne le conjugue de c1                                   */
/*-------------------------------------------------------------*/
 t_complexe i_compl_conj(t_complexe *c1)
{
 t_complexe c;
 i_compl_cmplx(&c, c1->reel, - (c1->imag) );
 return c;
}

/*----------------i_compl_pdiv---------------------------------*/
/* division de deux complexes : c1/=c2                         */
/*-------------------------------------------------------------*/
 void i_compl_pdiv(t_complexe * const c1,const t_complexe * const c2)
{
	double ratio, den;
	double abr, abi, cr;
    
	if( (abr = c2->reel) < 0.)
		abr = - abr;
	if( (abi = c2->imag) < 0.)
		abi = - abi;
	if( abr <= abi )
		{
		ratio = (double)c2->reel / c2->imag ;
		den = c2->imag * (1 + ratio*ratio);
		cr = (c1->reel*ratio + c1->imag) / den;
		c1->imag = (c1->imag*ratio - c1->reel) / den;
		}

	else
		{
		ratio = (double)c2->imag / c2->reel ;
		den = c2->reel * (1 + ratio*ratio);
		cr = (c1->reel + c1->imag*ratio) / den;
		c1->imag = (c1->imag - c1->reel*ratio) / den;
		}
	c1->reel = cr;
}

/*----------------i_compl_div2d--------------------------------*/
/* division d'un reel par un complexe c1/(c2_r+i*c2_i)         */
/*-------------------------------------------------------------*/
t_complexe i_compl_div2d(const double c1, 
                                       const double c2_r, 
                                       const double c2_i)
{
	register double ratio, den;
	register double abr, abi;
    t_complexe zres;
    
	if( (abr = c2_r) < 0.)
		abr = - abr;
	if( (abi = c2_i) < 0.)
		abi = - abi;
	if( abr <= abi )
		{
		ratio = (double)c2_r / c2_i ;
		den = c2_i * (1 + ratio*ratio);
		zres.reel = (c1*ratio) / den;
		zres.imag =  - c1 / den;
		}

	else
		{
		ratio = (double)c2_i / c2_r ;
		den = c2_r * (1 + ratio*ratio);
		zres.reel = c1 / den;
		zres.imag =  (- c1*ratio) / den;
		}
    return zres;
}

/*----------------i_compl_div4d--------------------------------*/
/* division de deux complexes : (c1_r+i*c1_i)/(c2_r+i*c2_i)    */
/*-------------------------------------------------------------*/
t_complexe i_compl_div4d(register const double c1_r,register const double c1_i,
                                       register const double  c2_r,register const double c2_i)
{
	register double ratio, den;
	register double abr, abi;
    t_complexe zres;
    
	if( (abr = c2_r) < 0.)
		abr = - abr;
	if( (abi = c2_i) < 0.)
		abi = - abi;
	if( abr <= abi )
		{
		ratio = (double)c2_r / c2_i ;
		den = c2_i * (1 + ratio*ratio);
		zres.reel = (c1_r*ratio + c1_i) / den;
		zres.imag = (c1_i*ratio - c1_r) / den;
		}

	else
		{
		ratio = (double)c2_i / c2_r ;
		den = c2_r * (1 + ratio*ratio);
		zres.reel = (c1_r + c1_i*ratio) / den;
		zres.imag = (c1_i - c1_r*ratio) / den;
		}
 return zres; 
}

/*----------------i_compl_div----------------------------------*/
/* division de deux complexes : c1/c2                          */
/*-------------------------------------------------------------*/
t_complexe i_compl_div(const t_complexe c1,const t_complexe c2)
{
 t_complexe c=c1;
 i_compl_pdiv(&c,&c2);
 return c;
}


/*----------------i_compl_pow2d--------------------------------*/
/*retourne la puissance (p_r+i*p_i)**b                         */
/*-------------------------------------------------------------*/
t_complexe i_compl_pow2d(const double p_r, const double p_i,int n)
{
	register unsigned long u;
	register double t;
	register double q_r=1,q_i=0,x_r,x_i;
    t_complexe zq;
    
	if(n == 0)
		goto done;
	if(n < 0)
		{
		n = -n;
		zq = i_compl_div2d(1,p_r,p_i);
		x_r = zq.reel;
		x_i = zq.imag;
		}
	else
		{
		x_r = p_r;
		x_i = p_i;
		}

	for(u = n; ; )
		{
		if(u & 01)
			{
			t = q_r * x_r - q_i * x_i;
			q_i = q_r * x_i + q_i * x_r;
			q_r = t;
			}
		if(u >>= 1)
			{
			t = x_r * x_r - x_i * x_i;
			x_i = 2 * x_r * x_i;
			x_r = t;
			}
		else
			break;
		}
 done:
	zq.reel=q_r;
	zq.imag=q_i;
	return zq;
}


/*----------------i_compl_powreel------------------------------*/
/*retourne la puissance a**b avec b reel                       */
/*-------------------------------------------------------------*/
/* v0.97 M. GASTINEAU 22/03/99: ajout de i_compl_powreel       */
t_complexe i_compl_powreel(const t_complexe z,double p_dK)
{
 double dMod = pow(i_compl_module(z),p_dK);
 double dAngle = p_dK*i_compl_angle(z);
 t_complexe zRes;
 zRes.reel=dMod*cos(dAngle);
 zRes.imag=dMod*sin(dAngle);
 return zRes;
}

/*----------------i_compl_pow----------------------------------*/
/*retourne la puissance a**b                                   */
/*-------------------------------------------------------------*/
 t_complexe i_compl_pow(const t_complexe a,int n)
{
#if 1
	unsigned long u;
	double t;
	t_complexe q={1.0, 0.0}, x;

	if(n == 0)
		goto done;
	if(n < 0)
		{
		n = -n;
		x.reel = 1.0;
		x.imag=0.0;
		i_compl_pdiv(&x, &a);
		}
	else
		{
		x.reel = a.reel;
		x.imag = a.imag;
		}

	for(u = n; ; )
		{
		if(u & 01)
			{
			t = q.reel * x.reel - q.imag * x.imag;
			q.imag = q.reel * x.imag + q.imag * x.reel;
			q.reel = t;
			}
		if(u >>= 1)
			{
			t = x.reel * x.reel - x.imag * x.imag;
			x.imag = 2 * x.reel * x.imag;
			x.reel = t;
			}
		else
			break;
		}
 done:
	return q;
#else
	t_complexe q;

	if(n == 0)
	{
	 q.reel=1;
	 q.imag=0;
	}
	else
	{
	 double dpow=pow(i_compl_module(a),n);
	 double anglea=i_compl_angle(a);
	 q.reel=dpow*cos(anglea*n);
	 q.imag=dpow*sin(anglea*n);
	}
	return q;
#endif /*0*/
}

/*----------------expcomplexe-----------------------------------*/
/* calcule et retourne exp(c1)                                  */
/*--------------------------------------------------------------*/
/* v0.96 M. GASTINEAU 04/09/98 : ajout */
t_complexe i_compl_exp(t_complexe c1)
{
 t_complexe r;
 double expx;
 expx = exp(c1.reel);
 r.reel = expx * cos(c1.imag);
 r.imag = expx * sin(c1.imag);
 return r;
}


/*----------------i_compl_psub----------------------------------*/
/* soustrait de deux nombres complexes c1 et c2 : c1-=c2        */
/*--------------------------------------------------------------*/
/*v0.96 M. GASTINEAU 01/12/98 : ajout */
void i_compl_psub(t_complexe *c1,t_complexe *c2) 
{
 c1->reel -= c2->reel;
 c1->imag -= c2->imag;
}

/*----------------i_compl_sub----------------------------------*/
/* soustrait de deux nombres complexes c1 et c2 : c1-c2        */
/*-------------------------------------------------------------*/
/*v0.96 M. GASTINEAU 01/12/98 : ajout */
t_complexe i_compl_sub(t_complexe c1,t_complexe c2)
{
 t_complexe c;
 c.reel = c1.reel - c2.reel;
 c.imag = c1.imag - c2.imag;
 return c;
}


/*----------------i_compl_cos----------------------------------*/
/* retourne le cosinus de c1 (cos x  cosh y  -  i sin x sinh y)*/ 
/*-------------------------------------------------------------*/
/* v0.97 M. GASTINEAU 15/02/99: ajout */
t_complexe i_compl_cos(t_complexe c1)
{
 t_complexe c;
 c.reel = cos(c1.reel)*cosh(c1.imag);
 c.imag = -sin(c1.reel)*sinh(c1.imag);
 return c;
}

/*----------------i_compl_sin----------------------------------*/
/* retourne le sinus de c1 (= sin x  cosh y  +  i cos x sinh y)*/
/*-------------------------------------------------------------*/
/* v0.97 M. GASTINEAU 15/02/99: ajout */
t_complexe i_compl_sin(t_complexe c1)
{
 t_complexe c;
 c.reel = sin(c1.reel)*cosh(c1.imag);
 c.imag = cos(c1.reel)*sinh(c1.imag);
 return c;
}

/*----------------i_compl_cosh---------------------------------*/
/* retourne le cosinus hyperbolique de c1                      */
/* (= cosh x  cos y  +  i sinh x sin y)                        */
/*-------------------------------------------------------------*/
/* v0.97 M. GASTINEAU 15/02/99: ajout */
t_complexe i_compl_cosh(t_complexe c1)
{
 t_complexe c;
 c.reel = cosh(c1.reel)*cos(c1.imag);
 c.imag = sinh(c1.reel)*sin(c1.imag);
 return c;
}

/*----------------i_compl_sinh---------------------------------*/
/* retourne le sinus hyperbolique de c1                        */
/* (= sinh x  cos y  +  i cosh x sin y)                        */
/*-------------------------------------------------------------*/
/* v0.97 M. GASTINEAU 15/02/99: ajout */
t_complexe i_compl_sinh(t_complexe c1)
{
 t_complexe c;
 c.reel = sinh(c1.reel)*cos(c1.imag);
 c.imag = cosh(c1.reel)*sin(c1.imag);
 return c;
}

/*----------------i_compl_tan----------------------------------*/
/* retourne la tangente de c1                                  */
/*-------------------------------------------------------------*/
/* v0.97 M. GASTINEAU 15/02/99: ajout */
t_complexe i_compl_tan(t_complexe c1)
{
 t_complexe c;
 double u2=2.E0*c1.reel;
 double v2=2.E0*c1.imag;
 double denom=cos(u2)+cosh(v2);
 c.reel = sin(u2)/denom;
 c.imag = sinh(v2)/denom;
 return c;
}

/*----------------i_compl_tanh---------------------------------*/
/* retourne la tangente hyperbolique de c1                     */
/*-------------------------------------------------------------*/
/* v0.97 M. GASTINEAU 15/02/99: ajout */
t_complexe i_compl_tanh(t_complexe c1)
{
 t_complexe c;
 double u2=2.E0*c1.reel;
 double v2=2.E0*c1.imag;
 double denom=cos(u2)+cosh(v2);
 c.reel = sinh(u2)/denom;
 c.imag = sin(v2)/denom;
 return c;
}

/*----------------i_compl_log----------------------------------*/
/* retourne le logarithme de c1                                */
/*-------------------------------------------------------------*/
/* v0.97 M. GASTINEAU 15/02/99: ajout */
t_complexe i_compl_log(t_complexe c1)
{
 t_complexe c;
 c.reel = log(i_compl_module(c1));
 c.imag = i_compl_angle(c1);
 return c;
}

/*----------------i_compl_log10--------------------------------*/
/* retourne le logarithme  base 10 de c1                       */
/*-------------------------------------------------------------*/
/* v0.97 M. GASTINEAU 15/02/99: ajout */
t_complexe i_compl_log10(t_complexe c1)
{
 t_complexe c;
 double norm=log(10);
 c.reel = log(i_compl_module(c1))/norm;
 c.imag = i_compl_angle(c1)/norm;
 return c;
}

#endif
