/* Interpolate.c 
   Accelerator Toolbox 
   Created: 14/11/08
   Z.Martí
 
  
  
 Numerical recipies functions for 2D interpolation adapted to Matlab 
 * leanguagge
 
*/

#include "at.h"
#include <math.h>

/*x1a is the direction of the colums: y and x1a the direction of the rows:x*/

/****************************************************************************/
/****************************Bilinear interpolation*****************************/
/****************************************************************************/

void linint(double *x1a, double *x2a, double *ya, int m, int n, double x1, double x2, double *y)
{
    int klo,khi,k,ilo,ihi,i;
    double u,t,y1,y2,y3,y4,f;

    if ((x1<=x1a[m-1])&&(x1>=x1a[0])&&(x2<=x2a[n-1])&&(x2>=x2a[0])) {
        ilo=0;
        ihi=m-1;	
        klo=0;
        khi=n-1;
	
        while (ihi-ilo > 1) {
            i=(ihi+ilo) >> 1;
            if (x1a[i] > x1) ihi=i;
            else ilo=i;
        }
	
        while (khi-klo > 1) {
            k=(khi+klo) >> 1;
            if (x2a[k] > x2) khi=k;
            else klo=k;
        }
        y1=ya[ilo+klo*m];
        y2=ya[ilo+khi*m];
        y3=ya[ihi+khi*m];
        y4=ya[ihi+klo*m];
/*      Result will be NaN...
        if (x1a[ihi]==x1a[ilo]||x2a[khi]==x2a[klo])
            mexErrMsgTxt("Bad xa input to routine linint");
*/
        u=(x1-x1a[ilo])/(x1a[ihi]-x1a[ilo]);
        t=(x2-x2a[klo])/(x2a[khi]-x2a[klo]);
        f=(1-t)*(1-u)*y1+t*(1-u)*y2+t*u*y3+(1-t)*u*y4;
    }
    else {  /*This is redundant...The limits have been taking as apertures previously*/
        /*mexPrintf("Out of Id data range\n");*/
        f=sqrt(-1);
    }
    *y=f; 
}


/****************************************************************************/
/****************************Cubic interpolation*****************************/
/****************************************************************************/



void splin2(double *x1a, double *x2a, double *ya, double *y2a, int m, int n, double x1, double x2, double *y)
{
	void spline(double *x, double *y, int n, double yp1, double ypn, double *y2);
    void splint(double xa[], double ya[], double y2a[], int n, double x, double *y);
    int j,i;
	double *ytmp,*yytmp,*yaj,*y2aj;
      
	ytmp=(double*)mxCalloc(m,sizeof(double));
	yytmp=(double*)mxCalloc(m,sizeof(double));   
	yaj=(double*)mxCalloc(n,sizeof(double));
	y2aj=(double*)mxCalloc(n,sizeof(double)); 
    

	for (j=0;j<m;j++)
    {       
        /*Take one row of ya and y2a*/    
        for (i=0;i<n;i++)
        {
            yaj[i]=ya[j+i*m];
            y2aj[i]=y2a[j+i*m];
        }
        splint(x2a,yaj,y2aj,n,x2,&yytmp[j]);

    }
    mxFree(yaj);    
    mxFree(y2aj);

	spline(x1a,yytmp,m,1.0e30,1.0e30,ytmp);
	splint(x1a,yytmp,ytmp,m,x1,y);
    mxFree(ytmp);
    mxFree(yytmp);
    
}


void splie2(double *x1a, double *x2a, double *ya, int m, int n, double *y2a)
{
	void spline(double *x, double *y, int n, double yp1, double ypn, double *y2);
    int j,i;
    double *yaj,*y2aj;
        

    yaj=(double*)mxCalloc(n,sizeof(double));
	y2aj=(double*)mxCalloc(n,sizeof(double)); 
    
      
    for (j=0;j<m;j++)
    {       
        /*Take one row of ya and y2a*/    
        for (i=0;i<n;i++)
        {
            yaj[i]=ya[j+i*m];
            y2aj[i]=y2a[j+i*m];
        }
        spline(x2a,yaj,n,1.0e30,1.0e30,y2aj);
    }    
    mxFree(yaj);  
    mxFree(y2aj);
}



void spline(double *x, double *y, int n, double yp1, double ypn, double *y2)
{
	int i,k;
	double p,qn,sig,un,*u;  
     
	u=(double*)mxCalloc(n-1,sizeof(double));
    
	if (yp1 > 0.99e30)
		y2[0]=u[0]=0.0;
	else {
		y2[0] = -0.5;
		u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
	}

      
	for (i=1;i<=n-2;i++) {
		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}

	if (ypn > 0.99e30)
		qn=un=0.0;
	else {
		qn=0.5;
		un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
	}
	
      
    y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);

    
	for (k=n-2;k>=0;k--)
		y2[k]=y2[k]*y2[k+1]+u[k];
	mxFree(u);

}


void splint(double xa[], double ya[], double y2a[], int n, double x, double *y)
{
	int klo,khi,k;
	double h,b,a;
    


	klo=1;
	khi=n;
	while (khi-klo > 1) {
		k=(khi+klo) >> 1;
		if (xa[k] > x) khi=k;
		else klo=k;
	}
	h=xa[khi]-xa[klo];
/*	if (h == 0.0) mexErrMsgTxt("Bad xa input to routine splint");*/ /* Result will bw NaN... */
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;
	*y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}
