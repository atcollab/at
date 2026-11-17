#ifndef ATLALIB_H
#define ATLALIB_H
#include "atelem.c"
/*   File: atlalib.c
     Matrix and Vector algebra operations for Accelerator Toolbox
     A.Terebilo   8/14/00


1. Use mxCalloc, mxMalloc , mxFree for memory allocation/deallocation
2. Vector and matrix input arguments are (*double) obtained with mxGetDoubles
    Note: MATLAB internally represents matrixes as column-by-column
   1-dimentional arrays. For example  
                  A B C
                  D E F         is represented as    A D G B E H C F I
                  G H I
3. All matrixes are 6-by-6
4. All vectors  are 1-by-6 or 6-by-1
*/

static void ATmultmv(double *r, const double* A)
/*	multiplies 6-component column vector r by 6x6 matrix R: as in A*r 
  The result is store in the memory area of r !!!
*/

{   int i,j;
    double temp[6];

	for(i=0;i<6;i++)
     	{	temp[i]=0;
			for(j=0;j<6;j++)
        		temp[i]+=A[i+j*6]*r[j];
		}
  	for(i=0;i<6;i++)
	r[i]=temp[i];
} 

static void ATmultmv55(double *r, const double* A)
/*	multiplies 5-component column vector r by 5x5 matrix R: as in A*r 
    The result is store in the memory area of r !!!
*/ 
{   int i,j;
    double temp[5];

	for(i=0;i<5;i++)
     	{	temp[i]=0;
			for(j=0;j<5;j++)
        		temp[i]+=A[i+j*5]*r[j];
		}
  	for(i=0;i<5;i++)
	r[i]=temp[i];
} 

static void ATaddvv(double *r, const double *dr)
{	/*	Add two 6-component vectors vectors.
		The result is store in the memory area of r !!!
    */
	int i;
	for(i=0;i<6;i++)
		r[i]+=dr[i];
}     	

static void ATChangePRef(double* r, double scaling)
{
    r[1] /= scaling;
    r[3] /= scaling;
    r[4] = (r[4] + 1.0 - scaling) / scaling;
}

static void ATtranspm(double *M)
{	/* Transpose matrix M	
	   The result replaces the original matrix 
	*/
	int i,j;
	double temp;
	for(i=1;i<6;i++)
		for(j=0;j<i;j++)
			{	temp = M[i+j*6];
				M[i+j*6] = M[j+i*6];
				M[j+i*6] = temp;
			}
}

static void ATmultmm(const double *M2 , double *M1)
{	/* Mutrix multiplication M2*M1, the result is stored in the M1 memory area */
	int i,j,k;
	double column_temp[6];


	for (i=0;i<6;i++) 
		{	for (j=0;j<6;j++)
				{	column_temp[j] = 0;
					for (k=0;k<6;k++)
						column_temp[j] += M2[j+6*k]*M1[k+6*i];
				}
			for (j=0;j<6;j++)
					M1[j+6*i] = 	column_temp[j];
		}
}

static void ATmultmm55(const double *M2 , double *M1)
{	/* Mutrix multiplication M2*M1, the result is stored in the M1 memory area */
	int i,j,k;
	double column_temp[5];


	for (i=0;i<5;i++) 
		{	for (j=0;j<5;j++)
				{	column_temp[j] = 0;
					for (k=0;k<5;k++)
						column_temp[j] += M2[j+5*k]*M1[k+5*i];
				}
			for (j=0;j<5;j++)
					M1[j+5*i] = 	column_temp[j];
		}
}

static void ATsandwichmmt(const double *M ,double *B)
/* calculates the matrix product M*B*M' (M' = M transposed)
   The result is stored in B memory area
*/
{	int i,j,k;
	double row_temp[6];

	ATmultmm(M,B);
	
	for (i=0;i<6;i++) 
		{	for (j=0;j<6;j++)
				{	row_temp[j] = 0;
					for (k=0;k<6;k++)
						row_temp[j] += B[i+6*k]*M[j+6*k];
				}
			for (j=0;j<6;j++)
					B[i+6*j] = 	row_temp[j];
		}
}

static void ATaddmm(const double *M2 , double *M1)
/* adds two 6-by-6  matrixes  M1, M2 element-by-element
   The result is stored in M1 memory area
*/
{	int i,j;
	for (i=0;i<6;i++) 
		for (j=0;j<6;j++)
			M1[i+6*j] += M2[i+6*j];

}

static void markaslost(double *r6,int idx)
{
    r6[idx] = atGetInf();
}

static void checkiflostRectangularAp(double *r6, const double *limits)
{
	/* check limits for X position */
    if (r6[0]<limits[0] || r6[0]>limits[1])      markaslost(r6,5);
	/* check limits for Y position */
    else if (r6[2]<limits[2] || r6[2]>limits[3]) markaslost(r6,5);
}

static void checkiflostLongtRectangularAp(double *r6, const double *limits)
{
        /* check limits for DP position */
    if (r6[4]<limits[0] || r6[4]>limits[1]) markaslost(r6,5);
        /* check limits for CT position */
    else if (r6[5]<limits[2] || r6[5]>limits[3]) markaslost(r6,5);
}

static void checkiflostEllipticalAp(double *r6, const double *axesptr)
{
	double xnorm = r6[0]/axesptr[0];
	double znorm = r6[2]/axesptr[1];
	if ((xnorm*xnorm + znorm*znorm) >= 1) markaslost(r6,5);
}

/* Cholesky factorisation S = L L^T of a symmetric 6x6 matrix, row-major.
   Returns 0 if S is not positive definite */
static int cholesky6(const double *S, double *L)
{
    int i, j, k;
    for (i = 0; i < 36; i++) L[i] = 0.0;
    for (j = 0; j < 6; j++) {
        double d = S[6*j+j];
        for (k = 0; k < j; k++) d -= L[6*j+k]*L[6*j+k];
        if (!(d > 0.0)) return 0;
        L[6*j+j] = sqrt(d);
        for (i = j+1; i < 6; i++) {
            double v = S[6*i+j];
            for (k = 0; k < j; k++) v -= L[6*i+k]*L[6*j+k];
            L[6*i+j] = v/L[6*j+j];
        }
    }
    return 1;
}

/* One-sided Jacobi (Hestenes) SVD of a 6x6 matrix A, row-major, modified in
   place: on exit the columns of A are U*s, and A_in V = A_out */
static void hestenes6(double *A, double *V, double *s)
{
    int i, p, q, sweep;
    for (i = 0; i < 36; i++) V[i] = (i%7 == 0) ? 1.0 : 0.0;
    for (sweep = 0; sweep < 60; sweep++) {
        int rotated = 0;
        for (p = 0; p < 5; p++) {
            for (q = p+1; q < 6; q++) {
                double alpha = 0.0, beta = 0.0, gamma = 0.0;
                for (i = 0; i < 6; i++) {
                    alpha += A[6*i+p]*A[6*i+p];
                    beta += A[6*i+q]*A[6*i+q];
                    gamma += A[6*i+p]*A[6*i+q];
                }
                if (fabs(gamma) > 1.0e-15*sqrt(alpha*beta)) {
                    double zeta = (beta - alpha)/(2.0*gamma);
                    double t = (zeta >= 0.0 ? 1.0 : -1.0)/(fabs(zeta) + sqrt(1.0 + zeta*zeta));
                    double c = 1.0/sqrt(1.0 + t*t), sn = c*t;
                    rotated = 1;
                    for (i = 0; i < 6; i++) {
                        double ap = A[6*i+p], aq = A[6*i+q];
                        double vp = V[6*i+p], vq = V[6*i+q];
                        A[6*i+p] = c*ap - sn*aq;
                        A[6*i+q] = sn*ap + c*aq;
                        V[6*i+p] = c*vp - sn*vq;
                        V[6*i+q] = sn*vp + c*vq;
                    }
                }
            }
        }
        if (!rotated) break;
    }
    for (p = 0; p < 6; p++) {
        double nrm = 0.0;
        for (i = 0; i < 6; i++) nrm += A[6*i+p]*A[6*i+p];
        s[p] = sqrt(nrm);
    }
}

/* Eigen emittances of a 6x6 sigma matrix (row-major), ordered by plane.
   The singular values of K = L^T J L, with sigma = L L^T, are the
   emittances, each twice. Each mode is assigned to the plane holding the
   largest share of its action. Returns 0 if sigma is not positive definite */
static int eigen_emittances(const double *sigma, double *emit)
{
    /* Symplectic unit matrix with the longitudinal block transposed for (delta, ct) */
    static const double jdiag[6] = {1.0, -1.0, 1.0, -1.0, -1.0, 1.0};
    double L[36], K[36], V[36], s[6], share[9];
    int idx[6], i, j, k, p, best[3], perm[6][3] = {{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
    double bestsum = -1.0;
    if (!cholesky6(sigma, L)) return 0;
    /* K = L^T J L: row i of J L is jdiag[i] times row i^1 of L */
    for (i = 0; i < 6; i++) {
        for (j = 0; j < 6; j++) {
            double v = 0.0;
            for (k = 0; k < 6; k++) v += L[6*k+i]*jdiag[k]*L[6*(k^1)+j];
            K[6*i+j] = v;
        }
    }
    hestenes6(K, V, s);
    /* Sort the singular values: consecutive pairs belong to the same mode */
    for (i = 0; i < 6; i++) idx[i] = i;
    for (i = 1; i < 6; i++) {
        int t = idx[i];
        for (j = i; j > 0 && s[idx[j-1]] > s[t]; j--) idx[j] = idx[j-1];
        idx[j] = t;
    }
    for (k = 0; k < 3; k++) {
        double a[6], b[6], tot = 0.0;
        for (i = 0; i < 6; i++) {
            a[i] = b[i] = 0.0;
            for (j = 0; j < 6; j++) {
                a[i] += L[6*i+j]*V[6*j+idx[2*k]];
                b[i] += L[6*i+j]*V[6*j+idx[2*k+1]];
            }
        }
        for (p = 0; p < 3; p++) {
            share[3*p+k] = fabs(a[2*p]*b[2*p+1] - a[2*p+1]*b[2*p]);
            tot += share[3*p+k];
        }
        for (p = 0; p < 3; p++) share[3*p+k] /= tot;
    }
    for (i = 0; i < 6; i++) {
        double sum = share[3*perm[i][0]+0] + share[3*perm[i][1]+1] + share[3*perm[i][2]+2];
        if (sum > bestsum) {
            bestsum = sum;
            for (k = 0; k < 3; k++) best[k] = perm[i][k];
        }
    }
    for (k = 0; k < 3; k++) emit[best[k]] = 0.5*(s[idx[2*k]] + s[idx[2*k+1]]);
    return 1;
}

#endif /*ATLALIB_H*/

