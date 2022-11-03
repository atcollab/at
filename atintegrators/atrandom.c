#include <math.h>

void init_seed(seed){
    static int initseed=1;
    if(initseed)
    {
        srand(seed);
        initseed=0;
    }
}

double atdrand(void)   /* uniform distribution, (0..1] */
{
    return (rand()+1.0)/(RAND_MAX+1.0);
}

double atrandn() /* gaussian distribution, sigma=1, mean=0 */
{
	static bool hasSpare = false;
	static double spare;
	static double u, v, s;

	if (hasSpare) {
		hasSpare = false;
		return spare;
	}

	hasSpare = true;	
	do {
		u = (rand() / (double) RAND_MAX) * 2.0 - 1.0;
		v = (rand() / (double) RAND_MAX) * 2.0 - 1.0;
		s = u * u + v * v;
	}
	while( (s >= 1.0) || (s == 0.0) );
	s = sqrt(-2.0 * log(s) / s);
	spare = v * s;
	return u * s;
}

int atrandp(double lamb) /* poisson distribution */
{
    
    int pk,k;
    double l,p,u,u1,u2,twopi;
    
    l = -lamb;
    k = 0;
    p = 0.0;
    twopi = 4.0*asin(1.0);
    pk=0;
    
    if(lamb<11){
        
        do{
            k = k + 1;
            u = atdrand();
            p = p + log(u);
        }while(p>l);
        
        pk = k-1;
        
    }else{
        u1 = atdrand();
        u2 = atdrand();
        if(u1==0.0){
            u1 = 1e-18;
        };
        pk = (int)floor((sqrt(-2.0*log(u1))*cos(twopi*u2))*sqrt(lamb)+lamb);
    }
    
    return pk;    
}
