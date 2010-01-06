#include "mex.h"
#include "elempass.h"
#include "atlalib.c"
#include <stdlib.h>
#include <math.h>

/******************************************************************************/
/* PHYSICS SECTION ************************************************************/

void bl11wiggler (double *r)
{	/* K - is the quadrupole strength defined as
	   (e/Eo)(dBz/dx) [1/m^2] 
	   another notation: g0 [DESY paper]
	*/
	double x, px, y, py, ct, dp;

	
	x  = r[0];
	px = r[1];
	y  = r[2];
	py = r[3];
    dp = r[4];
    ct = r[5];
	
r[0]  = -7.0901531e-9*dp + 2.1285135e-8*pow(dp,2) - 4.2582822e-8*pow(dp,3) + 
   7.0986884e-8*pow(dp,4) - 1.0650226e-7*pow(dp,5) + 0.17500251*px - 
   0.17500917*dp*px + 0.17501998000000002*pow(dp,2)*px - 
   0.17503495*pow(dp,3)*px + 0.17505407*pow(dp,4)*px + 
   0.00046615681999999994*pow(px,2) - 0.0013991186*dp*pow(px,2) + 
   0.0028002125*pow(dp,2)*pow(px,2) - 
   0.0046714454*pow(dp,3)*pow(px,2) + 
   0.00009415890799999999*pow(px,3) - 0.00051283918*dp*pow(px,3) + 
   0.0016227193*pow(dp,2)*pow(px,3) + 0.0026325028*pow(px,4) - 
   0.013178106*dp*pow(px,4) + 0.0017739152000000001*pow(px,5) - 
   0.0035151788*pow(py,2) + 0.010545058000000001*dp*pow(py,2) - 
   0.021089225*pow(dp,2)*pow(py,2) + 
   0.035147337*pow(dp,3)*pow(py,2) - 0.00017496961*px*pow(py,2) + 
   0.00099137193*dp*px*pow(py,2) - 
   0.0032073346*pow(dp,2)*px*pow(py,2) - 
   0.013165211000000001*pow(px,2)*pow(py,2) + 
   0.065902871*dp*pow(px,2)*pow(py,2) - 
   0.015448497*pow(px,3)*pow(py,2) + 0.0028251989*pow(py,4) - 
   0.014140680000000001*dp*pow(py,4) + 0.0053600538*px*pow(py,4) + 
   1.0000463*x - 0.00009252638599999999*dp*x + 0.00013881863*pow(dp,2)*x - 
   0.00018515365*pow(dp,3)*x + 0.00023154902*pow(dp,4)*x + 
   0.010655106000000001*px*x - 0.021325115000000002*dp*px*x + 
   0.032025542*pow(dp,2)*px*x - 0.042771913*pow(dp,3)*px*x + 
   0.0027780014999999997*pow(px,2)*x - 
   0.012957015999999998*dp*pow(px,2)*x + 
   0.035163931*pow(dp,2)*pow(px,2)*x + 0.27363095*pow(px,3)*x - 
   1.0952354*dp*pow(px,3)*x + 0.076715563*pow(px,4)*x - 
   0.0022329564999999997*pow(py,2)*x + 0.010693406*dp*pow(py,2)*x - 
   0.029377954*pow(dp,2)*pow(py,2)*x - 0.75876294*px*pow(py,2)*x + 
   3.0366792*dp*px*pow(py,2)*x - 0.39326769*pow(px,2)*pow(py,2)*x + 
   0.047007114*pow(py,4)*x + 6.4772815e-7*pow(x,2) - 
   0.000090698616*dp*pow(x,2) + 0.00027017976*pow(dp,2)*pow(x,2) - 
   0.0005391435999999999*pow(dp,3)*pow(x,2) + 0.031424381*px*pow(x,2) - 
   0.11535896999999999*dp*px*pow(x,2) + 
   0.25184769*pow(dp,2)*px*pow(x,2) + 6.0046889*pow(px,2)*pow(x,2) - 
   18.026116000000002*dp*pow(px,2)*pow(x,2) + 
   1.4565705*pow(px,3)*pow(x,2) - 5.735196*pow(py,2)*pow(x,2) + 
   17.214564000000003*dp*pow(py,2)*pow(x,2) - 
   3.7239512*px*pow(py,2)*pow(x,2) + 0.19265772*pow(x,3) - 
   0.38537434*dp*pow(x,3) + 0.57825596*pow(dp,2)*pow(x,3) + 
   45.749988*px*pow(x,3) - 91.590025*dp*px*pow(x,3) + 
   14.753778*pow(px,2)*pow(x,3) - 12.919905*pow(py,2)*pow(x,3) + 
   0.0014447703999999999*pow(x,4) - 0.25241305999999997*dp*pow(x,4) + 
   81.52341*px*pow(x,4) + 323.09572000000003*pow(x,5) - 
   0.045491972*py*y + 0.09096779299999999*dp*py*y - 
   0.13642689*pow(dp,2)*py*y + 0.1818687*pow(dp,3)*py*y - 
   0.0023456997*px*py*y + 0.012905401*dp*px*py*y - 
   0.037550616*pow(dp,2)*px*py*y - 0.77470123*pow(px,2)*py*y + 
   3.1003242*dp*pow(px,2)*py*y - 0.26313788*pow(px,3)*py*y + 
   0.22718526*pow(py,3)*y - 0.9093159099999999*dp*pow(py,3)*y + 
   0.17637334*px*pow(py,3)*y - 0.056483421*py*x*y + 0.21064962*dp*py*x*y - 
   0.46253504999999995*pow(dp,2)*py*x*y - 22.939951999999998*px*py*x*y + 
   68.853719*dp*px*py*x*y - 7.2708048*pow(px,2)*py*x*y + 
   1.735616*pow(py,3)*x*y - 134.14494*py*pow(x,2)*y + 
   268.48339*dp*py*pow(x,2)*y - 76.281661*px*py*pow(x,2)*y - 
   323.92665*py*pow(x,3)*y + 0.000084597284*pow(y,2) - 
   0.00017130694*dp*pow(y,2) + 
   0.00026010924000000003*pow(dp,2)*pow(y,2) - 
   0.00035096928*pow(dp,3)*pow(y,2) - 0.012019481*px*pow(y,2) + 
   0.05618572*dp*px*pow(y,2) - 
   0.13251253999999998*pow(dp,2)*px*pow(y,2) - 
   5.7342371*pow(px,2)*pow(y,2) + 17.209136*dp*pow(px,2)*pow(y,2) - 
   1.1354968*pow(px,3)*pow(y,2) + 5.0787701*pow(py,2)*pow(y,2) - 
   15.244260999999998*dp*pow(py,2)*pow(y,2) + 
   2.3110962*px*pow(py,2)*pow(y,2) - 0.53704976*x*pow(y,2) + 
   1.0740506*dp*x*pow(y,2) - 1.6111701*pow(dp,2)*x*pow(y,2) - 
   131.06749*px*x*pow(y,2) + 262.28827*dp*px*x*pow(y,2) - 
   34.633437*pow(px,2)*x*pow(y,2) + 25.753281*pow(py,2)*x*pow(y,2) + 
   0.27035345*pow(x,2)*pow(y,2) + 
   0.7088587900000001*dp*pow(x,2)*pow(y,2) - 
   421.93696*px*pow(x,2)*pow(y,2) - 3299.4019*pow(x,3)*pow(y,2) + 
   41.188846999999996*py*pow(y,3) - 82.423315*dp*py*pow(y,3) + 
   13.108201*px*py*pow(y,3) + 202.72644999999997*py*x*pow(y,3) - 
   0.059987627*pow(y,4) - 0.018133504*dp*pow(y,4) + 
   28.305623999999998*px*pow(y,4) + 1290.3802*x*pow(y,4);

r[1] = -8.1082608e-8*dp + 1.6223047e-7*pow(dp,2) - 2.4340624e-7*pow(dp,3) + 
   3.2462451e-7*pow(dp,4) - 4.0589979e-7*pow(dp,5) + 1.0000463*px - 
   0.000092525523*dp*px + 0.00013881517*pow(dp,2)*px - 
   0.00018514501000000002*pow(dp,3)*px + 
   0.00023153174000000001*pow(dp,4)*px - 
   0.0053276862999999995*pow(px,2) + 0.010647557*dp*pow(px,2) - 
   0.015951968*pow(dp,2)*pow(px,2) + 
   0.021233266999999997*pow(dp,3)*pow(px,2) + 0.0031789209*pow(px,3) - 
   0.010872446000000001*dp*pow(px,3) + 
   0.024417363*pow(dp,2)*pow(px,3) - 0.59382684*pow(px,4) + 
   2.3751982*dp*pow(px,4) + 0.1367672*pow(px,5) + 
   0.0053249784*pow(py,2) - 0.010637325*dp*pow(py,2) + 
   0.015927265*pow(dp,2)*pow(py,2) - 
   0.021185021000000002*pow(dp,3)*pow(py,2) - 
   0.008851347399999999*px*pow(py,2) + 
   0.029955970999999998*dp*px*pow(py,2) - 
   0.066716897*pow(dp,2)*px*pow(py,2) + 
   3.4840881*pow(px,2)*pow(py,2) - 13.933735*dp*pow(px,2)*pow(py,2) - 
   1.3931154*pow(px,3)*pow(py,2) - 0.56738827*pow(py,4) + 
   2.2687176*dp*pow(py,4) + 0.5408207199999999*px*pow(py,4) + 
   0.00052868308*x - 0.00052877001*dp*x + 0.00052905238*pow(dp,2)*x - 
   0.00052953024*pow(dp,3)*x + 0.00053020368*pow(dp,4)*x - 
   4.3376602999999996e-6*px*x - 0.00016450052*dp*px*x + 
   0.0005065690400000001*pow(dp,2)*px*x - 
   0.0010219708000000001*pow(dp,3)*px*x + 0.069722739*pow(px,2)*x - 
   0.18811301*dp*pow(px,2)*x + 0.35521111*pow(dp,2)*pow(px,2)*x - 
   12.009954*pow(px,3)*x + 36.024276*dp*pow(px,3)*x + 
   4.24886*pow(px,4)*x - 0.06611057299999999*pow(py,2)*x + 
   0.17726217*dp*pow(py,2)*x - 0.33347997*pow(dp,2)*pow(py,2)*x + 
   35.482701*px*pow(py,2)*x - 106.40817*dp*px*pow(py,2)*x - 
   26.303508*pow(px,2)*pow(py,2)*x + 3.5025797*pow(py,4)*x - 
   0.000024786842*pow(x,2) - 0.0009647889799999999*dp*pow(x,2) + 
   0.0019546764*pow(dp,2)*pow(x,2) - 
   0.0029451521*pow(dp,3)*pow(x,2) + 0.57799368*px*pow(x,2) - 
   1.1561609*dp*px*pow(x,2) + 1.7347961*pow(dp,2)*px*pow(x,2) - 
   68.636259*pow(px,2)*pow(x,2) + 137.16966*dp*pow(px,2)*pow(x,2) + 
   56.152606*pow(px,3)*pow(x,2) + 
   68.60583799999999*pow(py,2)*pow(x,2) - 
   137.03676000000002*dp*pow(py,2)*pow(x,2) - 
   176.14078*px*pow(py,2)*pow(x,2) + 2.2018812*pow(x,3) - 
   2.2025445*dp*pow(x,3) + 2.2043319*pow(dp,2)*pow(x,3) - 
   0.090243942*px*pow(x,3) - 0.72161976*dp*px*pow(x,3) + 
   402.41303999999997*pow(px,2)*pow(x,3) - 
   426.74719999999996*pow(py,2)*pow(x,3) - 
   0.25576783000000003*pow(x,4) - 2.3237497*dp*pow(x,4) + 
   1615.6956*px*pow(x,4) + 3692.9987*pow(x,5) - 0.000043323744*py*y + 
   0.00029489547*dp*py*y - 0.00075470007*pow(dp,2)*py*y + 
   0.0014227437*pow(dp,3)*py*y - 0.13145644*px*py*y + 
   0.3531415*dp*px*py*y - 0.66506515*pow(dp,2)*px*py*y + 
   34.939597*pow(px,2)*py*y - 104.765*dp*pow(px,2)*py*y - 
   19.081098*pow(px,3)*py*y - 11.459814999999999*pow(py,3)*y + 
   34.351481*dp*pow(py,3)*y + 13.993312*px*pow(py,3)*y - 
   1.0739636*py*x*y + 2.147873*dp*py*x*y - 3.2222425*pow(dp,2)*py*x*y + 
   268.19402*px*py*x*y - 535.51935*dp*px*py*x*y - 
   359.86807999999996*pow(px,2)*py*x*y + 
   91.20047299999999*pow(py,3)*x*y - 0.48587206*py*pow(x,2)*y + 
   4.4974555*dp*py*pow(x,2)*y - 2492.1132*px*py*pow(x,2)*y - 
   6597.9197*py*pow(x,3)*y - 0.000058737291999999996*pow(y,2) + 
   0.0010598295*dp*pow(y,2) - 
   0.0020606928999999997*pow(dp,2)*pow(y,2) + 
   0.0030613761*pow(dp,3)*pow(y,2) - 0.53690123*px*pow(y,2) + 
   1.0735984*dp*px*pow(y,2) - 1.6102522*pow(dp,2)*px*pow(y,2) + 
   65.518783*pow(px,2)*pow(y,2) - 130.81834*dp*pow(px,2)*pow(y,2) - 
   67.776979*pow(px,3)*pow(y,2) - 65.44291*pow(py,2)*pow(y,2) + 
   130.57113*dp*pow(py,2)*pow(y,2) + 
   145.52013*px*pow(py,2)*pow(y,2) - 6.1370529*x*pow(y,2) + 
   6.1363592*dp*x*pow(y,2) - 6.1381035*pow(dp,2)*x*pow(y,2) - 
   0.3828955*px*x*pow(y,2) + 4.1880728*dp*px*x*pow(y,2) - 
   1309.7892*pow(px,2)*x*pow(y,2) + 986.3035*pow(py,2)*x*pow(y,2) - 
   0.76002119*pow(x,2)*pow(y,2) + 18.889076*dp*pow(x,2)*pow(y,2) - 
   9895.699799999999*px*pow(x,2)*pow(y,2) - 
   37703.386*pow(x,3)*pow(y,2) + 0.29617093*py*pow(y,3) - 
   1.6805354*dp*py*pow(y,3) + 700.04992*px*py*pow(y,3) + 
   5159.2203*py*x*pow(y,3) - 0.29108761*pow(y,4) - 
   1.6811974*dp*pow(y,4) + 1289.4713*px*pow(y,4) + 
   14743.345000000001*x*pow(y,4);

r[2]  = 0.17497891999999998*py - 0.17493511*dp*py + 0.17486859*pow(dp,2)*py - 
   0.17477934*pow(dp,3)*py + 0.17466737999999998*pow(dp,4)*py - 
   0.00093220474*px*py + 0.0027976051000000004*dp*px*py - 
   0.0055983496*pow(dp,2)*px*py + 0.009337744799999999*pow(dp,3)*px*py - 
   0.00026269109*pow(px,2)*py + 0.0014444269999999999*dp*pow(px,2)*py - 
   0.0045954409*pow(dp,2)*pow(px,2)*py - 0.013323805*pow(px,3)*py + 
   0.066673783*dp*pow(px,3)*py - 0.0084400798*pow(px,4)*py - 
   0.00011204161*pow(py,3) + 0.0005196842400000001*dp*pow(py,3) - 
   0.0014778735000000002*pow(dp,2)*pow(py,3) + 
   0.017670001*px*pow(py,3) - 0.088403375*dp*px*pow(py,3) + 
   0.012871718*pow(px,2)*pow(py,3) - 0.0037092953*pow(py,5) - 
   0.010653862*py*x + 0.021319136*dp*py*x - 0.032009046*pow(dp,2)*py*x + 
   0.04273682*pow(dp,3)*py*x - 0.0049344296*px*py*x + 
   0.023426268*dp*px*py*x - 0.064103649*pow(dp,2)*px*py*x - 
   0.8527300200000001*pow(px,2)*py*x + 3.4127259*dp*pow(px,2)*py*x - 
   0.28913160000000004*pow(px,3)*py*x + 0.26354119*pow(py,3)*x - 
   1.0546669*dp*pow(py,3)*x + 0.22765327000000002*px*pow(py,3)*x - 
   0.027872518*py*pow(x,2) + 0.10469499*dp*py*pow(x,2) - 
   0.23049649*pow(dp,2)*py*pow(x,2) - 12.0078*px*py*pow(x,2) + 
   36.043349*dp*px*py*pow(x,2) - 4.0788048*pow(px,2)*py*pow(x,2) + 
   1.1012564*pow(py,3)*pow(x,2) - 45.744058*py*pow(x,3) + 
   91.56251400000001*dp*py*pow(x,3) - 27.340848*px*py*pow(x,3) - 
   75.332887*py*pow(x,4) + 0.99965487*y + 0.00069022618*dp*y - 
   0.0010352849*pow(dp,2)*y + 0.0013803022*pow(dp,3)*y - 
   0.0017252738*pow(dp,4)*y - 0.045503999999999996*px*y + 
   0.091015903*dp*px*y - 0.13654716*pow(dp,2)*px*y + 
   0.18210923*pow(dp,3)*px*y - 0.0027340192*pow(px,2)*y + 
   0.012779815000000002*dp*pow(px,2)*y - 
   0.034716912999999995*pow(dp,2)*pow(px,2)*y - 
   0.28954486*pow(px,3)*y + 1.1587379*dp*pow(px,3)*y - 
   0.07502729999999999*pow(px,4)*y - 0.0030851985*pow(py,2)*y + 
   0.010577695*dp*pow(py,2)*y - 0.023795324*pow(dp,2)*pow(py,2)*y + 
   0.9112131800000001*px*pow(py,2)*y - 3.6464184*dp*px*pow(py,2)*y + 
   0.31924628*pow(px,2)*pow(py,2)*y - 0.14545475000000002*pow(py,4)*y - 
   1.3039179e-6*x*y + 0.00016885406*dp*x*y - 0.00050266142*pow(dp,2)*x*y + 
   0.0010027663*pow(dp,3)*x*y - 0.056483681*px*x*y + 0.21065643*dp*px*x*y - 
   0.4625665*pow(dp,2)*px*x*y - 12.007387000000001*pow(px,2)*x*y + 
   36.041650000000004*dp*pow(px,2)*x*y - 2.7582384*pow(px,3)*x*y + 
   11.468661*pow(py,2)*x*y - 34.419985*dp*pow(py,2)*x*y + 
   6.2751134*px*pow(py,2)*x*y - 0.53704669*pow(x,2)*y + 
   1.0741246*dp*pow(x,2)*y - 1.6114845*pow(dp,2)*pow(x,2)*y - 
   134.16289*px*pow(x,2)*y + 268.55516*dp*px*pow(x,2)*y - 
   42.213002*pow(px,2)*pow(x,2)*y + 33.432414*pow(py,2)*pow(x,2)*y - 
   0.0070129182*pow(x,3)*y + 1.0343714*dp*pow(x,3)*y - 
   323.92934*px*pow(x,3)*y - 1649.6827*pow(x,4)*y - 
   0.036644354*py*pow(y,2) + 0.08978636699999999*dp*py*pow(y,2) - 
   0.15938364*pow(dp,2)*py*pow(y,2) + 11.468684*px*py*pow(y,2) - 
   34.419019*dp*px*py*pow(y,2) + 2.715799*pow(px,2)*py*pow(y,2) - 
   2.4393897*pow(py,3)*pow(y,2) + 131.05073000000002*py*x*pow(y,2) - 
   262.21781000000004*dp*py*x*pow(y,2) + 
   59.850387999999995*px*py*x*pow(y,2) + 
   368.13113*py*pow(x,2)*pow(y,2) - 0.17956436*pow(y,3) + 
   0.3590617*dp*pow(y,3) - 0.5384129200000001*pow(dp,2)*pow(y,3) + 
   41.179999*px*pow(y,3) - 82.387929*dp*px*pow(y,3) + 
   7.4329976*pow(px,2)*pow(y,3) - 21.927567*pow(py,2)*pow(y,3) - 
   0.1987844*x*pow(y,3) - 0.19618978*dp*x*pow(y,3) + 
   202.71693000000002*px*x*pow(y,3) + 2580.616*pow(x,2)*pow(y,3) - 
   116.16086999999999*py*pow(y,4) - 449.55066999999997*pow(y,5);

r[3] = 0.99965487*py + 0.00069022249*dp*py - 0.0010352702*pow(dp,2)*py + 
   0.0013802653*pow(dp,3)*py - 0.0017252*pow(dp,4)*py + 
   0.045505877*px*py - 0.090994545*dp*px*py + 0.13645047*pow(dp,2)*px*py - 
   0.1818581*pow(dp,3)*px*py - 0.0095389879*pow(px,2)*py + 
   0.032789357*dp*pow(px,2)*py - 0.073924605*pow(dp,2)*pow(px,2)*py + 
   2.2971365*pow(px,3)*py - 9.187553*dp*pow(px,3)*py - 
   0.7111575099999999*pow(px,4)*py - 0.0021715683*pow(py,3) + 
   0.010446210000000001*dp*pow(py,3) - 
   0.028752935*pow(dp,2)*pow(py,3) - 2.2337413*px*pow(py,3) + 
   8.9326549*dp*px*pow(py,3) + 1.2348332*pow(px,2)*pow(py,3) - 
   0.1664766*pow(py,5) + 0.000022753218999999996*py*x + 
   0.000096685312*dp*py*x - 0.00035832168000000003*pow(dp,2)*py*x + 
   0.00076218377*pow(dp,3)*py*x - 0.13148671*px*py*x + 
   0.35325976*dp*px*py*x - 0.66535371*pow(dp,2)*px*py*x + 
   34.953711*pow(px,2)*py*x - 104.83670000000001*dp*pow(px,2)*py*x - 
   17.792740000000002*pow(px,3)*py*x - 11.467856000000001*pow(py,3)*x + 
   34.390508*dp*pow(py,3)*x + 15.320924*px*pow(py,3)*x - 
   0.5370676600000001*py*pow(x,2) + 1.0741296*dp*py*pow(x,2) - 
   1.6113788*pow(dp,2)*py*pow(x,2) + 134.20481*px*py*pow(x,2) - 
   268.19083*dp*px*py*pow(x,2) - 176.2727*pow(px,2)*py*pow(x,2) + 
   49.43295499999999*pow(py,3)*pow(x,2) + 0.16411103*py*pow(x,3) + 
   0.52105131*dp*py*pow(x,3) - 830.94073*px*py*pow(x,3) - 
   1649.9142*py*pow(x,4) - 0.0039441358*y + 0.0039433574*dp*y - 
   0.0039423435*pow(dp,2)*y + 0.0039410943*pow(dp,3)*y - 
   0.0039396097*pow(dp,4)*y - 0.000010290338000000001*px*y + 
   0.00019580724*dp*px*y - 0.0005565474799999999*pow(dp,2)*px*y + 
   0.001092529*pow(dp,3)*px*y - 0.081955682*pow(px,2)*y + 
   0.22572643*dp*pow(px,2)*y - 0.43131816*pow(dp,2)*pow(px,2)*y + 
   11.470407*pow(px,3)*y - 34.397608000000005*dp*pow(px,3)*y - 
   5.1032118*pow(px,4)*y - 0.057622416999999995*pow(py,2)*y + 
   0.19294975*dp*pow(py,2)*y - 0.4058942*pow(dp,2)*pow(py,2)*y - 
   33.077385*px*pow(py,2)*y + 99.16680299999999*dp*px*pow(py,2)*y + 
   24.337009*pow(px,2)*pow(py,2)*y - 5.4433578*pow(py,4)*y - 
   0.00011760492*x*y + 0.0021201803*dp*x*y - 
   0.004122689099999999*pow(dp,2)*x*y + 0.0061253586*pow(dp,3)*x*y - 
   1.0739939*px*x*y + 2.1478846*dp*px*x*y - 3.2221072*pow(dp,2)*px*x*y + 
   131.08196999999998*pow(px,2)*x*y - 261.8142*dp*pow(px,2)*x*y - 
   128.22902*pow(px,3)*x*y - 131.02506*pow(py,2)*x*y + 
   261.6991*dp*pow(py,2)*x*y + 298.59988999999996*px*pow(py,2)*x*y - 
   6.1370527*pow(x,2)*y + 6.1363588*dp*pow(x,2)*y - 
   6.1381029*pow(dp,2)*pow(x,2)*y - 0.12925899000000002*px*pow(x,2)*y + 
   3.4277395*dp*px*pow(x,2)*y - 1310.0878*pow(px,2)*pow(x,2)*y + 
   986.49396*pow(py,2)*pow(x,2)*y - 0.5073468*pow(x,3)*y + 
   12.595381000000001*dp*pow(x,3)*y - 6598.2427*px*pow(x,3)*y - 
   18851.692*pow(x,4)*y - 0.5386129900000001*py*pow(y,2) + 
   1.0767589*dp*py*pow(y,2) - 1.6140145*pow(dp,2)*py*pow(y,2) - 
   123.44973999999999*px*py*pow(y,2) + 
   246.38680000000002*dp*px*py*pow(y,2) + 
   167.77439999999999*pow(px,2)*py*pow(y,2) - 
   75.480215*pow(py,3)*pow(y,2) + 0.092986515*py*x*pow(y,2) - 
   2.6562255*dp*py*x*pow(y,2) + 2100.7015*px*py*x*pow(y,2) + 
   7740.153600000001*py*pow(x,2)*pow(y,2) - 2.0518361*pow(y,3) + 
   2.0502438*dp*pow(y,3) - 2.0472498*pow(dp,2)*pow(y,3) + 
   0.15401796*px*pow(y,3) - 1.2541835*dp*px*pow(y,3) + 
   394.76833999999997*pow(px,2)*pow(y,3) - 
   554.03261*pow(py,2)*pow(y,3) - 1.1628486*x*pow(y,3) - 
   6.7307959*dp*x*pow(y,3) + 5159.2078*px*x*pow(y,3) + 
   29486.688000000002*pow(x,2)*pow(y,3) - 2246.4465*py*pow(y,4) - 
   5135.579199999999*pow(y,5);

r[5] = 1.*ct + 2.6274183e-6*dp - 3.9414962e-6*pow(dp,2) + 
   5.255346099999999e-6*pow(dp,3) - 6.569132299999999e-6*pow(dp,4) + 
   7.8827085e-6*pow(dp,5) + 7.0991885e-9*px - 2.8398793000000002e-8*dp*px + 
   7.1005537e-8*pow(dp,2)*px - 1.4203513e-7*pow(dp,3)*px + 
   2.4861574e-7*pow(dp,4)*px - 0.087500538*pow(px,2) + 
   0.17500378*dp*pow(px,2) - 0.26251191*pow(dp,2)*pow(px,2) + 
   0.35002708*pow(dp,3)*pow(px,2) - 0.00046588626*pow(px,3) + 
   0.0018634463*dp*pow(px,3) - 0.0046585726*pow(dp,2)*pow(px,3) - 
   0.000068537895*pow(px,4) + 0.00045322792999999997*dp*pow(px,4) - 
   0.0026287146*pow(px,5) - 0.087497755*pow(py,2) + 
   0.17498939*dp*pow(py,2) - 0.26247102*pow(dp,2)*pow(py,2) + 
   0.34993879*pow(dp,3)*pow(py,2) + 0.004447171799999999*px*pow(py,2) - 
   0.017788253*dp*px*pow(py,2) + 
   0.044470036000000004*pow(dp,2)*px*pow(py,2) + 
   0.00028792687*pow(px,2)*pow(py,2) - 
   0.0019747266*dp*pow(px,2)*pow(py,2) + 
   0.026466668000000002*pow(px,3)*pow(py,2) - 
   0.000032735306*pow(py,4) + 0.00023687572000000002*dp*pow(py,4) - 
   0.020486848000000002*px*pow(py,4) + 8.1691474e-8*x - 3.2488387e-7*dp*x + 
   7.3058138e-7*pow(dp,2)*x - 1.2987792e-6*pow(dp,3)*x + 
   2.0296544e-6*pow(dp,4)*x + 6.278515900000001e-9*px*x - 
   5.7121893e-8*dp*px*x + 2.2283613e-7*pow(dp,2)*px*x - 
   6.057685499999999e-7*pow(dp,3)*px*x - 
   0.010647044999999999*pow(px,2)*x + 0.031932439*dp*pow(px,2)*x - 
   0.063847455*pow(dp,2)*pow(px,2)*x - 0.0014425215*pow(px,3)*x + 
   0.0090555348*dp*pow(px,3)*x - 0.27344884999999997*pow(px,4)*x + 
   0.010646765*pow(py,2)*x - 0.031930595*dp*pow(py,2)*x + 
   0.063840982*pow(dp,2)*pow(py,2)*x + 0.0037615684*px*pow(py,2)*x - 
   0.024250787*dp*px*pow(py,2)*x + 1.6105983*pow(px,2)*pow(py,2)*x - 
   0.26349095*pow(py,4)*x + 0.00026437278*pow(x,2) - 
   0.00052900349*dp*pow(x,2) + 0.00079417303*pow(dp,2)*pow(x,2) - 
   0.0010601626*pow(dp,3)*pow(x,2) + 0.000083770635*px*pow(x,2) - 
   0.0003351236*dp*px*pow(x,2) + 0.00083796904*pow(dp,2)*px*pow(x,2) - 
   0.007094425099999999*pow(px,2)*pow(x,2) + 
   0.049428858*dp*pow(px,2)*pow(x,2) - 6.0009871*pow(px,3)*pow(x,2) + 
   0.0053408891*pow(py,2)*pow(x,2) - 
   0.042401124*dp*pow(py,2)*pow(x,2) + 
   17.732113000000002*px*pow(py,2)*pow(x,2) + 0.00032159699*pow(x,3) - 
   0.0013030584*dp*pow(x,3) + 0.0029448491*pow(dp,2)*pow(x,3) + 
   0.000035780088*px*pow(x,3) - 0.00030294937000000004*dp*px*pow(x,3) - 
   45.711375*pow(px,2)*pow(x,3) + 
   45.702132999999996*pow(py,2)*pow(x,3) + 0.55053427*pow(x,4) - 
   1.1017585*dp*pow(x,4) + 0.20187993999999998*px*pow(x,4) + 
   0.4647604*pow(x,5) - 2.0037392000000002e-8*py*y + 1.2884328e-7*dp*py*y - 
   4.4325671e-7*pow(dp,2)*py*y + 1.1286815e-6*pow(dp,3)*py*y + 
   0.09098868700000001*px*py*y - 0.27292937*dp*px*py*y + 
   0.54577507*pow(dp,2)*px*py*y + 0.0009044242000000001*pow(px,2)*py*y - 
   0.009880756*dp*pow(px,2)*py*y + 1.0636488*pow(px,3)*py*y - 
   0.0008466107900000001*pow(py,3)*y + 0.0059932824*dp*pow(py,3)*y - 
   1.1378777*px*pow(py,3)*y - 0.00018421122*py*x*y + 
   0.0007369218*dp*py*x*y - 0.0018426446999999999*pow(dp,2)*py*x*y + 
   0.022640595*px*py*x*y - 0.17280176*dp*px*py*x*y + 
   34.924864*pow(px,2)*py*x*y - 11.465421*pow(py,3)*x*y - 
   0.00020294750000000001*py*pow(x,2)*y + 0.0015193151*dp*py*pow(x,2)*y + 
   268.03727*px*py*pow(x,2)*y - 1.0950493*py*pow(x,3)*y - 
   0.0019723584999999997*pow(y,2) + 0.0039450655*dp*pow(y,2) - 
   0.005918447300000001*pow(dp,2)*pow(y,2) + 
   0.0078928312*pow(dp,3)*pow(y,2) - 0.00017735312*px*pow(y,2) + 
   0.00070943234*dp*px*pow(y,2) - 0.0017736968*pow(dp,2)*px*pow(y,2) - 
   0.018908593*pow(px,2)*pow(y,2) + 
   0.055548143*dp*pow(px,2)*pow(y,2) + 5.7304808*pow(px,3)*pow(y,2) - 
   0.0022227324*pow(py,2)*pow(y,2) + 
   0.028995814*dp*pow(py,2)*pow(y,2) - 
   16.537766*px*pow(py,2)*pow(y,2) - 0.0010600452*x*pow(y,2) + 
   0.004123027600000001*dp*x*pow(y,2) - 
   0.00919018*pow(dp,2)*x*pow(y,2) - 0.00017866282*px*x*pow(y,2) + 
   0.0014840844*dp*px*x*pow(y,2) + 130.93131*pow(px,2)*x*pow(y,2) - 
   130.98424*pow(py,2)*x*pow(y,2) - 3.0700139*pow(x,2)*pow(y,2) + 
   6.1454405*dp*pow(x,2)*pow(y,2) - 1.9233416*px*pow(x,2)*pow(y,2) - 
   6.2972657*pow(x,3)*pow(y,2) + 0.000012639294999999999*py*pow(y,3) - 
   0.00018420748*dp*py*pow(y,3) - 82.30040000000001*px*py*pow(y,3) + 
   0.78783182*py*x*pow(y,3) - 0.51326924*pow(y,4) + 
   1.0264577*dp*pow(y,4) + 0.20718673999999998*px*pow(y,4) + 
   1.6828279*x*pow(y,4);  
    
    
}

void BL11WolskiPass(double *r, int num_particles)
{	int c;
	double *r6;
	
	for(c = 0;c<num_particles;c++)
		{	r6 = r+c*6;
            
		    if(!mxIsNaN(r6[0]) & mxIsFinite(r6[4]))

		    {   
                bl11wiggler(r6);
			
			}
			
		}		
}

/********** END PHYSICS SECTION ***********************************************/
/******************************************************************************/

/********** WINDOWS DLL GATEWAY SECTION ***************************************/


ExportMode int* passFunction(const mxArray *ElemData, int *FieldNumbers,
								double *r_in, int num_particles, int mode)

#define NUM_FIELDS_2_REMEMBER 6

{	

    switch(mode)
		{	case NO_LOCAL_COPY:	/* Obsolete in AT1.3 et fields by names from MATLAB workspace  */
		        {	
				
				}	break;	
				
			case MAKE_LOCAL_COPY: 	/* Find field numbers first
										Save a list of field number in an array
										and make returnptr point to that array
								    */
				{						
					
				}	break;

			case	USE_LOCAL_COPY:	/* Get fields from MATLAB using field numbers
										The second argument ponter to the array of field 
										numbers is previously created with 
										QuadLinPass( ..., MAKE_LOCAL_COPY)
								    */
											
				{	

		
				}	break;

		}

	BL11WolskiPass (r_in, num_particles);
	return(NULL);	
}


/********** END WINDOWS DLL GATEWAY SECTION ***************************************/
/********** MATLAB GATEWAY  ***************************************/

void mexFunction(	int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{	
	int m, n;
	double *r_in;   

	

	if(nrhs)
	{
	
	/* ALLOCATE memory for the output array of the same size as the input */
	m = mxGetM(prhs[1]);
	n = mxGetN(prhs[1]);
	if(m!=6) 
		{mexErrMsgTxt("Second argument must be a 6 x N matrix");}	
	
	


	plhs[0] = mxDuplicateArray(prhs[1]);
	r_in = mxGetPr(plhs[0]);
	BL11WolskiPass (r_in, n);
	}
	else
	{   /* return list of required fields */
	    plhs[0] = mxCreateCellMatrix(0,0);
	   
	    
	    if(nlhs>1) /* Required and optional fields */ 
	    {   plhs[1] = mxCreateCellMatrix(0,0);
	        
	    }
	}
}


