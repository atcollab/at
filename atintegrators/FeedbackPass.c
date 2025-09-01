/* IdentityPass.c 
   Accelerator Toolbox
   Revision 7/16/03
   A.Terebilo
*/

#include "atelem.c"
#include "atlalib.c"

struct elem 
{
  double GX;
  double GY;
  double GZ;
  double syncZ;
};

void FeedbackPass(double *r_in, double gx, double gy, double gz, double syncZ, int num_particles)
{	
    double *r6;
    int c;
    double mx, my, mz, npart;
    mx=0.0;
    my=0.0;
    mz=0.0;
    npart = 0.0;
    
    for (c = 0; c<num_particles; c++) {	/*Loop over particles  */
        r6 = r_in+c*6;
        if (!atIsNaN(r6[0])) {
            mx += r6[0];
            my += r6[2];
            mz += r6[5]; 
            npart += 1; 
        }
    }
    
    if (npart>0.0){
        mx /= npart;
        my /= npart;
        mz /= npart;    
        for (c = 0; c<num_particles; c++) {	/*Loop over particles  */
            r6 = r_in+c*6;
            if (!atIsNaN(r6[0])) {
                r6[0] -= mx*gx;
                r6[2] -= my*gy;
                r6[5] -= mz*gz;    
            }
        }
    }
}

#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
        double *r_in, int num_particles, struct parameters *Param)
{
    if (!Elem) {
        double GX, GY, GZ;
        GX=atGetDouble(ElemData,"GX"); check_error();
        GY=atGetDouble(ElemData,"GY"); check_error();
        GZ=atGetDouble(ElemData,"GZ"); check_error();
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->GX=GX;
        Elem->GY=GY;
        Elem->GZ=GZ;
    }
    FeedbackPass(r_in,Elem->GX,Elem->GY,Elem->GZ,num_particles);
    return Elem;
}

MODULE_DEF(IdentityPass)        /* Dummy module initialisation */

#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/
