#include "atconstants.h"
#include "atelem.c"
#include <math.h>
#include <float.h>
#include <complex.h>
#ifdef MPI
#include <mpi.h>
#include <mpi4py/mpi4py.h>
#endif

static void IIR_init(double cutoff, double *IIRcoef, double *IIRout, double T1, int every, double Vc){
    
    if(cutoff==0){
        IIRcoef[0] =  1.0;
    }else{
        double omega = TWOPI * cutoff;
        double T = T1 * every;
        double alpha = cos(omega * T) - 1;
        double tmp = alpha*alpha - 2*alpha;
        if(tmp > 0){
            IIRcoef[0] = alpha + sqrt(tmp);
        }else{
            IIRcoef[0] = T * cutoff * TWOPI;
        }
    IIRout[0] = Vc;
    }
}

static void IIR(double complex input, double *IIRcoef, double *IIRout){
    /*

    """Return IIR filter output."""
    self.IIRout = (1 - self.IIRcoef) * self.IIRout + self.IIRcoef * input
    return self.IIRout
    
    */

    
    IIRout[0] = (1 - IIRcoef[0]) * IIRout[0] + IIRcoef[0] * creal(input);
    IIRout[1] = (1 - IIRcoef[0]) * IIRout[1] + IIRcoef[0] * cimag(input);
    


}
        

static double complex Vg2Ig_real(double vgen, double thetag, double psi, double RL){
    /*
    Return Ig from Vg (assuming constant Vg).

    Eq.25 of ref [2] assuming the dVg/dt = 0.
    */
    double complex vgen_phasor = vgen * cexp(_Complex_I * (thetag + TWOPI/4)); // phase shift needed for vgen def
    double complex Ig = (vgen_phasor / RL) * (1 - _Complex_I * tan(psi));
    return creal(Ig);
}

static double complex Vg2Ig_imag(double vgen, double thetag, double psi, double RL){
    /*
    Return Ig from Vg (assuming constant Vg).

    Eq.25 of ref [2] assuming the dVg/dt = 0.
    */
    double complex vgen_phasor = vgen * cexp(_Complex_I * (thetag + TWOPI/4));  // phase shift needed for vgen def
    
    double complex Ig = (vgen_phasor / RL) * (1 - _Complex_I * tan(psi));
    return cimag(Ig);
}


        
        
static void init_Ig2Vg_matrix(int ring_harmn, double *Ig2Vg_vec_real, double *Ig2Vg_vec_imag, double *Ig2Vg_tmp_real, double *Ig2Vg_tmp_imag, double filling_time, double psi, double T1, double *Ig2Vg_mat_real, double *Ig2Vg_mat_imag){
    /*
    Initialize matrix for Ig2Vg_matrix.

    Shoud be called before first use of Ig2Vg_matrix and after each cavity
    parameter change.
    */

    /*
        k = np.arange(0, self.ring.h)
        self.Ig2Vg_vec = np.exp(-1 / self.cav_res.filling_time *
                                (1 - 1j * np.tan(self.cav_res.psi)) *
                                self.ring.T1 * (k+1))
        tempV = np.exp(-1 / self.cav_res.filling_time * self.ring.T1 * k *
                       (1 - 1j * np.tan(self.cav_res.psi)))
        for idx in np.arange(self.ring.h):
           self.Ig2Vg_mat[idx:, idx] = tempV[:self.ring.h - idx]
    */                   

    int idx=0;

    for(idx=0;idx<ring_harmn;idx++){
        Ig2Vg_vec_real[idx] = 0.0;
        Ig2Vg_vec_imag[idx] = 0.0;
        Ig2Vg_tmp_real[idx] = 0.0;
        Ig2Vg_tmp_imag[idx] = 0.0;
    }
    
    for(idx=0;idx<ring_harmn*ring_harmn;idx++){
        Ig2Vg_mat_real[idx] = 0.0;
        Ig2Vg_mat_imag[idx] = 0.0;            
    }
    
    
    int ibucket=0;
    for(ibucket=0;ibucket<ring_harmn;ibucket++){
        double complex tmp1 = cexp(-1.0 / filling_time * (1.0 - _Complex_I * tan(psi)) * T1 * ( (double) ibucket + 1.0));
        Ig2Vg_vec_real[ibucket] = creal(tmp1);
        Ig2Vg_vec_imag[ibucket] = cimag(tmp1);
        
        double complex tmp2 = cexp((-1.0 / filling_time) * (T1 * (double) ibucket * (1.0 - _Complex_I * tan(psi))));
        Ig2Vg_tmp_real[ibucket] = creal(tmp2);
        Ig2Vg_tmp_imag[ibucket] = cimag(tmp2);
    }
   

    
    for(idx=0;idx<ring_harmn;idx++){
        int tempV_start_ind = 0;
        int tempV_end_ind = ring_harmn - idx;
                
        int idtmpv = 0;
        for(idtmpv=tempV_start_ind;idtmpv<tempV_end_ind;idtmpv++){
        
            Ig2Vg_mat_real[idx + idtmpv + ring_harmn*idx] = Ig2Vg_tmp_real[idtmpv];
            Ig2Vg_mat_imag[idx + idtmpv + ring_harmn*idx] = Ig2Vg_tmp_imag[idtmpv];
            
        }
    }        
}

static void init_FFconst(bool FF, double *ig_phasor_real, double *ig_phasor_imag, int ring_harmn, double *FFconst){
    //Initialize feedforward constant
    double FFconst_real=0.0;
    double FFconst_imag=0.0;
    
    int idx=0;
    if(FF){
        for(idx=0;idx<ring_harmn;idx++){
            FFconst_real += ig_phasor_real[idx];
            FFconst_imag += ig_phasor_imag[idx];
        }
        FFconst_real /= ring_harmn;
        FFconst_imag /= ring_harmn;
    }
    FFconst[0] = FFconst_real;
    FFconst[1] = FFconst_imag;
}

static void init_phasor_arrays(double vgen, double thetag, double *ig_phasor_real, double *ig_phasor_imag, double *ig_phasor_record_real, double *ig_phasor_record_imag, int ring_harmn, double RL, double psi, double *generator_phasor_record_real, double *generator_phasor_record_imag){

    double ig_real = Vg2Ig_real(vgen, thetag, psi, RL);
    double ig_imag = Vg2Ig_imag(vgen, thetag, psi, RL);
    
    double Vg_real = -vgen*sin(thetag);
    double Vg_imag = vgen*cos(thetag);
    int idx=0;
    for(idx=0;idx<ring_harmn;idx++){
        ig_phasor_real[idx]=ig_real;
        ig_phasor_imag[idx]=ig_imag;
        ig_phasor_record_real[idx]=ig_real;
        ig_phasor_record_imag[idx]=ig_imag;
        generator_phasor_record_real[idx]=Vg_real;
        generator_phasor_record_imag[idx]=Vg_imag;
    }
    
}


static void set_cavity_phasor(double vgen, double thetag, double *vbeam_phasor, double *vcav_phasor){
    // BE AWARE THAT MBTRACK2 INCLUDES LOSS FACTOR IN ADDITION TO u0. THIS MAY RESULT IN NUMBERS BEING SLIGHTLY DIFFERENT
    double complex generator_phasor = vgen*cexp(_Complex_I*(thetag+TWOPI/4));
    double complex beam_phasor = vbeam_phasor[0]*cexp(_Complex_I*vbeam_phasor[1]);
    double complex cavity_phasor = generator_phasor + beam_phasor;
    
    vcav_phasor[0] = creal(cavity_phasor);
    vcav_phasor[1] = cimag(cavity_phasor);

}


static void mat_dot(double *Ig2Vg_mat_real, double *Ig2Vg_mat_imag, double *ig_phasor_record_real, double *ig_phasor_record_imag, double *dot_output_real, double *dot_output_imag, int ring_harmn){	
    /*
    For whatever stupud reason I have mIg2V to be the inverse matrix in the memory. So the logic
    here is kind of inverted with i and j
    
    The problem here is the fucking complex numbers. 
    */
    int i,j;
    double val_real=0.0;
    double val_imag=0.0;
    double complex tmp=0;
    for (i=0;i<ring_harmn;i++){ 
        val_real = 0.0;
        val_imag = 0.0;
        for(j=0;j<ring_harmn;j++){
            tmp = (Ig2Vg_mat_real[ring_harmn*j + i] + _Complex_I*Ig2Vg_mat_imag[ring_harmn*j + i]) * (ig_phasor_record_real[j] + _Complex_I*ig_phasor_record_imag[j]);
            val_real += creal(tmp);
            val_imag += cimag(tmp);
        }
        dot_output_real[i] = val_real;
        dot_output_imag[i] = val_imag;
        
    }   

}

static void init_vc_previous(double *vc_previous_real, double *vc_previous_imag, int samplenum, double *vcav_phasor){
    /*
    self.vc_previous = np.ones(
    self.sample_num) * self.cav_res.cavity_phasor
    */
    int idx=0;
    for(idx=0;idx<samplenum;idx++){
        vc_previous_real[idx] = vcav_phasor[0];
        vc_previous_imag[idx] = vcav_phasor[1];
    }

}

static void concat_vc_list(double *vc_previous_real, double *vc_previous_imag, double *vc_list_real, double *vc_list_imag, double *cavity_phasor_record_real, double *cavity_phasor_record_imag, int ring_harmn, int samplenum){

    int idx=0;
    for(idx=0;idx<samplenum;idx++){
        vc_list_real[idx] = vc_previous_real[idx];
        vc_list_imag[idx] = vc_previous_imag[idx];
    }
    for(idx=0;idx<ring_harmn;idx++){
        vc_list_real[idx + samplenum] = cavity_phasor_record_real[idx];
        vc_list_imag[idx + samplenum] = cavity_phasor_record_imag[idx];

    }    

} 




static void init_sample_list(int *sample_list, int ring_harmn, int every){
    /*
    self.sample_list = range(0, self.ring.h, self.every)
    
    but then later

    self.sample_list = range(index + self.every - self.ring.h, self.ring.h,
                         self.every)
    */
    int tt=0;
    int idx=0;
    for(idx=0;idx<ring_harmn;idx=idx+every){
        sample_list[tt] = idx;
        tt += 1;
    }
}

void update_sample_list(int *sample_list, int index, int every, int ring_harmn){
    int idx=0;
    int tt = 0;
    for(idx=index+every-ring_harmn;idx<ring_harmn;idx=idx+every){
        sample_list[tt] = idx;
        tt += 1;
    }   
}


static void Ig2Vg_matrix(double *generator_phasor_record_real, double *generator_phasor_record_imag,
                         double *Ig2Vg_vec_real, double *Ig2Vg_vec_imag,
                         double *Ig2Vg_mat_real, double *Ig2Vg_mat_imag,
                         double *ig_phasor_record_real, double *ig_phasor_record_imag,
                         double *dot_output_real, double *dot_output_imag,
                         double kloss, double T1, int ring_harmn){
    /*
    Return Vg from Ig using matrix formalism.
    Warning: self.init_Ig2Vg should be called after each CavityResonator
    parameter change.


    generator_phasor_record = (
        self.Ig2Vg_vec * self.cav_res.generator_phasor_record[-1] +
        self.Ig2Vg_mat.dot(self.ig_phasor_record) *
        kloss * T1)
        
    */
    int idx;
    double complex tmp=0;
    mat_dot(Ig2Vg_mat_real, Ig2Vg_mat_imag, ig_phasor_record_real, ig_phasor_record_imag, dot_output_real, dot_output_imag, ring_harmn);
    
    for(idx=0;idx<ring_harmn;idx++){
        tmp =  (Ig2Vg_vec_real[idx] + _Complex_I*Ig2Vg_vec_imag[idx])*(generator_phasor_record_real[ring_harmn-1] + _Complex_I * generator_phasor_record_imag[ring_harmn-1]);
        generator_phasor_record_real[idx] = creal(tmp) +
                                            dot_output_real[idx] * kloss * T1;       
        generator_phasor_record_imag[idx] = cimag(tmp) +
                                            dot_output_imag[idx] * kloss * T1;     
    }        
}


static void Ig2Vg(double *generator_phasor_record_real, double *generator_phasor_record_imag,
                 double *Ig2Vg_vec_real, double *Ig2Vg_vec_imag,
                 double *Ig2Vg_mat_real, double *Ig2Vg_mat_imag,
                 double *ig_phasor_record_real, double *ig_phasor_record_imag,
                 double *dot_output_real, double *dot_output_imag,
                 double kloss, double T1, int ring_harmn, double *vgen_arr){
    /*
    def Ig2Vg(self):
        """
        Go from Ig to Vg.

        Apply new values to cav_res.generator_phasor_record, cav_res.Vg and
        cav_res.theta_g from ig_phasor_record.
        """
        self.cav_res.generator_phasor_record = self.Ig2Vg_matrix()
        self.cav_res.Vg = np.mean(np.abs(self.cav_res.generator_phasor_record))
        self.cav_res.theta_g = np.mean(
            np.angle(self.cav_res.generator_phasor_record))
    */
    
    Ig2Vg_matrix(generator_phasor_record_real, generator_phasor_record_imag,
                 Ig2Vg_vec_real, Ig2Vg_vec_imag,
                 Ig2Vg_mat_real, Ig2Vg_mat_imag,
                 ig_phasor_record_real, ig_phasor_record_imag,
                 dot_output_real, dot_output_imag,
                 kloss, T1, ring_harmn);
                 
    int idx;
    double mean_vg = 0.0;
    double mean_thetag=0.0;
    for(idx=0;idx<ring_harmn;idx++){
        mean_vg += sqrt(generator_phasor_record_real[idx]*generator_phasor_record_real[idx] + 
                        generator_phasor_record_imag[idx]*generator_phasor_record_imag[idx]);
        mean_thetag += -atan2(generator_phasor_record_real[idx], generator_phasor_record_imag[idx]);
    }

    mean_vg /= ring_harmn;
    mean_thetag /= ring_harmn;
    vgen_arr[0] = mean_vg;
    vgen_arr[1] = mean_thetag;
}



static void track_PIL(double *vc_previous_real, double *vc_previous_imag,
                      double *cavity_phasor_record_real, double *cavity_phasor_record_imag,
                      double *ig_phasor_real, double *ig_phasor_imag,
                      int *sample_list, int samplenum, int record_size, int samplelist_length,
                      double *diff_record_real, double *diff_record_imag,
                      double *FFconst, double *gain, double *I_record,
                      double frf,
                      double Vc, double theta,
                      double *generator_phasor_record_real, double *generator_phasor_record_imag,
                      double *Ig2Vg_vec_real, double *Ig2Vg_vec_imag,
                      double *Ig2Vg_mat_real, double *Ig2Vg_mat_imag,
                      double *ig_phasor_record_real, double *ig_phasor_record_imag,
                      double *dot_output_real, double *dot_output_imag,
                      double kloss, double T1, int ring_harmn, double *vgen_arr,
                      double *IIRout, double *IIRcoef,
                      double *vc_list_real, double *vc_list_imag,
                      int every,
                      double psi, double rshunt
                      ){
    /*
    def track(self, apply_changes: bool = True):
        """
        Tracking method for the Cavity PI control feedback.

        Returns
        -------
        None.

        """
        vc_list = np.concatenate([
            self.vc_previous, self.cav_res.cavity_phasor_record
        ])  #This line is slowing down the process.
        self.ig_phasor.fill(self.ig_phasor[-1])

        for index in self.sample_list:
            # 2) updating Ig using last item of the list
            diff = self.diff_record[-1] - self.FFconst
            self.I_record += diff / self.ring.f1
            fb_value = self.gain[0] * diff + self.gain[1] * self.I_record
            self.ig_phasor[index:] = self.Vg2Ig(fb_value) + self.FFconst
            # Shift the record
            self.diff_record = np.roll(self.diff_record, 1)
            # 1) recording diff as a first item of the list
            mean_vc = np.mean(vc_list[index:self.sample_num + index]) * np.exp(
                -1j * self.cav_res.theta)
            self.diff_record[0] = self.cav_res.Vc - self.IIR(mean_vc)
        # update sample_list for next turn
        self.sample_list = range(index + self.every - self.ring.h, self.ring.h,
                                 self.every)
        # update vc_previous for next turn
        self.vc_previous = self.cav_res.cavity_phasor_record[-self.sample_num:]

        self.ig_phasor = self.ig_phasor + self.ig_modulation_signal
        self.ig_phasor_record = self.ig_phasor

        if apply_changes:
            self.Ig2Vg()
            
    */            

    int idx=0;
    int index=0;
    int index2=0;
    int sample=0;
    double diff_real, diff_imag, fb_value_real, fb_value_imag;
    double fb_amp, fb_phase, fb_real, fbr, fbi;
    double mean_vc_arr[] = {0.0, 0.0};
    double complex mean_vc = 0.0 + _Complex_I * 0.0;
    concat_vc_list(vc_previous_real, vc_previous_imag, vc_list_real, vc_list_imag, cavity_phasor_record_real, cavity_phasor_record_imag, ring_harmn, samplenum);
    
    for(idx=0;idx<samplenum;idx++){
        index = sample_list[idx];
        
        diff_real = diff_record_real[record_size-1] - FFconst[0];
        diff_imag = diff_record_imag[record_size-1] - FFconst[1];
        I_record[0] += diff_real *T1; // f1;
        I_record[1] += diff_imag *T1; // f1;
        fb_value_real = gain[0] * diff_real + gain[1] * I_record[0];
        fb_value_imag = gain[0] * diff_imag + gain[1] * I_record[1];
        fb_amp = sqrt(fb_value_real*fb_value_real + fb_value_imag*fb_value_imag);
        fb_phase = -atan2(fb_value_real, fb_value_imag);

        fbr = Vg2Ig_real(fb_amp, fb_phase, psi, rshunt);
        fbi = Vg2Ig_imag(fb_amp, fb_phase, psi, rshunt);
        for(index2=index;index2<ring_harmn;index2++){
            ig_phasor_real[index2] = fbr + FFconst[0];
            ig_phasor_imag[index2] = fbi + FFconst[1];
            }

        roll_array(diff_record_real, record_size);
        roll_array(diff_record_imag, record_size);


        compute_mean_vc(vc_list_real, vc_list_imag, mean_vc_arr, idx, samplenum);
        
        
        mean_vc = (mean_vc_arr[0] + _Complex_I * mean_vc_arr[1])*cexp(-_Complex_I * (theta+TWOPI/4));
        
        
        IIR(mean_vc, IIRcoef, IIRout);

        diff_record_real[0] = Vc - IIRout[0];
        diff_record_imag[0] =  IIRout[1];
        };

        update_sample_list(sample_list, index, every, ring_harmn);
        
        
        update_vc_previous(vc_previous_real, vc_previous_imag,
                           samplenum, ring_harmn, 
                           cavity_phasor_record_real, cavity_phasor_record_imag);
                           
        update_ig_phasor(ig_phasor_real, ig_phasor_imag, ig_phasor_record_real, ig_phasor_record_imag, ring_harmn);
                           
        int apply_changes=1;
        //if(apply_changes==1){                           
        Ig2Vg(generator_phasor_record_real, generator_phasor_record_imag,
                     Ig2Vg_vec_real, Ig2Vg_vec_imag,
                     Ig2Vg_mat_real, Ig2Vg_mat_imag,
                     ig_phasor_record_real, ig_phasor_record_imag,
                     dot_output_real, dot_output_imag,
                     kloss, T1, ring_harmn, vgen_arr);
                         
        //};    
    };

void update_ig_phasor(double *ig_phasor_real, double *ig_phasor_imag, double *ig_phasor_record_real, double *ig_phasor_record_imag, int ring_harmn){
    int idx=0;
    for(idx=0;idx<ring_harmn;idx++){
        ig_phasor_record_real[idx] = ig_phasor_real[idx];
        ig_phasor_record_imag[idx] = ig_phasor_imag[idx];
    }    

}
void update_vc_previous(double *vc_previous_real, double *vc_previous_imag, int samplenum, int ring_harmn, double *cavity_phasor_record_real, double *cavity_phasor_record_imag){
    /*
    
    */
    int idx=0;
    for(idx=0;idx<samplenum;idx++){
        vc_previous_real[idx] = cavity_phasor_record_real[ring_harmn-samplenum+idx];
        vc_previous_imag[idx] = cavity_phasor_record_imag[ring_harmn-samplenum+idx];
    }
    
}


void compute_mean_vc(double *vc_list_real, double *vc_list_imag, double *vc_mean, int index, int samplenum){
    int idx=0;
    for(idx=index;idx<index+samplenum;idx++){
        vc_mean[0] += vc_list_real[idx]/samplenum;
        vc_mean[1] += vc_list_imag[idx]/samplenum;
    }
}




void roll_array(double *arr, int arr_len){


    if(arr_len>1){
        int idx = 0;
        double tmp=0.0;
        double tmp2=0.0;    
        tmp = arr[arr_len-1];
       
        for(idx=0;idx<arr_len;idx++){
            tmp2 = arr[idx];
            arr[idx] = tmp;
            tmp = tmp2;
        }
    }
    
}


