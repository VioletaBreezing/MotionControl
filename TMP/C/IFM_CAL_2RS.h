#ifndef _IFM_CAL_2RS_H_
#define _IFM_CAL_2RS_H_

int model_init_interface(double* phi, double* psi, 
                        double *a_x, double *a_y, double *b, 
                        double lambda,
                        int len, int max_iter);
int model_run_interface(double *h_meas, double *pos_est, int len);

#endif
