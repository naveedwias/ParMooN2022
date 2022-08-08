#ifndef _ANDERSON_ACC_
#define _ANDERSON_ACC_
extern "C" {
void anderson_acceleration(int ndim, 
                           int kdim,
                           double *sol_old, 
                           double *fpast,
                           double *delta,
                           double *work,
                           int lwork,
                           int *info_error,
			   double *alphas_x_i) ;
}
#endif
