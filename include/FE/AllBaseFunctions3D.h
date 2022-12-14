#include <string.h> // memset

#include <BF_C_T_P00_3D.h>
#include <BF_C_T_P0_3D.h>
#include <BF_C_T_P1_3D.h>
#include <BF_C_T_P2_3D.h>
#include <BF_C_T_P3_3D.h>

#include <BF_C_T_B2_3D.h>

#include <BF_C_H_Q00_3D.h>
#include <BF_C_H_Q0_3D.h>
#include <BF_C_H_Q1_3D.h>
#include <BF_C_H_Q2_3D.h>
#include <BF_C_H_Q3_3D.h>
#include <BF_C_H_Q4_3D.h>

#include <BF_N_H_Q1_3D.h>
#include <BF_N_T_P1_3D.h>

#include <BF_D_T_P1_3D.h>
#include <BF_D_T_P2_3D.h>
#include <BF_D_T_P3_3D.h>

#include <BF_D_H_P1_3D.h>
#include <BF_D_H_P2_3D.h>
#include <BF_D_H_P3_3D.h>

#include <BF_D_H_Q1_3D.h>
#include <BF_D_H_Q2_3D.h>

#include <BF_B_H_IB2_3D.h>

#include <BF_N_T_P2_3D.h>

#include <BF_N_H_Q2_3D.h>
#include <BF_N_H_Q3_3D.h>
#include <BF_N_H_Q4_3D.h>
#include <BF_C_H_UL1_3D.h>
#include <BF_C_H_UL2_3D.h>
#include <BF_C_H_UL3_3D.h>

// Raviart-Thomas (RT) basis functions, vector valued basis functions
#include <BF_N_T_RT0_3D.h>
#include <BF_N_T_RT1_3D.h>
#include <BF_N_T_RT2_3D.h>
#include <BF_N_T_RT3_3D.h>
#include <BF_N_H_RT0_3D.h>
#include <BF_N_H_RT1_3D.h>
#include <BF_N_H_RT2_3D.h>

// Brezzi-Douglas-Marini (BDM) basis functions, vector valued basis functions
#include <BF_N_T_BDDF1_3D.h>
#include <BF_N_T_BDDF2_3D.h>
#include <BF_N_T_BDDF3_3D.h>
#include <BF_N_H_BDDF1_3D.h>
#include <BF_N_H_BDDF2_3D.h>
#include <BF_N_H_BDDF3_3D.h>
