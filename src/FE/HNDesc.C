#include "HNDesc.h"
#include "MooNMD_Io.h"

THNDesc::THNDesc(HNDesc t)
 : type(t), coeff()
{
  switch(type)
  {
    case HNDesc::HN_C_P1_2D_0:
      coeff = { 0.5, 0.5 };
      break;
    case HNDesc::HN_C_P2_2D_0:
      coeff = { -0.125, 0.75, 0.375 };
      break;
    case HNDesc::HN_C_P2_2D_1:
      coeff = { 0.375, 0.75, -0.125 };
      break;
    case HNDesc::HN_C_P3_2D_0:
      coeff = {  0.0625, -0.3125,  0.9375,  0.3125 };
      break;
    case HNDesc::HN_C_P3_2D_1:
      coeff = { -0.0625,  0.5625,  0.5625, -0.0625 };
      break;
    case HNDesc::HN_C_P3_2D_2:
      coeff = {  0.3125,  0.9375, -0.3125,  0.0625 };
      break;
    case HNDesc::HN_C_P4_2D_0:
      coeff = {  -0.0390625, 0.21875, -0.546875, 1.09375, 0.2734375 };
      break;
    case HNDesc::HN_C_P4_2D_1:
      coeff = { 0.0234375, -0.15625, 0.703125, 0.46875, -0.0390625 };
      break;
    case HNDesc::HN_C_P4_2D_2:
      coeff = { -0.0390625, 0.46875, 0.703125, -0.15625, 0.0234375 };
      break;
    case HNDesc::HN_C_P4_2D_3:
      coeff = { 0.2734375, 1.09375, -0.546875, 0.21875, -0.0390625 };
      break;
    case HNDesc::HN_C_P5_2D_0:
      coeff = { 7./256., -45./256., 63./128., -105./128., 315./256., 63./256. };
      break;
    case HNDesc::HN_C_P5_2D_1:
      coeff = { -3./256., 21./256., -35./128., 105./128., 105./256., -7./256. };
      break;
    case HNDesc::HN_C_P5_2D_2:
      coeff = { 3./256., -25./256., 75./128., 75./128., -25./256., 3./256. };
      break;
    case HNDesc::HN_C_P5_2D_3:
      coeff = { -7./256., 105./256., 105./128., -35./128., 21./256., -3./256. };
      break;
    case HNDesc::HN_C_P5_2D_4:
      coeff = { 63./256., 315./256., -105./128., 63./128., -45./256., 7./256. };
      break;
    case HNDesc::HN_C_P6_2D_0:
    case HNDesc::HN_C_P6_2D_1:
    case HNDesc::HN_C_P6_2D_2:
    case HNDesc::HN_C_P6_2D_3:
    case HNDesc::HN_C_P6_2D_4:
    case HNDesc::HN_C_P6_2D_5:
    case HNDesc::HN_C_P7_2D_0:
    case HNDesc::HN_C_P7_2D_1:
    case HNDesc::HN_C_P7_2D_2:
    case HNDesc::HN_C_P7_2D_3:
    case HNDesc::HN_C_P7_2D_4:
    case HNDesc::HN_C_P7_2D_5:
    case HNDesc::HN_C_P7_2D_6:
    case HNDesc::HN_C_P8_2D_0:
    case HNDesc::HN_C_P8_2D_1:
    case HNDesc::HN_C_P8_2D_2:
    case HNDesc::HN_C_P8_2D_3:
    case HNDesc::HN_C_P8_2D_4:
    case HNDesc::HN_C_P8_2D_5:
    case HNDesc::HN_C_P8_2D_6:
    case HNDesc::HN_C_P8_2D_7:
    case HNDesc::HN_C_P9_2D_0:
    case HNDesc::HN_C_P9_2D_1:
    case HNDesc::HN_C_P9_2D_2:
    case HNDesc::HN_C_P9_2D_3:
    case HNDesc::HN_C_P9_2D_4:
    case HNDesc::HN_C_P9_2D_5:
    case HNDesc::HN_C_P9_2D_6:
    case HNDesc::HN_C_P9_2D_7:
    case HNDesc::HN_C_P9_2D_8:
      ErrThrow("Hanging node descriptor of type ", type,
               " is not yet implemented");
      break;
    case HNDesc::HN_N_P1_2D_0:
      coeff = { 0.5, 0.5 };
      break;
    case HNDesc::HN_N_P2_2D_0:
      coeff = { 0.75, -0.25, -0.75, -0.25 };
      break;
    case HNDesc::HN_N_P3_2D_0:
      coeff = { 0, -0.625, 0.125, 0,  0.625, 0.125 };
      break;
    case HNDesc::HN_N_P4_2D_0:
      coeff = { -0.4375, -0.4375,  0.4375, -0.0625,
                 0.4375, -0.4375, -0.4375, -0.0625 };
      break;
    case HNDesc::HN_N_P5_2D_0:
      coeff = { 0,  0.1875, 0.5625, -0.28125, 0.03125,
                0, -0.1875, 0.5625,  0.28125, 0.03125 };
      break;
    case HNDesc::HN_C_Q1_3D_E:
      coeff = { 0.5, 0.5 };
      break;
    case HNDesc::HN_C_Q1_3D_F:
      coeff = { 0.25, 0.25, 0.25, 0.25 };
      break;
    case HNDesc::HN_C_Q2_3D_E:
      coeff = { 0.375, 0.75, -0.125 };
      break;
    case HNDesc::HN_C_Q2_3D_F:
      coeff = { 0.140625, 0.28125, -0.046875, 0.28125, 0.56250, -0.093750,
               -0.046875, -0.093750, 0.015625 };
      break;
    case HNDesc::HN_C_Q3_3D_1:
      coeff = { 0.3125, 0.9375, -0.3125, 0.0625 };
      break;
    case HNDesc::HN_C_Q3_3D_2:
      coeff = { -0.0625, 0.5625, 0.5625, -0.0625 };
      break;
    case HNDesc::HN_C_Q3_3D_3:
      coeff = { -.01953125000, .1757812500, .1757812500, -.01953125000,
                -.05859375000, .5273437500, .5273437500, -.05859375000,
                 .01953125000, -.1757812500, -.1757812500, .01953125000,
                -.003906250000, .03515625000, .03515625000, -.003906250000 };
      break;
    case HNDesc::HN_C_Q3_3D_4:
      coeff = { .09765625000, .2929687500, -.09765625000, .01953125000,
                .2929687500, .8789062500, -.2929687500, .05859375000,
               -.09765625000, -.2929687500, .09765625000, -.01953125000,
                .01953125000, .05859375000, -.01953125000, .003906250000 };
      break;
    case HNDesc::HN_C_Q3_3D_5:
      coeff = { .003906250000, -.03515625000, -.03515625000, .003906250000,
               -.03515625000, .3164062500, .3164062500, -.03515625000,
               -.03515625000, .3164062500, .3164062500, -.03515625000,
                .003906250000, -.03515625000, -.03515625000, .003906250000 };
      break;
    case HNDesc::HN_C_P1_3D_E:
      coeff = { 0.5, 0.5 };
      break;
    case HNDesc::HN_C_P2_3D_E:
      coeff = { 0.375, 0.75, -0.125 };
      break;
    case HNDesc::HN_C_P2_3D_F:
      coeff = { 0.5, -0.125, 0.5, 0.25, -0.125 };
      break;
    case HNDesc::HN_C_P3_3D_E:
      coeff = { 0.3125, 0.9375, -0.3125, 0.0625 };
      break;
    case HNDesc::HN_C_P3_3D_M:
      coeff = { -0.0625, 0.5625, 0.5625, -0.0625 };
      break;
    case HNDesc::HN_C_P3_3D_F:
      coeff = { 0.5, 0.5, -0.25, 0.5, -0.25, 0.0625, -0.0625, -0.0625, 0.0625 };
      break;
    case HNDesc::HN_C_P3_3D_G:
      coeff = { -0.0625, 0.375, 0, 0, 0.1875, 0.75, 0, -0.1875, -0.125, 0.0625};
      break;
    case HNDesc::HN_N_P1_3D_E:
      coeff = { 0.25, 0.25, 0.25, 0.25 };
      break;
    case HNDesc::HN_N_P2_3D_0:
      coeff = { 0.25,  0.125, 0.125, 0,     0,     0.125,
                0,     0.125, 0,     0.125, 0,     0.125 };
      break;
    case HNDesc::HN_N_P2_3D_1:
      coeff = { 0,     0,     0.125, 0,     0.125, 0,
                0.25,  0.125, 0.125, 0,     0.125, 0.125 };
      break;
    case HNDesc::HN_N_P2_3D_2:
      coeff = { 0,     0.125, 0,     0.25,  0.125, 0.125,
                0,     0,     0.125, 0.125, 0.125, 0 };
      break;
    default:
      ErrThrow("unknown Hanging node descriptor type ", type);
  }
}
