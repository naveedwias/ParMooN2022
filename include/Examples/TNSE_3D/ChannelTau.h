// Time Navier-Stokes 3D problem for turbulent channel flow used for the paper:
// An assessment of two classes of variational multiscale methods for
// the simulation of incompressible turbulent flows, N. Ahmed and V. John, 2020
// (to be run with ChannelTauRoutines)
// In contrast to the literature, the coordinates y and z are interchanged!!!
// The geometry can be obtained with the Gmsh file ChannelTau.geo,
// do NOT use the parameter refinement_n_initial_steps, otherwise the mesh
// progression in the Z-direction will not be repected. Instead generate a new
// finer mesh with the Gmsh file ChannelTau.geo.
// The expected dimensions are:
//   - for RE = 180: X in [-2*M_PI ; 2*M_PI],
//                   Y in [-2*M_PI/3 ; 2*M_PI/3],
//                   Z in [0 ; 2]
//   - for RE = 395 or 590: X in [-M_PI ; M_PI],
//                          Y in [-M_PI/2 ; M_PI/2],
//                          Z in [0 ; 2]
//
// the initial setup is as in V. Gravemeier, J. Comput. Phys. (2006)
// with 10% random noise (positive and negative)

#ifndef __CHANNEL_TAU__
#define __CHANNEL_TAU__

#include <ChannelTauRoutines.h>


namespace channel_tau
{
  double RE; // reynolds number
  double INTERNAL_BULK_EXPECTED;
  double INTERNAL_BULK_SIMULATION;

  const double Z_MIN = 0.;
  const double Z_MAX = 2.;
  const double tol = 1e-6;
  const double noise = 0.1;

  // ===========================================================================
  // example file
  // ===========================================================================
  void ExampleFile(const ParameterDatabase& user_input_parameter_db)
  {
    if(user_input_parameter_db.contains("refinement_n_initial_steps"))
    {
      if(!user_input_parameter_db["refinement_n_initial_steps"].is(0))
      {
        ErrThrow("Example ChannelTau.h: do not use the parameter "
                 "'refinement_n_initial_steps' but provide a refined mesh, "
                 "otherwise the mesh progression in the Z-direction will not "
                 "be repected");
      }
    }

    double x_t, y_t;
    if(RE==180)
    {
      x_t = 2*M_PI;
      y_t = 2*M_PI/3;
      INTERNAL_BULK_EXPECTED = 15.6803;
      INTERNAL_BULK_SIMULATION = 15.6803;
    }
    else if(RE==395)
    {
      x_t = M_PI;
      y_t = M_PI/2;
      INTERNAL_BULK_EXPECTED = 17.5452;
      INTERNAL_BULK_SIMULATION = 17.5452;
    }
    else if(RE==590)
    {
      x_t = M_PI;
      y_t = M_PI/2;
      INTERNAL_BULK_EXPECTED = 18.6544005283276;
      INTERNAL_BULK_SIMULATION = 18.6544005283276;
    }
    else
    {
      ErrThrow("Example ChannelTau.h is only defined for Re = 180, 395 or 590");
    }

    // export values to the class TNS_ChannelTau
    TNS_ChannelTau::set_bulkVelocityExpected(INTERNAL_BULK_EXPECTED);
    TNS_ChannelTau::set_bulkVelocity_sim(INTERNAL_BULK_SIMULATION);

    int rank = 0;
#ifdef _MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    if(rank==0)
    {
      Output::print("Example: ChannelTau.h, Re_tau = ", RE ,
                    ", the expected domain dimensions are:"
                    "\n  X_min = ", -x_t, ", X_max = ", x_t,
                    "\n  Y_min = ", -y_t, ", Y_max = ", y_t,
                    "\n  Z_min = 0, Z_max = 2,\nwith periodic boundaries in "
                    "X and Y-direction.");
    }
  }

  // ===========================================================================
  // exact solution RE = 180
  // ===========================================================================
  double DNS_profile_180(double zz)
  {
    double val;
    const int nb_values = 65;
    double value = 0;

    double z[nb_values] =
    {
      0.0000e-00, 3.0118e-04, 1.2045e-03, 2.7095e-03, 4.8153e-03, 7.5205e-03,
      1.0823e-02, 1.4722e-02, 1.9215e-02, 2.4298e-02, 2.9969e-02, 3.6224e-02,
      4.3060e-02, 5.0472e-02, 5.8456e-02, 6.7007e-02, 7.6120e-02, 8.5790e-02,
      9.6011e-02, 1.0678e-01, 1.1808e-01, 1.2991e-01, 1.4227e-01, 1.5515e-01,
      1.6853e-01, 1.8242e-01, 1.9679e-01, 2.1165e-01, 2.2699e-01, 2.4279e-01,
      2.5905e-01, 2.7575e-01, 2.9289e-01, 3.1046e-01, 3.2844e-01, 3.4683e-01,
      3.6561e-01, 3.8477e-01, 4.0430e-01, 4.2419e-01, 4.4443e-01, 4.6500e-01,
      4.8590e-01, 5.0710e-01, 5.2860e-01, 5.5039e-01, 5.7244e-01, 5.9476e-01,
      6.1732e-01, 6.4010e-01, 6.6311e-01, 6.8632e-01, 7.0972e-01, 7.3329e-01,
      7.5702e-01, 7.8090e-01, 8.0491e-01, 8.2904e-01, 8.5327e-01, 8.7759e-01,
      9.0198e-01, 9.2644e-01, 9.5093e-01, 9.7546e-01, 1.0000e-00
    };

    double  Umean[nb_values] =
    {
      0.0000e+00, 5.3639e-02, 2.1443e-01, 4.8197e-01, 8.5555e-01, 1.3339e+00,
      1.9148e+00, 2.5939e+00, 3.3632e+00, 4.2095e+00, 5.1133e+00, 6.0493e+00,
      6.9892e+00, 7.9052e+00, 8.7741e+00, 9.5790e+00, 1.0311e+01, 1.0967e+01,
      1.1550e+01, 1.2066e+01, 1.2520e+01, 1.2921e+01, 1.3276e+01, 1.3590e+01,
      1.3870e+01, 1.4121e+01, 1.4349e+01, 1.4557e+01, 1.4750e+01, 1.4931e+01,
      1.5101e+01, 1.5264e+01, 1.5419e+01, 1.5569e+01, 1.5714e+01, 1.5855e+01,
      1.5993e+01, 1.6128e+01, 1.6260e+01, 1.6389e+01, 1.6515e+01, 1.6637e+01,
      1.6756e+01, 1.6872e+01, 1.6985e+01, 1.7094e+01, 1.7200e+01, 1.7302e+01,
      1.7400e+01, 1.7494e+01, 1.7585e+01, 1.7672e+01, 1.7756e+01, 1.7835e+01,
      1.7911e+01, 1.7981e+01, 1.8045e+01, 1.8103e+01, 1.8154e+01, 1.8198e+01,
      1.8235e+01, 1.8264e+01, 1.8285e+01, 1.8297e+01, 1.8301e+01
    };

    if(zz<=1)
    {
      val = zz;
    }
    else
    {
      val = 2-zz;
    }
    for(int i=0 ; i<nb_values-1 ; i++)
    {
      if(val>=z[i] && val<=z[i+1])
      {
        value = val*(Umean[i+1]-Umean[i])/(z[i+1]-z[i])
              + (Umean[i]*z[i+1]-Umean[i+1]*z[i])/(z[i+1]-z[i]);
      }
    }
    return(value);
  }

  // ===========================================================================
  // exact solution RE = 395
  // ===========================================================================
  double DNS_profile_395(double zz)
  {
    double val;
    const int nb_values = 129;
    double value = 0;

    double z[nb_values] =
    {
      0.0000e-00 ,7.5298e-05 ,3.0118e-04 ,6.7762e-04 ,1.2045e-03 ,1.8819e-03,
      2.7095e-03 ,3.6874e-03 ,4.8153e-03 ,6.0930e-03 ,7.5205e-03 ,9.0974e-03,
      1.0823e-02 ,1.2699e-02 ,1.4722e-02 ,1.6895e-02 ,1.9215e-02 ,2.1683e-02,
      2.4298e-02 ,2.7060e-02 ,2.9969e-02 ,3.3024e-02 ,3.6224e-02 ,3.9569e-02,
      4.3060e-02 ,4.6694e-02 ,5.0472e-02 ,5.4393e-02 ,5.8456e-02 ,6.2661e-02,
      6.7007e-02 ,7.1494e-02 ,7.6120e-02 ,8.0886e-02 ,8.5790e-02 ,9.0832e-02,
      9.6011e-02 ,1.0133e-01 ,1.0678e-01 ,1.1236e-01 ,1.1808e-01 ,1.2393e-01,
      1.2991e-01 ,1.3603e-01 ,1.4227e-01 ,1.4864e-01 ,1.5515e-01 ,1.6178e-01,
      1.6853e-01 ,1.7541e-01 ,1.8242e-01 ,1.8954e-01 ,1.9679e-01 ,2.0416e-01,
      2.1165e-01 ,2.1926e-01 ,2.2699e-01 ,2.3483e-01 ,2.4279e-01 ,2.5086e-01,
      2.5905e-01 ,2.6735e-01 ,2.7575e-01 ,2.8427e-01 ,2.9289e-01 ,3.0162e-01,
      3.1046e-01 ,3.1940e-01 ,3.2844e-01 ,3.3758e-01 ,3.4683e-01 ,3.5617e-01,
      3.6561e-01 ,3.7514e-01 ,3.8477e-01 ,3.9449e-01 ,4.0430e-01 ,4.1420e-01,
      4.2419e-01 ,4.3427e-01 ,4.4443e-01 ,4.5467e-01 ,4.6500e-01 ,4.7541e-01,
      4.8590e-01 ,4.9646e-01 ,5.0710e-01 ,5.1782e-01 ,5.2860e-01 ,5.3946e-01,
      5.5039e-01 ,5.6138e-01 ,5.7244e-01 ,5.8357e-01 ,5.9476e-01 ,6.0601e-01,
      6.1732e-01 ,6.2868e-01 ,6.4010e-01 ,6.5158e-01 ,6.6311e-01 ,6.7469e-01,
      6.8632e-01 ,6.9799e-01 ,7.0972e-01 ,7.2148e-01 ,7.3329e-01 ,7.4513e-01,
      7.5702e-01 ,7.6894e-01 ,7.8090e-01 ,7.9289e-01 ,8.0491e-01 ,8.1696e-01,
      8.2904e-01 ,8.4114e-01 ,8.5327e-01 ,8.6542e-01 ,8.7759e-01 ,8.8978e-01,
      9.0198e-01 ,9.1420e-01 ,9.2644e-01 ,9.3868e-01 ,9.5093e-01 ,9.6319e-01,
      9.7546e-01 ,9.8773e-01 ,1.0000e-00
    };

    double  Umean[nb_values] =
    {
      0.0000e+00 ,2.9538e-02 ,1.1811e-01 ,2.6562e-01 ,4.7198e-01 ,7.3701e-01,
      1.0605e+00 ,1.4417e+00 ,1.8799e+00 ,2.3729e+00 ,2.9177e+00 ,3.5093e+00,
      4.1409e+00 ,4.8032e+00 ,5.4854e+00 ,6.1754e+00 ,6.8611e+00 ,7.5309e+00,
      8.1754e+00 ,8.7870e+00 ,9.3607e+00 ,9.8937e+00 ,1.0385e+01 ,1.0836e+01,
      1.1248e+01 ,1.1624e+01 ,1.1966e+01 ,1.2278e+01 ,1.2563e+01 ,1.2822e+01,
      1.3060e+01 ,1.3278e+01 ,1.3479e+01 ,1.3664e+01 ,1.3837e+01 ,1.3998e+01,
      1.4148e+01 ,1.4290e+01 ,1.4425e+01 ,1.4552e+01 ,1.4673e+01 ,1.4790e+01,
      1.4902e+01 ,1.5011e+01 ,1.5117e+01 ,1.5221e+01 ,1.5322e+01 ,1.5421e+01,
      1.5518e+01 ,1.5614e+01 ,1.5707e+01 ,1.5799e+01 ,1.5890e+01 ,1.5979e+01,
      1.6067e+01 ,1.6153e+01 ,1.6239e+01 ,1.6324e+01 ,1.6409e+01 ,1.6493e+01,
      1.6576e+01 ,1.6659e+01 ,1.6741e+01 ,1.6823e+01 ,1.6903e+01 ,1.6984e+01,
      1.7063e+01 ,1.7141e+01 ,1.7218e+01 ,1.7294e+01 ,1.7369e+01 ,1.7443e+01,
      1.7517e+01 ,1.7590e+01 ,1.7664e+01 ,1.7738e+01 ,1.7812e+01 ,1.7886e+01,
      1.7960e+01 ,1.8034e+01 ,1.8108e+01 ,1.8182e+01 ,1.8254e+01 ,1.8326e+01,
      1.8396e+01 ,1.8466e+01 ,1.8535e+01 ,1.8603e+01 ,1.8669e+01 ,1.8734e+01,
      1.8797e+01 ,1.8859e+01 ,1.8919e+01 ,1.8978e+01 ,1.9035e+01 ,1.9092e+01,
      1.9148e+01 ,1.9202e+01 ,1.9256e+01 ,1.9308e+01 ,1.9359e+01 ,1.9408e+01,
      1.9456e+01 ,1.9503e+01 ,1.9548e+01 ,1.9593e+01 ,1.9636e+01 ,1.9678e+01,
      1.9719e+01 ,1.9758e+01 ,1.9796e+01 ,1.9832e+01 ,1.9865e+01 ,1.9897e+01,
      1.9927e+01 ,1.9955e+01 ,1.9981e+01 ,2.0004e+01 ,2.0026e+01 ,2.0046e+01,
      2.0064e+01 ,2.0080e+01 ,2.0094e+01 ,2.0106e+01 ,2.0116e+01 ,2.0123e+01,
      2.0129e+01 ,2.0132e+01 ,2.0133e+01
    };

    if(zz<=1)
    {
      val = zz;
    }
    else
    {
      val = 2-zz;
    }
    for(int i=0 ; i<nb_values-1 ; i++)
    {
      if(val>=z[i] && val<=z[i+1])
      {
        value = val*(Umean[i+1]-Umean[i])/(z[i+1]-z[i])
              + (Umean[i]*z[i+1]-Umean[i+1]*z[i])/(z[i+1]-z[i]);
      }
    }
    return(value);
  }

  // ===========================================================================
  // exact solution RE = 590
  // ===========================================================================
  double DNS_profile_590(double zz)
  {
    double val;
    const int nb_values = 129;
    double value = 0;

    double z[nb_values] =
    {
      0.0000e+00, 7.5298e-05, 3.0118e-04, 6.7762e-04, 1.2045e-03, 1.8819e-03,
      2.7095e-03, 3.6874e-03, 4.8153e-03, 6.0930e-03, 7.5205e-03, 9.0974e-03,
      1.0823e-02, 1.2699e-02, 1.4722e-02, 1.6895e-02, 1.9215e-02, 2.1683e-02,
      2.4298e-02, 2.7060e-02, 2.9969e-02, 3.3024e-02, 3.6224e-02, 3.9569e-02,
      4.3060e-02, 4.6694e-02, 5.0472e-02, 5.4393e-02, 5.8456e-02, 6.2661e-02,
      6.7007e-02, 7.1494e-02, 7.6120e-02, 8.0886e-02, 8.5790e-02, 9.0832e-02,
      9.6011e-02, 1.0133e-01, 1.0678e-01, 1.1236e-01, 1.1808e-01, 1.2393e-01,
      1.2991e-01, 1.3603e-01, 1.4227e-01, 1.4864e-01, 1.5515e-01, 1.6178e-01,
      1.6853e-01, 1.7541e-01, 1.8242e-01, 1.8954e-01, 1.9679e-01, 2.0416e-01,
      2.1165e-01, 2.1926e-01, 2.2699e-01, 2.3483e-01, 2.4279e-01, 2.5086e-01,
      2.5905e-01, 2.6735e-01, 2.7575e-01, 2.8427e-01, 2.9289e-01, 3.0162e-01,
      3.1046e-01, 3.1940e-01, 3.2844e-01, 3.3758e-01, 3.4683e-01, 3.5617e-01,
      3.6561e-01, 3.7514e-01, 3.8477e-01, 3.9449e-01, 4.0430e-01, 4.1420e-01,
      4.2419e-01, 4.3427e-01, 4.4443e-01, 4.5467e-01, 4.6500e-01, 4.7541e-01,
      4.8590e-01, 4.9646e-01, 5.0710e-01, 5.1782e-01, 5.2860e-01, 5.3946e-01,
      5.5039e-01, 5.6138e-01, 5.7244e-01, 5.8357e-01, 5.9476e-01, 6.0601e-01,
      6.1732e-01, 6.2868e-01, 6.4010e-01, 6.5158e-01, 6.6311e-01, 6.7469e-01,
      6.8632e-01, 6.9799e-01, 7.0972e-01, 7.2148e-01, 7.3329e-01, 7.4513e-01,
      7.5702e-01, 7.6894e-01, 7.8090e-01, 7.9289e-01, 8.0491e-01, 8.1696e-01,
      8.2904e-01, 8.4114e-01, 8.5327e-01, 8.6542e-01, 8.7759e-01, 8.8978e-01,
      9.0198e-01, 9.1420e-01, 9.2644e-01, 9.3868e-01, 9.5093e-01, 9.6319e-01,
      9.7546e-01, 9.8773e-01, 1.0000e-00
    };

    double Umean[nb_values] =
    {
      0.0000e+00, 4.4231e-02, 1.7699e-01, 3.9816e-01, 7.0750e-01, 1.1046e+00,
      1.5886e+00, 2.1573e+00, 2.8061e+00, 3.5272e+00, 4.3078e+00, 5.1307e+00,
      5.9748e+00, 6.8174e+00, 7.6375e+00, 8.4177e+00, 9.1455e+00, 9.8142e+00,
      1.0421e+01, 1.0968e+01, 1.1458e+01, 1.1895e+01, 1.2287e+01, 1.2636e+01,
      1.2950e+01, 1.3233e+01, 1.3489e+01, 1.3722e+01, 1.3935e+01, 1.4131e+01,
      1.4313e+01, 1.4482e+01, 1.4640e+01, 1.4789e+01, 1.4931e+01, 1.5066e+01,
      1.5196e+01, 1.5321e+01, 1.5442e+01, 1.5560e+01, 1.5674e+01, 1.5786e+01,
      1.5896e+01, 1.6003e+01, 1.6110e+01, 1.6214e+01, 1.6317e+01, 1.6418e+01,
      1.6518e+01, 1.6616e+01, 1.6713e+01, 1.6808e+01, 1.6902e+01, 1.6995e+01,
      1.7087e+01, 1.7179e+01, 1.7269e+01, 1.7358e+01, 1.7446e+01, 1.7533e+01,
      1.7620e+01, 1.7705e+01, 1.7789e+01, 1.7873e+01, 1.7956e+01, 1.8038e+01,
      1.8120e+01, 1.8202e+01, 1.8282e+01, 1.8362e+01, 1.8441e+01, 1.8520e+01,
      1.8598e+01, 1.8676e+01, 1.8754e+01, 1.8831e+01, 1.8907e+01, 1.8982e+01,
      1.9056e+01, 1.9128e+01, 1.9199e+01, 1.9269e+01, 1.9338e+01, 1.9406e+01,
      1.9473e+01, 1.9539e+01, 1.9604e+01, 1.9668e+01, 1.9731e+01, 1.9794e+01,
      1.9855e+01, 1.9916e+01, 1.9976e+01, 2.0035e+01, 2.0093e+01, 2.0149e+01,
      2.0204e+01, 2.0258e+01, 2.0311e+01, 2.0364e+01, 2.0415e+01, 2.0464e+01,
      2.0513e+01, 2.0561e+01, 2.0608e+01, 2.0653e+01, 2.0698e+01, 2.0741e+01,
      2.0784e+01, 2.0826e+01, 2.0866e+01, 2.0905e+01, 2.0943e+01, 2.0979e+01,
      2.1013e+01, 2.1046e+01, 2.1076e+01, 2.1105e+01, 2.1131e+01, 2.1156e+01,
      2.1178e+01, 2.1198e+01, 2.1215e+01, 2.1230e+01, 2.1242e+01, 2.1251e+01,
      2.1258e+01, 2.1262e+01, 2.1263e+01
    };

    if(zz<=1)
    {
      val = zz;
    }
    else
    {
      val = 2-zz;
    }
    for(int i=0 ; i<nb_values-1 ; i++)
    {
      if(val>=z[i] && val<=z[i+1])
      {
        value = val*(Umean[i+1]-Umean[i])/(z[i+1]-z[i])
              + (Umean[i]*z[i+1]-Umean[i+1]*z[i])/(z[i+1]-z[i]);
      }
    }
    return(value);
  }

  void ExactU1(double, double, double, double* values)
  {
    values[0] = 0;
    values[1] = 0;
    values[2] = 0;
    values[3] = 0;
    values[4] = 0;
  }

  void ExactU2(double, double, double, double* values)
  {
    values[0] = 0;
    values[1] = 0;
    values[2] = 0;
    values[3] = 0;
    values[4] = 0;
  }

  void ExactU3(double, double, double, double* values)
  {
    values[0] = 0;
    values[1] = 0;
    values[2] = 0;
    values[3] = 0;
    values[4] = 0;
  }

  void ExactP(double, double, double, double* values)
  {
    values[0] = 0;
    values[1] = 0;
    values[2] = 0;
    values[3] = 0;
    values[4] = 0;
  }

  // ===========================================================================
  // initial solution
  // ===========================================================================
  void InitialU1(double, double, double z, double* values)
  {
    if(RE==180)
    {
      values[0] = DNS_profile_180(z)
                + noise*INTERNAL_BULK_SIMULATION*(2*(double)rand()/RAND_MAX-1);
    }
    else if(RE==395)
    {
      values[0] = DNS_profile_395(z)
                + noise*INTERNAL_BULK_SIMULATION*(2*(double)rand()/RAND_MAX-1);
    }
    else if(RE==590)
    {
      values[0] = DNS_profile_395(z)
                + noise*INTERNAL_BULK_SIMULATION*(2*(double)rand()/RAND_MAX-1);
    }
    else
    {
      ErrThrow("wrong reynolds number ");
    }
    if((std::abs(z-Z_MIN)<tol)||(std::abs(Z_MAX-z)<tol))
    {
      values[0] = 0;
    }
  }

  void InitialU2(double, double, double z, double* values)
  {
    values[0] = noise*INTERNAL_BULK_SIMULATION*(2*(double)rand()/RAND_MAX-1);
    if((std::abs(z-Z_MIN)<tol)||(std::abs(Z_MAX-z)<tol))
    {
      values[0] = 0;
    }
  }

  void InitialU3(double, double, double z, double* values)
  {
    values[0] = noise*INTERNAL_BULK_SIMULATION*(2*(double)rand()/RAND_MAX-1);
    if((std::abs(z-Z_MIN)<tol)||(std::abs(Z_MAX-z)<tol))
    {
      values[0] = 0;
    }
  }

  void InitialP(double, double, double, double *values)
  {
    values[0] = 0;
    TDatabase::ParamDB->INTERNAL_PERIODIC_IDENTITY = 4;
  }

  // ===========================================================================
  // boundary conditions
  // ===========================================================================
  void BoundCondition(int ref, double, double, double, BoundCond& cond)
  {
    if(ref == 1 || ref == 2 || ref == 3)
    {
      cond = PERIODIC;
    }
    else if(ref == 4) // walls
    {
      cond = DIRICHLET;
    }
  }

  void U1BoundValue(int, double, double, double, double& value)
  {
    value = 0;
  }

  void U2BoundValue(int, double, double, double, double& value)
  {
    value = 0;
  }

  void U3BoundValue(int, double, double, double, double& value)
  {
    value = 0;
  }

  // ===========================================================================
  // coefficients for Stokes form: A, B1, B2, f1, f2
  // ===========================================================================
  void LinCoeffs(int n_points, const double*, const double*, const double*,
                 const double*const*, double** coeffs)
  {
    double* coeff;
    double eps = 1./RE;
    double u1 = INTERNAL_BULK_EXPECTED;
    double u2 = TNS_ChannelTau::get_bulkVelocity_sim();
    double dt = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
    //Output::print("bulk velocity: mean ", u1, " sim ", u2);

    for(int i=0 ; i<n_points ; i++)
    {
      coeff = coeffs[i];

      coeff[0] = eps;
      coeff[1] = 1 + (u1 - u2)/dt;
      coeff[2] = 0;
      coeff[3] = 0;
      coeff[4] = 0.;
      coeff[5] = 0.;
      coeff[6] = 0.;
    }
  }

} // namespace channel_tau
#endif
