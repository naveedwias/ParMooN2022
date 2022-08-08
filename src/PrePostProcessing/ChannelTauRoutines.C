#include <ChannelTauRoutines.h>

#include <BaseCell.h>
#include <Collection.h>
#include <Database.h>
#include <MooNMD_Io.h>
#ifdef _MPI
#include <ParFECommunicator3D.h>
#endif
#include <PeriodicJoint.h>

#include <algorithm>    // std::stable_sort
#include <functional>
#include <numeric>


// bulkVelocity_expected and bulkVelocity_sim are initialized in ChannelTau.h
double TNS_ChannelTau::bulkVelocity_expected = 0.;
double TNS_ChannelTau::bulkVelocity_sim = 0.;


/* ************************************************************************** */
TNS_ChannelTau::TNS_ChannelTau(const TDomain&           domain,
                               const ParameterDatabase& param_db)
  : TimeNavierStokes<3>(domain, param_db)
{
  reynolds_number = db["reynolds_number"];

  check_and_set_parameters();

  GetCoordinatesOfDof();

  set_up_memory();

  // compute the summation of layers dof for averaging
  count_dofs_per_layer(sum_layer_dofs);
}


/* ************************************************************************** */
void TNS_ChannelTau::check_and_set_parameters()
{
  int rank = 0;
#ifdef _MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  if(TDatabase::ParamDB->CELL_MEASURE==0)
  {
    TDatabase::ParamDB->CELL_MEASURE = 2;
    if(rank==0)
    {
      Output::print("CELL_MEASURE set to: ",
                    TDatabase::ParamDB->CELL_MEASURE);
    }
  }

  TDatabase::ParamDB->INTERNAL_QUAD_RULE = 2;
  if(rank==0)
  {
    Output::print("INTERNAL_QUAD_RULE set to: ",
                  TDatabase::ParamDB->INTERNAL_QUAD_RULE);
  }
}


/* ************************************************************************** */
void TNS_ChannelTau::GetCoordinatesOfDof()
{
#ifdef _MPI
  int rank;
  int size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  std::shared_ptr<const TFESpace3D> feSpace = this->get_velocity_space();

  size_t nDofs = feSpace->get_n_dof();
  xDofs.resize(nDofs); yDofs.resize(nDofs); zDofs.resize(nDofs);
  for(size_t i=0 ; i<nDofs ; i++)
  {
    feSpace->GetDOFPosition(i, xDofs[i], yDofs[i], zDofs[i]);
  }

//   // some checks on correct elements
//   size_t nCells = feSpace->GetCollection()->GetN_Cells();
//   for(size_t i=0 ; i<nCells ; i++)
//   {
//     FE_type cE = feSpace->get_fe_type(i);
// 
//     if((cE != C_Q2_3D_H_A && cE != C_Q2_3D_H_M)
//     && (cE != C_Q1_3D_H_A && cE != C_Q1_3D_H_M))
//     {
//       ErrThrow("coordinates of dofs are not tested for the finite element ",
//                cE, " yet");
//     }
//   }

  // get the different Z-positions of the layers
  // (with periodic boundaries, some dofs are identified)
  for(size_t i=0 ; i<nDofs ; i++)
  {
#ifdef _MPI
    const int* masters = this->get_velocity_space()
                             ->get_communicator().GetMaster();
    if(masters[i] != rank)
    {
      continue;
    }
#endif
    if(!is_containing(zLayers, zDofs[i]))
    {
      zLayers.push_back(zDofs[i]);
      nZLayers++;
    }
  }
#ifdef _MPI
  int nZLayers_glob = 0;
  std::vector<int> receivebuf(size);
  // gather the numbers nzlayers of layers to root
  MPI_Gather(&nZLayers, 1, MPI_INT,         // send
             receivebuf.data(), 1, MPI_INT, // receive
             0, MPI_COMM_WORLD);            // control

  // count the global number of zlayers on root
  std::vector<int> displs(size);
  if(rank==0)
  {
    for(int i=0 ; i<size ; i++)
    {
      displs[i] = nZLayers_glob;
      nZLayers_glob += receivebuf[i];
    }
  }

  double* sendzlayer = zLayers.data();
  std::vector<double> recievezlayer(nZLayers_glob);
  // gather the zlayers from all processes to root
  MPI_Gatherv(sendzlayer, nZLayers, MPI_DOUBLE, // send
              recievezlayer.data(), receivebuf.data(), displs.data(), // receive
              MPI_DOUBLE, 0, MPI_COMM_WORLD);  // control

  if(rank==0)
  {
    std::sort(recievezlayer.begin(), recievezlayer.end());
    auto new_end = unique_double(recievezlayer.begin(), recievezlayer.end());
    recievezlayer.erase(new_end, recievezlayer.end());
    nZLayers_glob = recievezlayer.size();
    zLayers.resize(nZLayers_glob);
    std::copy(recievezlayer.begin(), recievezlayer.end(), zLayers.begin());
  }

  // give the information back to all processes
  MPI_Bcast(&nZLayers_glob, 1, MPI_INT, 0, MPI_COMM_WORLD);
  nZLayers = nZLayers_glob;
  if(rank != 0)
  {
    zLayers.resize(nZLayers);
  }
  MPI_Bcast(&zLayers[0], nZLayers, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#else
  // sort the zlayers
  std::sort(zLayers.begin(), zLayers.end());
#endif
}


/* ************************************************************************** */
void TNS_ChannelTau::set_up_memory()
{
  std::vector<double> temp(nZLayers, 0.);
  for(int i=0 ; i<3 ; i++)
  {
    MeanVelocity.push_back(temp);
  }

  for(int i=0 ; i<6 ; i++)
  {
    MeanReynoldsStress.push_back(temp);
  }

  DerivMeanVelocity.resize(nZLayers, 0.);

  sum_layer_dofs.resize(nZLayers, 0);
}


/* ************************************************************************** */
void TNS_ChannelTau::count_dofs_per_layer(std::vector<int>& summ)
{
  for(size_t i=0 ; i<zDofs.size() ; i++)
  {
#ifdef _MPI
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::shared_ptr<const TFESpace3D> space = this->get_velocity_space();
    if(space->get_communicator().GetMaster()[i] != rank)
    {
      continue;
    }
#endif
    for(size_t j=0 ; j<nZLayers ; j++)
    {
      if(std::abs(zDofs[i]-zLayers[j]) < tolerance)
      {
        summ[j]++;
        break;
      }
    }
  }
#ifdef _MPI
  int n_elems = nZLayers;
  MPI_Allreduce(MPI_IN_PLACE, &summ[0], n_elems,
                MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif
}


/* ************************************************************************** */
void TNS_ChannelTau::output()
{
  this->TimeNavierStokes<3>::output();

  computeMeanValues();
}


/* ************************************************************************** */
void TNS_ChannelTau::computeMeanValues()
{
  // spatial mean
  // compute the summation of all values per layer...
  std::vector<double> temp(nZLayers, 0);
  std::deque<std::vector<double>> spatialVelocityMean;
  for(size_t i=0 ; i<3 ; i++)
  {
    spatialVelocityMean.push_back(temp);
  }
  std::deque<std::vector<double>> spatialUiUjMean;
  for(size_t i=0 ; i<6 ; i++)
  {
    spatialUiUjMean.push_back(temp);
  }
  compute_sum_u_uiuj_on_layers(spatialVelocityMean, spatialUiUjMean);

  // ... and then the average
  for(auto& v : spatialVelocityMean)
  {
    std::transform(v.begin(), v.begin()+nZLayers, sum_layer_dofs.begin(),
                   v.begin(), std::divides<double>());
  }
  for(auto& v : spatialUiUjMean)
  {
    std::transform(v.begin(), v.begin()+nZLayers, sum_layer_dofs.begin(),
                   v.begin(), std::divides<double>());
  }

  // compute the bulk velocity (streamwise velocity)
  compute_bulkVelocity(spatialVelocityMean[0]);

  // temporal mean
  // define effective starting time for averaging and initialize
  if(t0_avg<0.) // time averaging not started yet
  {
    MeanVelocity[0] = std::move(spatialVelocityMean[0]);
    MeanVelocity[1] = std::move(spatialVelocityMean[1]);
    MeanVelocity[2] = std::move(spatialVelocityMean[2]);

    MeanReynoldsStress[0] = std::move(spatialUiUjMean[0]);
    MeanReynoldsStress[1] = std::move(spatialUiUjMean[1]);
    MeanReynoldsStress[2] = std::move(spatialUiUjMean[2]);
    MeanReynoldsStress[3] = std::move(spatialUiUjMean[3]);
    MeanReynoldsStress[4] = std::move(spatialUiUjMean[4]);
    MeanReynoldsStress[5] = std::move(spatialUiUjMean[5]);

    double t = this->get_time_stepping_scheme().current_time_;
    double start_t_avg = this->get_db()["start_time_averaging_at"];
    if(t>=start_t_avg) // start the time averaging
    {
      t0_avg = t;
    }
  }
  else // time averaging is already started
  {
    for(size_t i=0 ; i<nZLayers ; i++)
    {
      temporalMean(spatialVelocityMean[0], MeanVelocity[0]);
      temporalMean(spatialVelocityMean[1], MeanVelocity[1]);
      temporalMean(spatialVelocityMean[2], MeanVelocity[2]);

      temporalMean(spatialUiUjMean[0], MeanReynoldsStress[0]);
      temporalMean(spatialUiUjMean[1], MeanReynoldsStress[1]);
      temporalMean(spatialUiUjMean[2], MeanReynoldsStress[2]);
      temporalMean(spatialUiUjMean[3], MeanReynoldsStress[3]);
      temporalMean(spatialUiUjMean[4], MeanReynoldsStress[4]);
      temporalMean(spatialUiUjMean[5], MeanReynoldsStress[5]);
    }
  }

  /// compute mean velocity z-derivative
  compute_derivative(MeanVelocity[0], DerivMeanVelocity);

  print_quantity_of_interest();

}


/* ************************************************************************** */
void TNS_ChannelTau::compute_sum_u_uiuj_on_layers(
                                  std::deque<std::vector<double>>& velocity_sum,
                                  std::deque<std::vector<double>>&  uiuj_mean)
{
  std::shared_ptr<const TFESpace3D> space = this->get_velocity_space();
  const TFEVectFunct3D& U = this->get_velocity();
  size_t nuDofs = U.GetComponent(0)->GetLength();
  // component of solution vector
  const std::vector<double> u1(U.GetComponent(0)->GetValues(),
                               U.GetComponent(0)->GetValues()+nuDofs);
  const std::vector<double> u2(U.GetComponent(1)->GetValues(),
                               U.GetComponent(1)->GetValues()+nuDofs);
  const std::vector<double> u3(U.GetComponent(2)->GetValues(),
                               U.GetComponent(2)->GetValues()+nuDofs);
  std::deque<std::vector<double>> u;
  u.push_back(u1);  u.push_back(u2);  u.push_back(u3);

  for(size_t i=0 ; i<nuDofs ; i++)
  {
#ifdef _MPI
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const int* masters = this->get_velocity_space()
                             ->get_communicator().GetMaster();
    if(masters[i] != rank)
    {
      continue;
    }
#endif
    for(size_t j=0 ; j<nZLayers ; j++)
    {
      if(std::abs(zDofs[i]-zLayers[j]) < tolerance)
      {
        velocity_sum[0][j] = velocity_sum[0][j] + u[0][i];
        velocity_sum[1][j] = velocity_sum[1][j] + u[1][i];
        velocity_sum[2][j] = velocity_sum[2][j] + u[2][i];

        uiuj_mean[0][j] = uiuj_mean[0][j] + u[0][i]*u[0][i];
        uiuj_mean[1][j] = uiuj_mean[1][j] + u[1][i]*u[1][i];
        uiuj_mean[2][j] = uiuj_mean[2][j] + u[2][i]*u[2][i];
        uiuj_mean[3][j] = uiuj_mean[3][j] + u[0][i]*u[1][i];
        uiuj_mean[4][j] = uiuj_mean[4][j] + u[0][i]*u[2][i];
        uiuj_mean[5][j] = uiuj_mean[5][j] + u[1][i]*u[2][i];
        break;
      }
    }
  }
  // communicate the meanvelocity and sum up
#ifdef _MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int n_elems = nZLayers;
  for(size_t comp=0 ; comp<velocity_sum.size() ; comp++)
  {
    MPI_Allreduce(MPI_IN_PLACE, &velocity_sum[comp][0],
                  n_elems, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }
  for(size_t comp=0 ; comp<uiuj_mean.size() ; comp++)
  {
    MPI_Allreduce(MPI_IN_PLACE, &uiuj_mean[comp][0],
                  n_elems, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }
#endif
}


/* ************************************************************************** */
void TNS_ChannelTau::compute_bulkVelocity(
                                        const std::vector<double>& spatialVMean)
{
  bulkVelocity_sim = 0;
  for(size_t i=0 ; i<nZLayers-1 ; i++)
  {
    bulkVelocity_sim += 0.5*(spatialVMean[i+1]+spatialVMean[i])
                           *(zLayers[i+1]-zLayers[i]);
  }
  bulkVelocity_sim /= (zLayers[nZLayers-1] - zLayers[0]);
}


/* ************************************************************************** */
void TNS_ChannelTau::temporalMean(std::vector<double> spatial_mean,
                                  std::vector<double>& temporal_mean)
{
  double step_length = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double currentTime = TDatabase::TimeDB->CURRENTTIME;
  double factor;

  if(currentTime - t0_avg <= step_length + tolerance)
  {
    factor = 0.5; // first averaging step
  }
  else
  {
    factor = step_length / (currentTime - t0_avg);
  }
  for(size_t i=0 ; i<spatial_mean.size() ; i++)
  {
    temporal_mean[i] = temporal_mean[i]
                        + factor*(spatial_mean[i] - temporal_mean[i]);
  }
}


/* ************************************************************************** */
void TNS_ChannelTau::compute_derivative(const std::vector<double>& values,
                                        std::vector<double>&       derivative)
{
  derivative[0] = values[1] / zLayers[1];
  derivative[nZLayers-1] = (values[nZLayers-2] - values[nZLayers-1])
                            / (zLayers[nZLayers-2] - zLayers[nZLayers-1]);

  for(size_t i=1 ; i<nZLayers-1 ; i++)
  {
    double D = (zLayers[i]-zLayers[i-1])
             * (zLayers[i+1]-zLayers[i-1])
             * (zLayers[i+1]-zLayers[i]);

    double b = ( (values[i]-values[i-1])
                 * (zLayers[i+1]*zLayers[i+1] - zLayers[i-1]*zLayers[i-1])
               - (values[i+1]-values[i-1])
                 * (zLayers[i]*zLayers[i] - zLayers[i-1]*zLayers[i-1]) )
             / D;

    double a = ( (zLayers[i]-zLayers[i-1]) * (values[i+1]-values[i-1])
               - (zLayers[i+1]-zLayers[i-1]) * (values[i]-values[i-1]) )
             / D;

    derivative[i] = 2*a*zLayers[i] + b;
  }
}


/* ************************************************************************** */
void TNS_ChannelTau::print_quantity_of_interest()
{
  int rank=0;
#ifdef _MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  if(rank==0)
  {
    double t = TDatabase::TimeDB->CURRENTTIME;
    Output::print("at t: ", t, "\t bulk velocity: ", setw(8), bulkVelocity_sim,
                  " tau: ", (bulkVelocity_expected-bulkVelocity_sim)+1);

    if(t0_avg < 0)
    {
      return;
    }
    else
    {
      // friction velocity u_tau and derivative of mean velocity dmu(streamwise)
      double u_tau = get_FrictionVelocity();
      for(size_t i=0 ; i<nZLayers ; i++)
      {
        Output::print("at t: ", std::scientific, t,
                      " z: ", setw(8), zLayers[i],
                      " z+: ", setw(8), reynolds_number
                                       *(1-std::abs(1-zLayers[i])),
                      "\t u_m: ", setw(8), MeanVelocity[0][i],
                      " du_m: ", setw(8), DerivMeanVelocity[i],
                      " v_m: ", setw(8), MeanVelocity[2][i],
                      " w_m: ", setw(8), MeanVelocity[1][i]);
      }
      // root mean square velocities rms
      double rmsu, rmsv, rmsw;
      double R12, R13, R23;
      double R12_abs=0., R13_abs=0., R23_abs=0.;
      std::deque<std::vector<double>> rms = get_rms();
      for(size_t i=0 ; i<nZLayers ; i++)
      {
        // formula 11 Volker's paper
        rmsu = std::sqrt(std::abs(rms[0][i]));
        rmsv = std::sqrt(std::abs(rms[1][i]));
        rmsw = std::sqrt(std::abs(rms[2][i]));
        // formula 10 Volker's paper
        R12 = MeanReynoldsStress[3][i]
            - MeanVelocity[0][i] * MeanVelocity[1][i];
        R12_abs += std::abs(R12);
        R13 = MeanReynoldsStress[4][i]
            - MeanVelocity[0][i] * MeanVelocity[2][i];
        R13_abs += std::abs(R13);
        R23 = MeanReynoldsStress[5][i]
            - MeanVelocity[1][i] * MeanVelocity[2][i];
        R23_abs += std::abs(R23);

        Output::print("at t: ", std::scientific, t,
                      " z: ", setw(8), zLayers[i],
                      " z+: ", setw(8), reynolds_number
                                       *(1-std::abs(1-zLayers[i])),
                      "\t rms_u: ", setw(8), rmsu,
                      " rms_v: ", setw(8), rmsw,
                      " rms_w: ", setw(8), rmsv,
                      " R_uv: ", setw(8), R13,
                      " R_uw: ", setw(8), R12,
                      " R_vw: ", setw(8), R23);
      }
      Output::print("u_tau ", u_tau, " zero statistics ", setw(8),
                    R12_abs/nZLayers, " ", setw(8), R13_abs/nZLayers, " ",
                    setw(8), R23_abs/nZLayers);
    }
  }
}


/* ************************************************************************** */
double TNS_ChannelTau::get_FrictionVelocity()
{
  double u_tau = 1./reynolds_number
               * 0.5*(DerivMeanVelocity[0]-DerivMeanVelocity[nZLayers-1]);
  return u_tau;
}


/* ************************************************************************** */
std::deque<std::vector< double>> TNS_ChannelTau::get_rms()
{
  std::deque<std::vector<double>> rms(3);
  for(size_t i=0 ; i<nZLayers ; i++)
  {
    double rms_u = 2./3.*( MeanReynoldsStress[0][i]
                         - std::pow(MeanVelocity[0][i],2))
                 - 1./3.*( MeanReynoldsStress[1][i]
                         - std::pow(MeanVelocity[1][i],2))
                 - 1./3.*( MeanReynoldsStress[2][i]
                         - std::pow(MeanVelocity[2][i],2));
    rms[0].push_back(rms_u);

    double rms_v = 2./3.*( MeanReynoldsStress[1][i]
                         - std::pow(MeanVelocity[1][i],2))
                 - 1./3.*( MeanReynoldsStress[0][i]
                         - std::pow(MeanVelocity[0][i],2))
                 - 1./3.*( MeanReynoldsStress[2][i]
                         - std::pow(MeanVelocity[2][i],2));
    rms[1].push_back(rms_v);

    double rms_w = 2./3.*( MeanReynoldsStress[2][i]
                         - std::pow(MeanVelocity[2][i],2))
                 - 1./3.*( MeanReynoldsStress[0][i]
                         - std::pow(MeanVelocity[0][i],2))
                 - 1./3.*( MeanReynoldsStress[1][i]
                         - std::pow(MeanVelocity[1][i],2));
    rms[2].push_back(rms_w);
  }
  return rms;
}


/* ************************************************************************** */
bool is_containing(const std::vector<double>& V,
                   const double x,
                   double tolerance)
{
  auto it = std::find_if(V.begin(), V.end(), [x, tolerance](double V_i)
            { return std::abs(x - V_i) <= tolerance; });
  if(it!=V.end())
  {
    return true;
  }
  else
  {
    return false;
  }
}


/* ************************************************************************** */
std::vector<double>::iterator unique_double(std::vector<double>::iterator first,
                                            std::vector<double>::iterator last,
                                            double tolerance)
{
  if(first==last)
  {
    return last;
  }

  std::vector<double>::iterator new_end = first;
  while(++first!=last)
  {
    if((std::abs(*new_end-*first) > tolerance) && ++new_end != first)
    {
      *new_end = std::move(*first);
    }
  }
  return ++new_end;
}
