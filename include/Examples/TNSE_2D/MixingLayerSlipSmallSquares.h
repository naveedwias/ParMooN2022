#ifndef MIXINGLAYERSLIPSMALLSQUARES_H
#define MIXINGLAYERSLIPSMALLSQUARES_H

double DIMENSIONLESS_VISCOSITY;

#define U_INFTY 1

void ExampleFile()
{
  Output::print("Example: MixingLayerSlipSmallSquares.h");
  TDatabase::ParamDB->INTERNAL_QUAD_RULE = 97;
}

// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double x, double y, double *values)
{
  double rho_0, z1, z2, cn;

  rho_0 = 1./28.;
  z1 = (4*y-2)/rho_0;
  z2 = (y-0.5)/rho_0;
  cn = 1.e-3;
  values[0] = U_INFTY * (std::exp(z1)-1)/(std::exp(z1)+1);
  
  values[0] -= cn * U_INFTY * std::exp(-z2*z2)*2*z2/rho_0*(std::cos(8.*M_PI*x) + std::cos(20.*M_PI*x));
}

void InitialU2(double x, double y, double *values)
{
  double rho_0, z2, cn;

  rho_0 = 1./28.;
  z2 = (y-0.5)/rho_0;
  cn = 1.e-3;
  values[0] = cn*U_INFTY * std::exp(-z2*z2)*(std::sin(8.*M_PI*x)*8.*M_PI + std::sin(20.*M_PI*x) * 20.*M_PI);
}

void InitialP(double, double, double *values)
{
  values[0] = 0;
}


// ========================================================================
// exact solution
// ========================================================================
auto& ExactU1 = unknown_solution_2d;
auto& ExactU2 = unknown_solution_2d;
auto& ExactP = unknown_solution_2d;

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int, double, BoundCond &cond)
{
  cond = SLIP_FRICTION_PENETRATION_RESISTANCE;
  TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION = 1;
}

void U1BoundValue(int BdComp, double, double &value)
{
  if(BdComp>3)
  {
    ErrThrow( "ERROR in file " , __FILE__ , ", line: ",  __LINE__ ,
              ": wrong boundary part number: " , BdComp);
  }
  value =0.;
}

void U2BoundValue(int, double, double &value)
{
  value = 0;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, const double *, const double *,
               const double *const*, double **coeffs)
{
  for(int i=0;i<n_points;i++)
  {
    coeffs[i][0] = DIMENSIONLESS_VISCOSITY;
    coeffs[i][1] = 0;
    coeffs[i][2] = 0;
    coeffs[i][3] = 0;
  }
}

void ComputeVorticiyThickness(const TFEFunction2D *Vorticity, double &thickness)
{
  int i,j,k,l,found,max_lines,N_Cells,N_Edges,N_Vort;
  double x0,x1,x,y0,y1,val[3],max;
  double *global_y_coord, *aver_vort;

  auto vorticity_space=Vorticity->GetFESpace2D();
  N_Vort = vorticity_space->get_n_dof();
  global_y_coord = new double[N_Vort];
  aver_vort =  new double[N_Vort];

  // get pointer to set of mesh cells which define the fe space
  auto Coll = vorticity_space->GetCollection();
  // get number of mesh cells
  N_Cells = Coll->GetN_Cells();

  max_lines = 0;
  // loop over all mesh cells
  for(i=0;i<N_Cells;i++)
  {
    // get current mesh cell
    auto cell = Coll->GetCell(i);
    // get number of edges
    N_Edges = cell->GetN_Edges();
    // for all vertices
    for (j=0;j<N_Edges;j++)
    {
      // compute coordinate of vertex
      cell->GetVertex(j)->GetCoords(x0, y0);
      // no vertices on lower and upper boundary
      if ((std::abs(y0-1)< 1e-6)||(std::abs(y0)< 1e-6))
        continue;
      // for all edges which has a larger number
      for (k=j+1; k < N_Edges;k++)
      {
        // compute coordinate of vertex
        cell->GetVertex(k)->GetCoords(x1, y1);
        // vertices do not posses the same y coordinate, continue
        if (std::abs(y0-y1)> 1e-6)
          continue;
        // compute integral by edge midpoint rule
        x = (x0+x1)/2;
        // compute value of fe function in midpoint of edge
        Vorticity->FindGradientLocal(cell,i,x,y0,val);
        //Output::print(x, " ", y0, " ", val[0]);
        // add value to entry in array
        found = 0;
        for (l=0;l<max_lines;l++)
        {
          // not the correct entry
          if (std::abs(y0-global_y_coord[l])> 1e-6)
            continue;
          // correct entry found
          found = 1;
          aver_vort[l] += std::abs(x0-x1)*val[0];
        }
        // new entry has to be created
        if (!found)
        {
          global_y_coord[max_lines] = y0;
          aver_vort[max_lines] = std::abs(x0-x1)*val[0];
          max_lines++;
          if (max_lines > N_Vort)
          {
            ErrThrow("Increase N_Vort in ComputeVorticiyThickness !!!");
          }
        }
      }
    }
  }

  // compute averages and vorticity thickness
  max = -1;
  for (i=0;i<max_lines;i++)
  {
    // each integral is computed twice
    aver_vort[i]/=2.0;
    // Output::print(i, " ", aver_vort[i], " ", global_y_coord[i]);
    if (std::abs(aver_vort[i]) > max)
      max = std::abs(aver_vort[i]);
  }

  thickness =  2.*U_INFTY/max; // (2*U_INFTY)/(max/2) (divide max by length of
                               // the integration domain)
  delete global_y_coord;
  delete aver_vort;
  return;
}

void EvaluateSolution(const TimeNavierStokes<2> &tnse2d, double &zero_vort)
{
  const TFEFunction2D& vorticity(tnse2d.get_vorticity_funct());
  double thickness;
  ComputeVorticiyThickness(&vorticity, thickness);

  double t=TDatabase::TimeDB->CURRENTTIME;
  const TFESpace2D *sp = vorticity.GetFESpace2D().get();
  MultiIndex2D allderiv[3]= { MultiIndex2D::D00, MultiIndex2D::D10,
                              MultiIndex2D::D01 };
  TAuxParam2D aux;
  double locerr[8];
  vorticity.GetErrors(unknown_solution_2d, 3, allderiv, 2, L2H1Errors, nullptr, &aux, 1, 
                              &sp, locerr);
  
  Output::print(t, " " , "enstrophy ", (locerr[0]*locerr[0])/2, " ", (locerr[1]*locerr[1])/2);
  
  if(zero_vort < 0)
    zero_vort = thickness;
  Output::print( t, " ", "vorticity thickness: ", thickness, " ", thickness/zero_vort);
  
}
#endif // MIXINGLAYERSLIPSMALLSQUARES_H
