// ==========================================================================
// instationary problem
// ==========================================================================

//===========================================================================
// example file
// =========================================================================
#define __SIN4__

void ExampleFile()
{
  Output::print<1>("Example: Sin4.h");
}

// exact solution
void Exact(double x, double y, double z, double *values)
{
  double t;

    t = TDatabase::TimeDB->CURRENTTIME;
    values[0] = std::sin(t)*(std::sin(2*Pi*x)*std::sin(2*Pi*y)*std::sin(2*Pi*z)+1);
    values[1] = std::sin(t)*2*Pi*(std::cos(2*Pi*x)*std::sin(2*Pi*y)*std::sin(2*Pi*z));
    values[2] = std::sin(t)*2*Pi*(std::sin(2*Pi*x)*std::cos(2*Pi*y)*std::sin(2*Pi*z));
    values[3] = std::sin(t)*2*Pi*(std::sin(2*Pi*x)*std::sin(2*Pi*y)*std::cos(2*Pi*z));
    values[4] = std::sin(t)*4*Pi*Pi*(-std::sin(2*Pi*x)*std::sin(2*Pi*y)*std::sin(2*Pi*z)
				-std::sin(2*Pi*x)*std::sin(2*Pi*y)*std::sin(2*Pi*z)
				-std::sin(2*Pi*x)*std::sin(2*Pi*y)*std::sin(2*Pi*z));
}

// initial conditon
void InitialCondition(double x, double y, double z, double *values)
{
    values[0] = 0;
}

// kind of boundary condition
void BoundCondition(int, double x, double y, double z, BoundCond &cond)
{
  if (TDatabase::ParamDB->SOLD_PARAMETER_TYPE== FEM_FCT)
    cond = NEUMANN;
  else
    cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(int, double x, double y, double z, double &value)
{
  double t;

  t = TDatabase::TimeDB->CURRENTTIME;
  value = std::sin(t);
}

void BilinearCoeffs(int n_points, double *X, double *Y, double *Z,
        double **parameters, double **coeffs)
{
  double eps = 1.0/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff;                                  // *param;
  double x, y, z, c, a[3], b[3], s[3], h;
  double t = TDatabase::TimeDB->CURRENTTIME;
  
  b[0] = 1;
  b[1] = 1;
  b[2] = 1;  
  c = 0;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    // param = parameters[i];

    x = X[i];
    y = Y[i];
    z = Z[i];

    // diffusion
    coeff[0] = eps;
    // convection in x direction
    coeff[1] = b[0];
    // convection in y direction
    coeff[2] = b[1];
    // convection in z direction
    coeff[3] = b[2];
    // reaction
    coeff[4] = c;
     // rhs
    coeff[5] = std::cos(t)*(std::sin(2*Pi*x)*std::sin(2*Pi*y)*std::sin(2*Pi*z)+1)
	-coeff[0] * (std::sin(t)*4*Pi*Pi*(-std::sin(2*Pi*x)*std::sin(2*Pi*y)*std::sin(2*Pi*z)
			-std::sin(2*Pi*x)*std::sin(2*Pi*y)*std::sin(2*Pi*z)
			-std::sin(2*Pi*x)*std::sin(2*Pi*y)*std::sin(2*Pi*z)))
	+coeff[1] * std::sin(t)*2*Pi*(std::cos(2*Pi*x)*std::sin(2*Pi*y)*std::sin(2*Pi*z))
	+coeff[2] * std::sin(t)*2*Pi*(std::sin(2*Pi*x)*std::cos(2*Pi*y)*std::sin(2*Pi*z))
	+coeff[3] * std::sin(t)*2*Pi*(std::sin(2*Pi*x)*std::sin(2*Pi*y)*std::cos(2*Pi*z))
	+coeff[4] * std::sin(t)*(std::sin(2*Pi*x)*std::sin(2*Pi*y)*std::sin(2*Pi*z)+1);
    // rhs from previous time step
    coeff[6] = 0;
  }
}

/****************************************************************/
//
// for FEM_TVD
//
/****************************************************************/
void CheckWrongNeumannNodes(TCollection *Coll, TFESpace3D *fespace,
int &N_neum_to_diri, int* &neum_to_diri,
			    double* &neum_to_diri_x, double* &neum_to_diri_y, double* &neum_to_diri_z)
{
   const int max_entries = 50000;  
  int i, j, N_, min_val;
  int N_Cells, N_V, diri_counter = 0, found, diri_counter_1 = 0;
  int *dof;
  int boundary_vertices[8], tmp_diri[max_entries];
  double x[8], y[8], z[8], eps = 1e-6, tmp_x[max_entries], tmp_y[max_entries], tmp_z[max_entries];
  TBaseCell *cell;
  TVertex *vertex;
  FE3D CurrentElement;
  
  // number of mesh cells
  N_Cells = Coll->GetN_Cells();

  diri_counter = 0;
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_V = cell->GetN_Vertices();
    found = 0;
    for (j=0;j<N_V;j++)
    {
      // read coordinates of the mesh cell
      boundary_vertices[j] = 0;
      vertex = cell->GetVertex(j);
      vertex->GetCoords(x[j], y[j], z[j]);
      // vertex on the upper lid
      if ((std::abs(x[j])<eps)||(std::abs(y[j])<eps)||(std::abs(x[j]-1)<eps)||(std::abs(y[j]-1)<eps)
	  ||(std::abs(z[j])<eps)||(std::abs(z[j]-1)<eps))
      {
	  // Dirichlet boundary
	  boundary_vertices[j] = 1;
	  found++;
      }
    }
    
    // no cell with face with vertex on the boundary
    if (found<3) 
	continue;
    // finite element on the mesh cell
    CurrentElement = fespace->get_fe_type(i);
    // number of basis functions (= number of d.o.f.)
    N_ = fespace.get_n_local_dof(i);
    // the array which gives the mapping of the local to the global d.o.f.
    dof = fespace->GetGlobalDOF(i);
    switch(CurrentElement)
    {
	// P_1, Q_1
	case C_P1_3D_T_A:
	case C_Q1_3D_H_A:
	case C_Q1_3D_H_M:
	    for (j=0;j<N_V;j++)
	    {
		// vertex on the boundary
		if (boundary_vertices[j])
		{
		    if (CurrentElement==C_P1_3D_T_A)
			tmp_diri[diri_counter] = dof[j];
		    else
		    {
			if ((j==0)||(j==1)||(j==4)||(j==5))
			{
			    tmp_diri[diri_counter] = dof[j];
			}
			else
			{
			    if (j==2)
				tmp_diri[diri_counter] = dof[3];
			    if (j==3)
				tmp_diri[diri_counter] = dof[2];
			    if (j==6)
				tmp_diri[diri_counter] = dof[7];
			    if (j==7)
				tmp_diri[diri_counter] = dof[6];
			}
		    }
		    if (diri_counter > max_entries)
		    {
			ErrThrow("tmp_diri too short !!!");
		    }
		    tmp_x[diri_counter] = x[j];
		    tmp_y[diri_counter] = y[j];
		    tmp_z[diri_counter] = z[j];
		    diri_counter++;
        Output::print(j, " ", tmp_x[diri_counter-1], " ", x[j]);
		}
	    }
	    break;
	default:
	    ErrThrow("CheckNeumannNodesForVelocity not implemented for element ",
               CurrentElement);
	    break;
    }	    
  }
    
  // condense
  for (i=0;i<diri_counter;i++)
  {
      if (tmp_diri[i] == -1)
	  continue;
      diri_counter_1++;
      for (j=i+1;j<diri_counter;j++)
      {
	  if (tmp_diri[i] == tmp_diri[j])
	  {
	      tmp_diri[j] = -1;
	  }
      }
  }
  
  Output::print("CheckNeumannNodesForVelocity: N_neum_to_diri ", diri_counter_1);
  N_neum_to_diri = diri_counter_1;
  // allocate array for the indices
  neum_to_diri = new int[diri_counter_1];
  // allocate array for the corresponding boundary coordinates
  neum_to_diri_x = new double[diri_counter_1];
  neum_to_diri_y = new double[diri_counter_1];
  neum_to_diri_z = new double[diri_counter_1];
  // fill array and sort
  for (i=0;i<diri_counter_1;i++)
  {
      min_val = tmp_diri[0];
      found = 0;
      for (j=1;j<diri_counter;j++)
      {
	  if ((tmp_diri[j]>0) && ((tmp_diri[j] < min_val) || 
				  (min_val == -1)))
	  {
	      min_val =  tmp_diri[j];
	      found = j;
	  }
      }
      neum_to_diri[i] = tmp_diri[found];
      neum_to_diri_x[i] = tmp_x[found];
      neum_to_diri_y[i] = tmp_y[found];
      neum_to_diri_z[i] = tmp_z[found];
      tmp_diri[found] = -1;
  }
  
  for (i=0;i<diri_counter_1;i++)
  {
    Output::print(i, " ", neum_to_diri[i], " ", neum_to_diri_x[i], " ",
                  neum_to_diri_y[i], " ", neum_to_diri_z[i]);
  }
}

  
