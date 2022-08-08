// ======================================================================
// Example from P.W. Hemker
// ======================================================================
#define __HEMKER1996__

double DIFFUSION;
// solution is 1 (false) or piecewise linear (true) at the cylinder,
bool modified_Hemker = false;

void ExampleFile()
{
  Output::print("Example: Hemker1996.h with diffusion ", DIFFUSION, 
          modified_Hemker ? " and modified boundary data at the cylinder": "");
}

// exact solution
auto& Exact = unknown_solution_2d;

void BoundCondition(int BdComp, double, BoundCond &cond)
{
    switch(BdComp)
    {
	case 0:
	case 1:
	case 2:
	    cond = NEUMANN;
	    break;
	default:
	    cond = DIRICHLET;
    }
}


// value of boundary condition
void BoundValue(int BdComp, double param, double &value)
{
  switch(BdComp)
  {
    case 4:
      value = 1;
      if(modified_Hemker)
      {
        const double p2 = 2 * M_PI;
        if(param >= 0 && param < 0.25) // top left quarter of the circle
          value = std::sin(p2*param); // = y
        else if(param >= 0.25 && param < 0.5) // top right quarter
          value = 1.;
        else if(param >= 0.5 && param < 0.75) // lower right quarter
          value = 1 + std::sin(p2*param); // = 1+y
        else // lower left quarter
          value = 0.;
      }
      break;
    default:
      value = 0;
  }
}

// initial conditon
void InitialCondition(double, double, double *values)
{
  values[0] = 0;
}

void BoundConditionAdjoint(int BdComp, double, BoundCond &cond)
{
    switch(BdComp)
    {
	case 0:
	case 2:
	case 3:
	    cond = NEUMANN;
	    break;
	default:
	    cond = DIRICHLET;
    }
}


// value of boundary condition
void BoundValueAdjoint(int, double, double &value)
{
    value = 0;
}

void BilinearCoeffs(int n_points, const double *, const double *,
                    const double *const*, double **coeffs)
{
//   double eps=1/TDatabase::ParamDB->RE_NR;
  double eps=DIFFUSION;
  double angle = 0, v1, v2;
  int i;
  double *coeff;

  v1 = std::cos(angle);
  v2 = std::sin(angle);

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;
    coeff[1] = v1;
    coeff[2] = v2;
    coeff[3] = 0;

    coeff[4] = 0;
  }
}
 
/** compute curve of the outflow boundary */
void ComputeOutflowBoundary(int , TFEFunction2D *ufct)
{
  double h, x=4,values[3],y;
  int i, bound_points = 401;
  h = 6.0/(bound_points-1);
  for (i=0;i<bound_points; i++)
  {
      y = -3+i*h;
      ufct->FindGradient(x,y,values);
      Output::print("cutline ", x, " ", y, " ", values[0]);
  }
}

// computation of some global errors, only for P1 or Q1 !!!
//
// values[0] : absolute value of largest negative undershoot in 
//             a circle around the cylinder
// values[1] : difference of largest positive value in a circle
//             around the cylinder and 1
// values[2] : absolute value of largest negative undershoot 
//             for x > 2
// values[3] : difference of largest positive value for x>2 and 1
// 
void ComputeLocalExtrema(TFEFunction2D *ufct, double *values)
{
  int N_BaseFunct;
  double xi, eta;
  double *uorig, *uxorig, *uyorig, *uref, *uxiref, *uetaref, u;
  double *Values, x, y, val;
  int N_Cells, N_Edges;
  int i, j, k;
  double extr[4];

  extr[0] = -1;
  extr[1] = -1;
  extr[2] = -1;
  extr[3] = 0;

  auto FESpace2D = ufct->GetFESpace2D();
  Values = ufct->GetValues();  

  auto Coll = FESpace2D->GetCollection();
  N_Cells = Coll->GetN_Cells();

  // loop over all edges
  for(i=0;i<N_Cells;i++)
  {
    auto cell = Coll->GetCell(i);
    N_Edges=cell->GetN_Edges();
    //double diam = cell->GetDiameter();
    
    auto FE_Obj = FESpace2D->get_fe(i);
    auto RefTrans = FE_Obj.GetRefTransID();

    // get base function object
    auto bf = FE_Obj.GetBaseFunct();
    N_BaseFunct = bf->GetDimension();
    
    uorig = new double[N_BaseFunct];
    uxorig = new double[N_BaseFunct];
    uyorig = new double[N_BaseFunct];
    
    uref = new double[N_BaseFunct];
    uxiref = new double[N_BaseFunct];
    uetaref = new double[N_BaseFunct];
    
    // set cell for reference transformation
    FEDatabase::SetCellForRefTrans(cell, RefTrans);
    for (j=0;j<N_Edges;j++)
    {
      // compute coordinates
      x = cell->GetVertex(j)->GetX();
      y = cell->GetVertex(j)->GetY();
      if (x<-1.5)
	  continue;
      // find local coordinates of the given point
      //cout << " x: " << x << endl;
      //cout << " y: " << y << endl;
      FEDatabase::GetRefFromOrig(RefTrans, x, y, xi, eta);
      //cout << " xi: " << xi << endl;
      //cout << "eta: " << eta << endl;

      bf->GetDerivatives(MultiIndex2D::D00, xi, eta, uref);
      bf->GetDerivatives(MultiIndex2D::D10, xi, eta, uxiref);
      bf->GetDerivatives(MultiIndex2D::D01, xi, eta, uetaref);
      
      FEDatabase::GetOrigValues(RefTrans, xi, eta, bf, Coll, (TGridCell *)cell,
                                uref, uxiref, uetaref, uorig, uxorig, uyorig);
      // compute value of fe function at (x,y)
      u = 0;
      auto Numbers = FESpace2D->GetGlobalDOF(i);
      for(k=0;k<N_BaseFunct;k++)
      {
        val = Values[Numbers[k]];
        u += uorig[k]*val;
      }
      //Output::print(x, " ", y, " ", u);
      // strip with the circle
      if ((x>-1.5)&&(x<1.5))
      {
	  if ((u<=0)&&(std::abs(u)>extr[0]))
	      extr[0] = std::abs(u);
	  if ((u>=1)&&(std::abs(u-1)>extr[1]))
	      extr[1] = std::abs(u-1);
      }
      if (x>2)
      {
	  if ((u<=0)&&(std::abs(u)>extr[2]))
	      extr[2] = std::abs(u);
	  if ((u>=1)&&(std::abs(u-1)>extr[3]))
	      extr[3] = std::abs(u-1);
      }
    } // endfor (j) N_Edges
  
    delete uorig;
    delete uxorig;
    delete uyorig;
    delete uref;
    delete uxiref;
    delete uetaref;
  } // endfor
    
  values[0] = extr[0];
  values[1] = extr[1];
  values[2] = extr[2];
  values[3] = extr[3];
 }

/****************************************************************/
//
// for FEM_TVD
//
/****************************************************************/

void CheckWrongNeumannNodes(const TCollection *Coll, TFESpace2D *fespace,
			    int &N_neum_to_diri, int* &neum_to_diri,
			    int* &neum_to_diri_bdry, 
			    double* &neum_to_diri_param)
{
  const int max_entries = 50000;  
  int i, j, min_val, type;
  int N_Cells, N_V, diri_counter = 0, found, diri_counter_1 = 0;
  int boundary_vertices[4], tmp_diri[max_entries], tmp_bdry[max_entries];
  double x[4], y[4], eps = 1e-6, tmp_param[max_entries];
  TBaseCell *cell;
  TVertex *vertex;
  FE_type CurrentElement;

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
      vertex->GetCoords(x[j], y[j]);
      if ((std::abs(x[j]+3)<eps)//||(std::abs(y[j]+3)<eps)||(std::abs(y[j]-3)<eps)
	  || (std::abs(std::sqrt(x[j]*x[j]+y[j]*y[j])-1)<eps))
      {
	   boundary_vertices[j] = 1;
	   found++;
      }
    }
    // no cell with edge with vertex on the boundary
    if (found<2) 
	continue;
       // finite element on the mesh cell
    CurrentElement = fespace->get_fe_type(i);
    // the array which gives the mapping of the local to the global d.o.f.
    auto dof = fespace->GetGlobalDOF(i);
    switch(CurrentElement)
    {
	// P_1, Q_1
	case C_P1_2D_T_A:
	case C_Q1_2D_Q_A:
	case C_Q1_2D_Q_M:
	    for (j=0;j<N_V;j++)
	    {
		// vertex on the boundary
		if (boundary_vertices[j])
		{
		    if (CurrentElement==C_P1_2D_T_A)
			tmp_diri[diri_counter] = dof[j];
		    else
		    {
			if (j<2){
			    tmp_diri[diri_counter] = dof[j];
			}
			else
			{
			    if (j==2)
				tmp_diri[diri_counter] = dof[3];
			    else
				tmp_diri[diri_counter] = dof[2];
			}
		    }
		    if (diri_counter > max_entries)
		    {
			ErrThrow("tmp_diri too short !!!");
		    }
		    // inflow x = -3
		    if (std::abs(x[j]+3)<eps) 
		    {
			tmp_bdry[diri_counter] = 3;
			tmp_param[diri_counter] = (-y[j]+3)/6.0;
		    }
		    // circle
		    if (std::abs(std::sqrt(x[j]*x[j]+y[j]*y[j])-1)<eps) 
		    {
			tmp_bdry[diri_counter] = 4;
			// parameter does not matter, since b.c. equal to 1
			tmp_param[diri_counter] =  0;
		    }
		    diri_counter++;
		}
	    }
	    break;
	// P_2, Q_2
	case C_P2_2D_T_A:
	case C_Q2_2D_Q_A:
	case C_Q2_2D_Q_M:
            // loop over the edges
 	    for (j=0;j<N_V;j++)
	    {
              // check of edge j is on boundary  
              if (boundary_vertices[j] && boundary_vertices[(j+1)%N_V])
              {
		// check if this is a boundary edge
		type = cell->GetJoint(j)->GetType();
		if (!((type == BoundaryEdge)||(type == IsoBoundEdge)))
		  continue;
	        switch(j)
                {
                   case 0:
                     tmp_diri[diri_counter] = dof[0];
                     tmp_diri[diri_counter+1] = dof[1];
                     tmp_diri[diri_counter+2] = dof[2];
                   break;
                  case 1:
                     if (N_V==3)
                     {
                       tmp_diri[diri_counter] = dof[2];
                       tmp_diri[diri_counter+1] = dof[4];
                       tmp_diri[diri_counter+2] = dof[5];
                     }
                     else
                     {
                       tmp_diri[diri_counter] = dof[2];
                       tmp_diri[diri_counter+1] = dof[5];
                       tmp_diri[diri_counter+2] = dof[8];
                     }
                   break;
                  case 2:
                     if (N_V==3)
                     {
                       tmp_diri[diri_counter] = dof[5];
                       tmp_diri[diri_counter+1] = dof[3];
                       tmp_diri[diri_counter+2] = dof[0];
                     }
                     else
                     {
                       tmp_diri[diri_counter] = dof[8];
                       tmp_diri[diri_counter+1] = dof[7];
                       tmp_diri[diri_counter+2] = dof[6];
                     }
                   break;
                   case 3:
                     tmp_diri[diri_counter] = dof[6];
                     tmp_diri[diri_counter+1] = dof[3];
                     tmp_diri[diri_counter+2] = dof[0];
                   break;

                }
              
		if (diri_counter+2 > max_entries)
		{
			ErrThrow("tmp_diri too short !!!");
		}

		// inflow x = -3
		if ((std::abs(x[j]+3)<eps)&&(std::abs(x[(j+1)%N_V]+3)<eps)) 
		{
		    tmp_bdry[diri_counter] = 3;
		    tmp_bdry[diri_counter+1] = 3;
		    tmp_bdry[diri_counter+2] = 3;
		    tmp_param[diri_counter] = (-y[j]+3)/6.0;
		    tmp_param[diri_counter+2] = (-y[(j+1)%N_V]+3)/6.0;
		    tmp_param[diri_counter+1] = (tmp_param[diri_counter] +  tmp_param[diri_counter+2])/2.0;
		}
		// circle
		if ((std::abs(std::sqrt(x[j]*x[j]+y[j]*y[j])-1)<eps) && 
		    (std::abs(std::sqrt(x[(j+1)%N_V]*x[(j+1)%N_V]+y[(j+1)%N_V]*y[(j+1)%N_V])-1)<eps))
		{
		    tmp_bdry[diri_counter] = 4;
		    tmp_bdry[diri_counter+1] = 4;
		    tmp_bdry[diri_counter+2] = 4;
		    // parameter does not matter, since b.c. equal to 1
		    tmp_param[diri_counter] = 0;
		    tmp_param[diri_counter+1] = 0;
		    tmp_param[diri_counter+2] = 0;
		}
		diri_counter +=3;
	      }
	    }
	    break;
	// P_3, Q_3
	case C_P3_2D_T_A:
	case C_Q3_2D_Q_A:
	case C_Q3_2D_Q_M:
            // loop over the edges
 	    for (j=0;j<N_V;j++)
	    {
              // check of edge j is on boundary  
              if (boundary_vertices[j] && boundary_vertices[(j+1)%N_V])
              {
		// check if this is a boundary edge
		type = cell->GetJoint(j)->GetType();
		if (!((type == BoundaryEdge)||(type == IsoBoundEdge)))
		  continue;

               // P3: local dof 0, 1, 2, 3 are on the boundary
               // Q3: local dof 0, 1, 2, 3 are on the boundary
	        switch(j)
                {
                   case 0:
                     tmp_diri[diri_counter] = dof[0];
                     tmp_diri[diri_counter+1] = dof[1];
                     tmp_diri[diri_counter+2] = dof[2];
		     tmp_diri[diri_counter+3] = dof[3];
                   break;
                  case 1:
                     if (N_V==3)
                     {
                       tmp_diri[diri_counter] = dof[3];
                       tmp_diri[diri_counter+1] = dof[6];
                       tmp_diri[diri_counter+2] = dof[8];
		       tmp_diri[diri_counter+3] = dof[9];
                     }
                     else
                     {
                       tmp_diri[diri_counter] = dof[3];
                       tmp_diri[diri_counter+1] = dof[7];
                       tmp_diri[diri_counter+2] = dof[11];
		       tmp_diri[diri_counter+3] = dof[15];
                     }
                   break;
                  case 2:
                     if (N_V==3)
                     {
                       tmp_diri[diri_counter] = dof[9];
                       tmp_diri[diri_counter+1] = dof[7];
                       tmp_diri[diri_counter+2] = dof[4];
                       tmp_diri[diri_counter+3] = dof[0];
		     }
                     else
                     {
                       tmp_diri[diri_counter] = dof[15];
                       tmp_diri[diri_counter+1] = dof[14];
                       tmp_diri[diri_counter+2] = dof[13];
			tmp_diri[diri_counter+3] = dof[12];
                     }
                   break;
                   case 3:
                     tmp_diri[diri_counter] = dof[12];
                     tmp_diri[diri_counter+1] = dof[8];
                     tmp_diri[diri_counter+2] = dof[4];
		     tmp_diri[diri_counter+3] = dof[0];
                   break;
                }
              
		if (diri_counter+3 > max_entries)
		{
			ErrThrow("tmp_diri too short !!!");
		}

		// inflow x = -3
		if ((std::abs(x[j]+3)<eps)&&(std::abs(x[(j+1)%N_V]+3)<eps)) 
		{
		    tmp_bdry[diri_counter] = 3;
		    tmp_bdry[diri_counter+1] = 3;
		    tmp_bdry[diri_counter+2] = 3;
		    tmp_bdry[diri_counter+3] = 3;
		    tmp_param[diri_counter] = (-y[j]+3)/6.0;
		    tmp_param[diri_counter+3] = (-y[(j+1)%N_V]+3)/6.0;
		    tmp_param[diri_counter+1] = (2*tmp_param[diri_counter] +  tmp_param[diri_counter+3])/3.0;
		    tmp_param[diri_counter+2] = (tmp_param[diri_counter] +  2*tmp_param[diri_counter+3])/2.0;
		}
		// circle
		if ((std::abs(std::sqrt(x[j]*x[j]+y[j]*y[j])-1)<eps) && 
		    (std::abs(std::sqrt(x[(j+1)%N_V]*x[(j+1)%N_V]+y[(j+1)%N_V]*y[(j+1)%N_V])-1)<eps))
		{
		    tmp_bdry[diri_counter] = 4;
		    tmp_bdry[diri_counter+1] = 4;
		    tmp_bdry[diri_counter+2] = 4;
		    tmp_bdry[diri_counter+3] = 4;
		    // parameter does not matter, since b.c. equal to 1
		    tmp_param[diri_counter] = 0;
		    tmp_param[diri_counter+1] = 0;
		    tmp_param[diri_counter+2] = 0;
		    tmp_param[diri_counter+3] = 0;
		}
		diri_counter +=4;
	      }
	    }
	    break;
	default:
	    ErrThrow("CheckNeumannNodesForVelocity not implemented for element ",
               CurrentElement);
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

  //Output::print("CheckNeumannNodesForVelocity: N_neum_to_diri ", diri_counter_1);
  N_neum_to_diri = diri_counter_1;
  // allocate array for the indices
  neum_to_diri = new int[diri_counter_1];
  // allocate array for the corresponding boundary numbers
  neum_to_diri_bdry = new int[diri_counter_1];
  // allocate array for the corresponding boundary parameters
  neum_to_diri_param = new double[diri_counter_1];
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
      neum_to_diri_bdry[i] = tmp_bdry[found];
      neum_to_diri_param[i] = tmp_param[found];
      tmp_diri[found] = -1;
  }
/*
  for (i=0;i<diri_counter_1;i++)
  {
      Output::print(i, " ", neum_to_diri[i], " ", neum_to_diri_bdry[i],
	     " ", neum_to_diri_param[i]);
  }
*/
}

void ComputeDifferenceToCoarseLevel(const TCollection *Coll_fine,
				    const TCollection *Coll_coarse,
				    TFEFunction2D *u_fine, 
				    TFEFunction2D *u_coarse)
{
    int i, j, k, N_Cells, N_Edges, coarse_no;
    double x, y, x_c, y_c, val_fine[4], val_coarse[4], c1err = -1, c1err_coarse = -1;
    double x_err, y_err, x_err_c, y_err_c;
    const TBaseCell *cell, *parent;
    
    // number of cells
    N_Cells = Coll_fine->GetN_Cells();
    
    // loop over all edges
    for(i=0;i<N_Cells;i++)
    {
	// cell
	cell = Coll_fine->GetCell(i);
	// parent cell
	parent = cell->GetParent();
	coarse_no = Coll_coarse->get_cell_index(parent);
	//Output::print(coarse_no);
	// number of edges
	N_Edges=cell->GetN_Edges();
	for (j=0;j<N_Edges;j++)
	{
	    cell->GetVertex(j)->GetCoords(x, y);
	    u_fine->FindGradientLocal(cell, i, x, y, val_fine);
	    u_coarse->FindGradientLocal(parent, coarse_no, x, y, val_coarse);
	    if (std::abs(val_fine[0] - val_coarse[0]) > c1err)
	    {
		c1err = std::abs(val_fine[0] - val_coarse[0]);
		x_err = x;
		y_err = y;
	    }
	    for (k=0;k<N_Edges;k++)
	    {
		parent->GetVertex(k)->GetCoords(x_c, y_c);
		if ((std::abs(x_c -x ) < 1e-6) && (std::abs(y_c -y ) < 1e-6))
		{
		    if (std::abs(val_fine[0] - val_coarse[0]) > c1err_coarse)
		    {
			c1err_coarse = std::abs(val_fine[0] - val_coarse[0]);
			x_err_c = x;
			y_err_c = y;
		    }
		}
	    }
	}
    } 

    //Output::print("C1 error f ", c1err, " \\& ", x_err, ",", y_err);
    //Output::print(" C1 error c ", c1err_coarse, " \\& ", x_err_c, ",", y_err_c);

    Output::print("C1 error f & ", c1err, " &  ( ", x_err, ",", y_err, ")", "\\\\\\hline");
    Output::print("C1 error c ", " & ", c1err_coarse, " &  ( ", x_err_c, ",", y_err_c, ")", "\\\\\\hline");
}
void ComputeCutLines_X(const TCollection *Coll, const TFEFunction2D *ufct, int level)
{
    double h, val[3], x, y, *cutvalues, y0=0;
  int i, j, N_Cells, bound_points = 10001;
  
  h = 3.0/(bound_points-1);

  cutvalues = new double[6*bound_points];
  memset(cutvalues , 0 , 6*bound_points*sizeof(double));

  N_Cells = Coll->GetN_Cells();
  for(i=0;i<N_Cells;i++)
  {   
    auto cell = Coll->GetCell(i);
    x = -1;
    for (j=0;j<bound_points;j++)
    {
	y = y0 + h * j;
	cutvalues[j] = y;
	if(cell->PointInCell(parmoon::Point(x,y)))
	{
	    ufct->FindGradientLocal(cell, i, x, y, val);
	    cutvalues[j+bound_points] = val[0];
	}
    }
    x = 0;
    for (j=0;j<bound_points;j++)
    {
	y = y0 + h * j;
	if(cell->PointInCell(parmoon::Point(x,y)))
	{
	    ufct->FindGradientLocal(cell, i, x, y, val);
	    cutvalues[j+2*bound_points] = val[0];
	}
    }
    x = 1;
    for (j=0;j<bound_points;j++)
    {
	y = y0 + h * j;
	if(cell->PointInCell(parmoon::Point(x,y)))
	{
	    ufct->FindGradientLocal(cell, i, x, y, val);
	    cutvalues[j+3*bound_points] = val[0];
	}
    }
    x = 2;
    for (j=0;j<bound_points;j++)
    {
	y = y0 + h * j;
	if(cell->PointInCell(parmoon::Point(x,y)))
	{
	    ufct->FindGradientLocal(cell, i, x, y, val);
	    cutvalues[j+4*bound_points] = val[0];
	}
    }
    x = 4;
    for (j=0;j<bound_points;j++)
    {
	y = y0 + h * j;
	if(cell->PointInCell(parmoon::Point(x,y)))
	{
	    ufct->FindGradientLocal(cell, i, x, y, val);
	    cutvalues[j+5*bound_points] = val[0];
	}
    }
  }

  for (j=0;j<bound_points;j++)
  {
      Output::print("cutx ", level, " ", cutvalues[j], " ", cutvalues[j+bound_points], " ", cutvalues[j+2*bound_points], " ", cutvalues[j+3*bound_points], " ", cutvalues[j+4*bound_points], " ", cutvalues[j+5*bound_points]);
  }
}

void ComputeCutLines_Y(const TCollection *Coll, const TFEFunction2D *ufct, int level)
{
    double h, val[3], x, y, *cutvalues;
  int i, j, N_Cells, bound_points = 20001;
    TBaseCell *cell;
  
  h = 10.0/(bound_points-1);

  cutvalues = new double[3*bound_points];
  memset(cutvalues , 0 , 3*bound_points*sizeof(double));
  
  N_Cells = Coll->GetN_Cells();
  for(i=0;i<N_Cells;i++)
  {   
    cell = Coll->GetCell(i);
    y = 0;
    for (j=0;j<bound_points;j++)
    {
	x = -2 + h * j;
	cutvalues[j] = x;
	if(cell->PointInCell(parmoon::Point(x,y)))
	{
	    ufct->FindGradientLocal(cell, i, x, y, val);
	    cutvalues[j+bound_points] = val[0];
	}
    }
    y = 1;
    for (j=0;j<bound_points;j++)
    {
	x = -2 + h * j;
	if(cell->PointInCell(parmoon::Point(x,y)))
	{
	    ufct->FindGradientLocal(cell, i, x, y, val);
	    /*if (cutvalues[j+2*bound_points]!=4711)
	    {
		ErrThrow("belegt");
	    }
	    else*/
	    cutvalues[j+2*bound_points] = val[0];
	}
    }
  }

  for (j=0;j<bound_points;j++)
  {
      Output::print("cuty ", level, " ", cutvalues[j], " ",
                    cutvalues[j+bound_points], " ", cutvalues[j+2*bound_points]);
  }
}

void ComputeCutLines_epsY(const TCollection *Coll, TFEFunction2D *ufct, int level)
{
    double h, val[3], x, y, *cutvalues;
  int i, j, N_Cells, bound_points = 20001;
    TBaseCell *cell;
    double eps=DIFFUSION;
  
  h = 10.0/(bound_points-1);

  cutvalues = new double[3*bound_points];
  memset(cutvalues , 0 , 3*bound_points*sizeof(double));

  N_Cells = Coll->GetN_Cells();
  for(i=0;i<N_Cells;i++)
  {   
    cell = Coll->GetCell(i);
    y = 1-std::sqrt(eps);
    for (j=0;j<bound_points;j++)
    {
	x = -2 + h * j;
	cutvalues[j] = x;
	if(cell->PointInCell(parmoon::Point(x,y)))
	{
	    ufct->FindGradientLocal(cell, i, x, y, val);
	    cutvalues[j+bound_points] = val[0];
	}
    }
    y = 1+std::sqrt(eps);
    for (j=0;j<bound_points;j++)
    {
	x = -2 + h * j;
	if(cell->PointInCell(parmoon::Point(x,y)))
	{
	    ufct->FindGradientLocal(cell, i, x, y, val);
	    cutvalues[j+2*bound_points] = val[0];
	}
    }
  }

  for (j=0;j<bound_points;j++)
  {
      Output::print("cutyeps", level, " ", cutvalues[j], " ",
                    cutvalues[j+bound_points], " ", cutvalues[j+2*bound_points]);
  }
}
void ComputeCutLines_eps_radial(const TCollection *Coll, TFEFunction2D *ufct, int level)
{
    double h, val[3], x, y, *cutvalues, tmp, r;
  int i, j, N_Cells, bound_points = 10001;
    TBaseCell *cell;
    double eps=DIFFUSION;
  
  h = 2*M_PI/(bound_points-1);

  cutvalues = new double[10*bound_points];
  memset(cutvalues , 0 , 10*bound_points*sizeof(double));

  N_Cells = Coll->GetN_Cells();
  for(i=0;i<N_Cells;i++)
  {   
    cell = Coll->GetCell(i);
    r = (1+eps)*(1+eps);
    for (j=0;j<bound_points;j++)
    {
	x = r*std::cos(h * j-M_PI);
	y = r*std::sin(h*j-M_PI);
	cutvalues[j] = h*j -M_PI;
	cutvalues[j+bound_points] = x;
	cutvalues[j+2*bound_points] = y;
	if(cell->PointInCell(parmoon::Point(x,y)))
	{
	    ufct->FindGradientLocal(cell, i, x, y, val);
	    cutvalues[j+3*bound_points] = val[0];
	}
    }
    tmp=std::pow(eps,2.0/3.0);
    r = (1+tmp)*(1+tmp);
    for (j=0;j<bound_points;j++)
    {
	x = r*std::cos(h * j-M_PI);
	y = r*std::sin(h*j-M_PI);
	cutvalues[j+4*bound_points] = x;
	cutvalues[j+5*bound_points] = y;
	if(cell->PointInCell(parmoon::Point(x,y)))
	{
	    ufct->FindGradientLocal(cell, i, x, y, val);
	    cutvalues[j+6*bound_points] = val[0];
	}
    }
    tmp=std::sqrt(eps);
    r = (1+tmp)*(1+tmp);
    for (j=0;j<bound_points;j++)
    {
	x = r*std::cos(h * j-M_PI);
	y = r*std::sin(h*j-M_PI);
	cutvalues[j+7*bound_points] = x;
	cutvalues[j+8*bound_points] = y;
	if(cell->PointInCell(parmoon::Point(x,y)))
	{
	    ufct->FindGradientLocal(cell, i, x, y, val);
	    cutvalues[j+9*bound_points] = val[0];
	}
    }

  }

  for (j=0;j<bound_points;j++)
  {
      Output::print("cutradeps", level, " ", cutvalues[j], " ",
            cutvalues[j+bound_points], " ", cutvalues[j+2*bound_points], " ",
            cutvalues[j+3*bound_points], " ", cutvalues[j+4*bound_points], " ",
            cutvalues[j+5*bound_points], " ", cutvalues[j+6*bound_points], " ",
            cutvalues[j+7*bound_points], " ", cutvalues[j+8*bound_points], " ",
            cutvalues[j+9*bound_points]);
  }
}

double ComputeBdLayer(const TFEFunction2D& ufct, double xcut)
{
  // Consider a line starting at (xcut, y0) going in the positive y-direction.
  // This should be behind the cylinder where the solution is one and then 
  // (going up) zero.
  // We will go along that line to find the first point where ufct is below 
  // valmax, then we go further along that line to find the first point where 
  // ufct is below valmin. Remembering where these two transitions appeared, we
  // are able to compute the "layer width". This is described also in 
  // Augustin, Caiazzo, Fiebach, Fuhrmann, John, Linke, Umla: "An assessment of 
  // discretizations for convection-dominated convection–diffusion equations",
  // 2011.
  // We avoid having to search all cells for every point, which would be very
  // expensive, by only searching the neighboring cells from where the previous
  // point was found.
  constexpr int bound_points = 100001;
  constexpr double valmin = 0.1;
  constexpr double valmax = 0.9;
  constexpr double y0 = 0.;
  constexpr double Width = 6;

  double ystart = y0;
  double yend = Width;
  double h = 0.5 * Width / (bound_points-1);
  parmoon::Point p(xcut, ystart);

  // First, we have to determine the cell in which the starting point (xcut, 0)
  // lies. Therefore, we have to do a loop over all cells and determine the
  // cell. This cell is called starting_cell_nr.
  auto collection = ufct.GetFESpace()->GetCollection();
  auto n_cells = collection->GetN_Cells();
  int starting_cell_nr = 0;
  for(int cell_i = 0; cell_i < n_cells; ++cell_i)
  {
    auto current_cell = collection->GetCell(cell_i);
    bool is_point_in_cell = current_cell->PointInCell(p);
    if (is_point_in_cell)
    {
      starting_cell_nr = cell_i;
      break;
    }
  }
  
  // Now a loop over all y values can be performed. In each iteration we have to
  // determine the cell in which the current point lies. For this purpose the
  // current cell number and the current cell are determined. In the beginning
  // the current cell number is of course the previously computed
  // starting_cell_nr.
  int current_cell_nr = starting_cell_nr;
  int jstart = 0;
  for(int j = 0; j < bound_points; j++)
  {
    p.y = y0 + h * j;
    collection->find_cell(p, current_cell_nr);
    double val;
    auto current_cell = collection->GetCell(current_cell_nr);
    ufct.FindValueLocal(current_cell, current_cell_nr, p.x, p.y, &val);
    if(val < valmax)
    {
      ystart = p.y;
      jstart = j;
      break;
    }
  }

  for(int j=jstart+1;j<bound_points;j++)
  {
    p.y = y0 + h * j;
    collection->find_cell(p, current_cell_nr);
    double val;
    auto current_cell = collection->GetCell(current_cell_nr);
    ufct.FindValueLocal(current_cell, current_cell_nr, p.x, p.y, &val);
    if(val < valmin)
    {
      yend = p.y;
      break;
    }
  }
  auto lw = std::abs(yend - ystart);

  if (lw >= std::abs(Width - ystart))
  {
    Output::warn("Hemker1996::compute_layer_width", "Failed to compute layer ",
        "width at x=",xcut,". I did not find an end point.");
  }
  if (lw <= h)
  {
    Output::warn("Hemker1996::compute_layer_width", "Layer width equals step ",
        "width at x=",xcut,"!");
  }
  return lw;
 }

/*******************************************************************************/
// computes the data for the evaluation of the solution 
/*******************************************************************************/
void hemker_postprocessing(ConvectionDiffusion<2> & cd2d)
{
  auto & uh = cd2d.get_function();

  // Compute cutlines
  std::array<double, 3> xcut{1.1, 4.0, 9.0}; // add more as desired
  for(auto x : xcut)
  {
    auto lw = ComputeBdLayer(uh, x);
    Output::print("layer width at x=", x, ": ", std::setprecision(12), lw);
  }

  // Compute (an approximation of) the global minimal and maximal value of uh
  uh.PrintMinMax("solution", false);
  // Compute the global minimal and maximal value of the FE-vector that is
  // cell wise transformed to Pk elements of the same order (if not Pk elements
  // are used already)
  uh.PrintMinMax("Pk solution vector", true);

  // Compute mean oscillation using cell wise minimum and maximum of uh
  constexpr double umin = 0.;
  constexpr double umax = 1.;
  double osc_mean = uh.compute_mean_oscillation(umin, umax, false);
  Output::print("mean oscillation ", std::setprecision(12), osc_mean);
  // Compute mean oscillation using cell wise minimum and maximum of the
  // FE-vector transformed to Pk elements of the same order
  osc_mean = uh.compute_mean_oscillation(umin, umax, true);
  Output::print("mean oscillation using Pk nodal functionals ",
      std::setprecision(12), osc_mean);
  // Compute mean oscillations of the FE-vector itself
  if (!cd2d.get_space()->is_discontinuous())
  {
    auto& entries = cd2d.get_solution().get_entries_vector();
    double osc_mean_vec = 0.;
    for(double e : entries)
    {
      osc_mean_vec += std::max(0., e-umax) + std::max(0., umin-e);
    }
    osc_mean_vec /= entries.size();
    Output::print("mean vector oscillation ", std::setprecision(12),
        osc_mean_vec);
  }
}
