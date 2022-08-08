// ======================================================================
// sharp characteristic interior layer
// Knopp, Lube, Rapin, CMAME 2002
// ======================================================================
#define __THREE_INTERIOR_LAYERS_INLET_M__

#include <IsoBoundEdge.h>
#include <BoundComp.h>

double DIFFUSION;

void ExampleFile()
{
  Output::print<1>("Example: ThreeInteriorLayers_inlet_M.h");
}
// exact solution
auto& Exact = unknown_solution_2d;

// kind of boundary condition (for FE space needed)
void BoundCondition(int i, double, BoundCond &cond)
{
    if (i==3)
	   cond = NEUMANN;
    else
	   cond = DIRICHLET;
    // just to fix the point (0,0), (1,0) to Dirichlet
    if (i==0)
       cond = DIRICHLET;
    if (i==2)
       cond = DIRICHLET;	    
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
    double regu_param = 1e-3;
    value = 0;
  
  if (BdComp==0)
  {
      if ((Param>0.375-regu_param)&& (Param<=0.375))
         value = (Param-0.375)/regu_param + 1;
      if ((Param>0.375)&& (Param<=0.5))
         value = -0.75*(Param-0.5)/0.125 + 0.25;
      if ((Param>0.5)&& (Param<=0.625))
         value = 0.25*(Param-0.625)/0.125 + 0.5;
      if ((Param>0.625)&& (Param<=0.625 + regu_param))
         value = -0.5*(Param-0.625)/regu_param + 0.5;
  }
}

void BilinearCoeffs(int n_points, const double *X, const double *Y,
                    const double *const*, double **coeffs)
{
  double eps = DIFFUSION;

  for(int i = 0; i < n_points; i++)
  {
    double x = X[i];
    double y = Y[i];
    coeffs[i][0] = eps;
    coeffs[i][1] = -y;
    coeffs[i][2] = x;
    coeffs[i][3] = 0;
    coeffs[i][4] = 0;
    coeffs[i][5] = std::sqrt(x*x + y*y);
  }
}

/*****************************************************************
 * POST PROCESSING
 ***************************************************************** */
/** compute curve of the outflow boundary */
void ComputeOutflowBoundary(const TFEFunction2D *ufct)
{
  const int max_bound_points = 100001;
  int i,j,k, N_Cells;
  TBaseCell *cell;
  const TCollection *Coll;
  double val[5];
  TJoint *joint;
  const TBoundComp *BoundComp;
  const TBoundEdge *boundedge;

  int comp, found, N_Edges;
  double x, y;
  double value, eps=1e-10;
  double y_coord[max_bound_points], uval[max_bound_points], min;
  int bound_points, index;

  Coll = ufct->GetFESpace2D()->GetCollection();
  N_Cells = Coll->GetN_Cells();

  bound_points = 0;
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_Edges=cell->GetN_Edges();
    found = 0;
    for(j=0;j<N_Edges;j++)              // loop over all edges of cell
    {                                   // find edges on boundary part 3  
      joint=cell->GetJoint(j);          // this is x=0   
      if ((joint->GetType() == BoundaryEdge)||
          (joint->GetType() == IsoBoundEdge)) // boundary edge 
      {
        
        boundedge = (TBoundEdge *)joint;  
        BoundComp = boundedge->GetBoundComp();  // get boundary component
        comp=BoundComp->GetID();              // boundary id 
        if (comp==3)
        {
          found = 1;
          break;
        }
      }
    }
    // no edge on outflow boundary found
    if (!found) continue;   
    
    // loop over the edges of the mesh cell 
    for (j=0;j<N_Edges;j++)
    {
      // check all vertices if they are on the boundary 
      x = cell->GetVertex(j)->GetX();
      // point not on boundary
      if (fabs(x)>eps)
         continue;    
      y = cell->GetVertex(j)->GetY();
      //check if this node is already in the array of boundpoints
      found = 0;
      for (k=bound_points-1;k>=0; k--)
      {
         if (fabs(y- y_coord[k]) < eps)
         {
            found = 1;
            break;
         }
      }
      if (found)
         continue;
      // new node
      y_coord[bound_points] = y;
      bound_points++;      
      if ( bound_points > max_bound_points)
      {
         Output::print("ThreeInteriorLayers_inlet_M.h: maximal number of boundary points reached !!!");
         exit(4711);
      }
      ufct->FindGradientLocal(cell, i, x, y, val);
      uval[bound_points-1] = val[0];
    } // endfor 
  } // endfor
  // order the arrays

  for (i=0;i<bound_points; i++)
  {
     min = 1e6;
     for(j=i;j<bound_points; j++)
     {
        if (y_coord[j]< min)
        {
           min = y_coord[j];
           index = j;
        }
     }
     // change the entries
     y_coord[index] = y_coord[i];
     y_coord[i] = min;
     value = uval[i];
     uval[i] = uval[index];
     uval[index] = value;
  }
  
  double loc_ext[20];
  loc_ext[0] = loc_ext[1] = loc_ext[2] = -1.0;    // coordinates 
  loc_ext[3] = loc_ext[5] = -4711.0; // values
  loc_ext[4] = 4711.0;
  
  for (i=0;i<bound_points; i++)
  {
     //OutPut("outflow " << level << " "  <<   y_coord[i] << 
     //       " " <<  uval[i] << endl);
     if (i>0 && i < bound_points -1)
     {
       if (y_coord[i] > 0.33 && y_coord[i] < 0.44 && uval[i] > loc_ext[3] && uval[i]>= uval[i-1] &&  uval[i]>= uval[i+1] )
       {
	 loc_ext[0] = y_coord[i];
	 loc_ext[3] = uval[i];
       }
       if ( y_coord[i] > 0.45 && y_coord[i] < 0.55 && uval[i] < loc_ext[4] && uval[i]<= uval[i-1] &&  uval[i]<= uval[i+1] )
       {
	 loc_ext[1] = y_coord[i];
	 loc_ext[4] = uval[i];
       }
       if ( y_coord[i] > 0.56 && y_coord[i] < 0.7 && uval[i] > loc_ext[5] && uval[i]>= uval[i-1] &&  uval[i]>= uval[i+1])
       {
	 loc_ext[2] = y_coord[i];
	 loc_ext[5] = uval[i];
       }
      }
  }
 
  Output::print("dof ",ufct->GetLength(), " extr1 ",loc_ext[0], " ", loc_ext[3], 
     " extr2 ", loc_ext[1], " ", loc_ext[4],
     " extr3 ", loc_ext[2], " ", loc_ext[5]);
  
   int no_extr = 0;

   for (i=0;i<bound_points; i++)
   {
     //OutPut("outflow " << level << " "  <<   y_coord[i] << 
     //       " " <<  uval[i] << endl);
    if (uval[i] < 0.2)
      continue;

      // local maximum 
      if (uval[i]>= uval[i-1] &&  uval[i]>= uval[i+1] )
       {
	 loc_ext[no_extr] = y_coord[i];
	 loc_ext[no_extr+10] = uval[i];
	 no_extr++;
       }

      // local minimum
      if (uval[i]<= uval[i-1] &&  uval[i]<= uval[i+1] )
       {
	 loc_ext[no_extr] = y_coord[i];
	 loc_ext[no_extr+10] = uval[i];
	 no_extr++;
       }
       if (no_extr == 10)
       {
	 Output::print("maximal number of extrema reached, stop further calculations ", no_extr);
	 break;
       }
  }

  Output::print("dof " ,ufct->GetLength() ," no extr " ,no_extr );
  for (i=0; i<no_extr; i++)
    Output::print(" extr[" ,i ,"] " ,loc_ext[i] ," " ,loc_ext[i+10]);

/*
  // data output for reference curve
  bound_points = 100001;
  int ii, last_cell = -1;
  double values[3];
  double h = 1.0/(bound_points-1);
  double h0 = h/10.0;

  x=0.0;

  for (i=0;i<bound_points; i++)
  {
    y = i*h;    
    if (last_cell >=0)
    {
      cell = Coll->GetCell(last_cell);
      // new point is in the same cell as former point
      if (cell->PointInCell(parmoon::Point(0.0,y)))
      {
	ufct->FindGradientLocal(cell, last_cell, x, y, values);
	OutPut("outflow " << level << " "  <<   y << " " <<  values[0] << endl);
	continue;
      }
    }
    // loop over the cells
    for(ii=0;ii<N_Cells;ii++)
    {
      cell = Coll->GetCell(ii);
      // this approach assumes that the fe function is continuous
      if (cell->PointInCell(parmoon::Point(0.0,y)))
      {
	ufct->FindGradientLocal(cell, ii, x, y, values);
	OutPut("outflow " << level << " "  <<   y << " " <<  values[0] << endl);
	last_cell = ii;
	// loop over the cells
	break;
      }
    }
  }
*/
  //OutPut(lay[0] << " " <<  lay[1] << " " << lay[2] << " " << lay[3] << endl);
 /* 
  double h, values[3];
  x = 0;
  bound_points = 100001;
  h = 1.0/(bound_points-1);
  for (i=0;i<bound_points; i++)
  {
      y = i*h;
      ufct->FindGradient(x,y,values);
      OutPut("outflow " << level << " " <<  y << 
            " " <<  values[0] << endl);
  }
  */
}

void ComputeExtremalValues(int N, const double *sol, double  *values)
{
   int i;
   double max, min;

   min = 1e10;
   max = -1e10;
   
   for(i=0;i<N;i++)
   {
      if(sol[i]-1 > max)
         max = sol[i]-1;
      if(sol[i] < min)
         min = sol[i];
   }

   values[0] = min;
   values[1] = max;
}

void ComputeDataForEvaluationOfSolution(const TFEFunction2D *u)
//					ofstream &errfile,ofstream &layerfile)
{
  double errors[5];
  
  auto n_dof = u->GetLength();
  const double * sol = u->GetValues();

  // over and undershoots
  ComputeExtremalValues(n_dof,sol,errors);
  Output::print(setprecision(8), "undershoots ", errors[0], " overshoots ",
                errors[1]);
  ComputeOutflowBoundary(u);
  return;
}
 
void three_interior_layers_inlet_M_postprocessing(ConvectionDiffusion<2> & cd2d)
{
  auto & u = cd2d.get_function();
  ComputeDataForEvaluationOfSolution(&u);
}
