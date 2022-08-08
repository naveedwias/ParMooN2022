// =======================================================================
// @(#)MainUtilities.C        1.11 07/03/00
//
// Purpose:     main utilities
//
// Authors:     Gunar Matthies  18.10.99
//
// =======================================================================

#ifdef _MPI
#include <mpi.h>
#include "ParFECommunicator3D.h"
#endif

#include <MainUtilities.h>
#include <Database.h>
#include <Convolution.h>
#include <LinAlg.h>
#include <ConvDiff.h>
#include "BaseCell.h"

#include "FEDatabase.h"
#ifdef __2D__
#include <FEVectFunct2D.h>
#include <NSE2DSUPG.h>
#include "SquareMatrix2D.h"
#endif

#ifdef __3D__
#include <FEVectFunct2D.h>
#include <FEVectFunct3D.h>
#include <CommonRoutineTNSE3D.h>
#include <BoundFace.h>
#include <QuadratureFormulaDatabase.h>
#include <RefTrans3D.h>
#include "SquareMatrix3D.h"
#endif

#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
// #include <malloc.h>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <stdio.h>
#include <cmath>

#include <set>
#include <unordered_map>

#ifdef __MAC64__
#include <malloc/malloc.h> // defines _MALLOC_MALLOC_H_
#else
#include <malloc.h>
#endif

double GetTime()
{
  struct rusage usage;
  double ret;
    
  if(getrusage(RUSAGE_SELF, &usage) == -1) 
    cerr << "Error in GetTime!" << endl;

  ret = ((double) usage.ru_utime.tv_usec)/1000000;

  ret += usage.ru_utime.tv_sec;

  return ret;
}

int GetMemory()
{
#ifdef _MALLOC_MALLOC_H_
 struct mstats info;
 info = mstats();
 return info.bytes_free+info.chunks_free;
#else  
  struct mallinfo MALLINFO;
  MALLINFO = mallinfo();
  return MALLINFO.usmblks+MALLINFO.uordblks;
#endif
}

//Print all memory information available through mallinfo.
void display_mallinfo(const std::string& program_part)
{
#ifdef _MALLOC_MALLOC_H_
  Output::print("Memory usage info called in program part: ", program_part);
  Output::print(GetMemory());  
#else
   
  struct mallinfo mi;
  
  mi = mallinfo();
  Output::print("Memory usage info called in program part: ", program_part);
  Output::print("Total non-mmapped bytes (arena):       ",mi.arena);
  Output::print("# of free chunks (ordblks):            ",mi.ordblks);
  Output::print("# of free fastbin blocks (smblks):     ",mi.smblks);
  Output::print("# of mapped regions (hblks):           ",mi.hblks);
  Output::print("Bytes in mapped regions (hblkhd):      ",mi.hblkhd);
  Output::print("Max. total allocated space (usmblks):  ",mi.usmblks);
  Output::print("Free bytes held in fastbins (fsmblks): ",mi.fsmblks);
  Output::print("Total allocated space (uordblks):      ",mi.uordblks);
  Output::print("Total free space (fordblks):           ",mi.fordblks);
  Output::print("Topmost releasable block (keepcost):   ",mi.keepcost);
#endif
}

#ifdef __2D__
int VertexNumber(int IEH, int NVE)
{
  int ret;

  if(NVE==4)
  {
    switch(IEH)
    {
      case 2: ret=3; break;
      case 3: ret=2; break;
      default: ret=IEH;
    }
  }
  else
  {
    ret = IEH;
  }

  return ret;
}

void StreamFunction(const TFESpace2D *velo, double *u1, double *u2,
                    const TFESpace2D *stream, double *psi)
{
  int IEH,IEL,k,l, N_Cells, NVE, NVE2, NVT, m, l1;
  int ListPointer, ListInput;
  int *KVIND;
  int IH, JVE, IHV, IND, INDH, IVTH, IVT, IELH;
  double px1, px2, py1, py2, dn1, dn2;
  double v1, v2;
  int *JointDOFs;

  auto Coll = velo->GetCollection();
  N_Cells = Coll->GetN_Cells();
  NVT = stream->get_n_dof();
  KVIND = new int [NVT];
  memset(KVIND, 0, sizeof(int)*NVT);
  memset(psi, 0, sizeof(double)*NVT);
  const TBaseCell ** CellList = new const TBaseCell* [N_Cells];
  memset(CellList, 0, sizeof(const TBaseCell*)*N_Cells);

  for(IEH=0;IEH<N_Cells;IEH++)
    Coll->GetCell(IEH)->SetClipBoard(IEH);

  l = stream->GetGlobalDOF(0)[0];
  psi[l] = 0.0;
  KVIND[l] = 1;

  const TBaseCell * cell = Coll->GetCell(0);
  CellList[0] = cell;
  l = -(10+0);
  cell->SetClipBoard(l);
  ListPointer = 0; 
  ListInput = 1;

  while(ListPointer<ListInput)
  {
    cell = CellList[ListPointer];

    l = cell->GetClipBoard();
    IEL = IEH = -(10+l);
    
    // Output::print("cell: ", IEL);
    IH = 0;
    NVE2 = cell->GetN_Edges();

    for(k=0;k<NVE2;k++)
    {
      JVE = stream->GetGlobalDOF(IEL)[VertexNumber(k, NVE2)];
      IHV = KVIND[JVE];
      IH = IH+IHV;
      if(IHV >= 1) IND = k;
    } // endfor k

    if((IH >= NVE2) || (IH == 0))
    {
      // already all value computed in this cell
      cell->SetClipBoard(-1);

      for(k=0;k<NVE2;k++)
      {
        auto neigh = cell->GetJoint(k)->GetNeighbour(cell);
   
        if(neigh)
        {
          if(neigh->GetClipBoard()>=0)
          { 
            // neigh was not handled
            IELH = neigh->GetClipBoard();
            IH = 0;
            NVE = neigh->GetN_Edges();
  
            for(l=0;l<NVE;l++)
            {
              JVE = stream->GetGlobalDOF(IELH)[VertexNumber(l, NVE)];
              IHV = KVIND[JVE];
              IH = IH+IHV;
              if(IHV>=1) INDH = l;
            } // endfor l
    
            if((IH<NVE) && (IH>0))
            {
              // add cell to list since not all values were calculated
              CellList[ListInput] = neigh;
              l = -(10+neigh->GetClipBoard());
              if(l<0)
                neigh->SetClipBoard(l);
              ListInput++;
            }
          } // >=0
        } // if neigh
      } // endfor k;
    }
    else
    {
      // at least one value must be calculated
      // cout << "cell number: " << IEL << endl;
      for(k=0;k<NVE2-1;k++)
      {
        INDH = (IND+1) % NVE2; 
        IVTH = stream->GetGlobalDOF(IEL)[VertexNumber(INDH, NVE2)];

        auto element = velo->get_fe(IEL);
        auto basefunct = element.GetBaseFunct_ID();
        auto desc = element.GetFEDesc();
        // int N_JointDOFs = desc->GetN_JointDOF();
        if(KVIND[IVTH]<1)
        {
          KVIND[IVTH]=1;

          l = VertexNumber(IND, NVE2);
          IVT = stream->GetGlobalDOF(IEL)[l];
  
          // do calculation
          cell->GetVertex(IND)->GetCoords(px1, py1);
          cell->GetVertex(INDH)->GetCoords(px2, py2);
          dn1 = py2-py1;
          dn2 = px1-px2;
  
          l1 = 0;
          while(stream->GetGlobalDOF(IEL)[VertexNumber(l1,NVE2)] !=IVT) l1++;
          // cout << "l1= " << l1 << " IVT: " << IVT;
          // cout << " IVTH: " << IVTH << endl;

          switch(basefunct)
          {
            case BF_N_T_P1_2D:
            case BF_N_Q_Q1_2D:
                m = velo->GetGlobalDOF(IEL)[IND];
                psi[IVTH] = psi[IVT] + u1[m]*dn1 + u2[m]*dn2;
              break;
            case BF_C_T_P0_2D:
            case BF_C_Q_Q0_2D:
                m=velo->GetGlobalDOF(IEL)[0];
                psi[IVTH] = psi[IVT] + u1[m]*dn1 + u2[m]*dn2;
              break;
            case BF_C_T_P1_2D:
            case BF_C_Q_Q1_2D:
                JointDOFs = desc->GetJointDOF(l1);
                m=velo->GetGlobalDOF(IEL)[JointDOFs[0]];
                v1 = u1[m];
                v2 = u2[m];
                m=velo->GetGlobalDOF(IEL)[JointDOFs[1]];
                v1 += u1[m];
                v2 += u2[m];
                v1 /= 2;
                v2 /= 2;
  
                psi[IVTH] = psi[IVT] + v1*dn1 + v2*dn2;
              break;
            case BF_C_T_P2_2D:
            case BF_C_Q_Q2_2D:
            case BF_C_T_B2_2D:
                JointDOFs = desc->GetJointDOF(l1);
                m=velo->GetGlobalDOF(IEL)[JointDOFs[0]];
                v1 = u1[m];
                v2 = u2[m];
                // cout << m << "  " << u1[m] << endl;
                m=velo->GetGlobalDOF(IEL)[JointDOFs[1]];
                v1 += 4*u1[m];
                v2 += 4*u2[m];
                // cout << m << "  " << u1[m] << endl;
                m=velo->GetGlobalDOF(IEL)[JointDOFs[2]];
                v1 += u1[m];
                v2 += u2[m];
                // cout << m << "  " << u1[m] << endl;
                v1 /= 6;
                v2 /= 6;
  
                psi[IVTH] = psi[IVT] + v1*dn1 + v2*dn2;
                // cout << "IVT: " << IVT << endl;
                // cout << "IVTH: " << IVTH << endl;
                // cout << "v1: " << v1 << "     " << "v2: " << v2 << endl;
                // cout << "delta psi: " << v1*dn1 + v2*dn2 << endl;
                // cout << endl;
              break;
            case BF_C_T_P3_2D:
            case BF_C_Q_Q3_2D:
            case BF_C_T_B3_2D:
                JointDOFs = desc->GetJointDOF(l1);
                m=velo->GetGlobalDOF(IEL)[JointDOFs[0]];
                v1 = u1[m];
                v2 = u2[m];
                m=velo->GetGlobalDOF(IEL)[JointDOFs[1]];
                v1 += 3*u1[m];
                v2 += 3*u2[m];
                m=velo->GetGlobalDOF(IEL)[JointDOFs[2]];
                v1 += 3*u1[m];
                v2 += 3*u2[m];
                m=velo->GetGlobalDOF(IEL)[JointDOFs[3]];
                v1 += u1[m];
                v2 += u2[m];
                v1 /= 8;
                v2 /= 8;
  
                psi[IVTH] = psi[IVT] + v1*dn1 + v2*dn2;
              break;
            case BF_C_T_P4_2D:
            case BF_C_Q_Q4_2D:
                JointDOFs = desc->GetJointDOF(l1);
                m=velo->GetGlobalDOF(IEL)[JointDOFs[0]];
                v1 = 7*u1[m];
                v2 = 7*u2[m];
                m=velo->GetGlobalDOF(IEL)[JointDOFs[1]];
                v1 += 32*u1[m];
                v2 += 32*u2[m];
                m=velo->GetGlobalDOF(IEL)[JointDOFs[2]];
                v1 += 12*u1[m];
                v2 += 12*u2[m];
                m=velo->GetGlobalDOF(IEL)[JointDOFs[3]];
                v1 += 32*u1[m];
                v2 += 32*u2[m];
                m=velo->GetGlobalDOF(IEL)[JointDOFs[4]];
                v1 += 7*u1[m];
                v2 += 7*u2[m];
                v1 /= 45;
                v2 /= 45;
  
                psi[IVTH] = psi[IVT] + v1*dn1 + v2*dn2;
              break;
            case BF_C_T_P5_2D:
            case BF_C_Q_Q5_2D:
                JointDOFs = desc->GetJointDOF(l1);
                m=velo->GetGlobalDOF(IEL)[JointDOFs[0]];
                v1 = 19*u1[m];
                v2 = 19*u2[m];
                m=velo->GetGlobalDOF(IEL)[JointDOFs[1]];
                v1 += 75*u1[m];
                v2 += 75*u2[m];
                m=velo->GetGlobalDOF(IEL)[JointDOFs[2]];
                v1 += 50*u1[m];
                v2 += 50*u2[m];
                m=velo->GetGlobalDOF(IEL)[JointDOFs[3]];
                v1 += 50*u1[m];
                v2 += 50*u2[m];
                m=velo->GetGlobalDOF(IEL)[JointDOFs[4]];
                v1 += 75*u1[m];
                v2 += 75*u2[m];
                m=velo->GetGlobalDOF(IEL)[JointDOFs[5]];
                v1 += 19*u1[m];
                v2 += 19*u2[m];
                v1 /= 144;
                v2 /= 144;
  
                psi[IVTH] = psi[IVT] + v1*dn1 + v2*dn2;
              break;
            case BF_C_T_P6_2D:
            case BF_C_Q_Q6_2D:
                JointDOFs = desc->GetJointDOF(l1);
                m=velo->GetGlobalDOF(IEL)[JointDOFs[0]];
                v1 = 41*u1[m];
                v2 = 41*u2[m];
                m=velo->GetGlobalDOF(IEL)[JointDOFs[1]];
                v1 += 216*u1[m];
                v2 += 216*u2[m];
                m=velo->GetGlobalDOF(IEL)[JointDOFs[2]];
                v1 += 27*u1[m];
                v2 += 27*u2[m];
                m=velo->GetGlobalDOF(IEL)[JointDOFs[3]];
                v1 += 272*u1[m];
                v2 += 272*u2[m];
                m=velo->GetGlobalDOF(IEL)[JointDOFs[4]];
                v1 += 27*u1[m];
                v2 += 27*u2[m];
                m=velo->GetGlobalDOF(IEL)[JointDOFs[5]];
                v1 += 216*u1[m];
                v2 += 216*u2[m];
                m=velo->GetGlobalDOF(IEL)[JointDOFs[6]];
                v1 += 41*u1[m];
                v2 += 41*u2[m];
                v1 /= 420;
                v2 /= 420;
  
                psi[IVTH] = psi[IVT] + v1*dn1 + v2*dn2;
              break;
            case BF_C_T_P7_2D:
            case BF_C_Q_Q7_2D:
                JointDOFs = desc->GetJointDOF(l1);
                m=velo->GetGlobalDOF(IEL)[JointDOFs[0]];
                v1 = 751*u1[m];
                v2 = 751*u2[m];
                m=velo->GetGlobalDOF(IEL)[JointDOFs[1]];
                v1 += 3577*u1[m];
                v2 += 3577*u2[m];
                m=velo->GetGlobalDOF(IEL)[JointDOFs[2]];
                v1 += 1323*u1[m];
                v2 += 1323*u2[m];
                m=velo->GetGlobalDOF(IEL)[JointDOFs[3]];
                v1 += 2989*u1[m];
                v2 += 2989*u2[m];
                m=velo->GetGlobalDOF(IEL)[JointDOFs[4]];
                v1 += 2989*u1[m];
                v2 += 2989*u2[m];
                m=velo->GetGlobalDOF(IEL)[JointDOFs[5]];
                v1 += 1323*u1[m];
                v2 += 1323*u2[m];
                m=velo->GetGlobalDOF(IEL)[JointDOFs[6]];
                v1 += 3577*u1[m];
                v2 += 3577*u2[m];
                m=velo->GetGlobalDOF(IEL)[JointDOFs[7]];
                v1 += 751*u1[m];
                v2 += 751*u2[m];
                v1 /= 8640;
                v2 /= 8640;
  
                psi[IVTH] = psi[IVT] + v1*dn1 + v2*dn2;
              break;
            default:
                cerr << "no stream function calculation possible for";
                cerr << " the base function set: " << basefunct << endl;
          }
        } // endif

        IND = INDH;
      } // endfor k

      // cout << "NVE2: " << NVE2 << endl;
      for(k=0;k<NVE2;k++)
      {
        auto neigh = cell->GetJoint(k)->GetNeighbour(cell);
   
        if(neigh)
        {
          if(neigh->GetClipBoard()>=0)
          { 
            // neigh was not handled
            IELH = neigh->GetClipBoard();
            IH = 0;
            NVE = neigh->GetN_Edges();
  
            for(l=0;l<NVE;l++)
            {
              JVE = stream->GetGlobalDOF(IELH)[VertexNumber(l, NVE)];
              IHV = KVIND[JVE];
              IH = IH+IHV;
              if(IHV>=1) INDH = l;
            } // endfor l
    
            if((IH<NVE) && (IH>0))
            {
              // add cell to list since not all values were calculated
              CellList[ListInput] = neigh;
              l = -(10+neigh->GetClipBoard());
              if(l<0)
                neigh->SetClipBoard(l);
              ListInput++;
            }
          } // !=-1
        } // if neigh
      } // endfor k;
      cell->SetClipBoard(-2);
    }
    ListPointer++;
  } // endwhile

/*
  for(IEH=0;IEH<NVT;IEH++)
    if(KVIND[IEH]==0)
      cerr << IEH << "   " << KVIND[IEH] << endl;
 */

  delete [] KVIND;

  delete [] CellList;
}

// determine L2 error, 2D
void L2Error(int N_Points, std::array<const double*, 2>,
                const double *AbsDetjk, const double *Weights, double,
                const double *const* Der, const double *const* Exact,
                const double *const*, double *LocError)
{
  LocError[0] = 0.0;

  for (int i = 0; i < N_Points; i++)
  {
    const double *deriv = Der[i];
    const double *exactval = Exact[i];

    double w = Weights[i] * AbsDetjk[i];

    double t = deriv[0] - exactval[0];
    LocError[0] += w * t * t;
  } // endfor i
}

// determine L2 and H1 error, 2D
void L2H1Errors(int N_Points, std::array<const double*, 2>,
                const double *AbsDetjk, const double *Weights, double,
                const double *const* Der, const double *const* Exact,
                const double *const*, double *LocError)
{
  LocError[0] = 0.0;
  LocError[1] = 0.0;

  for(int i=0;i<N_Points;i++)
  {
    const double *deriv = Der[i];
    const double *exactval = Exact[i];
    double w = Weights[i]*AbsDetjk[i];

    double t = deriv[0]-exactval[0];
    LocError[0] += w*t*t;

    t = deriv[1]-exactval[1];
    LocError[1] += w*t*t;
      
    t = deriv[2]-exactval[2];
    LocError[1] += w*t*t;
  } // endfor i
//   cout << "LocError[0]: " << LocError[0] << endl;
  // cout << "LocError[1]: " << LocError[1] << endl;
}

// determine L2-error, divergence error and H1 error, 2D
void L2DivH1Errors(int N_Points, std::array<const double*, 2>,
                   const double *AbsDetjk, const double *Weights, double,
                   const double *const* Der, const double *const* Exact,
                   const double *const*, double *LocError)
{
  int i;
  const double *deriv, *exactval;
  double w, t;

  LocError[0] = 0.0;
  LocError[1] = 0.0;
  LocError[2] = 0.0;
 // cout << endl;
  for(i=0;i<N_Points;i++)
  {
    deriv = Der[i];
    exactval = Exact[i];
    w = Weights[i]*AbsDetjk[i];
    
    t = deriv[0]-exactval[0]; // x-component
    LocError[0] += w*t*t;
    t = deriv[3]-exactval[4]; // y-component
    LocError[0] += w*t*t;

    // L2-error of divergence
    t  = deriv[1]-exactval[1]; // x-derivative of x-component
    t += deriv[5]-exactval[6]; // y-derivative of y-component
    LocError[1] += w*t*t;
    
    // H1 semi norm
    t = deriv[1]-exactval[1]; // x-derivative of x-component
    LocError[2] += w*t*t;
    t = deriv[2]-exactval[2]; // y-derivative of x-component
    LocError[2] += w*t*t;
    t = deriv[4]-exactval[5]; // x-derivative of y-component
    LocError[2] += w*t*t;
    t = deriv[5]-exactval[6]; // y-derivative of y-component
    LocError[2] += w*t*t;
  } // endfor i
  // cout << "LocError[0]: " << LocError[0] << endl;
  // cout << "LocError[1]: " << LocError[1] << endl;
}



void SDFEMErrors(int N_Points, std::array<const double*, 2>,
                 const double *AbsDetjk,  const double *Weights, double hK,
                 const double *const* Der, const double *const* Exact,
                 const double *const* coeffs, double *LocError)
{
  int i;
  const double *deriv, *exactval;
  double w;
  const double *coeff;
  double c0, c1, c2, c3, c5;
  double e0, e1, e2, e3;
  double loc0, loc1, loc2, loc3, loc4;

  //static double delta0 = TDatabase::ParamDB->DELTA0;
  //static double delta1 = TDatabase::ParamDB->DELTA1;
  //double time_step = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double delta;

  loc0 = 0.0;
  loc1 = 0.0;
  loc2 = 0.0;
  loc3 = 0.0;
  loc4 = 0.0;

  for(i=0;i<N_Points;i++)
  {
    coeff = coeffs[i];
    c0 = coeff[0];
    c1 = coeff[1];
    c2 = coeff[2];
    c3 = coeff[3];
    c5 = std::max(std::abs(c1),std::abs(c2));

    // if (X[i]>TDatabase::ParamDB->P6) continue;
    // if (Y[i]>TDatabase::ParamDB->P6) continue;

    // SUPG parameter
    delta = Compute_SDFEM_delta<2>(hK, c0, {{c1, c2}}, c3, c5);
    
    deriv = Der[i];
    exactval = Exact[i];
    w = Weights[i]*AbsDetjk[i];

    // error in solution
    e0 = deriv[0]-exactval[0];
    if (std::abs(e0)>loc4)
      loc4 = std::abs(e0);
    loc0 += w*e0*e0;
    if (std::abs(e0) > loc3)
	     loc3 = std::abs(e0);
    
    // error in derivative of solution
    e1 = deriv[1]-exactval[1];
    loc1 += w*e1*e1;
    e2 = deriv[2]-exactval[2];
    loc1 += w*e2*e2;
    
    // sd error
    // THIS IS ONLY CORRECT IF DIV(b) = 0
    e3 = c1*e1+c2*e2;
    loc2 += w*(c0*(e1*e1+e2*e2) + c3*e0*e0 + delta*e3*e3);
    
  } // endfor i
  LocError[0] = loc0;
  LocError[1] = loc1;
  LocError[2] = loc2;
  LocError[3] = loc3;

  //cout << "LocError[3]: " << LocError[3] << endl;
  // cout << "LocError[1]: " << LocError[1] << endl;
}

// determine L2, H1 and SDFEM error, in (0,P6)^2
void SDFEMErrorsSmooth(int N_Points, double *X, double *Y, double *AbsDetjk, 
                 const double *Weights, double hK, double **Der, double **Exact,
                 double **coeffs, double *LocError)
{
  int i, sd_type = TDatabase::ParamDB->SDFEM_TYPE;
  double *deriv, *exactval, w;
  double *coeff, c0, c1, c2, c5, alpha;
  double e0, e1, e2, e3;
  double loc0, loc1, loc2;

  static double delta0 = TDatabase::ParamDB->DELTA0;
  static double delta1 = TDatabase::ParamDB->DELTA1;
  double delta;

  loc0 = 0.0;
  loc1 = 0.0;
  loc2 = 0.0;

  for(i=0;i<N_Points;i++)
  {
    if (X[i]>TDatabase::ParamDB->P6) continue;
    if (Y[i]>TDatabase::ParamDB->P6) continue;

    coeff = coeffs[i];
    c0 = coeff[0];
    c1 = coeff[1];
    c2 = coeff[2];
    c5 = std::max(std::abs(c1),std::abs(c2));

    // if (X[i]>TDatabase::ParamDB->P6) continue;
    // if (Y[i]>TDatabase::ParamDB->P6) continue;

    if (sd_type==0)
    {
       if(c0 < hK*c5)
          delta = delta0 * hK/c5;
       else
          delta = delta1 *hK*hK/c0 ;
    }
    else
    {
       if (c5>0)
       {
          alpha = c5*hK/(2*c0);
          delta = hK*(1/std::tanh(alpha) - 1/alpha)/(2*c5);
       }
       else
          delta = 0;
    }
    deriv = Der[i];
    exactval = Exact[i];
    w = Weights[i]*AbsDetjk[i];

    e0 = deriv[0]-exactval[0];
    loc0 += w*e0*e0;

    e1 = deriv[1]-exactval[1];
    loc1 += w*e1*e1;
    e2 = deriv[2]-exactval[2];
    loc1 += w*e2*e2;

    e3 = c1*e1+c2*e2;
    loc2 += w*(c0*(e1*e1+e2*e2) + e0*e0 + delta*e3*e3);

  } // endfor i

  LocError[0] = loc0;
  LocError[1] = loc1;
  LocError[2] = loc2;

  // cout << "LocError[0]: " << LocError[0] << endl;
  // cout << "LocError[1]: " << LocError[1] << endl;
}

// determine L2, H1 and SDFEM error for smooth region in the
// example JohnMaubachTobiska1997 (x-0.5)^2+(y-0.5)^2 > r^2 
void SDFEMErrorsSmooth_JohnMaubachTobiska1997
(int N_Points, double *X, double *Y, double *AbsDetjk, 
 const double *Weights, double hK, double **Der, double **Exact,
 double **coeffs, double *LocError)
{
    int i, sd_type = TDatabase::ParamDB->SDFEM_TYPE;
    double *deriv, *exactval, w;
    double *coeff, c0, c1, c2, c3, c5, alpha;
    double e0, e1, e2, e3;
    double loc0, loc1, loc2;
    
    static double delta0 = TDatabase::ParamDB->DELTA0;
    static double delta1 = TDatabase::ParamDB->DELTA1;
    double time_step = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
    double delta, norm_b, x, y, r2=0.3*0.3;
    
    loc0 = 0.0;
    loc1 = 0.0;
    loc2 = 0.0;
    
    for(i=0;i<N_Points;i++)
    {
	x = X[i]-0.5;
	y = Y[i]-0.5;   
	if (x*x+y*y<=r2)
	  continue;
	coeff = coeffs[i];
	c0 = coeff[0];
	c1 = coeff[1];
	c2 = coeff[2];
	c3 = coeff[3];
	c5 = std::max(std::abs(c1),std::abs(c2));
		
	// SUPG parameter
	switch(sd_type)
	{
	    case 0:
		if(c0 < hK*c5)
		    delta = delta0 * hK/c5;
		else
		    delta = delta1 *hK*hK/c0;
		break;
	    case 1:
		norm_b = std::sqrt(c1*c1+c2*c2);
		//norm_b = linfb;
		if (norm_b > 0)
	    {
		alpha = norm_b*hK/(2*c0);
		delta = hK*(1/std::tanh(alpha) - 1/alpha)/(2*norm_b);
	    }
		else
		    delta = 0;
		break;
	    case 9:
		if (c0 <= hK)
		{
		    delta = delta0*time_step;	
                    //delta = delta0*time_step*time_step;
		}
		else
		{
		    delta = delta0*time_step;
                    //delta = delta0*time_step*time_step;
		}
		break;
	    case 10: 
		// second estimate in paper with Julia Novo
		// get the unscaled parameters
		if(c0 <= hK)
		    delta = delta0 * hK;
		else
		    delta = delta0 *hK*hK/c0 ;
		break;
	    case 11:
		// for estimate in paper with Julia Novo
		norm_b = std::sqrt(c1*c1+c2*c2);
		delta = delta0 * hK * std::sqrt(time_step)/norm_b;
		break;
	    default:
	      //Output::print("CHECK IF CORRECT DELTA IN ERROR COMPUTATION !!!");
		if(c0 < hK*c5)
		    delta = delta0 * hK/c5;
		else
		    delta = delta1 *hK*hK/c0;
		break;
	}
	deriv = Der[i];
	exactval = Exact[i];
	w = Weights[i]*AbsDetjk[i];
	
	// error in solution
	e0 = deriv[0]-exactval[0];
	loc0 += w*e0*e0;
	
	// error in derivative of solution
	e1 = deriv[1]-exactval[1];
	loc1 += w*e1*e1;
	e2 = deriv[2]-exactval[2];
	loc1 += w*e2*e2;
	
	// sd error
	// THIS IS ONLY CORRECT IF DIV(b) = 0
	e3 = c1*e1+c2*e2;
	loc2 += w*(c0*(e1*e1+e2*e2) + c3*e0*e0 + delta*e3*e3);
    } // endfor i
    LocError[0] = loc0;
    LocError[1] = loc1;
    LocError[2] = loc2;
}

// determine errors to interpolant
// paper with Julia Novo
void SDFEMErrorsInterpolant(int N_Points, double *, double *, double *AbsDetjk, 
                 const double *Weights, double hK, double **Der, double **Exact,
                 double **coeffs, double *LocError)
{
  int i;
  double *deriv, *exactval, w;
  double *coeff, c0, c1, c2, c5;
  double e0, e1, e2, e3, e4;
  double loc0, loc1, loc2, loc3;
  double delta, c_inv;

  loc0 = 0.0;
  loc1 = 0.0;
  loc2 = 0.0;
  loc3 = 0.0;

  switch (TDatabase::ParamDB->INTERNAL_POLYNOMIAL_DEGREE)
  {
    case 1:  c_inv = 1.0;
    break;
    case 2:
        // triangle
        c_inv = std::sqrt(48.0);
        // quad
        c_inv = std::sqrt(24.0);
        break;
    case 3:
        // triangle
        c_inv = std::sqrt((435+std::sqrt(26025.0))/4.0);
        // quad
        c_inv = std::sqrt((244+std::sqrt(9136.0))/3.0);
        break;
    default:
      ErrThrow("c_inv not defined ");
      break;
  }

  for(i=0;i<N_Points;i++)
  {
    // compute SUPG parameter
    coeff = coeffs[i];
    c0 = coeff[0];
    c1 = coeff[1];
    c2 = coeff[2];
    // double c3 = coeff[3];
    c5 = std::max(std::abs(c1),std::abs(c2));
    delta = Compute_SDFEM_delta<2>(hK, coeff[0], {{coeff[1], coeff[2]}}, coeff[3], c5);
    if (TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY == 120814)
      delta = TDatabase::ParamDB->INTERNAL_P1_Array[TDatabase::ParamDB->INTERNAL_LEVEL];

    // NOTE: order in derivatives in Derivatives_SD
    deriv = Der[i];
    exactval = Exact[i];
    w = Weights[i]*AbsDetjk[i];

    // error in solution
    e0 = deriv[2]-exactval[0];
    loc0 += w*e0*e0/delta;
    
    // error in streamline derivative of solution
    e1 = c1*(deriv[0]-exactval[1]);
    loc1 += w*delta*e1*e1;
    e2 = c2*(deriv[1]-exactval[2]);
    loc1 += w*delta*e2*e2;

    // first additional term, with gradient
    e3 = c0 * c_inv * (deriv[0]-exactval[1])/ hK;
    loc2 += w * 16.0 * delta * e3 * e3;
    e3 = c0 * c_inv * (deriv[1]-exactval[2])/ hK;
    loc2 += w * 16.0 * delta * e3 * e3;

    // second additional term, with Laplacian
     e4 = c0 * (deriv[3] + deriv[4] - exactval[3]);
    loc3 += w * 8.0 * delta * e4 * e4;
  } // endfor i
  LocError[0] = loc0;
  LocError[1] = loc1;
  LocError[2] = loc2;
  LocError[3] = loc3;
 //cout << "LocError[3]: " << LocError[3] << endl;
  // cout << "LocError[1]: " << LocError[1] << endl;
}

// determine L1 error, 2D
void L1Error(int N_Points, std::array<const double*, 2> xy,
             const double *AbsDetjk, const double *Weights, double hK,
             const double *const* Der, const double *const* Exact,
             const double *const*, double *LocError)
{
  int i, index;
  const double *deriv, *exactval;
  double w, t, area, v[3], va[3], xa[3], ya[3];
  double a[2], b[2], fac, areaa;
  LocError[0] = 0.0;

  // conforming linear finite elemens
  if (hK==-4711)
  {
      // area of the mesh cell
      area = AbsDetjk[0]/2;
      // compute values of the fe function in the vertices of the mesh cell
      // giveen are the values in the midpoints of the edges
      // Der[0][0] - edge between P0 and P1
      // Der[1][0] - edge between P1 and P2
      // Der[2][0] - edge between P2 and P0
      v[0] = Der[0][0] - Der[1][0] +   Der[2][0];
      v[1] = Der[1][0] - Der[2][0] +   Der[0][0];
      v[2] = Der[2][0] - Der[0][0] +   Der[1][0];
      //for (i=0;i<3;i++)
      //          Output::print(X[i], " ", Y[i], " ", v[i]);
      // four different cases
      // case one and two: no change of sign
      if ((v[0] >= 0)&& (v[1] >= 0) && (v[2] >= 0))
      {
          LocError[0] = area*(v[0]+v[1]+v[2])/3;
          return;
      }
      if ((v[0] <= 0)&& (v[1] <= 0) && (v[2] <= 0))
      {
          LocError[0] = -area*(v[0]+v[1]+v[2])/3;
          return;
      }
      // change of sign, one value has different sign than the other two ones
      // this value gets index 2
      t = v[0]*v[1]*v[2];
      // two negative, one positive values
      if (t>0)
      {
          for (i=0;i<3;i++)
          {
              if (v[i] > 0)
              {
                  index = i;
                  break;
              }
          }
      }
      else
      // two positive, one negative values
      {
          for (i=0;i<3;i++)
          {
              if (v[i] < 0)
              {
                  index = i;
                  break;
              }
          }
      }
      for (i=0;i<3;i++)
      {
          va[i] = v[(index+1+i)%3]; 
          xa[i] = xy[0][(index+1+i)%3]; 
          ya[i] = xy[1][(index+1+i)%3]; 
      }     
      //for (i=0;i<3;i++)
      //          Output::print(xa[i], " ", ya[i], " ", va[i]);
      // compute area of triangle with different sign of the fe function
      a[0] = xa[0] - xa[2];
      a[1] = ya[0] - ya[2];
      fac = -va[2]/(va[0]-va[2]);
      a[0] *= fac;
      a[1] *= fac;

      b[0] = xa[1] - xa[2];
      b[1] = ya[1] - ya[2];
      fac = -va[2]/(va[1]-va[2]);
      b[0] *= fac;
      b[1] *= fac;

      areaa = std::abs(a[0]*b[1]-b[0]*a[1])/2;

      //Output::print(a[0], " ", a[1], " ", b[0], " ", b[1], " ", areaa);
      LocError[0] = area*(v[0]+v[1]+v[2])/3 - 2*areaa*va[2]/3;
      LocError[0] = std::abs(LocError[0]);
  }
  else
  // use quadrature rule for other types of finite elements
  {
      for(i=0;i<N_Points;i++)
      {
          deriv = Der[i];
          exactval = Exact[i];
          w = Weights[i]*AbsDetjk[i];
          
          t = deriv[0]-exactval[0];
          LocError[0] += w*std::abs(t);
      } // endfor i
  }
  //Output::print("LocError(L1): ", LocError[0]);
}

// determine deformation tensor error
void DeformationTensorError(int N_Points, double *, double *,
                double *AbsDetjk, const double *Weights, double,
                double **Der, double **Exact,
                double **, double *LocError)
{
  int i;
  double *deriv_x, *exactval, *deriv_y, *exactval1, w, t;

  LocError[0] = 0.0;

  for(i=0;i<N_Points;i++)
  {
    // first component
    deriv_x = Der[i];
    exactval = Exact[i];
    // second component
    deriv_y = Der[i+N_Points];
    exactval1 = Exact[i+N_Points];
    // weight
    w = Weights[i]*AbsDetjk[i];

    // left upper term
    t = deriv_x[1]-exactval[1];
    LocError[0] += w*t*t;
    // right upper and left lower term
    t = ((deriv_x[2]-exactval[2])+(deriv_y[1]-exactval1[1]));
    LocError[0] += w*t*t/2.0;
    // right lower term
    t = deriv_y[2]-exactval1[2];
    LocError[0] += w*t*t;

  } // endfor i

  // cout << "LocError[0]: " << LocError[0] << endl;
  // cout << "LocError[1]: " << LocError[1] << endl;
}

// determine L2 and H1 error, 2D
void H1Norm(int N_Points, double *, double *, double *AbsDetjk, 
            const double *Weights, double,
            double **Der, double **Exact,
            double **, double *LocError)
{
  int i;
  double *deriv, *exactval, w, t;

  LocError[0] = 0.0;

  for(i=0;i<N_Points;i++)
  {
    deriv = Der[i];
    exactval = Exact[i];
    w = Weights[i]*AbsDetjk[i];

    t = deriv[0]-exactval[0];
    LocError[0] += w*t*t;

    t = deriv[1]-exactval[1];
    LocError[0] += w*t*t;
    t = deriv[2]-exactval[2];
    LocError[0] += w*t*t;
  } // endfor i

  // cout << "LocError[0]: " << LocError[0] << endl;
  // cout << "LocError[1]: " << LocError[1] << endl;
}
// compute the error in the divergence
void DivergenceError(int N_Points, double *, double *,
		     double *AbsDetjk, const double *Weights, double,
		     double **Der, double **,
		     double **, double *LocError)
{
  int i;
  double *deriv_x, *deriv_y, w, t;

  LocError[0] = 0.0;
  LocError[1] = 0.0;
 
  for(i=0;i<N_Points;i++)
  {
    // first component
    deriv_x = Der[i];
    // second component
    deriv_y = Der[i+N_Points];
    // weight
    w = Weights[i]*AbsDetjk[i];

     // u_x + v_y
    t = std::abs(deriv_x[1] + deriv_y[2]);
    LocError[0] += w*t;
    LocError[1] += w*t*t;    
  } // endfor i

  //cout << "LocError[0]: " << LocError[0] << endl;
}

// compute the error in the grad-div term for Oseen
void DivergenceErrorGradDivOseen(int N_Points, double *, double *,
         double *AbsDetjk, const double *Weights, double hK, 
         double **Der, double **,
         double **coeffs, double *LocError)
{
  int i;
  double *deriv_x, *deriv_y, w, t, nu, b1, b2, mu, *coeff;

  LocError[0] = 0.0;
 
  for(i=0;i<N_Points;i++)
  {
    coeff = coeffs[i];
    nu = coeff[0];
    b1 = coeff[3];
    b2 = coeff[4];

   // get stabilization parameters
    mu = graddiv_parameterOseen(hK, nu, b1, b2);
   // first component
    deriv_x = Der[i];
    // second component
    deriv_y = Der[i+N_Points];
    // weight
    w = Weights[i]*AbsDetjk[i];
    
     // u_x + v_y
    t = std::abs(deriv_x[1] + deriv_y[2]);
    LocError[0] += w*mu*t*t;    
  } // endfor i

  //cout << "LocError[0]: " << LocError[0] << endl;
}


// mesh cell parameters for shock capturing scheme DC_CD
void Parameters_DC_CD(int N_Points, double *, double *, double *AbsDetjk, 
           const double *Weights, double,
           double **Der, double **Exact,
           double **coeffs, double *LocError)
{
  int i;
  double *deriv, *exactval, w, t, *coeff;
  double eps, b1, b2, c, f;

  LocError[0] = 0.0;
  LocError[1] = 0.0;

  for(i=0;i<N_Points;i++)
  {
    coeff = coeffs[i];
    eps = coeff[0];
    b1 = coeff[1];
    b2 = coeff[2];
    c = coeff[3];
    f = coeff[4];

    deriv = Der[i];
    exactval = Exact[i];
    w = Weights[i]*AbsDetjk[i];

    t = deriv[2]-exactval[2];
   // LocError[0] += w*t*t;

    t = deriv[0]-exactval[0];
    LocError[0] += w*t*t;
    t = deriv[1]-exactval[1];
    LocError[0] += w*t*t;
   
    //   t= -eps*deriv[3]+ b1*deriv[1]+b2*deriv[2] + c*deriv[0] - f;
    t= -eps*(deriv[3] + deriv[4]) + b1*deriv[0]+b2*deriv[1] + c*deriv[2] - f;
    LocError[1] += w*t*t;
    

  } // endfor i

  // cout << "LocError[0]: " << LocError[0] << endl;
  // cout << "LocError[1]: " << LocError[1] << endl;
}

// mesh cell values for gradient and residual 
void Parameters_Gradient_Residual(int N_Points, double *, double *, double *AbsDetjk,
           const double *Weights, double,
           double **Der, double **Exact,
           double **coeffs, double *LocError)
{
  int i;
  double *deriv, *exactval, w, t, *coeff;
  double eps, b1, b2, c, f;

  LocError[0] = 0.0;
  LocError[1] = 0.0;

  for(i=0;i<N_Points;i++)
  {
    coeff = coeffs[i];
    eps = coeff[0];
    b1 = coeff[1];
    b2 = coeff[2];
    c = coeff[3];
    f = coeff[4];

    deriv = Der[i];
    exactval = Exact[i];
    w = Weights[i]*AbsDetjk[i];

    // L^2 norm of gradient gradient
    t = deriv[0]-exactval[0];
    LocError[0] += w*t*t;
    t = deriv[1]-exactval[1];
    LocError[0] += w*t*t;

    // L^2 norm of residual
    t= -eps*(deriv[3] + deriv[4]) + b1*deriv[0]+b2*deriv[1] + c*deriv[2] - f;
    LocError[1] += w*t*t;
  } // endfor i
}

#endif // 2D

#ifdef __3D__

// determine L2 error, 3D
void L2Error(int N_Points, std::array<const double*, 3>,
                const double *AbsDetjk, 
                const double *Weights, double,
                const double *const* Der, const double *const* Exact,
                const double *const*, double *LocError)
{
  int i;
  const double *deriv, *exactval;
  double w, t;

  LocError[0] = 0.0;

  for (i = 0; i < N_Points; i++)
  {
    deriv = Der[i];
    exactval = Exact[i];
    w = Weights[i] * AbsDetjk[i];

    t = deriv[0] - exactval[0];
    LocError[0] += w * t * t;
  } // endfor i
}

// determine L2 and H1 error, 3D
void L2H1Errors(int N_Points, std::array<const double*, 3>,
                const double *AbsDetjk, 
                const double *Weights, double,
                const double *const* Der, const double *const* Exact,
                const double *const*, double *LocError)
{
  int i;
  const double *deriv, *exactval;
  double w, t;

  LocError[0] = 0.0;
  LocError[1] = 0.0;

  for(i=0;i<N_Points;i++)
  {
    deriv = Der[i];
    exactval = Exact[i];
    w = Weights[i]*AbsDetjk[i];

    t = deriv[0]-exactval[0];
    LocError[0] += w*t*t;

    t = deriv[1]-exactval[1];
    LocError[1] += w*t*t;
    t = deriv[2]-exactval[2];
    LocError[1] += w*t*t;
    t = deriv[3]-exactval[3];
    LocError[1] += w*t*t;
  } // endfor i

  // cout << "LocError[0]: " << LocError[0] << endl;
  // cout << "LocError[1]: " << LocError[1] << endl;
}

void L2H1ErrorsSmooth(int N_Points, double *X, double *Y, double *Z, 
                double *AbsDetjk, 
                const double *Weights, double,
                double **Der, double **Exact,
                double **, double *LocError)
{
  int i;
  double *deriv, *exactval, w, t;

  LocError[0] = 0.0;
  LocError[1] = 0.0;

  for(i=0;i<N_Points;i++)
  {
    if (X[i]>TDatabase::ParamDB->P6) continue;
    if (Y[i]>TDatabase::ParamDB->P6) continue;
    if (Z[i]>TDatabase::ParamDB->P6) continue;
    
    deriv = Der[i];
    exactval = Exact[i];
    w = Weights[i]*AbsDetjk[i];

    t = deriv[0]-exactval[0];
    LocError[0] += w*t*t;

    t = deriv[1]-exactval[1];
    LocError[1] += w*t*t;
    t = deriv[2]-exactval[2];
    LocError[1] += w*t*t;
    t = deriv[3]-exactval[3];
    LocError[1] += w*t*t;
  } // endfor i

  // cout << "LocError[0]: " << LocError[0] << endl;
  // cout << "LocError[1]: " << LocError[1] << endl;
}

void L2DivH1Errors(int N_Points, std::array<const double*, 3>,
                   const double *AbsDetjk, const double *Weights, double,
                   const double *const* Der, const double *const* Exact,
                   const double *const*, double *LocError)
{
  LocError[0] = 0.0;
  LocError[1] = 0.0;
  LocError[2] = 0.0;
  
  for(int i = 0; i < N_Points; i++)
  {
    double w = Weights[i]*AbsDetjk[i];
    
    // L2-error
    double t = Der[i][0] - Exact[i][0]; // x-component
    LocError[0] += w*t*t;
    t = Der[i][4] - Exact[i][5];        // y-component
    LocError[0] += w*t*t;
    t = Der[i][8] - Exact[i][10];       // z-component
    LocError[0] += w*t*t;
    
    // L2-error of divergence
    t  = Der[i][1] - Exact[i][1];  // x-derivative of x-component
    t += Der[i][6] - Exact[i][7];  // y-derivative of y-component
    t += Der[i][11]- Exact[i][13]; // z-derivative of z-component
    LocError[1] += w*t*t;
    
    // H1 semi norm
    t = Der[i][1] - Exact[i][1];  // x-derivative of x-component
    LocError[2] += w*t*t;
    t = Der[i][2] - Exact[i][2];  // y-derivative of x-component
    LocError[2] += w*t*t;
    t = Der[i][3] - Exact[i][3];  // z-derivative of x-component
    LocError[2] += w*t*t;
    t = Der[i][5] - Exact[i][6];  // x-derivative of y-component
    LocError[2] += w*t*t;
    t = Der[i][6] - Exact[i][7];  // y-derivative of y-component
    LocError[2] += w*t*t;
    t = Der[i][7] - Exact[i][8];  // z-derivative of y-component
    LocError[2] += w*t*t;
    t = Der[i][9] - Exact[i][11]; // x-derivative of z-component
    LocError[2] += w*t*t;
    t = Der[i][10]- Exact[i][12]; // y-derivative of z-component
    LocError[2] += w*t*t;
    t = Der[i][11]- Exact[i][13]; // z-derivative of z-component
    LocError[2] += w*t*t;
  }
}

// determine L1 error, 3D
void L1Error(int, std::array<const double*, 3>, const double *, const double *,
             double, const double *const*, const double *const*,
             const double *const*, double *LocError)
{
    Output::print("computation of L1-error not implemented !!!");
    LocError[0] = 0;
}

// determine deformation tensor error
void DeformationTensorError(int N_Points, double *, double *, double *,
                double *AbsDetjk, const double *Weights, double, 
                double **Der, double **Exact,
                double **, double *LocError)
{
  int i;
  double *deriv_x, *exactval, *deriv_y, *exactval1, *deriv_z, *exactval2, w, t;

  LocError[0] = 0.0;

  for(i=0;i<N_Points;i++)
  {
    // first component
    deriv_x = Der[i];
    exactval = Exact[i];
    // second component
    deriv_y = Der[i+N_Points];
    exactval1 = Exact[i+N_Points];
    // third component
    deriv_z = Der[i+2*N_Points];
    exactval2 = Exact[i+2*N_Points];
    // weight
    w = Weights[i]*AbsDetjk[i];

    // left upper term
    t = deriv_x[1]-exactval[1];
    LocError[0] += w*t*t;
    // 12 and 21 term 
    t = ((deriv_x[2]-exactval[2])+(deriv_y[1]-exactval1[1]));
    LocError[0] += w*t*t/2.0;
    // 22 term
    t = deriv_y[2]-exactval1[2];
    LocError[0] += w*t*t;
    // 13 and 31 term
    t = ((deriv_x[3]-exactval[3])+(deriv_z[1]-exactval2[1]));
    LocError[0] += w*t*t/2.0;
    // 23 and 32 term
    t = ((deriv_y[3]-exactval1[3])+(deriv_z[2]-exactval2[2]));
    LocError[0] += w*t*t/2.0;
    // 33 term
    t = deriv_z[3]-exactval2[3];
    LocError[0] += w*t*t;       
  } // endfor i

  // cout << "LocError[0]: " << LocError[0] << endl;
  // cout << "LocError[1]: " << LocError[1] << endl;
}

// compute the error in the divergence
void DivergenceError(int N_Points, double *, double *, double *,
		     double *AbsDetjk, const double *Weights, double,
		     double **Der, double **,
		     double **, double *LocError)
{
  int i;
  double *deriv_x, *deriv_y, *deriv_z, w, t;

  LocError[0] = 0.0;
  LocError[1] = 0.0;
 
  for(i=0;i<N_Points;i++)
  {
    // first component
    deriv_x = Der[i];
    // second component
    deriv_y = Der[i+N_Points];
    // third component
    deriv_z = Der[i+2*N_Points];
    // weight
    w = Weights[i]*AbsDetjk[i];

     // u_x + v_y + w_z
    t = std::abs(deriv_x[1] + deriv_y[2] + deriv_z[3]);
    LocError[0] += w*t;
    LocError[1] += w*t*t;    
  } // endfor i

  //cout << "LocError[0]: " << LocError[0] << endl;
}
// mesh cell parameters for shock capturing scheme DC_CD
void Parameters_DC_CD(int N_Points, double *, double *, double *,
                      double *AbsDetjk, 
                      const double *Weights, double,
                      double **Der, double **Exact,
                      double **coeffs, double *LocError)
{
  int i;
  double *deriv, *exactval, w, t, *coeff;
  double b1, b2, b3, c, f;  //eps;

  LocError[0] = 0.0;
  LocError[1] = 0.0;

  for(i=0;i<N_Points;i++)
  {
    coeff = coeffs[i];
//    eps = coeff[0];
    b1 = coeff[1];
    b2 = coeff[2];
    b3 = coeff[3];
    c = coeff[4];
    f = coeff[5];

    deriv = Der[i];
    exactval = Exact[i];
    w = Weights[i]*AbsDetjk[i];

    t = deriv[3]-exactval[3];
    LocError[0] += w*t*t;

    t = deriv[0]-exactval[0];
    LocError[0] += w*t*t;
    t = deriv[1]-exactval[1];
    LocError[0] += w*t*t;
    t = deriv[2]-exactval[2];
    LocError[0] += w*t*t;
    
    //t= -eps*(deriv[4] + deriv[5] + deriv[6]) + b1*deriv[1]+ b2*deriv[2] + b3*deriv[3] + c*deriv[0] - f;
    t= b1*deriv[0]+ b2*deriv[1] + b3*deriv[2] + c*deriv[3] - f;
    LocError[1] += w*t*t;
    

  } // endfor i

  // cout << "LocError[0]: " << LocError[0] << endl;
  // cout << "LocError[1]: " << LocError[1] << endl;
}

/*******************************************************************************/
//
// compute the Q criterion for a flow field
//
/*******************************************************************************/
void Q_criterion(TCollection *Coll,
TFEFunction3D *velocity1, TFEFunction3D *velocity2,
TFEFunction3D *velocity3, double *Qcrit)
{
  int i, j, N_V, N_Cells;
  double x[8],y[8],z[8],values[4], eps = 1e-6;
  double grad_velo_xx, grad_velo_xy, grad_velo_xz,  grad_velo_yx, grad_velo_yy,  grad_velo_yz;
  double grad_velo_zx, grad_velo_zy, grad_velo_zz; 
  double x_sp, y_sp, z_sp, val;
 
  TBaseCell *cell;

   // number of cells
  N_Cells = Coll->GetN_Cells();
  //loop over the mesh cells of the global grid
  for(i=0;i<N_Cells;i++)
  {
    // get cell
    cell = Coll->GetCell(i);
    // number of vertices per cell
    N_V = cell->GetN_Vertices();
    // compute barycenter
    x_sp=0.;
    y_sp=0.;
    z_sp=0.;
    for (j=0;j<N_V;j++)
    {
      // read coordinates of the mesh cell
      cell->GetVertex(j)->GetCoords(x[j], y[j], z[j]);
      x_sp += x[j];
      y_sp += y[j];
      z_sp += z[j];
     }
     x_sp /= N_V;
     y_sp /= N_V;
     z_sp /= N_V;     
      // compute gradient in barycenter
      velocity1->FindGradientLocal(i,x_sp,y_sp,z_sp,values);
      grad_velo_xx = values[1];
      grad_velo_xy = values[2];
      grad_velo_xz = values[3];
      velocity2->FindGradientLocal(i,x_sp,y_sp,z_sp,values);
      grad_velo_yx = values[1];
      grad_velo_yy = values[2];
      grad_velo_yz = values[3];
      velocity3->FindGradientLocal(i,x_sp,y_sp,z_sp,values);
      grad_velo_zx = values[1];
      grad_velo_zy = values[2];
      grad_velo_zz = values[3];
      
     // Q criterion : 0.5 * (Q:Q-S:S)
     val=-0.5*((grad_velo_xx*grad_velo_xx)
                 +(grad_velo_yy*grad_velo_yy)
                 +(grad_velo_zz*grad_velo_zz)
                 +(grad_velo_xy*grad_velo_yx)
                 +(grad_velo_xz*grad_velo_zx)
                 +(grad_velo_zy*grad_velo_yz)
                 +(grad_velo_yx*grad_velo_xy)
                 +(grad_velo_zx*grad_velo_xz)
                 +(grad_velo_yz*grad_velo_zy));

     Output::print(val);
     if (std::abs(val) < eps)
     {
  Qcrit[i]=0;
     }
     else
     {
       if (val > 0)
  Qcrit[i] = 1; 
  else 
    Qcrit[i]=-1;
     }
  }
}

struct WallShearStressFaceInfo
{
  public:
    // set of boundary DOF
    std::set<int> boundary_dof;

    // boundary dof to adjacent boundary faces
    std::unordered_map<int, std::vector<std::pair<int, double>>> faces;

    // boundary face to adjacent cell
    std::vector<int> cells;

    // boundary face to normal
    std::vector<double> normals;

    // boundary face to centroid
    std::vector<double> points;
};

struct WallShearStressInfo
{
  public:
    // set of boundary DOF
    std::set<int> boundary_dof;

    // boundary DOF to (accumulated) normal on the lowest adjacent boundary ref
    std::unordered_map<int, std::vector<double>> normals;

    // boundary DOF to vertex position
    std::unordered_map<int, std::vector<double>> points;

    // boundary DOF to list of adjacent cells
    std::unordered_map<int, std::vector<std::pair<int, double>>> cells;
};

std::unordered_map<const TFESpace*, WallShearStressFaceInfo> wss_face_info_cache;
std::unordered_map<const TFESpace*, WallShearStressInfo> wss_info_cache;

double GetTetraSolidAngle(const TBaseCell* cell, int i)
{
  // solid angle of a tetrahedron at vertex i:
  // sum of the dihedral angles at its edges minus pi
  // very inefficient but legible implementation...

  double Omega = -M_PI;

  double x_i, y_i, z_i;
  cell->GetVertex(i)->GetCoords(x_i, y_i, z_i);

  for (int j = 0; j < 4; j++)
  {
    if (i == j)
    {
      continue;
    }

    // find the other two vertices k and m

    int k = 0;
    while (i == k || j == k)
    {
      ++k;
    }

    int m = k + 1;
    while (i == m || j == m)
    {
      ++m;
    }

    double x_j, y_j, z_j;
    cell->GetVertex(j)->GetCoords(x_j, y_j, z_j);
    x_j -= x_i;
    y_j -= y_i;
    z_j -= z_i;

    double x_k, y_k, z_k;
    cell->GetVertex(k)->GetCoords(x_k, y_k, z_k);
    x_k -= x_i;
    y_k -= y_i;
    z_k -= z_i;

    double x_m, y_m, z_m;
    cell->GetVertex(m)->GetCoords(x_m, y_m, z_m);
    x_m -= x_i;
    y_m -= y_i;
    z_m -= z_i;

    // compute the angle between planes ijk and ijm by computing the dot
    // product of their normals. notice that normal orientation is
    // irrelevant here, so we needn't worry about the signs of the cross
    // products below

    // (j-i) x (k-i)
    double x_n1, y_n1, z_n1;
    x_n1 = y_j * z_k - z_j * y_k;
    y_n1 = z_j * x_k - x_j * z_k;
    z_n1 = x_j * y_k - y_j * x_k;

    // (j-i) x (m-i)
    double x_n2, y_n2, z_n2;
    x_n2 = y_j * z_m - z_j * y_m;
    y_n2 = z_j * x_m - x_j * z_m;
    z_n2 = x_j * y_m - y_j * x_m;

    double mag2_n1 = x_n1 * x_n1 + y_n1 * y_n1 + z_n1 * z_n1;
    double mag2_n2 = x_n2 * x_n2 + y_n2 * y_n2 + z_n2 * z_n2;

    Omega += std::acos((x_n1 * x_n2 + y_n1 * y_n2 + z_n1 * z_n2)
      / std::sqrt(mag2_n1 * mag2_n2));
  }

  return Omega;
}

const WallShearStressFaceInfo& GetWallShearStressFaceInfo(
  const TFEVectFunct3D &u,
  const TFEVectFunct3D &tau_w)
{
  // Caches or computes all the helper data needed for the wall shear stress
  // computation - adjacent cells, normals, etc.

  const auto velocity_space = u.GetFESpace();
  const auto coll = velocity_space->GetCollection();

  const auto wall_shear_stress_space = tau_w.GetFESpace();

  const auto* key = velocity_space.get();

  if (wss_face_info_cache.count(key) > 0)
  {
    return wss_face_info_cache[key];
  }

  auto& wss_face_info = wss_face_info_cache[key] = WallShearStressFaceInfo();

  auto& boundary_dof = wss_face_info.boundary_dof;
  auto& faces = wss_face_info.faces;
  auto& normals = wss_face_info.normals;
  auto& points = wss_face_info.points;
  auto& cells = wss_face_info.cells;

  int num_relevant_faces = 0;

  // notice that our loop includes halo cells since we need them for
  // computing the relevant averages
  int n_cells = coll->GetN_Cells();

  int min_ref = 0;
  bool any_ref = false;

  for (int i = 0; i < n_cells; i++)
  {
    const TBaseCell* cell = coll->GetCell(i);

    int n_joints = cell->GetN_Joints();

    // check all face joints
    for (int j = 0; j < n_joints; j++)
    {
      const auto joint = cell->GetJoint(j);

      if (joint->GetType() == BoundaryFace || joint->GetType() == IsoBoundFace)
      {
        const auto bd_face = reinterpret_cast<const TBoundFace*>(joint);
        const auto bd_comp = bd_face->GetBoundComp();

        int bc_id = bd_comp->get_physical_id();

        if (any_ref)
        {
          if (bc_id < min_ref)
          {
            min_ref = bc_id;
          }
        }
        else
        {
          min_ref = bc_id;
          any_ref = true;
        }
      }
    }
  }

#ifdef _MPI
  MPI_Allreduce(MPI_IN_PLACE, &min_ref, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
#endif

  for (int i = 0; i < n_cells; i++)
  {
    const TBaseCell* cell = coll->GetCell(i);

    int n_joints = cell->GetN_Joints();

    const int *wssDOF = wall_shear_stress_space->GetGlobalDOF(i);

    const TShapeDesc* shape_desc = cell->GetShapeDesc();

    const int* fv;
    const int* vpf;
    int max_vpf;

    // get face-vertex data
    shape_desc->GetFaceVertex(fv, vpf, max_vpf);

    // running offset into face-vertex array
    int fvo = 0;

    // check all face joints
    for (int j = 0; j < n_joints; j++)
    {
      const auto joint = cell->GetJoint(j);

      if (joint->GetType() == BoundaryFace || joint->GetType() == IsoBoundFace)
      {
        const auto bd_face = reinterpret_cast<const TBoundFace*>(joint);
        const auto bd_comp = bd_face->GetBoundComp();

        int bc_id = bd_comp->get_physical_id();

        if (bc_id != min_ref)
        {
          continue;
        }

        int num_face_verts = vpf[j];

        double cx = 0.0;
        double cy = 0.0;
        double cz = 0.0;

        // this is a boundary face! check out its vertices
        for (int r = 0; r < num_face_verts; r++)
        {
          int vertex = fv[fvo + r];
          int global_dof = wssDOF[vertex];

          boundary_dof.insert(global_dof);

          if (faces.count(global_dof) == 0)
          {
            faces.emplace(global_dof, 0);
          }

          int vertex_prev = fv[fvo + ((r + num_face_verts - 1) % num_face_verts)];
          int vertex_next = fv[fvo + ((r + 1) % num_face_verts)];

          double x, y, z;
          double x_prev, y_prev, z_prev;
          double x_next, y_next, z_next;

          cell->GetVertex(vertex)->GetCoords(x, y, z);

          cx += x;
          cy += y;
          cz += z;

          cell->GetVertex(vertex_prev)->GetCoords(x_prev, y_prev, z_prev);
          cell->GetVertex(vertex_next)->GetCoords(x_next, y_next, z_next);

          x_prev -= x;
          y_prev -= y;
          z_prev -= z;

          x_next -= x;
          y_next -= y;
          z_next -= z;

          double mag2_prev = x_prev * x_prev + y_prev * y_prev + z_prev * z_prev;
          double mag2_next = x_next * x_next + y_next * y_next + z_next * z_next;

          // angle beta at vertex b of triangle abc:
          // (a - b, c - b) = |a - b| |c - b| cos(beta)
          double face_angle = std::acos(
            (x_prev * x_next + y_prev * y_next + z_prev * z_next)
            / std::sqrt(mag2_prev * mag2_next));

          faces[global_dof].push_back(std::make_pair(num_relevant_faces, face_angle));
        }

        cx /= (double)num_face_verts;
        cy /= (double)num_face_verts;
        cz /= (double)num_face_verts;

        double nx, ny, nz;

        bd_comp->get_normal_vector(cx, cy, cz, nx, ny, nz);

        double nm = std::sqrt(nx * nx + ny * ny + nz * nz);

        if (nm > 0.0)
        {
          nx /= nm;
          ny /= nm;
          nz /= nm;
        }

        ++num_relevant_faces;
        cells.push_back(i);
        points.insert(points.end(), { cx, cy, cz });
        normals.insert(normals.end(), { nx, ny, nz });
      }

      // increment offset
      fvo += vpf[j];
    }
  }

  return wss_face_info;
}

const WallShearStressInfo& GetWallShearStressInfo(const TFEVectFunct3D &u,
  const TFEVectFunct3D &tau_w)
{
  // Caches or computes all the helper data needed for the wall shear stress
  // computation - adjacent cells, normals, etc.

  const auto velocity_space = u.GetFESpace();
  const auto coll = velocity_space->GetCollection();

  const auto wall_shear_stress_space = tau_w.GetFESpace();

  const auto* key = velocity_space.get();

  if (wss_info_cache.count(key) > 0)
  {
    return wss_info_cache[key];
  }

  auto& wss_info = wss_info_cache[key] = WallShearStressInfo();

  auto& boundary_dof = wss_info.boundary_dof;
  auto& normals = wss_info.normals;
  auto& points = wss_info.points;
  auto& cells = wss_info.cells;

  // boundary DOF to index of lowest adjacent boundary ref
  std::unordered_map<int, int> min_ref;

  // notice that our loop includes halo cells since we need them for
  // computing the relevant averages
  int n_cells = coll->GetN_Cells();

  double nx, ny, nz;
  double x, y, z;

  double x_prev, y_prev, z_prev;
  double x_next, y_next, z_next;

  // identify boundary DOF (in tau_w's space) and gather info
  for (int i = 0; i < n_cells; i++)
  {
    const TBaseCell* cell = coll->GetCell(i);
    bool is_tetra = cell->GetN_Vertices() == 4;

    const TShapeDesc* shape_desc = cell->GetShapeDesc();

    const int* fv;
    const int* vpf;
    int max_vpf;

    // get face-vertex data
    shape_desc->GetFaceVertex(fv, vpf, max_vpf);

    const int *wssDOF = wall_shear_stress_space->GetGlobalDOF(i);

    int n_joints = cell->GetN_Joints();

    // running offset into face-vertex array
    int fvo = 0;

    // check all face joints
    for (int j = 0; j < n_joints; j++)
    {
      const auto joint = cell->GetJoint(j);

      if (joint->GetType() == BoundaryFace || joint->GetType() == IsoBoundFace)
      {
        const auto bd_face = reinterpret_cast<const TBoundFace*>(joint);
        const auto bd_comp = bd_face->GetBoundComp();

        int bc_id = bd_comp->get_physical_id();
        int num_face_verts = vpf[j];

        // this is a boundary face! check out its vertices
        for (int r = 0; r < num_face_verts; r++)
        {
          int vertex = fv[fvo + r];
          int global_dof = wssDOF[vertex];

          boundary_dof.insert(global_dof);

          int vertex_prev = fv[fvo + ((r + num_face_verts - 1) % num_face_verts)];
          int vertex_next = fv[fvo + ((r + 1) % num_face_verts)];

          bool has_ref = min_ref.count(global_dof) > 0;
          int current_ref = has_ref ? min_ref[global_dof] : 0;

          cell->GetVertex(vertex)->GetCoords(x, y, z);
          bd_comp->get_normal_vector(x, y, z, nx, ny, nz);

          cell->GetVertex(vertex_prev)->GetCoords(x_prev, y_prev, z_prev);
          cell->GetVertex(vertex_next)->GetCoords(x_next, y_next, z_next);

          x_prev -= x;
          y_prev -= y;
          z_prev -= z;

          x_next -= x;
          y_next -= y;
          z_next -= z;

          double mag2_prev = x_prev * x_prev + y_prev * y_prev + z_prev * z_prev;
          double mag2_next = x_next * x_next + y_next * y_next + z_next * z_next;

          // angle beta at vertex b of triangle abc:
          // (a - b, c - b) = |a - b| |c - b| cos(beta)
          double face_angle = std::acos(
            (x_prev * x_next + y_prev * y_next + z_prev * z_next)
            / std::sqrt(mag2_prev * mag2_next));

          // apply weight to normal
          nx *= face_angle;
          ny *= face_angle;
          nz *= face_angle;

          if (points.count(global_dof) == 0)
          {
            // we haven't seen this vertex before, save its coordinates
            points.emplace(global_dof, 3);
            points[global_dof][0] = x;
            points[global_dof][1] = y;
            points[global_dof][2] = z;
          }

          if (!has_ref || bc_id < current_ref)
          {
            // we haven't seen this vertex before OR only on a later bc
            min_ref[global_dof] = bc_id;

            if (normals.count(global_dof) == 0)
            {
              // new vertex, allocate normal list
              normals.emplace(global_dof, 3);
            }

            // reset the accumulated normal to the one from this face
            normals[global_dof][0] = nx;
            normals[global_dof][1] = ny;
            normals[global_dof][2] = nz;
          }
          else if (bc_id == current_ref)
          {
            // we've memorized this bc for this vertex, accumulate normal
            normals[global_dof][0] += nx;
            normals[global_dof][1] += ny;
            normals[global_dof][2] += nz;
          }
          // else: we're on a later bc than we've memorized for this vertex

          if (cells.count(global_dof) == 0)
          {
            // we haven't seen this vertex before, allocate a cell list for it
            cells.emplace(global_dof, 0);
          }

          if (cells[global_dof].size() < 1 || cells[global_dof].back().first != i)
          {
            // add this cell to the vertex's cell list, if we haven't already
            // (note that since we're going through the cells in order, any
            // previous entries from this cell will be at the end)
            cells[global_dof].push_back(std::make_pair(i, is_tetra ? GetTetraSolidAngle(cell, vertex) : 1.0));
          }
        }
      }

      // increment offset
      fvo += vpf[j];
    }
  }

  // loop over the boundary vertices and normalize normal vectors at each one
  for (int global_dof: boundary_dof)
  {
    std::vector<double> &normal = normals[global_dof];

    double nx = normal[0];
    double ny = normal[1];
    double nz = normal[2];

    double inv_normal = 1.0 / std::sqrt(nx * nx + ny * ny + nz * nz);

    normal[0] *= inv_normal;
    normal[1] *= inv_normal;
    normal[2] *= inv_normal;
  }

  return wss_info;
}

void ComputeWallShearStress(const double nu, const TFEVectFunct3D &u,
  TFEVectFunct3D &tau_w, const ViscositySettings& settings)
{
  // Computes the (vector-valued) wall shear stress $\tau_w$ on each boundary
  // vertex, i.e. the normal derivative of the tangential velocity at that
  // point. The normal used is the averaged normal of each adjacent boundary
  // face belonging to the lowest adjacent boundary ref; the velocity gradient
  // is averaged from all adjacent cells. The normals are weighted by the
  // interior angle of each face at the vertex; the cells are weighted by the
  // interior *solid* angle at the vertex, if they are tetrahedra; otherwise
  // the weights are identically one.
  //
  // So if $n$ is the normal vector at the vertex, the vector computed is
  //
  // $$(\tau_w)_j = \nu n \cdot \nabla v_j$$
  //
  // where
  //
  // $$v = u - (u \cdot n) n.$$
  //
  // So:
  //
  // $$(\tau_w)_j = \nu \sum_k n_k (\partial_k u_j - n_j \sum_i n_i \partial_k u_i).$$

  bool newtonian = settings.mode == ViscosityMode::Newtonian;

  if (u.GetN_Components() != 3)
  {
    ErrThrow("Velocity function for wall shear stress computation"
      "has to have 3 components!");
  }

  if (tau_w.GetN_Components() != 3)
  {
    ErrThrow("Wall shear stress function has to have 3 components!");
  }

  const auto wall_shear_stress_space = tau_w.GetFESpace();

  for (auto fe: wall_shear_stress_space->GetUsedElements())
  {
    if (fe != C_P1_3D_T_A && fe != C_Q1_3D_H_M && fe != C_Q1_3D_H_A)
    {
      ErrThrow("Invalid wall shear stress function space!");
    }
  }

  const auto velocity_space = u.GetFESpace();
  int wss_n_global_dof = wall_shear_stress_space->get_n_dof();
  double* wall_shear_stress = tau_w.GetValues();

  const auto& wss_info = GetWallShearStressInfo(u, tau_w);

  const auto& boundary_dof = wss_info.boundary_dof;
  const auto& normals = wss_info.normals;
  const auto& points = wss_info.points;
  const auto& cells = wss_info.cells;

  // length is 4 since ...[0] is the value
  std::vector<double> u1_grad_tmp(4);
  std::vector<double> u2_grad_tmp(4);
  std::vector<double> u3_grad_tmp(4);

  std::vector<double> u1_grad(4);
  std::vector<double> u2_grad(4);
  std::vector<double> u3_grad(4);

  // loop over the boundary vertices and compute WSS at each one
  for (int global_dof: boundary_dof)
  {
    const std::vector<double> &normal = normals.at(global_dof);
    const std::vector<double> &point = points.at(global_dof);

    double nx = normal[0];
    double ny = normal[1];
    double nz = normal[2];

    double x = point[0];
    double y = point[1];
    double z = point[2];

    for (int s = 1; s < 4; s++)
    {
      u1_grad[s] = 0.0;
      u2_grad[s] = 0.0;
      u3_grad[s] = 0.0;
    }

    double local_cell_sum = 0.0;

    // accumulate velocity gradient values from adjacent cells
    for (const auto& pair: cells.at(global_dof))
    {
      int i = pair.first;
      double weight = pair.second;

      u.FindVectGradientLocal(i, x, y, z,
        u1_grad_tmp.data(), u2_grad_tmp.data(), u3_grad_tmp.data());

      for (int s = 1; s < 4; s++)
      {
        u1_grad[s] += weight * u1_grad_tmp[s];
        u2_grad[s] += weight * u2_grad_tmp[s];
        u3_grad[s] += weight * u3_grad_tmp[s];
      }

      local_cell_sum += weight;
    }

    // average gradient values
    double inv_local_cells = 1.0 / local_cell_sum;

    for (int s = 1; s < 4; s++)
    {
      u1_grad[s] *= inv_local_cells;
      u2_grad[s] *= inv_local_cells;
      u3_grad[s] *= inv_local_cells;
    }

    // $$(\tau_w)_j = \nu \sum_k n_k (\partial_k u_j - n_j \sum_i n_i \partial_k u_i)$$

    double w1 = nx * (u1_grad[1] - nx * (nx * u1_grad[1] + ny * u2_grad[1] + nz * u3_grad[1]))
              + ny * (u1_grad[2] - nx * (nx * u1_grad[2] + ny * u2_grad[2] + nz * u3_grad[2]))
              + nz * (u1_grad[3] - nx * (nx * u1_grad[3] + ny * u2_grad[3] + nz * u3_grad[3]));

    double w2 = nx * (u2_grad[1] - ny * (nx * u1_grad[1] + ny * u2_grad[1] + nz * u3_grad[1]))
              + ny * (u2_grad[2] - ny * (nx * u1_grad[2] + ny * u2_grad[2] + nz * u3_grad[2]))
              + nz * (u2_grad[3] - ny * (nx * u1_grad[3] + ny * u2_grad[3] + nz * u3_grad[3]));

    double w3 = nx * (u3_grad[1] - nz * (nx * u1_grad[1] + ny * u2_grad[1] + nz * u3_grad[1]))
              + ny * (u3_grad[2] - nz * (nx * u1_grad[2] + ny * u2_grad[2] + nz * u3_grad[2]))
              + nz * (u3_grad[3] - nz * (nx * u1_grad[3] + ny * u2_grad[3] + nz * u3_grad[3]));

    if (newtonian)
    {
      // multiply with viscosity
      wall_shear_stress[global_dof] = nu * w1;
      wall_shear_stress[global_dof + wss_n_global_dof] = nu * w2;
      wall_shear_stress[global_dof + 2 * wss_n_global_dof] = nu * w3;
    }
    else
    {
      double grad_u[9];

      for (int i = 0; i < 3; i++)
      {
        grad_u[0 + i] = u1_grad[i + 1];
        grad_u[3 + i] = u2_grad[i + 1];
        grad_u[6 + i] = u3_grad[i + 1];
      }

      double nu_eff = NonNewtonianViscosity<3>(grad_u, settings);

      // multiply with viscosity
      wall_shear_stress[global_dof] = nu_eff * w1;
      wall_shear_stress[global_dof + wss_n_global_dof] = nu_eff * w2;
      wall_shear_stress[global_dof + 2 * wss_n_global_dof] = nu_eff * w3;
    }
  }

#if _MPI
  // make level 3 consistent

  const TParFECommunicator3D& comm = ((const TFESpace3D*)wall_shear_stress_space.get())->get_communicator();

  comm.queue_consistency_update(wall_shear_stress, 3);
  comm.queue_consistency_update(wall_shear_stress + wss_n_global_dof, 3);
  comm.queue_consistency_update(wall_shear_stress + 2 * wss_n_global_dof, 3);

  TParFECommunicator3D::flush_consistency_updates();
#endif
}

void ComputeWallShearStressFace(const double nu, const TFEVectFunct3D &u,
  TFEVectFunct3D &tau_w, const ViscositySettings& settings)
{
  // Computes the (vector-valued) wall shear stress $\tau_w$ on each boundary
  // face, i.e. the normal derivative of the tangential velocity at the face's
  // centroid.
  //
  // So if $n$ is the face's normal vector, the vector computed is
  //
  // $$(\tau_w)_j = \nu n \cdot \nabla v_j$$
  //
  // where
  //
  // $$v = u - (u \cdot n) n.$$
  //
  // So:
  //
  // $$(\tau_w)_j = \nu \sum_k n_k (\partial_k u_j - n_j \sum_i n_i \partial_k u_i).$$
  //
  // These are then averaged (weighted by face angle) over all adjacent
  // lowest-ref boundary faces on each boundary vertex.

  bool newtonian = settings.mode == ViscosityMode::Newtonian;

  if (u.GetN_Components() != 3)
  {
    ErrThrow("Velocity function for wall shear stress computation"
      "has to have 3 components!");
  }

  if (tau_w.GetN_Components() != 3)
  {
    ErrThrow("Wall shear stress function has to have 3 components!");
  }

  const auto wall_shear_stress_space = tau_w.GetFESpace();

  for (auto fe: wall_shear_stress_space->GetUsedElements())
  {
    if (fe != C_P1_3D_T_A && fe != C_Q1_3D_H_M && fe != C_Q1_3D_H_A)
    {
      ErrThrow("Invalid wall shear stress function space!");
    }
  }

  const auto velocity_space = u.GetFESpace();
  int wss_n_global_dof = wall_shear_stress_space->get_n_dof();
  double* wall_shear_stress = tau_w.GetValues();

  const auto& wss_face_info = GetWallShearStressFaceInfo(u, tau_w);

  auto& boundary_dof = wss_face_info.boundary_dof;
  auto& faces = wss_face_info.faces;
  auto& normals = wss_face_info.normals;
  auto& points = wss_face_info.points;
  auto& cells = wss_face_info.cells;

  int num_faces = (int)cells.size();
  std::vector<double> wss_per_face(3 * num_faces);

  // length is 4 since ...[0] is the value
  std::vector<double> u1_grad(4);
  std::vector<double> u2_grad(4);
  std::vector<double> u3_grad(4);

  // loop over the boundary faces and compute WSS at each one
  for (int i = 0; i < num_faces; i++)
  {
    int m = 3 * i;

    double cx = points[m];
    double cy = points[m + 1];
    double cz = points[m + 2];

    double nx = normals[m];
    double ny = normals[m + 1];
    double nz = normals[m + 2];

    u.FindVectGradientLocal(cells[i], cx, cy, cz,
      u1_grad.data(), u2_grad.data(), u3_grad.data());

    // $$(\tau_w)_j = \nu \sum_k n_k (\partial_k u_j - n_j \sum_i n_i \partial_k u_i)$$

    double w1 = nx * (u1_grad[1] - nx * (nx * u1_grad[1] + ny * u2_grad[1] + nz * u3_grad[1]))
              + ny * (u1_grad[2] - nx * (nx * u1_grad[2] + ny * u2_grad[2] + nz * u3_grad[2]))
              + nz * (u1_grad[3] - nx * (nx * u1_grad[3] + ny * u2_grad[3] + nz * u3_grad[3]));

    double w2 = nx * (u2_grad[1] - ny * (nx * u1_grad[1] + ny * u2_grad[1] + nz * u3_grad[1]))
              + ny * (u2_grad[2] - ny * (nx * u1_grad[2] + ny * u2_grad[2] + nz * u3_grad[2]))
              + nz * (u2_grad[3] - ny * (nx * u1_grad[3] + ny * u2_grad[3] + nz * u3_grad[3]));

    double w3 = nx * (u3_grad[1] - nz * (nx * u1_grad[1] + ny * u2_grad[1] + nz * u3_grad[1]))
              + ny * (u3_grad[2] - nz * (nx * u1_grad[2] + ny * u2_grad[2] + nz * u3_grad[2]))
              + nz * (u3_grad[3] - nz * (nx * u1_grad[3] + ny * u2_grad[3] + nz * u3_grad[3]));

    wss_per_face[m] = nu * w1;
    wss_per_face[m + 1] = nu * w2;
    wss_per_face[m + 2] = nu * w3;

    if (newtonian)
    {
      wss_per_face[m] = nu * w1;
      wss_per_face[m + 1] = nu * w2;
      wss_per_face[m + 2] = nu * w3;
    }
    else
    {
      double grad_u[9];

      for (int i = 0; i < 3; i++)
      {
        grad_u[0 + i] = u1_grad[i + 1];
        grad_u[3 + i] = u2_grad[i + 1];
        grad_u[6 + i] = u3_grad[i + 1];
      }

      double nu_eff = NonNewtonianViscosity<3>(grad_u, settings);

      wss_per_face[m] = nu_eff * w1;
      wss_per_face[m + 1] = nu_eff * w2;
      wss_per_face[m + 2] = nu_eff * w3;
    }
  }

  // loop over boundary dofs and average face WSS
  for (int global_dof: boundary_dof)
  {
    double w1 = 0.0;
    double w2 = 0.0;
    double w3 = 0.0;
    double ws = 0.0;

    for (const auto& pair: faces.at(global_dof))
    {
      int m = 3 * pair.first;
      double weight = pair.second;

      ws += weight;
      w1 += weight * wss_per_face[m];
      w2 += weight * wss_per_face[m + 1];
      w3 += weight * wss_per_face[m + 2];
    }

    if (ws > 0.0)
    {
      wall_shear_stress[global_dof] = w1 / ws;
      wall_shear_stress[global_dof + wss_n_global_dof] = w2 / ws;
      wall_shear_stress[global_dof + 2 * wss_n_global_dof] = w3 / ws;
    }
  }

#if _MPI
  // make level 3 consistent

  const TParFECommunicator3D& comm = ((const TFESpace3D*)wall_shear_stress_space.get())->get_communicator();

  comm.queue_consistency_update(wall_shear_stress, 3);
  comm.queue_consistency_update(wall_shear_stress + wss_n_global_dof, 3);
  comm.queue_consistency_update(wall_shear_stress + 2 * wss_n_global_dof, 3);

  TParFECommunicator3D::flush_consistency_updates();
#endif
}

void TransformReference(std::shared_ptr<const TFESpace> space,
  int cell_index, int n_points,
  const double* xi, const double* eta, const double* zeta,
  double* x, double* y, double* z)
{
  const TBaseCell* cell = space->GetCollection()->GetCell(cell_index);
  const FiniteElement& element = space->get_fe(cell_index);

  ReferenceTransformation_type ref_trans = element.GetRefTransID();
  BFRefElements ref_element = element.GetBaseFunct()->GetRefElement();

  if (TDatabase::ParamDB->USE_ISOPARAMETRIC
    && cell->has_isoparametric_joint())
  {
    switch (ref_element)
    {
      case BFRefElements::BFUnitHexahedron:
        ref_trans = ReferenceTransformation_type::HexaIsoparametric;
        break;

      case BFRefElements::BFUnitTetrahedron:
        ref_trans = ReferenceTransformation_type::TetraIsoparametric;
        break;

      default:
        Output::root_warn("RefTrans", "Unknown isoparametric reference "
          "transformation for reference element ", ref_element);
        break;
    }
  }

  TRefTrans3D *F_K = FEDatabase::GetRefTrans3D(ref_trans);

  F_K->SetCell(cell);
  F_K->GetOrigFromRef(n_points,
    xi, eta, zeta, x, y, z);
}

double ComputeTurbulentKineticEnergySmagorinsky(const TFEVectFunct3D &u, double scale)
{
  if (u.GetN_Components() != 3)
  {
    ErrThrow("Velocity function for turbulent kinetic energy computation"
      "has to have 3 components!");
  }

  const auto velocity_space = u.GetFESpace();
  const auto coll = velocity_space->GetCollection();
  int n_cells = coll->GetN_Cells();

  double sgs_energy = 0.0;

  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;

  std::array<double, 4> u1_grad;
  std::array<double, 4> u2_grad;
  std::array<double, 4> u3_grad;

  std::array<double, 3> u_val;
  std::array<double, 9> u_grad;

  for (int c = 0; c < n_cells; c++)
  {
    const TBaseCell* cell = coll->GetCell(c);

    const FiniteElement& element = velocity_space->get_fe(c);
    int fe_degree = element.GetBaseFunct()->GetPolynomialDegree();

#ifdef _MPI
    if (cell->IsHaloCell())
    {
      continue;
    }
#endif

    const TQuadFormula* quad_formula = QuadratureFormulaDatabase::
      qf_from_degree(2 * fe_degree, element.GetBaseFunct()->GetRefElement());

    int n_points = quad_formula->GetN_QuadPoints();

    x.resize(n_points);
    y.resize(n_points);
    z.resize(n_points);

    const double* xi = quad_formula->get_xi();
    const double* eta = quad_formula->get_eta();
    const double* zeta = quad_formula->get_zeta();

    TransformReference(velocity_space, c, n_points,
      xi, eta, zeta, x.data(), y.data(), z.data());

    double hK = cell->Get_hK(TDatabase::ParamDB->CELL_MEASURE);
    double v = cell->GetMeasure();

    for (int i = 0; i < n_points; i++)
    {
      double weight = v * quad_formula->get_weight(i);

      for (int j = 0; j < 4; j++)
      {
        u1_grad[j] = 0.0;
        u2_grad[j] = 0.0;
        u3_grad[j] = 0.0;
      }

      u.FindVectGradientLocal(c, x[i], y[i], z[i],
        u1_grad.data(), u2_grad.data(), u3_grad.data());

      u_val[0] = u1_grad[0];
      u_val[1] = u2_grad[0];
      u_val[2] = u3_grad[0];

      for (int j = 0; j < 3; j++)
      {
        u_grad[j + 0] = u1_grad[j];
        u_grad[j + 3] = u2_grad[j];
        u_grad[j + 6] = u3_grad[j];
      }

      sgs_energy += weight * sgsKineticEnergy3D(hK,
          u_val.data(), u_grad.data(),
          nullptr,
          &(x[i]), &(y[i]), &(z[i]),
          -1.0, scale);
    }
  }

#ifdef _MPI
  MPI_Allreduce(MPI_IN_PLACE, &sgs_energy, 1, MPI_DOUBLE,
    MPI_SUM, MPI_COMM_WORLD);
#endif

  return sgs_energy;
}

double ComputeTurbulentKineticEnergyRBVMS(const TFEVectFunct3D &u,
  const TFEVectFunct3D &u_m1, const TFEFunction3D &p, double nu,
  RBVMS_Settings settings)
{
  if (u.GetN_Components() != 3)
  {
    ErrThrow("Velocity function for turbulent kinetic energy computation"
      "has to have 3 components!");
  }

  double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;

  const auto velocity_space = u.GetFESpace();
  const auto coll = velocity_space->GetCollection();
  int n_velocity_dofs = velocity_space->get_n_dof();
  int n_cells = coll->GetN_Cells();

  const std::array<const double*, 3> u_entries = { u.GetValues(),
    u.GetValues() + n_velocity_dofs,
    u.GetValues() + 2 * n_velocity_dofs };

  const std::array<const double*, 3> u_m1_entries = { u_m1.GetValues(),
    u_m1.GetValues() + n_velocity_dofs,
    u_m1.GetValues() + 2 * n_velocity_dofs };

  double sgs_energy = 0.0;

  std::array<double, 3> u_val;
  std::array<double, 3> ut_val;
  std::array<double, 3> u_x;
  std::array<double, 3> u_y;
  std::array<double, 3> u_z;
  std::array<double, 3> u_lap;
  std::array<double, 3> Dp;

  std::array<double, 9> J_inv;

  std::array<double, 9> G;
  std::array<double, 3> g;

  std::array<double, 3> res_m;

  for (int c = 0; c < n_cells; c++)
  {
    TBaseCell* cell = coll->GetCell(c);

    const FiniteElement& element = velocity_space->get_fe(c);
    const BaseFunctions* bf = element.GetBaseFunct();

    const int* cell_dof = velocity_space->GetGlobalDOF(c);
    const int n_local_dof = velocity_space->get_n_local_dof(c);

#ifdef _MPI
    if (cell->IsHaloCell())
    {
      continue;
    }
#endif

    TQuadFormula qf_ref(QuadratureFormula_type::BaryCenterTetra); // dummy type
    TQuadFormula qf_orig(qf_ref);

    bool tmp = true;

    auto ref_trans = FEDatabase::GetOrig({ &element }, coll, cell,
      &tmp, qf_ref, qf_orig);

    TRefTrans3D *F_K = FEDatabase::GetRefTrans3D(ref_trans);
    F_K->SetCell(cell);

    const double* xi = qf_orig.get_xi();
    const double* eta = qf_orig.get_eta();
    const double* zeta = qf_orig.get_zeta();

    double **v = FEDatabase::GetOrigElementValues(*bf, MultiIndex3D::D000);
    double **v_x = FEDatabase::GetOrigElementValues(*bf, MultiIndex3D::D100);
    double **v_y = FEDatabase::GetOrigElementValues(*bf, MultiIndex3D::D010);
    double **v_z = FEDatabase::GetOrigElementValues(*bf, MultiIndex3D::D001);
    double **v_xx = FEDatabase::GetOrigElementValues(*bf, MultiIndex3D::D200);
    double **v_yy = FEDatabase::GetOrigElementValues(*bf, MultiIndex3D::D020);
    double **v_zz = FEDatabase::GetOrigElementValues(*bf, MultiIndex3D::D002);

    int n_points = qf_orig.GetN_QuadPoints();

    double hK = cell->Get_hK(TDatabase::ParamDB->CELL_MEASURE);

    for (int i = 0; i < n_points; i++)
    {
      u_val.fill(0.0);
      ut_val.fill(0.0);
      u_x.fill(0.0);
      u_y.fill(0.0);
      u_z.fill(0.0);
      u_lap.fill(0.0);
      Dp.fill(0.0);

      auto point = qf_orig.get_point(i);

      p.FindGradientLocal(c, point.x, point.y, point.z, Dp.data());

      for (int j = 0; j < n_local_dof; j++)
      {
        int dof = cell_dof[j];

        for (int s = 0; s < 3; s++)
        {
          u_val[s] += v[i][j] * u_entries[s][dof];
          ut_val[s] += v[i][j] * u_m1_entries[s][dof];

          u_x[s] += v_x[i][j] * u_entries[s][dof];
          u_y[s] += v_y[i][j] * u_entries[s][dof];
          u_z[s] += v_z[i][j] * u_entries[s][dof];

          u_lap[s] += v_xx[i][j] * u_entries[s][dof];
          u_lap[s] += v_yy[i][j] * u_entries[s][dof];
          u_lap[s] += v_zz[i][j] * u_entries[s][dof];
        }
      }

      for (int s = 0; s < 3; s++)
      {
        ut_val[s] = (u_val[s] - ut_val[s]) / tau;
      }

      double tau_m = 1.0;

      switch (settings.mode)
      {
        case RBVMS_ParamMode::G:
        {
          F_K->GetInverseTransformationDerivatives(xi[i], eta[i], zeta[i],
            J_inv.data());

          for (int s = 0; s < 3; s++)
          {
            g[s] = 0.0;

            for (int j = 0; j < 3; j++)
            {
              g[s] += J_inv[3 * s + j]; // d x^_j / d x_i

              int m = s + 3 * j; // G_ij
              G[m] = 0.0;
              for (int k = 0; k < 3; k++)
              {
                // d x^_k / d x_i * d x^_k / d x_j
                G[m] += J_inv[3 * s + k] * J_inv[3 * j + k];
              }
            }
          }

          tau_m = RBVMS_Param_G<3>(u_val.data(), G.data(), g.data(),
            tau, nu, settings).first;

          break;
        }

        case RBVMS_ParamMode::H:
        {
          tau_m = RBVMS_Param_H<3>(hK, settings).first;
          break;
        }

        case RBVMS_ParamMode::Codina:
        {
          ErrThrow("Cannot use Codina parameters in static RBVMS!");
          break;
        }
      }

      // the subgrid velocity is modelled as \tau_m times the momentum residue
      // (so the SGS energy is estimated as half the square integral of that)

      for (int s = 0; s < 3; s++)
      {
        res_m[s] = ut_val[s] - nu * u_lap[s]
        + (u_val[0] * u_x[s] + u_val[1] * u_y[s] + u_val[2] * u_z[s])
        + Dp[s];
      }

      double k_i = tau_m * tau_m *
        (res_m[0] * res_m[0] + res_m[1] * res_m[1] + res_m[2] * res_m[2]);

      sgs_energy += qf_orig.get_weight(i) * k_i;
    }
  }

#ifdef _MPI
  MPI_Allreduce(MPI_IN_PLACE, &sgs_energy, 1, MPI_DOUBLE,
    MPI_SUM, MPI_COMM_WORLD);
#endif

  return 0.5 * sgs_energy;
}

double ComputeTurbulentKineticEnergyRBVMSTime(const TFEVectFunct3D &u,
  std::shared_ptr<PointwiseAssemblyData> persistent_data)
{
  if (persistent_data == nullptr)
  {
    ErrThrow("Tried to use a method that needs persistent pointwise data "
        "without providing storage!");
  }

  double sgs_energy = 0.0;

  const auto velocity_space = u.GetFESpace();
  const auto coll = velocity_space->GetCollection();
  int n_cells = coll->GetN_Cells();

  TQuadFormula qf_ref(QuadratureFormula_type::BaryCenterTetra); // dummy type
  TQuadFormula qf_orig(qf_ref);

  bool second_derivatives = false;

  for (int c = 0; c < n_cells; c++)
  {
    TBaseCell* cell = coll->GetCell(c);

    FEDatabase::GetOrig({ &(velocity_space->get_fe(c)) },
      coll, cell, &second_derivatives, qf_ref, qf_orig);

    int n_points = qf_orig.GetN_QuadPoints();

    persistent_data->SetCurrentCell(c, n_points, 2 * 3);

    for (int i = 0; i < n_points; i++)
    {
      persistent_data->SetCurrentPoint(i);

      double w = qf_orig.get_weight(i);

      const double* u_prime = PointwiseAssemblyData::GetCurrentData();

      sgs_energy += w * (u_prime[0] * u_prime[0]
                       + u_prime[1] * u_prime[1]
                       + u_prime[2] * u_prime[2]);
    }
  }

#ifdef _MPI
  MPI_Allreduce(MPI_IN_PLACE, &sgs_energy, 1, MPI_DOUBLE,
    MPI_SUM, MPI_COMM_WORLD);
#endif

  return 0.5 * sgs_energy;
}

void ComputeTurbulentKineticEnergyCell(const TFEVectFunct3D &u,
  TFEFunction3D &tke, double scale)
{
  if (u.GetN_Components() != 3)
  {
    ErrThrow("Velocity function for turbulent kinetic energy computation"
      "has to have 3 components!");
  }

  const auto velocity_space = u.GetFESpace();
  const auto tke_space = tke.GetFESpace();

  int tke_n_global_dof = tke_space->get_n_dof();

  double* turbulent_kinetic_energy = tke.GetValues();

  const auto coll = velocity_space->GetCollection();
  if (coll != tke_space->GetCollection())
  {
    ErrThrow("Velocity space and TKE space live on different grids!");
  }

  int n_cells = coll->GetN_Cells();

  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;

  std::vector<int> num_dof_cells(tke_n_global_dof, 0);
  std::vector<double> k_at_dof(tke_n_global_dof, 0.0);

  std::array<double, 4> u1_grad;
  std::array<double, 4> u2_grad;
  std::array<double, 4> u3_grad;

  std::array<double, 3> u_val;
  std::array<double, 9> u_grad;

  for (int c = 0; c < n_cells; c++)
  {
    const TBaseCell* cell = coll->GetCell(c);

    int nv = cell->GetN_Vertices();
    for (int i = 0; i < nv; i++)
    {
      TDatabase::ParamDB->INTERNAL_VERTEX_X[i] = cell->GetVertex(i)->GetX();
      TDatabase::ParamDB->INTERNAL_VERTEX_Y[i] = cell->GetVertex(i)->GetY();
      TDatabase::ParamDB->INTERNAL_VERTEX_Z[i] = cell->GetVertex(i)->GetZ();
    }

    const FiniteElement& element = velocity_space->get_fe(c);
    int fe_degree = element.GetBaseFunct()->GetPolynomialDegree();

    const TQuadFormula* quad_formula = QuadratureFormulaDatabase::
      qf_from_degree(2 * fe_degree, element.GetBaseFunct()->GetRefElement());

    int n_points = quad_formula->GetN_QuadPoints();

    x.resize(n_points);
    y.resize(n_points);
    z.resize(n_points);

    const double* xi = quad_formula->get_xi();
    const double* eta = quad_formula->get_eta();
    const double* zeta = quad_formula->get_zeta();

    TransformReference(velocity_space, c, n_points,
      xi, eta, zeta, x.data(), y.data(), z.data());

    double hK = cell->Get_hK(TDatabase::ParamDB->CELL_MEASURE);

    double k_loc = 0.0;
    double w_sum = 0.0;

    for (int i = 0; i < n_points; i++)
    {
      double weight = quad_formula->get_weight(i);
      w_sum += weight;

      for (int j = 0; j < 4; j++)
      {
        u1_grad[j] = 0.0;
        u2_grad[j] = 0.0;
        u3_grad[j] = 0.0;
      }

      u.FindVectGradientLocal(c, x[i], y[i], z[i],
        u1_grad.data(), u2_grad.data(), u3_grad.data());

      u_val[0] = u1_grad[0];
      u_val[1] = u2_grad[0];
      u_val[2] = u3_grad[0];

      for (int j = 0; j < 3; j++)
      {
        u_grad[j + 0] = u1_grad[j];
        u_grad[j + 3] = u2_grad[j];
        u_grad[j + 6] = u3_grad[j];
      }

      k_loc += weight * sgsKineticEnergy3D(hK,
          u_val.data(), u_grad.data(),
          nullptr,
          &(x[i]), &(y[i]), &(z[i]),
          -1.0, scale);
    }

    k_loc /= w_sum;

    const int* DOF = tke_space->GetGlobalDOF(c);
    int n_dof_local = tke_space->get_n_local_dof(c);

    for (int i = 0; i < n_dof_local; i++)
    {
      int dof = DOF[i];

      num_dof_cells[dof]++;
      k_at_dof[dof] += k_loc;
    }
  }

  for (int dof = 0; dof < tke_n_global_dof; dof++)
  {
    if (num_dof_cells[dof] > 0)
    {
      turbulent_kinetic_energy[dof] = k_at_dof[dof] / num_dof_cells[dof];
    }
    else
    {
      turbulent_kinetic_energy[dof] = 0.0;
    }
  }

#if _MPI
  // make level 3 consistent

  const TParFECommunicator3D& comm = ((const TFESpace3D*)tke_space.get())->get_communicator();

  comm.consistency_update(turbulent_kinetic_energy, 3);

  TParFECommunicator3D::flush_consistency_updates();
#endif
}

void ComputeTurbulentKineticEnergy(const TFEVectFunct3D &u,
  TFEFunction3D &tke, double scale)
{
  if (u.GetN_Components() != 3)
  {
    ErrThrow("Velocity function for turbulent kinetic energy computation"
      "has to have 3 components!");
  }

  const auto velocity_space = u.GetFESpace();
  const auto tke_space = tke.GetFESpace();

  int tke_n_global_dof = tke_space->get_n_dof();

  double* turbulent_kinetic_energy = tke.GetValues();

  const auto coll = velocity_space->GetCollection();
  if (coll != tke_space->GetCollection())
  {
    ErrThrow("Velocity space and TKE space live on different grids!");
  }

  int n_cells = coll->GetN_Cells();

  std::vector<double> x(MaxN_PointsForNodal3D);
  std::vector<double> y(MaxN_PointsForNodal3D);
  std::vector<double> z(MaxN_PointsForNodal3D);

  std::vector<int> num_dof_cells(tke_n_global_dof, 0);
  std::vector<double> k_at_dof(tke_n_global_dof, 0.0);

  std::array<double, 4> u1_grad;
  std::array<double, 4> u2_grad;
  std::array<double, 4> u3_grad;

  std::array<double, 3> u_val;
  std::array<double, 9> u_grad;

  const double* xi;
  const double* eta;
  const double* zeta;

  for (int c = 0; c < n_cells; c++)
  {
    const TBaseCell* cell = coll->GetCell(c);

    const int* DOF = tke_space->GetGlobalDOF(c);

    const FiniteElement& element = tke_space->get_fe(c);
    const NodalFunctional* nf = element.GetNodalFunctional();

    int n_points_local;
    nf->GetPointsForAll(n_points_local, xi, eta, zeta);

    int n_dof_local = element.GetN_DOF();

    TransformReference(tke_space, c, n_points_local,
      xi, eta, zeta,
      x.data(), y.data(), z.data());

    double hK = cell->Get_hK(TDatabase::ParamDB->CELL_MEASURE);

    for (int i = 0; i < n_dof_local; i++)
    {
      int dof = DOF[i];

      num_dof_cells[dof]++;

      for (int j = 0; j < 4; j++)
      {
        u1_grad[j] = 0.0;
        u2_grad[j] = 0.0;
        u3_grad[j] = 0.0;
      }

      u.FindVectGradientLocal(c, x[i], y[i], z[i],
        u1_grad.data(), u2_grad.data(), u3_grad.data());

      u_val[0] = u1_grad[0];
      u_val[1] = u2_grad[0];
      u_val[2] = u3_grad[0];

      for (int j = 0; j < 3; j++)
      {
        u_grad[j + 0] = u1_grad[j];
        u_grad[j + 3] = u2_grad[j];
        u_grad[j + 6] = u3_grad[j];
      }

      k_at_dof[dof] += sgsKineticEnergy3D(hK,
          u_val.data(), u_grad.data(),
          nullptr,
          &(x[i]), &(y[i]), &(z[i]),
          -1.0, scale);
    }

    for (int dof = 0; dof < tke_n_global_dof; dof++)
    {
      if (num_dof_cells[dof] > 0)
      {
        turbulent_kinetic_energy[dof] = k_at_dof[dof] / num_dof_cells[dof];
      }
      else
      {
        turbulent_kinetic_energy[dof] = 0.0;
      }
    }
  }

#if _MPI
  // make level 3 consistent

  const TParFECommunicator3D& comm = ((const TFESpace3D*)tke_space.get())->get_communicator();

  comm.consistency_update(turbulent_kinetic_energy, 3);

  TParFECommunicator3D::flush_consistency_updates();
#endif
}

void ComputeEffectiveViscosity(const TFEVectFunct3D &u,
  TFEFunction3D &nu_eff, const ViscositySettings& settings)
{
  if (u.GetN_Components() != 3)
  {
    ErrThrow("Velocity function for effective viscosity computation"
      "has to have 3 components!");
  }

  const auto velocity_space = u.GetFESpace();
  const auto nu_space = nu_eff.GetFESpace();

  int nu_n_global_dof = nu_space->get_n_dof();

  double* effective_viscosity = nu_eff.GetValues();

  const auto coll = velocity_space->GetCollection();
  if (coll != nu_space->GetCollection())
  {
    ErrThrow("Velocity space and viscosity space live on different grids!");
  }

  int n_cells = coll->GetN_Cells();

  std::vector<double> x(MaxN_PointsForNodal3D);
  std::vector<double> y(MaxN_PointsForNodal3D);
  std::vector<double> z(MaxN_PointsForNodal3D);

  std::vector<int> num_dof_cells(nu_n_global_dof, 0);
  std::vector<double> nu_at_dof(nu_n_global_dof, 0.0);

  std::array<double, 4> u1_grad;
  std::array<double, 4> u2_grad;
  std::array<double, 4> u3_grad;

  std::array<double, 3> u_val;
  std::array<double, 9> u_grad;

  const double* xi;
  const double* eta;
  const double* zeta;

  for (int c = 0; c < n_cells; c++)
  {
    const int* DOF = nu_space->GetGlobalDOF(c);

    const FiniteElement& element = nu_space->get_fe(c);
    const NodalFunctional* nf = element.GetNodalFunctional();

    int n_points_local;
    nf->GetPointsForAll(n_points_local, xi, eta, zeta);

    int n_dof_local = element.GetN_DOF();

    TransformReference(nu_space, c, n_points_local,
      xi, eta, zeta,
      x.data(), y.data(), z.data());

    for (int i = 0; i < n_dof_local; i++)
    {
      int dof = DOF[i];

      num_dof_cells[dof]++;

      for (int j = 0; j < 4; j++)
      {
        u1_grad[j] = 0.0;
        u2_grad[j] = 0.0;
        u3_grad[j] = 0.0;
      }

      u.FindVectGradientLocal(c, x[i], y[i], z[i],
        u1_grad.data(), u2_grad.data(), u3_grad.data());

      u_val[0] = u1_grad[0];
      u_val[1] = u2_grad[0];
      u_val[2] = u3_grad[0];

      for (int j = 0; j < 3; j++)
      {
        u_grad[j + 0] = u1_grad[j];
        u_grad[j + 3] = u2_grad[j];
        u_grad[j + 6] = u3_grad[j];
      }

      nu_at_dof[dof] += NonNewtonianViscosity<3>(u_grad.data(), settings);
    }

    for (int dof = 0; dof < nu_n_global_dof; dof++)
    {
      if (num_dof_cells[dof] > 0)
      {
        effective_viscosity[dof] = nu_at_dof[dof] / num_dof_cells[dof];
      }
      else
      {
        effective_viscosity[dof] = 0.0;
      }
    }
  }

#if _MPI
  // make level 3 consistent

  const TParFECommunicator3D& comm = ((const TFESpace3D*)nu_space.get())->get_communicator();

  comm.consistency_update(effective_viscosity, 3);

  TParFECommunicator3D::flush_consistency_updates();
#endif
}

#endif // 3D

#ifdef __2D__
// ========================================================================
// put l infinity norm of u in coeff5
// ========================================================================
void LInfU(int N_Points, double **, double **Params, TBaseCell *)
{
  int i;
  double max, u1, u2, u, *param;

  max = -1;

  for(i=0;i<N_Points;i++)
  {
    param = Params[i];
    u1 = param[0];
    u2 = param[1];

    u = std::max(std::abs(u1), std::abs(u2));

    if(u>max) max = u;
  }

  for(i=0;i<N_Points;i++)
    Params[i][2] = max;
}

void linfb(int N_Points, double **Coeffs, double **, TBaseCell *)
{
  int i;
  double max, *coeff, b1, b2, b;

  max = -1;

  for(i=0;i<N_Points;i++)
  {
    coeff = Coeffs[i];
    b1 = coeff[1];
    b2 = coeff[2];

    b = std::max(std::abs(b1), std::abs(b2));
    if(b>max) max = b;
  }

  for(i=0;i<N_Points;i++)
    Coeffs[i][5] = max;
}
void ave_l2b_quad_points(int N_Points, double **Coeffs, double **, TBaseCell *)
{
  int i;
  double max, *coeff, b1, b2, b;

  max = 0;

  for(i=0;i<N_Points;i++)
  {
    coeff = Coeffs[i];
    b1 = coeff[1];
    b2 = coeff[2];

    b= std::sqrt(b1*b1+b2*b2); 
    max += b;
  }

  max /= N_Points;

  for(i=0;i<N_Points;i++)
    Coeffs[i][5] = max;
}

#endif // 2D

#ifdef __3D__
void linfb(int N_Points, double **Coeffs, double **, TBaseCell *)
{
  int i;
  double max, *coeff, b1, b2, b3, b;

  max = -1;

  for(i=0;i<N_Points;i++)
  {
    coeff = Coeffs[i];
    b1 = coeff[1];
    b2 = coeff[2];
    b3 = coeff[3];

    b = std::max(std::abs(b1), std::abs(b2));
    b = std::max(std::abs(b3), b);

    if(b>max) max = b;
  }

  for(i=0;i<N_Points;i++)
    Coeffs[i][6] = max;
}

void ave_l2b_quad_points(int N_Points, double **Coeffs, double **, TBaseCell *)
{
  int i;
  double max, *coeff, b1, b2, b3, b;

  max = 0;

  for(i=0;i<N_Points;i++)
  {
    coeff = Coeffs[i];
    b1 = coeff[1];
    b2 = coeff[2];
    b3 = coeff[3];

    b =  std::sqrt(b1*b1+b2*b2+b3*b3);

    max += b;
  }

  max /= N_Points;

  for(i=0;i<N_Points;i++)
    Coeffs[i][6] = max;
}
#endif // 3D


void unknown_solution_3d(double, double, double, double *values)
{
  values[0] =0;
  values[1] =0;
  values[2] =0;
  values[3] =0;
  values[4] =0;
}

void unknown_solution_2d(double, double, double *values)
{
  values[0] =0;
  values[1] =0;
  values[2] =0;
  values[3] =0;
}

void BoundConditionVMM(int, double, BoundCond &cond)
{
   cond = NEUMANN;
}

void BoundConditionNoBoundCondition(int, double, BoundCond &cond)
{
   cond = NEUMANN;
}
void BoundConditionNoBoundCondition(int, double, double, double, BoundCond &cond)
{
   cond = NEUMANN;
}
void BoundaryValueHomogenous(int, double, double &value)
{
  value = 0;
}
void BoundaryValueHomogenous(int, double, double, double, double &value)
{
  value = 0;
}
void BoundaryValueNoBoundaryValue(int, double, double &value)
{
  value = 0;
}

void BoundConditionNSE(int, double, BoundCond &cond)
{
   cond = DIRICHLET;
}

void BoundaryConditionPressSep(int, double, BoundCond &cond)
{
   cond = NEUMANN;
}

void BoundaryValuePressSep(int, double, double &value)
{
  value = 0;
}
void BoundaryConditionPressSep3D(int, double, double, double, BoundCond &cond)
{
  cond = NEUMANN;
}

void BoundaryValuePressSep3D(int, double, double, double, double &value)
{
  value = 0;
}

void BoundaryConditionNewton(int, double, double, double, BoundCond &cond)
{
  cond = DIRICHLET;
}

void BoundaryValueNewton(int, double, double, double, double &value)
{
  value = 0;
}

// boundary condition for assembling in the FEM-FCT scheme
void BoundCondition_FEM_FCT(int, double, BoundCond &cond)
{
    cond = NEUMANN;
}

void BoundValue_FEM_FCT(int, double, double &value)
{
    value = 0;
}
void BoundCondition_FEM_FCT(int, double, double, double, BoundCond &cond)
{
    cond = NEUMANN;
}
void BoundValue_FEM_FCT(int, double, double, double, double &value)
{
  value = 0;
}

// ========================================================================
// boundary values for auxiliary problem in Galdi/Layton model
// ========================================================================

// void BoundConditionAuxProblem(int i, double t, BoundCond &cond)
// {
//   cond = NEUMANN;
//   //cond = DIRICHLET;
// }
// 
// void BoundValueAuxProblem(int BdComp, double Param, double &value)
// {
//   value = 0;
// }

// ========================================================================
// boundary values for higher order fe in VMS
// ========================================================================

// void ho_BoundCondition(int i, double t, BoundCond &cond)
// {
//   cond = DIRICHLET;
// }
// 
// void ho_BoundValue(int BdComp, double Param, double &value)
// {
//   value = 0;
// }


void SetPolynomialDegree()
{
    if ((TDatabase::ParamDB->ANSATZ_ORDER >0)
        && (TDatabase::ParamDB->ANSATZ_ORDER <10))
    {
        TDatabase::ParamDB->INTERNAL_POLYNOMIAL_DEGREE =
            TDatabase::ParamDB->ANSATZ_ORDER;
        return;
    }
    
    if ((TDatabase::ParamDB->ANSATZ_ORDER <0)
        && (TDatabase::ParamDB->ANSATZ_ORDER > -10))
    {
        TDatabase::ParamDB->INTERNAL_POLYNOMIAL_DEGREE =
            -TDatabase::ParamDB->ANSATZ_ORDER;
        return;
    }
    if (TDatabase::ParamDB->ANSATZ_ORDER == -101)
    {
        TDatabase::ParamDB->INTERNAL_POLYNOMIAL_DEGREE = 1;
        return;
    }

    // Discontious Galerkin Methods
    if ((TDatabase::ParamDB->ANSATZ_ORDER <-10)
        && (TDatabase::ParamDB->ANSATZ_ORDER > -15))
    {
        TDatabase::ParamDB->INTERNAL_POLYNOMIAL_DEGREE =
            -TDatabase::ParamDB->ANSATZ_ORDER - 10;
        return;
    }
    if ((TDatabase::ParamDB->ANSATZ_ORDER <-100)
        && (TDatabase::ParamDB->ANSATZ_ORDER > -150))
    {
        TDatabase::ParamDB->INTERNAL_POLYNOMIAL_DEGREE =
            (-TDatabase::ParamDB->ANSATZ_ORDER - 100)/10;
        return;
    }

    //==========LOCALPROJECTION=================
    if (TDatabase::ParamDB->ANSATZ_ORDER == 100)
    {
        TDatabase::ParamDB->INTERNAL_POLYNOMIAL_DEGREE = 1;
        return;
    }
    if (TDatabase::ParamDB->ANSATZ_ORDER == 201
      || TDatabase::ParamDB->ANSATZ_ORDER == 200
      || TDatabase::ParamDB->ANSATZ_ORDER == 211
      || TDatabase::ParamDB->ANSATZ_ORDER == 221
      || TDatabase::ParamDB->ANSATZ_ORDER == 222)
    {
        TDatabase::ParamDB->INTERNAL_POLYNOMIAL_DEGREE = 2;
        return;
    }
    if (TDatabase::ParamDB->ANSATZ_ORDER == 302
      || TDatabase::ParamDB->ANSATZ_ORDER == 301
      || TDatabase::ParamDB->ANSATZ_ORDER == 312
      || TDatabase::ParamDB->ANSATZ_ORDER == 322)
    {
        TDatabase::ParamDB->INTERNAL_POLYNOMIAL_DEGREE = 3;
        return;
    }
    if (TDatabase::ParamDB->ANSATZ_ORDER == 403
      || TDatabase::ParamDB->ANSATZ_ORDER == 402
      || TDatabase::ParamDB->ANSATZ_ORDER == 413
      || TDatabase::ParamDB->ANSATZ_ORDER == 423)
    {
        TDatabase::ParamDB->INTERNAL_POLYNOMIAL_DEGREE = 4;
        return;
    }
    if (TDatabase::ParamDB->ANSATZ_ORDER == 504
      || TDatabase::ParamDB->ANSATZ_ORDER == 503
      || TDatabase::ParamDB->ANSATZ_ORDER == 514
      || TDatabase::ParamDB->ANSATZ_ORDER == 524)
    {
        TDatabase::ParamDB->INTERNAL_POLYNOMIAL_DEGREE = 5;
        return;
    }
    if (TDatabase::ParamDB->ANSATZ_ORDER == 605
      || TDatabase::ParamDB->ANSATZ_ORDER == 604
      || TDatabase::ParamDB->ANSATZ_ORDER == 615
      || TDatabase::ParamDB->ANSATZ_ORDER == 625)
    {
        TDatabase::ParamDB->INTERNAL_POLYNOMIAL_DEGREE = 6;
        return;
    }
    if (TDatabase::ParamDB->ANSATZ_ORDER == 706
      || TDatabase::ParamDB->ANSATZ_ORDER == 705
      || TDatabase::ParamDB->ANSATZ_ORDER == 716
      || TDatabase::ParamDB->ANSATZ_ORDER == 726)
    {
        TDatabase::ParamDB->INTERNAL_POLYNOMIAL_DEGREE = 7;
        return;
    }
    if (TDatabase::ParamDB->ANSATZ_ORDER == 807
      || TDatabase::ParamDB->ANSATZ_ORDER == 806
      || TDatabase::ParamDB->ANSATZ_ORDER == 817
      || TDatabase::ParamDB->ANSATZ_ORDER == 827)
    {
        TDatabase::ParamDB->INTERNAL_POLYNOMIAL_DEGREE = 8;
        return;
    }
    if (TDatabase::ParamDB->ANSATZ_ORDER == 908
      || TDatabase::ParamDB->ANSATZ_ORDER == 907
      || TDatabase::ParamDB->ANSATZ_ORDER == 918
      || TDatabase::ParamDB->ANSATZ_ORDER == 928)
    {
        TDatabase::ParamDB->INTERNAL_POLYNOMIAL_DEGREE = 9;
        return;
    }
    //===========================================
    ErrThrow("INTERNAL_POLYNOMIAL_DEGREE not defined for ANSATZ_ORDER ",
             TDatabase::ParamDB->ANSATZ_ORDER);
}

void CheckMaximumPrinciple(TSquareMatrix *A, double *sol, int N_Active,
                           double *errors)
{
  int i,j,k,index,ii,jj,kk,found;
  double tol=1e-8, minimum, maximum, val, maxerror;
  const int *RowPtr, *KCol;
  int N_DOF, maxprinciple = 0; 

  RowPtr = A->get_row_ptr();
  KCol = A->get_vector_columns();
  //double *Entries = A->GetEntries();
  N_DOF = A->get_n_rows();
  maxerror = 0;

  j = RowPtr[0];
  for(i=0;i<N_DOF;i++)
  {
      k = RowPtr[i+1];
      val = sol[i];
      maximum = -1e10;
      minimum = 1e10;
      if (i<N_Active)
      {
          for(;j<k;j++)
          {
              index = KCol[j];
              if (index==i)
                  continue;
              if (sol[index]>maximum)
                  maximum = sol[index];
              if (sol[index]<minimum)
                  minimum = sol[index];
          }
      }
      else
          // Dirichlet nodes
      {
          jj = RowPtr[0];
          for (ii=0;ii<N_Active;ii++)
          {
              found = 0;
              kk = RowPtr[ii+1];
              for (;jj<kk;jj++)
              {
                  index = KCol[jj];
                  // off diagonal entry in row ii and column i
                  if (index==i)
                  {
                      found++;
                      if (sol[ii]>maximum)
                          maximum = sol[ii];
                      if (sol[ii]<minimum)
                          minimum = sol[ii];
                  }
                  // check for neighbour Dirichlet nodes
                  if (found)
                  {
                      jj = RowPtr[ii];
                      for (;jj<kk;jj++)
                      {
                          index = KCol[jj];
                          // off diagonal entry in row ii and column >= N_Active
                          if ((index>=N_Active)&&(index!=i))
                          {
                              if (sol[index]>maximum)
                                  maximum = sol[index];
                              if (sol[index]<minimum)
                                  minimum = sol[index];
                          }
                          
                      }
                  }
              }
          }
      }   
      // computation of maximum and minimum of neighbours done
      minimum = minimum - tol;
      maximum = maximum + tol;
      if ((minimum<=val) && (maximum>=val))
      {
          maxprinciple++;
          continue;
      }
      if ((maximum>-2) &&  (val-maximum > maxerror))
          maxerror = val-maximum;
      if ((minimum<2) && (minimum-val > maxerror))
          maxerror = minimum-val;
  } // endfor i
  Output::print("# dof ", N_DOF, " max principle fulfilled ", maxprinciple, " ",
                maxprinciple*100.0/N_DOF, "% maxerror ", maxerror);
  errors[0] = maxprinciple*100.0/N_DOF ;
  errors[1] = maxerror;
} // end CheckMaximumPrinciple

//save array of double pointers in a file
void SaveData(char *name, int N_Array, double **sol, int *N_Unknowns)
{
    int i;
    char OldString[] = "_old", MvString[] = "mv ", SpaceString[] = " ", tmp[200];
    
    // move the previous data set
    strcpy(tmp,MvString);
    strcat(tmp,name);
    strcat(tmp,SpaceString);
    strcat(tmp,name);
    strcat(tmp,OldString);
    Output::print(tmp);
    //system(tmp);

  std::ofstream dat(name);

  if(!dat)
  {
    Output::print("CANNOT OPEN FILE '", name, "' FOR SAVING DATA!");
    return;
  }
  for (i=0;i<N_Array;i++)
  {
      Output::print(N_Unknowns[i]);
   	  dat.write((char *)sol[i],sizeof(double)*N_Unknowns[i]);
  }  
  dat.close();
  
  Output::print("wrote output into file: ", name);
}
void ReadData(char *name, int N_Array, double **sol, int *N_Unknowns)
{
    int i;
  std::ifstream dat(name);

  if(!dat)
  {
    ErrThrow("CANNOT OPEN FILE '", name, "' FOR READING DATA!");
  }

  for (i=0;i<N_Array;i++)
  {
   	  dat.read((char *)sol[i],sizeof(double)*N_Unknowns[i]);
	  //SwapDoubleArray(sol[i], N_Unknowns[i]); 
  }  
  dat.close();
  
  Output::print("read input from file: ", name);
}

// save sol into a file
void SaveData(const std::string& basename, double *sol, int nDOF)
{
  std::string filename = basename + "init";
  std::ofstream ofile;
  ofile.open( filename.c_str() , std::ios::out | std::ios::trunc );
  if( !ofile.is_open() )
  {
    ErrThrow("CANNOT OPEN FILE '", filename, "' FOR SAVING DATA!");
  }
  ofile << setprecision( 12 );
  ofile.write((char *)sol,sizeof(double)*nDOF);
  ofile.close();
  Output::print("saving data into file: ", filename);
}

void ReadData(const std::string& filename, double *sol, int nDOF)
{
  std::ifstream ifile(filename.c_str());
  if(!ifile)
  {
    ErrThrow("CANNOT OPEN FILE '", filename, "' FOR READING DATA!");
  }
  ifile.read((char *)sol,sizeof(double)*nDOF);  
  ifile.close();
  Output::print("reading data from file: ", filename);
}


double graddiv_parameterOseen(double hK, double, double, double)
{
  double tau;
  
  switch(TDatabase::ParamDB->DIV_DIV_STAB_TYPE)
  {
    // constant
    case 0:
      tau = TDatabase::ParamDB->DIV_DIV_STAB_C1;
      break;
      // depending on mesh width 
    case 1:
      //  if (nu < hK)
        tau = TDatabase::ParamDB->DIV_DIV_STAB_C1 * std::pow(hK,TDatabase::ParamDB->DIV_DIV_STAB_C2);
      // else
      //  tau = TDatabase::ParamDB->DIV_DIV_STAB_C1;
      break;
    default:
      ErrThrow("DIV_DIV_STAB_TYPE ", TDatabase::ParamDB->DIV_DIV_STAB_TYPE,
               " not implemented");
  }
  return(tau);
}

#ifdef __2D__
void ComputeVorticityDivergence(TFEFunction2D *u1, TFEFunction2D *u2,
                                const TFESpace2D *vorticity_space, 
                                double *vort,  double *div)
{
  int N_loc_dofVort;
  const TBaseCell *cell;
  const double * xi_ref, *eta_ref;
  double X_orig[MaxN_PointsForNodal2D], Y_orig[MaxN_PointsForNodal2D];
  double val[3];
  double PointValuesDiv[MaxN_PointsForNodal2D];
  double FunctionalValuesDiv[MaxN_PointsForNodal2D];
  double PointValuesVort[MaxN_PointsForNodal2D];
  double FunctionalValuesVort[MaxN_PointsForNodal2D];
  int N_Vort = vorticity_space->get_n_dof();
  memset(vort,0,N_Vort*sizeof(double));
  memset(div,0,N_Vort*sizeof(double));
  int *N_Found = new int[N_Vort];
  memset(N_Found,0,N_Vort*sizeof(int));
  auto Coll = vorticity_space->GetCollection();

  for(int i=0;i<Coll->GetN_Cells();i++)
  {
    cell = Coll->GetCell(i);        
    auto Element = vorticity_space->get_fe(i);
    auto RefTrans = Element.GetRefTransID();
    FEDatabase::SetCellForRefTrans(cell,RefTrans);
    auto nf = Element.GetNodalFunctional();
    nf->GetPointsForAll(N_loc_dofVort, xi_ref, eta_ref);
    FEDatabase::GetOrigFromRef(RefTrans, N_loc_dofVort, xi_ref, eta_ref,
                               X_orig, Y_orig);
    auto DOF = vorticity_space->GetGlobalDOF(i);
    memset(PointValuesDiv, 0, MaxN_PointsForNodal2D*sizeof(double));
    memset(PointValuesVort, 0, MaxN_PointsForNodal2D*sizeof(double));
    for (int j=0;j<N_loc_dofVort;j++)
    {
      u1->FindGradientLocal(cell,i,X_orig[j],Y_orig[j],val);
      PointValuesVort[j] -= val[2];
      PointValuesDiv[j] += val[1];

      u2->FindGradientLocal(cell,i,X_orig[j],Y_orig[j],val); 
      PointValuesVort[j] += val[1];
      PointValuesDiv[j] += val[2];
    }    
    nf->GetAllFunctionals(Coll, cell, PointValuesVort,
                          FunctionalValuesVort);
    nf->GetAllFunctionals(Coll, cell, PointValuesDiv,
                          FunctionalValuesDiv);
    int N_LocalDOF = Element.GetN_DOF();
    for(int j=0;j<N_LocalDOF;j++)
    {
      int index = DOF[j];
      vort[index] += FunctionalValuesVort[j];
      div[index] += FunctionalValuesDiv[j];
      N_Found[index]++;
    }
  } 
  for (int i=0;i<N_Vort;i++)
  {
    vort[i]/=(double)N_Found[i];
    div[i]/=(double)N_Found[i];
  }
  delete [] N_Found;
}
#endif
