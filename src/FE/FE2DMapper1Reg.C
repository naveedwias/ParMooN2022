#include <FE2DMapper1Reg.h>
#include <HNDesc.h>
#include <MooNMD_Io.h>

FEMapper1Reg2D::FEMapper1Reg2D(FEMapper_type t)
 : FEMapper()
{
  this->FEMapper::type = t;
  switch(this->FEMapper::type)
  {
    case FEMapper_type::C0_2_C0_2_2D:
    case FEMapper_type::C1_2_C1_2_2D:
    case FEMapper_type::C2_2_C2_2_2D:
    case FEMapper_type::C3_2_C3_2_2D:
    case FEMapper_type::C4_2_C4_2_2D:
    case FEMapper_type::C5_2_C5_2_2D:
    case FEMapper_type::C6_2_C6_2_2D:
    case FEMapper_type::C7_2_C7_2_2D:
    case FEMapper_type::C8_2_C8_2_2D:
    case FEMapper_type::C9_2_C9_2_2D:
    case FEMapper_type::N1_2_N1_2_2D:
    case FEMapper_type::N2_2_N2_2_2D:
    case FEMapper_type::N3_2_N3_2_2D:
    case FEMapper_type::N4_2_N4_2_2D:
    case FEMapper_type::N5_2_N5_2_2D:
      ErrThrow("you tried to construct an FEMapper1Reg2D of type ", type,
               " which is used for regular grids in 2D only. For that you need "
               "to use the base class FEMapper");
      break;
    // ONE regular grid, same pattern on both sides
    case FEMapper_type::C0_2_C0_2_1Reg_2D:
      name = "C0_2_C0_2_1Reg";
      description = "conforming P0 or Q0 element, one regular grid";
      n_dof0 = 0;
      n_dof1 = 0;
      N_DOF2 = 0;
      //aux.resize(0);
      N_Mid = 0;
      n_pairs = 0;
      N_NoOpposite = 0;
      N_Hanging = 0;
      break;
    case FEMapper_type::C1_2_C1_2_1Reg_2D:
      name = "C1_2_C1_2_1Reg";
      description = "conforming P1 or Q1 element, one regular grid";
      n_dof0 = 2;
      n_dof1 = 2;
      N_DOF2 = 2;
      aux.resize(6);
      N_Mid = 1;
      Mid = {3,4};
      n_pairs = 2;
      pairs.push_back({ {0,5}, {1,2} });
      N_NoOpposite = 0;
      N_Hanging = 1;
      Hanging = { 3 };
      HangingTypes = { HNDesc::HN_C_P1_2D_0 };
      coupling = { {0, 1} };
      break;
    case FEMapper_type::C2_2_C2_2_1Reg_2D:
      name = "C2_2_C2_2_1Reg";
      description = "conforming P2 or Q2 element, one regular grid";
      n_dof0 = 3;
      n_dof1 = 3;
      N_DOF2 = 3;
      aux.resize(9);
      N_Mid = 1;
      Mid = {5,6};
      n_pairs = 3;
      pairs.push_back({ {0,8}, {1,5}, {2,3} });
      N_NoOpposite = 0;
      N_Hanging = 2;
      Hanging = { 4, 7 };
      HangingTypes = { HNDesc::HN_C_P2_2D_0, HNDesc::HN_C_P2_2D_1 };
      coupling = { { 0, 1, 2 }, { 0, 1, 2 } };
      break;
    case FEMapper_type::C3_2_C3_2_1Reg_2D:
      name = "C3_2_C3_2_1Reg";
      description = "conforming P3 or Q3 element, one regular grid";
      n_dof0 = 4;
      n_dof1 = 4;
      N_DOF2 = 4;
      aux.resize(12);
      N_Mid = 1;
      Mid = {7,8};
      n_pairs = 4;
      pairs.push_back({ {0,11}, {1,9}, {2,6}, {3,4} });
      N_NoOpposite = 0;
      N_Hanging = 3;
      Hanging = { 5, 7, 10 };
      HangingTypes = { HNDesc::HN_C_P3_2D_0, HNDesc::HN_C_P3_2D_1,
                       HNDesc::HN_C_P3_2D_2 };
      coupling = { { 0, 1, 2, 3 }, { 0, 1, 2, 3 }, { 0, 1, 2, 3 } };
      break;
    case FEMapper_type::C4_2_C4_2_1Reg_2D:
      name = "C4_2_C4_2_1Reg";
      description = "conforming P4 or Q4 element, one regular grid";
      n_dof0 = 5;
      n_dof1 = 5;
      N_DOF2 = 5;
      aux.resize(15);
      N_Mid = 1;
      Mid = {9,10};
      n_pairs = 5;
      pairs.push_back({ {0,14}, {1,12}, {2,9}, {3,7}, {4,5} });
      N_NoOpposite = 0;
      N_Hanging = 4;
      Hanging = { 6, 8, 11, 13 };
      HangingTypes = { HNDesc::HN_C_P4_2D_0, HNDesc::HN_C_P4_2D_1,
                       HNDesc::HN_C_P4_2D_2, HNDesc::HN_C_P4_2D_3 };
      coupling = { { 0, 1, 2, 3, 4 }, { 0, 1, 2, 3, 4 }, { 0, 1, 2, 3, 4 },
                   { 0, 1, 2, 3, 4 } };
      break;
    case FEMapper_type::C5_2_C5_2_1Reg_2D:
      name = "C5_2_C5_2_1Reg";
      description = "conforming P5 or Q5 element, one regular grid";
      n_dof0 = 6;
      n_dof1 = 6;
      N_DOF2 = 6;
      aux.resize(18);
      N_Mid = 1;
      Mid = {11,12};
      n_pairs = 6;
      pairs.push_back({ {0,17}, {1,15}, {2,13}, {3,10}, {4,8}, {5,6} });
      N_NoOpposite = 0;
      N_Hanging = 5;
      Hanging = { 7, 9, 11, 14, 16 };
      HangingTypes = { HNDesc::HN_C_P5_2D_0, HNDesc::HN_C_P5_2D_1,
                       HNDesc::HN_C_P5_2D_2, HNDesc::HN_C_P5_2D_3,
                       HNDesc::HN_C_P5_2D_4 };
      coupling = { { 0, 1, 2, 3, 4, 5 }, { 0, 1, 2, 3, 4, 5 },
                   { 0, 1, 2, 3, 4, 5 }, { 0, 1, 2, 3, 4, 5 },
                   { 0, 1, 2, 3, 4, 5 } };
      break;
    case FEMapper_type::N1_2_N1_2_1Reg_2D:
      name = "N1_2_N1_2_1Reg";
      description = "nonconforming P1 or Q1 element, one regular grid";
      n_dof0 = 1;
      n_dof1 = 1;
      N_DOF2 = 1;
      aux.resize(3);
      N_Mid = 0;
      n_pairs = 0;
      N_NoOpposite = 2;
      NoOpposite = { 1, 2 };
      N_Hanging = 1;
      Hanging = { 0 };
      HangingTypes = { HNDesc::HN_N_P1_2D_0 };
      coupling = { { 1, 2 } };
      break;
    case FEMapper_type::N2_2_N2_2_1Reg_2D:
      name = "N2_2_N2_2_1Reg";
      description = "nonconforming P2 or Q2 element, one regular grid";
      n_dof0 = 2;
      n_dof1 = 2;
      N_DOF2 = 2;
      aux.resize(6);
      N_Mid = 0;
      n_pairs = 0;
      N_NoOpposite = 4;
      NoOpposite = { 2, 3, 4, 5 };
      N_Hanging = 2;
      Hanging = { 0, 1 };
      HangingTypes = { HNDesc::HN_N_P1_2D_0, HNDesc::HN_N_P2_2D_0 };
      coupling = { { 2, 4, }, { 2, 3, 4, 5 } };
      break;
    case FEMapper_type::N3_2_N3_2_1Reg_2D:
      name = "N3_2_N3_2_1Reg";
      description = "nonconforming P3 or Q3 element, one regular grid";
      n_dof0 = 3;
      n_dof1 = 3;
      N_DOF2 = 3;
      aux.resize(9);
      N_Mid = 0;
      n_pairs = 0;
      N_NoOpposite = 6;
      NoOpposite = { 3, 4, 5, 6, 7, 8 };
      N_Hanging = 3;
      Hanging = { 0, 1, 2 };
      HangingTypes = { HNDesc::HN_N_P1_2D_0, HNDesc::HN_N_P2_2D_0,
                       HNDesc::HN_N_P3_2D_0 };
      coupling = { { 3, 6 }, { 3, 4, 6, 7 }, { 3, 4, 5, 6, 7, 8 } };
      break;
    case FEMapper_type::N4_2_N4_2_1Reg_2D:
      name = "N4_2_N4_2_1Reg";
      description = "nonconforming P4 or Q4 element, one regular grid";
      n_dof0 = 4;
      n_dof1 = 4;
      N_DOF2 = 4;
      aux.resize(12);
      N_Mid = 0;
      n_pairs = 0;
      N_NoOpposite = 8;
      NoOpposite = { 4, 5, 6, 7, 8, 9, 10, 11 };
      N_Hanging = 4;
      Hanging = { 0, 1, 2, 3 };
      HangingTypes = { HNDesc::HN_N_P1_2D_0, HNDesc::HN_N_P2_2D_0,
                       HNDesc::HN_N_P3_2D_0, HNDesc::HN_N_P4_2D_0 };
      coupling = { { 4, 8 }, { 4, 5, 8, 9 }, { 4, 5, 6, 8, 9, 10 },
                   { 4, 5, 6, 7, 8, 9, 10, 11 } };
      break;
    case FEMapper_type::N5_2_N5_2_1Reg_2D:
      name = "N5_2_N5_2_1Reg";
      description = "nonconforming P5 or Q5 element, one regular grid";
      n_dof0 = 5;
      n_dof1 = 5;
      N_DOF2 = 5;
      aux.resize(15);
      N_Mid = 0;
      n_pairs = 0;
      N_NoOpposite = 10;
      NoOpposite = { 5, 6, 7, 8, 9, 10, 11, 12, 13, 14 };
      N_Hanging = 5;
      Hanging = { 0, 1, 2, 3, 4  };
      HangingTypes = { HNDesc::HN_N_P1_2D_0, HNDesc::HN_N_P2_2D_0,
                       HNDesc::HN_N_P3_2D_0, HNDesc::HN_N_P4_2D_0,
                       HNDesc::HN_N_P5_2D_0 };
      coupling = { { 5, 10 }, { 5, 6, 10, 11 }, { 5, 6, 7, 10, 11, 12 },
                   { 5, 6, 7, 8, 10, 11, 12, 13 },
                   { 5, 6, 7, 8, 9, 10, 11, 12, 13, 14 } };
      break;
    // ONE regular grid, different pattern
    case FEMapper_type::N2_2_N1_2_1Reg_2D: // coarse: second order, fine: first order
      name = "N2_2_N1_2_1Reg";
      description = "nonconforming P2 or Q2 element, one regular grid";
      n_dof0 = 2;
      n_dof1 = 1;
      N_DOF2 = 1;
      aux.resize(4);
      N_Mid = 0;
      n_pairs = 0;
      N_NoOpposite = 3;
      NoOpposite = { 1, 2, 3};
      N_Hanging = 1;
      Hanging = { 0 };
      HangingTypes = { HNDesc::HN_N_P1_2D_0 };
      coupling = { { 2, 3 } };
      break;
    case FEMapper_type::N3_2_N2_2_1Reg_2D: // coarse: third order,  fine: second order
      name = "N3_2_N2_2_1Reg";
      description = "nonconforming P3 or Q3 element, one regular grid";
      n_dof0 = 3;
      n_dof1 = 2;
      N_DOF2 = 2;
      aux.resize(7);
      N_Mid = 0;
      n_pairs = 0;
      N_NoOpposite = 5;
      NoOpposite = { 2, 3, 4, 5, 6 };
      N_Hanging = 2;
      Hanging = { 0, 1 };
      HangingTypes = { HNDesc::HN_N_P1_2D_0, HNDesc::HN_N_P2_2D_0 };
      coupling = { { 3, 5 }, { 3, 4, 5, 6 } };
      break;
    case FEMapper_type::N4_2_N3_2_1Reg_2D: // coarse: fourth order, fine: third order
      name = "N4_2_N3_2_1Reg";
      description = "nonconforming P4 or Q4 element, one regular grid";
      n_dof0 = 4;
      n_dof1 = 3;
      N_DOF2 = 3;
      aux.resize(10);
      N_Mid = 0;
      n_pairs = 0;
      N_NoOpposite = 7;
      NoOpposite = { 3, 4, 5, 6, 7, 8, 9 };
      N_Hanging = 3;
      Hanging = { 0, 1, 2 };
      HangingTypes = { HNDesc::HN_N_P1_2D_0, HNDesc::HN_N_P2_2D_0,
                       HNDesc::HN_N_P3_2D_0 };
      coupling = { { 4, 7 }, { 4, 5, 7, 8 }, { 4, 5, 6, 7, 8, 9 } };
      break;
    case FEMapper_type::N5_2_N4_2_1Reg_2D: // coarse: fifth order,  fine: fourth order
      name = "N5_2_N5_2_1Reg";
      description = "nonconforming P5 or Q5 element, one regular grid";
      n_dof0 = 5;
      n_dof1 = 4;
      N_DOF2 = 4;
      aux.resize(13);
      N_Mid = 0;
      n_pairs = 0;
      N_NoOpposite = 9;
      NoOpposite = { 4, 5, 6, 7, 8, 9, 10, 11, 12 };
      N_Hanging = 4;
      Hanging = { 0, 1, 2, 3 };
      HangingTypes = { HNDesc::HN_N_P1_2D_0, HNDesc::HN_N_P2_2D_0,
                       HNDesc::HN_N_P3_2D_0, HNDesc::HN_N_P4_2D_0 };
      coupling = { { 5, 9 }, { 5, 6, 9, 10 }, { 5, 6, 7, 9, 10, 11 },
                   { 5, 6, 7, 8, 9, 10, 11, 12 } };
      break;
    // 3D mapper types
    case FEMapper_type::D_D_3D:
    case FEMapper_type::P1_P1_3D:
    case FEMapper_type::P2_P2_3D:
    case FEMapper_type::P3_P3_3D:
    case FEMapper_type::P2B_P2B_3D:
    case FEMapper_type::Q1_Q1_3D:
    case FEMapper_type::Q2_Q2_3D:
    case FEMapper_type::Q3_Q3_3D:
    case FEMapper_type::Q4_Q4_3D:
    case FEMapper_type::N1_N1_3D:
    case FEMapper_type::N2_N2_3D:
    case FEMapper_type::N3_N3_3D:
    case FEMapper_type::N4_N4_3D:
    case FEMapper_type::NP1_NP1_3D:
    case FEMapper_type::NP2_NP2_3D:
      ErrThrow("you tried to construct an FEMapper1Reg2D of type ", type,
               " which is used regular grids in 3D only. For that you need to "
               "use the base class FEMapper");
      break;
    // 1-regular grid, same pattern on both sides
    case FEMapper_type::P1_P1_1Reg_3D:
    case FEMapper_type::P2_P2_1Reg_3D:
    case FEMapper_type::P3_P3_1Reg_3D:
    case FEMapper_type::Q1_Q1_1Reg_3D:
    case FEMapper_type::Q2_Q2_1Reg_3D:
    case FEMapper_type::Q3_Q3_1Reg_3D:
    case FEMapper_type::NP1_NP1_1Reg_3D:
    case FEMapper_type::NP2_NP2_1Reg_3D:
      ErrThrow("you tried to construct an FEMapper1Reg2D of type ", type,
               " which is used for one-regular grids in 3D only. For that you "
               "need to use the derived class FEMapper1Reg3D");
      break;
    default:
      ErrThrow("unknown fe mapper type ", t);
  }
}

/** map the given local degrees of freedom,
    coarse cell has lower number,
    0 is coarser side 0; 1,2 are on finer side 1 */
void FEMapper1Reg2D::Map_1Reg(
  bool coarse_fine, int *Global, int I_K0, int I_K1, int I_K2,
  const int *Indices0, const int *Indices1, const int *Indices2,
  int &Counter, int LowerFirstChild, std::vector<THangingNode*>& vect,
  std::vector<int>& numbers, std::map<HNDesc, THNDesc *>& hn_descriptors) const
{
  int i, v, w;

  for(i=0;i<n_dof0;i++)
  {
    aux[i]=I_K0+Indices0[i];
  }

  for(i=0;i<n_dof1;i++)
  {
    aux[i+n_dof0]=I_K1+Indices1[i];
  }

  for(i=0;i<N_DOF2;i++)
  {
    aux[i+n_dof0+n_dof1]=I_K2+Indices2[i];
  }

  if(!pairs.empty())
  {
    for(auto pair : pairs[0])
    {
      if(pair.first!=-1)
      {
        if(coarse_fine)
          map_single_dof(Global, aux[pair.first], aux[pair.second], Counter);
        else
          map_single_dof(Global, aux[pair.second], aux[pair.first], Counter);
        // cout << pair.first << "  " << pair.second << endl;
      }
    }
  }

  for(i=0;i<N_Mid;i++)
  {
    if(LowerFirstChild)
      map_single_dof(Global, aux[Mid[i*2]], aux[Mid[i*2+1]], Counter);
    else
      map_single_dof(Global, aux[Mid[i*2+1]], aux[Mid[i*2]], Counter);
  }

  for(i=0;i<N_Hanging;i++)
  {
    w=aux[Hanging[i]];
    while( (v=Global[w]) > -1)
    {
      w=v;
    }
    Global[w]=HANGINGNODE;
    auto hn_desc = get_hn_descriptor(HangingTypes[i], hn_descriptors);
    std::vector<int> partner_dofs;
    partner_dofs.reserve(N_Hanging-1);
    for(int k = 0; k < N_Hanging; ++k)
    {
      if(k != i)
      {
        partner_dofs.push_back(aux[Hanging[k]]);
      }
    }
    auto hn = new THangingNode(hn_desc, aux.data(), coupling[i].data(),
                               partner_dofs);
    vect.push_back(hn);
    numbers.push_back(aux[Hanging[i]]);
  }

  for(i=0;i<N_NoOpposite;i++)
  {
    w=aux[NoOpposite[i]];
    while( (v=Global[w]) > -1 )
    {
      w=v;
    }

    if( Global[w] == -1 )
    {
      Counter--;
      Global[aux[NoOpposite[i]]]=Counter;
    }
  } // endfor i
}

THNDesc * FEMapper1Reg2D::get_hn_descriptor(
  HNDesc type, std::map<HNDesc, THNDesc *>& hn_descriptors) const
{
  auto it = hn_descriptors.find(type);
  if(it == hn_descriptors.end())
  {
    it = hn_descriptors.insert(std::pair<HNDesc, THNDesc *>(
                                   type, new THNDesc(type))).first;
  }
  return it->second;
}
