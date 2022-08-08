#include "FEMapper.h"
#include "MooNMD_Io.h"

FEMapper::FEMapper(FEMapper_type t)
 : type(t)
{
  switch(type)
  {
    case FEMapper_type::C0_2_C0_2_2D:
      name = "C0_2_C0_2";
      description = "conforming P0 or Q0 element";
      n_dof0 = 0;
      n_dof1 = 0;
      n_pairs = 0;
      break;
    case FEMapper_type::C1_2_C1_2_2D:
      name = "C1_2_C1_2";
      description = "conforming P1 or Q1 element";
      n_dof0 = 2;
      n_dof1 = 2;
      n_pairs = 2;
      pairs.push_back({ {0,3}, {1,2} });
      aux.resize(4);
      break;
    case FEMapper_type::C2_2_C2_2_2D:
      name = "C2_2_C2_2";
      description = "conforming P2 or Q2 element";
      n_dof0 = 3;
      n_dof1 = 3;
      n_pairs = 3;
      pairs.push_back({ {0,5}, {1,4}, {2,3} });
      aux.resize(6);
      break;
    case FEMapper_type::C3_2_C3_2_2D:
      name = "C3_2_C3_2";
      description = "conforming P3 or Q3 element";
      n_dof0 = 4;
      n_dof1 = 4;
      n_pairs = 4;
      pairs.push_back({ {0,7}, {1,6}, {2,5}, {3,4} });
      aux.resize(8);
      break;
    case FEMapper_type::C4_2_C4_2_2D:
      name = "C4_2_C4_2";
      description = "conforming P4 or Q4 element";
      n_dof0 = 5;
      n_dof1 = 5;
      n_pairs = 5;
      pairs.push_back({ {0,9}, {1,8}, {2,7}, {3,6}, {4,5} });
      aux.resize(10);
      break;
    case FEMapper_type::C5_2_C5_2_2D:
      name = "C5_2_C5_2";
      description = "conforming P5 or Q5 element";
      n_dof0 = 6;
      n_dof1 = 6;
      n_pairs = 6;
      pairs.push_back({ {0,11}, {1,10}, {2,9}, {3,8}, {4,7}, {5,6} });
      aux.resize(12);
      break;
    case FEMapper_type::C6_2_C6_2_2D:
      name = "C6_2_C6_2";
      description = "conforming P6 or Q6 element";
      n_dof0 = 7;
      n_dof1 = 7;
      n_pairs = 7;
      pairs.push_back({ {0,13}, {1,12}, {2,11}, {3,10}, {4,9}, {5,8}, {6,7} });
      aux.resize(14);
      break;
    case FEMapper_type::C7_2_C7_2_2D:
      name = "C7_2_C7_2";
      description = "conforming P7 or Q7 element";
      n_dof0 = 8;
      n_dof1 = 8;
      n_pairs = 8;
      pairs.push_back({ {0,15}, {1,14}, {2,13}, {3,12}, {4,11}, {5,10}, {6,9},
                        {7,8} });
      aux.resize(16);
      break;
    case FEMapper_type::C8_2_C8_2_2D:
      name = "C8_2_C8_2";
      description = "conforming P8 or Q8 element";
      n_dof0 = 9;
      n_dof1 = 9;
      n_pairs = 9;
      pairs.push_back({ {0,17}, {1,16}, {2,15}, {3,14}, {4,13}, {5,12}, {6,11},
                        {7,10}, {8,9} });
      aux.resize(18);
      break;
    case FEMapper_type::C9_2_C9_2_2D:
      name = "C9_2_C9_2";
      description = "conforming P9 or Q9 element";
      n_dof0 = 10;
      n_dof1 = 10;
      n_pairs = 10;
      pairs.push_back({ {0,19}, {1,18}, {2,17}, {3,16}, {4,15}, {5,14}, {6,13},
                        {7,12}, {8,11}, {9,10} });
      aux.resize(20);
      break;
    case FEMapper_type::N1_2_N1_2_2D:
      name = "N1_2_N1_2";
      description = "nonconforming P1 or Q1 element";
      n_dof0 = 1;
      n_dof1 = 1;
      n_pairs = 1;
      pairs.push_back({ {0,1} });
      aux.resize(2);
      break;
    case FEMapper_type::N2_2_N2_2_2D:
      name = "N2_2_N2_2";
      description = "nonconforming Q2 element";
      n_dof0 = 2;
      n_dof1 = 2;
      n_pairs = 2;
      pairs.push_back({ {0,2}, {1,3} });
      aux.resize(4);
      break;
    case FEMapper_type::N3_2_N3_2_2D:
      name = "N3_2_N3_2";
      description = "nonconforming Q3 element";
      n_dof0 = 3;
      n_dof1 = 3;
      n_pairs = 3;
      pairs.push_back({ {0,3}, {1,4}, {2,5} });
      aux.resize(6);
      break;
    case FEMapper_type::N4_2_N4_2_2D:
      name = "N4_2_N4_2";
      description = "nonconforming Q4 element";
      n_dof0 = 4;
      n_dof1 = 4;
      n_pairs = 4;
      pairs.push_back({ {0,4}, {1,5}, {2,6}, {3,7} });
      aux.resize(8);
      break;
    case FEMapper_type::N5_2_N5_2_2D:
      name = "N5_2_N5_2";
      description = "nonconforming Q5 element";
      n_dof0 = 5;
      n_dof1 = 5;
      n_pairs = 5;
      pairs.push_back({ {0,5}, {1,6}, {2,7}, {3,8}, {4,9} });
      aux.resize(10);
      break;
    // ONE regular grid, same pattern on both sides
    case FEMapper_type::C0_2_C0_2_1Reg_2D:
    case FEMapper_type::C1_2_C1_2_1Reg_2D:
    case FEMapper_type::C2_2_C2_2_1Reg_2D:
    case FEMapper_type::C3_2_C3_2_1Reg_2D:
    case FEMapper_type::C4_2_C4_2_1Reg_2D:
    case FEMapper_type::C5_2_C5_2_1Reg_2D:
    case FEMapper_type::N1_2_N1_2_1Reg_2D:
    case FEMapper_type::N2_2_N2_2_1Reg_2D:
    case FEMapper_type::N3_2_N3_2_1Reg_2D:
    case FEMapper_type::N4_2_N4_2_1Reg_2D:
    case FEMapper_type::N5_2_N5_2_1Reg_2D:
    // ONE regular grid, different pattern
    case FEMapper_type::N2_2_N1_2_1Reg_2D: // coarse: second order, fine: first order
    case FEMapper_type::N3_2_N2_2_1Reg_2D: // coarse: third order,  fine: second order
    case FEMapper_type::N4_2_N3_2_1Reg_2D: // coarse: fourth order, fine: third order
    case FEMapper_type::N5_2_N4_2_1Reg_2D: // coarse: fifth order,  fine: fourth order
      ErrThrow("you tried to construct an FEMapper of type ", type, " which is "
               "used for one-regular grids only. For that you need to use the "
               "derived class FEMapper1Reg2D");
      break;
    // 3D mapper types
    case FEMapper_type::D_D_3D:
      name = "D_D";
      description = "discontinuous elements";
      n_dof0 = 0;
      n_dof1 = 0;
      n_pairs = 0;
      break;
    case FEMapper_type::P1_P1_3D:
      name = "P1_P1";
      description = "conforming P1 element";
      n_dof0 = 3;
      n_dof1 = 3;
      n_pairs = 3;
      pairs.push_back({ {0,3}, {1,5}, {2,4} });
      pairs.push_back({ {0,4}, {1,3}, {2,5} });
      pairs.push_back({ {0,5}, {1,4}, {2,3} });
      aux.resize(6);
      break;
    case FEMapper_type::P2_P2_3D:
      name = "P2_P2";
      description = "conforming P2 element";
      n_dof0 = 6;
      n_dof1 = 6;
      n_pairs = 6;
      pairs.push_back({ {0,6}, {1,9}, {2,11}, {3,7}, {4,10}, {5,8} });
      pairs.push_back({ {0,8}, {1,7}, {2,6}, {3,10}, {4,9}, {5,11} });
      pairs.push_back({ {0,11}, {1,10}, {2,8}, {3,9}, {4,7}, {5,6} });
      aux.resize(12);
      break;
    case FEMapper_type::P3_P3_3D:
      name = "P3_P3";
      description = "conforming P3 element";
      n_dof0 = 10;
      n_dof1 = 10;
      n_pairs = 10;
      pairs.push_back({ {0,10}, {1,14}, {2,17}, {3,19}, {4,11}, {5,15}, {6,18},
                        {7,12}, {8,16}, {9,13} });
      pairs.push_back({ {0,13}, {1,12}, {2,11}, {3,10}, {4,16}, {5,15}, {6,14},
                        {7,18}, {8,17}, {9,19} });
      pairs.push_back({ {0,19}, {1,18}, {2,16}, {3,13}, {4,17}, {5,15}, {6,12},
                        {7,14}, {8,11}, {9,10} });
      aux.resize(20);
      break;
    case FEMapper_type::P2B_P2B_3D:
      name = "P2B_P2B";
      description = "conforming P2 element with face and cell bubbles";
      n_dof0 = 7;
      n_dof1 = 7;
      n_pairs = 7;
      pairs.push_back({ {0,7}, {1,10}, {2,12}, {3,8}, {4,11}, {5,9}, {6,13} });
      pairs.push_back({ {0,9}, {1,8}, {2,7}, {3,11}, {4,10}, {5,12}, {6,13} });
      pairs.push_back({ {0,12}, {1,11}, {2,9}, {3,10}, {4,8}, {5,7}, {6,13} });
      aux.resize(14);
      break;
    case FEMapper_type::Q1_Q1_3D:
      name = "Q1_Q1";
      description = "conforming Q1 element";
      n_dof0 = 4;
      n_dof1 = 4;
      n_pairs = 4;
      pairs.push_back({ {0,4}, {1,6}, {2,5}, {3,7} });
      pairs.push_back({ {0,5}, {1,4}, {2,7}, {3,6} });
      pairs.push_back({ {0,7}, {1,5}, {2,6}, {3,4} });
      pairs.push_back({ {0,6}, {1,7}, {2,4}, {3,5} });
      aux.resize(8);
      break;
    case FEMapper_type::Q2_Q2_3D:
      name = "Q2_Q2";
      description = "conforming Q2 element";
      n_dof0 = 9;
      n_dof1 = 9;
      n_pairs = 9;
      pairs.push_back({ {0,9}, {1,12}, {2,15}, {3,10}, {4,13}, {5,16}, {6,11},
                        {7,14}, {8,17} });
      pairs.push_back({ {0,11}, {1,10}, {2,9}, {3,14}, {4,13}, {5,12}, {6,17},
                        {7,16}, {8,15} });
      pairs.push_back({ {0,17}, {1,14}, {2,11}, {3,16}, {4,13}, {5,10}, {6,15},
                        {7,12}, {8, 9} });
      pairs.push_back({ {0,15}, {1,16}, {2,17}, {3,12}, {4,13}, {5,14}, {6,9},
                        {7,10}, {8,11} });
      aux.resize(18);
      break;
    case FEMapper_type::Q3_Q3_3D:
      name = "Q3_Q3";
      description = "conforming Q3 element";
      n_dof0 = 16;
      n_dof1 = 16;
      n_pairs = 16;
      pairs.push_back({ {0,16}, {1,20}, {2,24}, {3,28}, {4,17}, {5,21}, {6,25},
                        {7,29}, {8,18}, {9,22}, {10,26}, {11,30}, {12,19},
                        {13,23}, {14,27}, {15,31} });
      pairs.push_back({ {0,19}, {1,18}, {2,17}, {3,16}, {4,23}, {5,22}, {6,21},
                        {7,20}, {8,27}, {9,26}, {10,25}, {11,24}, {12,31},
                        {13,30}, {14,29}, {15,28} });
      pairs.push_back({ {0,31}, {1,27}, {2,23}, {3,19}, {4,30}, {5,26}, {6,22},
                        {7,18}, {8,29}, {9,25}, {10,21}, {11,17}, {12,28},
                        {13,24}, {14,20}, {15,16} });
      pairs.push_back({ {0,28}, {1,29}, {2,30}, {3,31}, {4,24}, {5,25}, {6,26},
                        {7,27}, {8,20}, {9,21}, {10,22}, {11,23}, {12,16},
                        {13,17}, {14,18}, {15,19} });
      aux.resize(32);
      break;
    case FEMapper_type::Q4_Q4_3D:
      name = "Q4_Q4";
      description = "conforming Q4 element";
      n_dof0 = 25;
      n_dof1 = 25;
      n_pairs = 25;
      pairs.push_back({ {0,25}, {1,30}, {2,35}, {3,40}, {4,45}, {5,26}, {6,31},
                        {7,36}, {8,41}, {9,46}, {10,27}, {11,32}, {12,37},
                        {13,42}, {14,47}, {15,28}, {16,33}, {17,38}, {18,43},
                        {19,48}, {20,29}, {21,34}, {22,39}, {23,44}, {24,49} });
      pairs.push_back({ {0,29}, {1,28}, {2,27}, {3,26}, {4,25}, {5,34}, {6,33},
                        {7,32}, {8,31}, {9,30}, {10,39}, {11,38}, {12,37},
                        {13,36}, {14,35}, {15,44}, {16,43}, {17,42}, {18,41},
                        {19,40}, {20,49}, {21,48}, {22,47}, {23,46}, {24,45} });
      pairs.push_back({ {0,49}, {1,44}, {2,39}, {3,34}, {4,29}, {5,48}, {6,43},
                        {7,38}, {8,33}, {9,28}, {10,47}, {11,42}, {12,37},
                        {13,32}, {14,27}, {15,46}, {16,41}, {17,36}, {18,31},
                        {19,26}, {20,45}, {21,40}, {22,35}, {23,30}, {24,25} });
      pairs.push_back({ {0,45}, {1,46}, {2,47}, {3,48}, {4,49}, {5,40}, {6,41},
                        {7,42}, {8,43}, {9,44}, {10,35}, {11,36}, {12,37},
                        {13,38}, {14,39}, {15,30}, {16,31}, {17,32}, {18,33},
                        {19,34}, {20,25}, {21,26}, {22,27}, {23,28}, {24,29} });
      aux.resize(50);
      break;
    case FEMapper_type::N1_N1_3D:
      name = "N1_N1";
      description = "nonconforming Q1Rot element";
      n_dof0 = 1;
      n_dof1 = 1;
      n_pairs = 1;
      pairs.push_back({ {0,1} });
      pairs.push_back({ {0,1} });
      pairs.push_back({ {0,1} });
      pairs.push_back({ {0,1} });
      aux.resize(2);
      break;
    case FEMapper_type::N2_N2_3D:
      name = "N2_N2";
      description = "nonconforming element of order 2";
      n_dof0 = 3;
      n_dof1 = 3;
      n_pairs = 3;
      pairs.push_back({ {0,3}, {1,5}, {2,4} });
      pairs.push_back({ {0,3}, {1,4}, {2,5} });
      pairs.push_back({ {0,3}, {1,5}, {2,4} });
      pairs.push_back({ {0,3}, {1,4}, {2,5} });
      aux.resize(6);
      break;
    case FEMapper_type::N3_N3_3D:
      name = "N3_N3";
      description = "nonconforming element of order 3";
      n_dof0 = 6;
      n_dof1 = 6;
      n_pairs = 6;
      pairs.push_back({ {0,6}, {1,8}, {2,7}, {3,11}, {4,10}, {5,9} });
      pairs.push_back({ {0,6}, {1,7}, {2,8}, {3,9}, {4,10}, {5,11} });
      pairs.push_back({ {0,6}, {1,8}, {2,7}, {3,11}, {4,10}, {5,9} });
      pairs.push_back({ {0,6}, {1,7}, {2,8}, {3,9}, {4,10}, {5,11} });
      aux.resize(12);
      break;
    case FEMapper_type::N4_N4_3D:
      name = "N4_N4";
      description = "nonconforming element of order 4";
      n_dof0 = 10;
      n_dof1 = 10;
      n_pairs = 10;
      pairs.push_back({ {0,10}, {1,12}, {2,11}, {3,15}, {4,14}, {5,13}, {6,19},
                        {7,18}, {8,17}, {9,16} });
      pairs.push_back({ {0,10}, {1,11}, {2,12}, {3,13}, {4,14}, {5,15}, {6,16},
                        {7,17}, {8,18}, {9,19} });
      pairs.push_back({ {0,10}, {1,12}, {2,11}, {3,15}, {4,14}, {5,13}, {6,19},
                        {7,18}, {8,17}, {9,16} });
      pairs.push_back({ {0,10}, {1,11}, {2,12}, {3,13}, {4,14}, {5,15}, {6,16},
                        {7,17}, {8,18}, {9,19} });
      aux.resize(20);
      break;
    case FEMapper_type::NP1_NP1_3D:
      name = "NP1_NP1";
      description = "nonconforming P1 element";
      n_dof0 = 1;
      n_dof1 = 1;
      n_pairs = 1;
      pairs.push_back({ {0,1} });
      pairs.push_back({ {0,1} });
      pairs.push_back({ {0,1} });
      aux.resize(2);
      break;
    case FEMapper_type::NP2_NP2_3D:
      name = "NP2_NP2";
      description = "nonconforming P2 element";
      n_dof0 = 3;
      n_dof1 = 3;
      n_pairs = 3;
      pairs.push_back({ {0,3}, {1,5}, {2,4} });
      pairs.push_back({ {0,4}, {1,3}, {2,5} });
      pairs.push_back({ {0,5}, {1,4}, {2,3} });
      aux.resize(6);
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
      ErrThrow("you tried to construct an FEMapper of type ", type, " which is "
               "used for one-regular grids only. For that you need to use the "
               "derived class FEMapper1Reg3D");
      break;
    default:
      ErrThrow("unknown fe mapper type ", t);
  }
}

void FEMapper::map(int *GlobalNumbers, int I_K0, int I_K1,
                   const int *Indices0, const int *Indices1,
                   int &counter, int map_type) const
{
  for(int i = 0; i < n_dof0; i++)
  {
    aux[i] = I_K0 + Indices0[i];
  }

  for(int i = 0; i < n_dof1; i++)
  {
    aux[i+n_dof0] = I_K1 + Indices1[i];
  }

  if(n_pairs != 0)
  {
    auto& current_pairs = pairs.at(map_type);
    for(auto pair : current_pairs)
      map_single_dof(GlobalNumbers, aux[pair.first], aux[pair.second], counter);
  }
}

void FEMapper::map_boundary(int *GlobalNumbers, int I_K, const int *Indices,
                            int &BoundCounter,
                            std::vector<int>& hanging_numbers) const
{
  // note that `GlobalNumbers` is `TFESpace::GlobalNumbers`
  for(int i = 0; i < n_dof0; i++)
  {
    // i is the local (wrt the joint) dof
    // Indices[i] is the local (wrt the cell) dof
    map_single_boundary_dof(GlobalNumbers, I_K + Indices[i], BoundCounter,
                            hanging_numbers);
  }
}

void FEMapper::map_boundary_edge(int N_EdgeDOF, int *GlobalNumbers, int I_K,
                                 const int *Indices, int &BoundCounter,
                                 std::vector<int>& hanging_numbers) const
{
  for(int i=0;i<N_EdgeDOF;i++)
  {
    // i is the local (wrt the edge) dof
    // Indices[i] is the local (wrt the cell) dof
    map_single_boundary_dof(GlobalNumbers, I_K + Indices[i], BoundCounter,
                            hanging_numbers);
  }
}

void FEMapper::map_boundary_vertex(int* GlobalNumbers, int I_K, int Index,
                                   int& BoundCounter,
                                   std::vector<int>& hanging_numbers) const
{
  map_single_boundary_dof(GlobalNumbers, I_K + Index, BoundCounter,
                          hanging_numbers);
}

int FEMapper::find_in_global_numbers(const int *GlobalNumbers, int w) const
{
  int v;
  // -1 is the initial value for all entries in `GlobalNumbers`
  // non-negative values are understood as indices into the same array 
  // `GlobalNumbers`. This means two dofs are identified in the resulting fe
  // space
  while( (v=GlobalNumbers[w]) > -1)
  {
    w=v;
  }
  return w;
}

void FEMapper::map_single_dof(int *GlobalNumbers, int dof0, int dof1,
                              int &Counter) const
{
  int w0 = find_in_global_numbers(GlobalNumbers, dof0);
  int w1 = find_in_global_numbers(GlobalNumbers, dof1);

  if( (GlobalNumbers[dof0] == -1) && (GlobalNumbers[dof1] == -1) )
  {
    // there are no information => lower index dof[01] gets right to define
    Counter--;
    if(dof0<dof1)
    {
      GlobalNumbers[dof0]=Counter; 
      GlobalNumbers[dof1]=dof0;
    }
    else
    {
      GlobalNumbers[dof1]=Counter;
      GlobalNumbers[dof0]=dof1;
    }
    // cout << "both" << endl;
  }
  else
  {
    if(GlobalNumbers[dof0] == -1)
    {
      // => dof0=w0
      // cout << "one 0" << endl;
      if(w1>dof0)
      {
        GlobalNumbers[dof0]=GlobalNumbers[w1];
        GlobalNumbers[dof1]=dof0;
        GlobalNumbers[w1]=dof0;
      }
      else
      {
        GlobalNumbers[dof0]=w1;
      }
    }
    else
    {
      if(GlobalNumbers[dof1] == -1)
      {
        // => dof1=w1
        // cout << "one 1" << endl;
        if(w0>dof1)
        {
          GlobalNumbers[dof1]=GlobalNumbers[w0];
          GlobalNumbers[dof0]=dof1;
          GlobalNumbers[w0]=dof1;
        }
        else
        {
          GlobalNumbers[dof1]=w0;
        }
      }
      else
      {
        if(GlobalNumbers[w0]!=GlobalNumbers[w1])
        {
          int e = std::max(GlobalNumbers[w0],GlobalNumbers[w1]);
          int w = std::min(w0,w1);

          GlobalNumbers[w0] = w;
          GlobalNumbers[w1] = w;
          GlobalNumbers[w] = e;
        }
        else
        {
          if(w0!=w1)
          {
            ErrThrow("This should not happen. The details are not fully "
                     "understood at the moment");
          }
        }
      }
    }
  }

  // some sanity checks
  if(GlobalNumbers[dof0]>=dof0)
  {
    ErrThrow("error with dof0 (", dof0, ")");
  }
  if(GlobalNumbers[w0]>=w0)
  {
    ErrThrow("error with w0");
  }
  if(GlobalNumbers[dof1]>=dof1)
  {
    ErrThrow("error with dof1 (", dof1, ")");
  }
  if(GlobalNumbers[w1]>=w1)
  {
    ErrThrow("error with w1");
  }
}

void FEMapper::map_single_boundary_dof(int *GlobalNumbers, int dof,
                                       int &BoundCounter,
                                       std::vector<int>& hanging_numbers) const
{
  int w = find_in_global_numbers(GlobalNumbers, dof);

  // `BoundCounter` is always negative and smaller or equal to -9
  // -1 is the initial value for all entries in `GlobalNumbers`
  if(GlobalNumbers[w] == -1) // there are no information
  {
    BoundCounter--;
    GlobalNumbers[dof]=BoundCounter;
  } // endif
  else
  {
    if(GlobalNumbers[w] < BoundCounter)
    {
      // dof is marked as an inner dof or a later boundary dof => modify mark
      // (later in the sense of the enumeration `BoundCond` in Constants.h.)
      BoundCounter--;
      GlobalNumbers[w]=BoundCounter;
      if(dof!=w) GlobalNumbers[dof]=w; // do not allow recursive linking
    } // endif
    else
    {
      if(GlobalNumbers[w] == HANGINGNODE)
      {
        /// @note we can only get here in 3D
        /// @note This could be wrong: maybe hanging nodes on boundary edges
        ///       should be left hanging, because otherwise the space might be
        ///       discontinuous (or otherwise non conforming). This depends also
        ///       on the boundary condition and data, for homogeneous Dirichlet
        ///       data this does work.
        // remove hanging node, mark it as boundary node
        int N_ = hanging_numbers.size();
        for(int j=0;j<N_;j++)
        {
          if(dof == hanging_numbers[j])
          {
            hanging_numbers[j] = -1;
            BoundCounter--;
            GlobalNumbers[w]=BoundCounter;
            if(dof!=w) GlobalNumbers[dof]=w; // do not allow recursive linking
          }
        } // endfor j
      }
    }
  } //endelse
}

std::ostream & operator<<(std::ostream& out, const FEMapper_type t)
{
  const char* s = 0;
#define PROCESS_VAL(p)                                                         \
  case(FEMapper_type::p):                                             \
    s = #p;                                                                    \
    break;
  switch(t)
  {
    PROCESS_VAL(C0_2_C0_2_2D);
    PROCESS_VAL(C1_2_C1_2_2D);
    PROCESS_VAL(C2_2_C2_2_2D);
    PROCESS_VAL(C3_2_C3_2_2D);
    PROCESS_VAL(C4_2_C4_2_2D);
    PROCESS_VAL(C5_2_C5_2_2D);
    PROCESS_VAL(C6_2_C6_2_2D);
    PROCESS_VAL(C7_2_C7_2_2D);
    PROCESS_VAL(C8_2_C8_2_2D);
    PROCESS_VAL(C9_2_C9_2_2D);
    PROCESS_VAL(N1_2_N1_2_2D);
    PROCESS_VAL(N2_2_N2_2_2D);
    PROCESS_VAL(N3_2_N3_2_2D);
    PROCESS_VAL(N4_2_N4_2_2D);
    PROCESS_VAL(N5_2_N5_2_2D);
    // ONE regular grid, same pattern on both sides
    PROCESS_VAL(C0_2_C0_2_1Reg_2D);
    PROCESS_VAL(C1_2_C1_2_1Reg_2D);
    PROCESS_VAL(C2_2_C2_2_1Reg_2D);
    PROCESS_VAL(C3_2_C3_2_1Reg_2D);
    PROCESS_VAL(C4_2_C4_2_1Reg_2D);
    PROCESS_VAL(C5_2_C5_2_1Reg_2D);
    PROCESS_VAL(N1_2_N1_2_1Reg_2D);
    PROCESS_VAL(N2_2_N2_2_1Reg_2D);
    PROCESS_VAL(N3_2_N3_2_1Reg_2D);
    PROCESS_VAL(N4_2_N4_2_1Reg_2D);
    PROCESS_VAL(N5_2_N5_2_1Reg_2D);
    // ONE regular grid, different pattern
    PROCESS_VAL(N2_2_N1_2_1Reg_2D); // coarse: second order, fine: first order
    PROCESS_VAL(N3_2_N2_2_1Reg_2D); // coarse: third order,  fine: second order
    PROCESS_VAL(N4_2_N3_2_1Reg_2D); // coarse: fourth order, fine: third order
    PROCESS_VAL(N5_2_N4_2_1Reg_2D); // coarse: fifth order,  fine: fourth order
    // 3D mapper types
    PROCESS_VAL(D_D_3D);
    PROCESS_VAL(P1_P1_3D);
    PROCESS_VAL(P2_P2_3D);
    PROCESS_VAL(P3_P3_3D);
    PROCESS_VAL(P2B_P2B_3D);
    PROCESS_VAL(Q1_Q1_3D);
    PROCESS_VAL(Q2_Q2_3D);
    PROCESS_VAL(Q3_Q3_3D);
    PROCESS_VAL(Q4_Q4_3D);
    PROCESS_VAL(N1_N1_3D);
    PROCESS_VAL(N2_N2_3D);
    PROCESS_VAL(N3_N3_3D);
    PROCESS_VAL(N4_N4_3D);
    PROCESS_VAL(NP1_NP1_3D);
    PROCESS_VAL(NP2_NP2_3D);
      break;
    // 1-regular grid, same pattern on both sides
    PROCESS_VAL(P1_P1_1Reg_3D);
    PROCESS_VAL(P2_P2_1Reg_3D);
    PROCESS_VAL(P3_P3_1Reg_3D);
    PROCESS_VAL(Q1_Q1_1Reg_3D);
    PROCESS_VAL(Q2_Q2_1Reg_3D);
    PROCESS_VAL(Q3_Q3_1Reg_3D);
    PROCESS_VAL(NP1_NP1_1Reg_3D);
    PROCESS_VAL(NP2_NP2_1Reg_3D);
      break;
    default:
      s = "unknown fe mapper type";
  }
#undef PROCESS_VAL
  return out << s;
}
