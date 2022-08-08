#include <FE3DMapper1Reg.h>
#include <HNDesc.h>
#include <MooNMD_Io.h>

FEMapper1Reg3D::FEMapper1Reg3D(FEMapper_type t)
 : FEMapper(), 
   twist_permutation{}, CurrentPairs{nullptr}, CurrentNoOp{nullptr},
   N_NoOpposite{0}, no_opposite{}, N_Hanging{0}, Hanging{}, HangingTypes{},
   coupling{}
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
      ErrThrow("you tried to construct an FEMapper1Reg3D of type ", type,
               " which is used for regular grids in 2D only. For that you need "
               "to use the base class FEMapper");
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
      ErrThrow("you tried to construct an FEMapper1Reg3D of type ", type,
               " which is used for one-regular grids in 2D only. For that you "
               "need to use the class FEMapper1Reg2D");
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
      ErrThrow("you tried to construct an FEMapper1Reg3D of type ", type,
               " which is used regular grids in 3D only. For that you need to "
               "use the base class FEMapper");
      break;
    // 1-regular grid, same pattern on both sides
    case FEMapper_type::P1_P1_1Reg_3D:
      name = "P1_P1_1Reg";
      description = "conforming P1 element, 1-regular";
      n_dof0 = 3;
      n_dof1 = 3;
      n_pairs = 9;
      pairs.push_back({ {0,12}, {1,5}, {2,7}, {3,14}, {4,8}, {5,9}, {6,13},
                        {7,11}, {8,10} });
      aux.resize(15);
      N_Hanging = 3;
      Hanging = { 1, 4, 2 };
      HangingTypes = { HNDesc::HN_C_P1_3D_E, HNDesc::HN_C_P1_3D_E,
                       HNDesc::HN_C_P1_3D_E };
      coupling = { { 0, 3 }, { 3, 6 }, { 0, 6 } };
      twist_permutation = { { 0, 1, 2 }, { 1, 2, 0 }, { 2, 0, 1 } };
      break;
    case FEMapper_type::P2_P2_1Reg_3D:
      name = "P2_P2_1Reg";
      description = "conforming P2 element, 1-regular";
      n_dof0 = 6;
      n_dof1 = 6;
      n_pairs = 15;
      pairs.push_back({ {0,24}, {2,11}, {4,21}, {5,14}, {6,29}, {8,17}, {10,19},
                        {11,18}, {12,26}, {14,23}, {16,22}, {17,20}, {18,27},
                        {20,28}, {23,25} });
      aux.resize(30);
      N_Hanging = 9;
      Hanging = { 1, 3, 4, 7, 9, 10, 13, 15, 16 };
      HangingTypes = { HNDesc::HN_C_P2_3D_E, HNDesc::HN_C_P2_3D_E,
                       HNDesc::HN_C_P2_3D_F, HNDesc::HN_C_P2_3D_E,
                       HNDesc::HN_C_P2_3D_E, HNDesc::HN_C_P2_3D_F,
                       HNDesc::HN_C_P2_3D_E, HNDesc::HN_C_P2_3D_E,
                       HNDesc::HN_C_P2_3D_F };
      coupling = { { 0, 2, 6 }, { 0, 5, 12 }, { 2, 6, 5, 8, 12 }, { 6, 8, 12 },
                   { 6, 2, 0 }, { 8, 12, 2, 5, 0 }, { 12, 5, 0 }, { 12, 8, 6 },
                   { 5, 0, 8, 2, 6 } };
      twist_permutation = { { 0, 1, 2, 3, 4, 5 }, { 2, 4, 5, 1, 3, 0 },
                            { 5, 3, 0, 4, 1, 2 } };
      break;
    case FEMapper_type::P3_P3_1Reg_3D:
      name = "P3_P3_1Reg";
      description = "conforming P3 element, 1-regular";
      n_dof0 = 10;
      n_dof1 = 10;
      n_pairs = 22;
      pairs.push_back({ {0,40}, {2,44}, {3,19}, {6,34}, {7,41}, {8,37}, {9,23},
                        {10,49}, {12,48}, {13,29}, {16,32}, {17,47}, {18,31},
                        {19,30}, {20,43}, {22,42}, {23,39}, {26,38}, {27,46},
                        {28,36}, {29,33}, {35,45} });
      aux.resize(50);
      N_Hanging = 18;
      Hanging = { 1, 4, 11, 14, 21, 24, 3, 13, 23, 5, 15, 25, 6, 8, 16, 18, 26,
                  28 };
      HangingTypes = { HNDesc::HN_C_P3_3D_E, HNDesc::HN_C_P3_3D_E,
                       HNDesc::HN_C_P3_3D_E, HNDesc::HN_C_P3_3D_E,
                       HNDesc::HN_C_P3_3D_E, HNDesc::HN_C_P3_3D_E,
                       HNDesc::HN_C_P3_3D_M, HNDesc::HN_C_P3_3D_M,
                       HNDesc::HN_C_P3_3D_M, HNDesc::HN_C_P3_3D_F,
                       HNDesc::HN_C_P3_3D_F, HNDesc::HN_C_P3_3D_F,
                       HNDesc::HN_C_P3_3D_G, HNDesc::HN_C_P3_3D_G,
                       HNDesc::HN_C_P3_3D_G, HNDesc::HN_C_P3_3D_G,
                       HNDesc::HN_C_P3_3D_G, HNDesc::HN_C_P3_3D_G };
      coupling = { { 0, 2, 17, 10 }, { 0, 7, 22, 20 }, { 10, 12, 27, 20 },
                   { 10, 17, 2, 0 }, { 20, 22, 7, 0 }, { 20, 27, 12, 10 },
                   { 0, 2, 17, 10 }, { 10, 12, 27, 20 }, { 0, 7, 22, 20 },
                   { 2, 7, 17, 35, 22, 10, 12, 27, 20 },
                   { 12, 17, 27, 35, 2, 20, 22, 7, 0 },
                   { 22, 27, 7, 35, 12, 0, 2, 17, 10 },
                   { 0, 2, 17, 10, 7, 35, 12, 22, 27, 20 },
                   { 0, 7, 22, 20, 2, 35, 27, 17, 12, 10 },
                   { 10, 12, 27, 20, 17, 35, 22, 2, 7, 0 },
                   { 10, 17, 2, 0, 12, 35, 7, 27, 22, 20 },
                   { 20, 22, 7, 0, 27, 35, 2, 12, 17, 10 },
                   { 20, 27, 12, 10, 22, 35, 17, 7, 2, 0 } };
      twist_permutation = { { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 },
                            { 3, 6, 8, 9, 2, 5, 7, 1, 4, 0 },
                            { 9, 7, 4, 0, 8, 5, 1, 6, 2, 3 }, };
      break;
    case FEMapper_type::Q1_Q1_1Reg_3D:
      name = "Q1_Q1_1Reg";
      description = "conforming Q1 element, 1-regular";
      n_dof0 = 4;
      n_dof1 = 4;
      n_pairs = 11;
      pairs.push_back({ {0,16}, { 1,6}, {2,13}, {3,7}, {4,18}, {5,10}, {7,11},
                        {8,19}, {9,14}, {11,15}, {12,17} });
      aux.resize(20);
      N_Hanging = 5;
      Hanging = { 1, 5, 9, 2, 3 };
      HangingTypes = { HNDesc::HN_C_Q1_3D_E, HNDesc::HN_C_Q1_3D_E,
                       HNDesc::HN_C_Q1_3D_E, HNDesc::HN_C_Q1_3D_E,
                       HNDesc::HN_C_Q1_3D_F };
      coupling = { { 0, 4 }, { 4, 8 }, { 8, 12 }, { 0, 12 }, { 0, 4, 8, 12 } };
      twist_permutation = { { 0, 1, 2, 3 }, { 1, 3, 0, 2 }, { 3, 2, 1, 0 },
                            { 2, 0, 3, 1 } };
      break;
    case FEMapper_type::Q2_Q2_1Reg_3D:
      name = "Q2_Q2_1Reg";
      description = "conforming Q2 element, 1-regular";
      n_dof0 = 9;
      n_dof1 = 9;
      n_pairs = 20;
      pairs.push_back({ {0,36}, {2,39}, {5,16}, {6,37}, {7,32}, {8,40}, {9,42},
                        {11,43}, {14,25}, {2,15}, {8,17}, {18,44}, {20,41},
                        {23,34}, {11,24}, {17,26}, {27,38}, {6,29}, {20,33},
                        {26,35} });
      aux.resize(45);
      N_Hanging = 16;
      Hanging = { 1, 3, 10, 12, 19, 21, 28, 30, 5, 14, 23, 7, 4, 13, 22, 31 };
      HangingTypes = { HNDesc::HN_C_Q2_3D_E, HNDesc::HN_C_Q2_3D_E,
                       HNDesc::HN_C_Q2_3D_E, HNDesc::HN_C_Q2_3D_E,
                       HNDesc::HN_C_Q2_3D_E, HNDesc::HN_C_Q2_3D_E,
                       HNDesc::HN_C_Q2_3D_E, HNDesc::HN_C_Q2_3D_E,
                       HNDesc::HN_C_Q2_3D_E, HNDesc::HN_C_Q2_3D_E,
                       HNDesc::HN_C_Q2_3D_E, HNDesc::HN_C_Q2_3D_E,
                       HNDesc::HN_C_Q2_3D_F, HNDesc::HN_C_Q2_3D_F,
                       HNDesc::HN_C_Q2_3D_F, HNDesc::HN_C_Q2_3D_F };
      coupling = { { 0, 2, 9 }, { 0, 6, 27 }, { 9, 11, 18 }, { 9, 2, 0 },
                   { 18, 20, 27 }, { 18, 11, 9 }, { 27, 29, 0 }, { 27, 20, 18 },
                   { 2, 8, 20 }, { 11, 8, 6 }, { 20, 8, 2 }, { 6, 8, 11 },
                   { 0, 6, 27, 2, 8, 20, 9, 11, 18 },
                   { 9, 2, 0, 11, 8, 6, 18, 20, 27 },
                   { 18, 11, 9, 20, 8, 2, 27, 6, 0 },
                   { 27, 20, 18, 6, 8, 11, 0, 2, 9 } };
      twist_permutation = { { 0, 1, 2, 3, 4, 5, 6, 7, 8 },
                            { 2, 5, 8, 1, 4, 7, 0, 3, 6 },
                            { 8, 7, 6, 5, 4, 3, 2, 1, 0 },
                            { 6, 3, 0, 7, 4, 1, 8, 5, 2 } };
      N_NoOpposite = 12;
      no_opposite.push_back({ 4, 13, 22, 31, 1, 3, 10, 12, 19, 21, 28, 30 });
      break;
    case FEMapper_type::Q3_Q3_1Reg_3D:
      name = "Q3_Q3_1Reg";
      description = "conforming Q3 element, 1-regular";
      n_dof0 = 16;
      n_dof1 = 16;
      n_pairs = 31;
      pairs.push_back({ {0,64}, {2,68}, {3,28}, {7,29}, {8,65}, {10,69},
                        {11,30}, {12,51}, {13,55}, {14,59}, {15,31}, {16,76},
                        {18,77}, {19,44}, {23,45}, {24,72}, {26,73}, {27,46},
                        {31,47}, {32,79}, {34,75}, {35,60}, {39,61}, {40,78},
                        {42,74}, {43,62}, {47,63}, {48,67}, {50,66}, {56,71},
                        {58,70} });
      aux.resize(80);
      N_Hanging = 33;
      Hanging = { 1, 20, 9, 22, 54, 41, 52, 33, 4, 49, 6, 57, 25, 38, 17, 36,
                  3, 11, 43, 35, 12, 14, 27, 19, 7, 23, 39, 13, 5, 21, 37, 53,
                  15 };
      HangingTypes = { HNDesc::HN_C_Q3_3D_1, HNDesc::HN_C_Q3_3D_1,
                       HNDesc::HN_C_Q3_3D_1, HNDesc::HN_C_Q3_3D_1,
                       HNDesc::HN_C_Q3_3D_1, HNDesc::HN_C_Q3_3D_1,
                       HNDesc::HN_C_Q3_3D_1, HNDesc::HN_C_Q3_3D_1,
                       HNDesc::HN_C_Q3_3D_1, HNDesc::HN_C_Q3_3D_1,
                       HNDesc::HN_C_Q3_3D_1, HNDesc::HN_C_Q3_3D_1,
                       HNDesc::HN_C_Q3_3D_1, HNDesc::HN_C_Q3_3D_1,
                       HNDesc::HN_C_Q3_3D_1, HNDesc::HN_C_Q3_3D_1,
                       HNDesc::HN_C_Q3_3D_2, HNDesc::HN_C_Q3_3D_2,
                       HNDesc::HN_C_Q3_3D_2, HNDesc::HN_C_Q3_3D_2,
                       HNDesc::HN_C_Q3_3D_2, HNDesc::HN_C_Q3_3D_2,
                       HNDesc::HN_C_Q3_3D_2, HNDesc::HN_C_Q3_3D_2,
                       HNDesc::HN_C_Q3_3D_3, HNDesc::HN_C_Q3_3D_3,
                       HNDesc::HN_C_Q3_3D_3, HNDesc::HN_C_Q3_3D_3,
                       HNDesc::HN_C_Q3_3D_4, HNDesc::HN_C_Q3_3D_4,
                       HNDesc::HN_C_Q3_3D_4, HNDesc::HN_C_Q3_3D_4,
                       HNDesc::HN_C_Q3_3D_5 };
      coupling = {{ 0, 2, 24, 16 }, { 16, 24, 2, 0 }, { 8, 10, 26, 18 },
                  { 18, 26, 10, 8 }, { 50, 58, 42, 40 }, { 40, 42, 58, 50 },
                  { 48, 56, 34, 32 }, { 32, 34, 56, 48 }, { 0, 8, 50, 48 },
                  { 48, 50, 8, 0 }, { 2, 10, 58, 56 }, { 56, 58, 10, 2 },
                  { 24, 26, 42, 34 }, { 34, 42, 26, 24 }, { 16, 18, 40, 32 },
                  { 32, 40, 18, 16 }, { 0, 2, 24, 16 }, { 8, 10, 26, 18 },
                  { 50, 58, 42, 40 }, { 48, 56, 34, 32 }, { 0, 8, 50, 48 },
                  { 2, 10, 58, 56 }, { 24, 26, 42, 34 }, { 16, 18, 40, 32 },
                  {0, 2, 24, 16, 8, 10, 26, 18, 50, 58, 42, 40, 48, 56, 34, 32},
                  {16, 18, 40, 32, 24, 26, 42, 34, 2, 10, 58, 56, 0, 8,  50,48},
                  {32, 34, 56, 48, 40, 42, 58, 50, 18, 26, 10, 8, 16, 24, 2, 0},
                  {48, 50, 8, 0, 56, 58, 10, 2, 34, 42, 26, 24, 32, 40, 18, 16},
                  {0, 2, 24, 16, 8, 10, 26, 18, 50, 58, 42, 40, 48, 56, 34, 32},
                  {16, 18, 40, 32, 24, 26, 42, 34, 2, 10, 58, 56, 0, 8,  50,48},
                  {32, 34, 56, 48, 40, 42, 58, 50, 18, 26, 10, 8, 16, 24, 2, 0},
                  {48, 50, 8, 0, 56, 58, 10, 2, 34, 42, 26, 24, 32, 40, 18, 16},
                  {0, 2, 24, 16, 8, 10, 26, 18, 50, 58, 42, 40, 48, 56, 34, 32}
                 };
      twist_permutation = {
                     {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 },
                     { 3, 7, 11, 15, 2, 6, 10, 14, 1, 5, 9, 13, 0, 4, 8, 12 },
                     { 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0 },
                     { 12, 8, 4, 0, 13, 9, 5, 1, 14, 10, 6, 2, 15, 11, 7, 3 } };
      break;
    case FEMapper_type::NP1_NP1_1Reg_3D:
      name = "NP1_NP1_1Reg";
      description = "conforming NP1 element, 1-regular";
      n_dof0 = 1;
      n_dof1 = 1;
      n_pairs = 0;
      aux.resize(5);
      N_Hanging = 1;
      Hanging = { 4 };
      HangingTypes = { HNDesc::HN_N_P1_3D_E };
      coupling = { { 0, 1, 2, 3 } };
      twist_permutation = { { 0 }, { 0 }, { 0 } };
      N_NoOpposite = 4;
      no_opposite.push_back({ 0, 1, 2, 3 });
      break;
    case FEMapper_type::NP2_NP2_1Reg_3D:
      name = "NP2_NP2_1Reg";
      description = "conforming NP2 element, 1-regular";
      n_dof0 = 3;
      n_dof1 = 3;
      n_pairs = 0;
      aux.resize(15);
      N_Hanging = 3;
      Hanging = { 12, 13, 14 };
      HangingTypes = { HNDesc::HN_N_P2_3D_0, HNDesc::HN_N_P2_3D_1,
                       HNDesc::HN_N_P2_3D_2 };
      coupling = { { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 },
                   { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 },
                   { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 } };
      twist_permutation = { { 0, 1, 2 }, { 1, 2, 0 }, { 2, 0, 1 } };
      N_NoOpposite = 12;
      no_opposite.push_back({ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 });
      break;
      ErrThrow("you tried to construct an FEMapper1Reg2D of type ", type,
               " which is used for one-regular grids in 3D only. For that you "
               "need to use the derived class FEMapper1Reg3D");
      break;
    default:
      ErrThrow("unknown fe mapper type ", t);
  }

  CurrentPairs = pairs[0].data();
  CurrentNoOp = N_NoOpposite ? no_opposite[0].data() : nullptr;
}

/** map the given local degrees of freedom,
    coarse cell has lower number,
    F0, F1, F2, F3 on finer side, C on coarser side */
void FEMapper1Reg3D::Map(int *Global,
             int I_KF0, int I_KF1, int I_KF2, int I_KF3,
             int I_KC,
             const int *IndicesF0, const int *IndicesF1, const int *IndicesF2,
             const int *IndicesF3, const int *IndicesC,
             int TwistIndexF0, int TwistIndexF1, int TwistIndexF2,
             int TwistIndexF3, int TwistIndexC,
             int DirichletBound,
             int &Counter,
             std::vector<THangingNode*>& vect,
             std::vector<int>& numbers,
             std::map<HNDesc, THNDesc *>& hn_descriptors) const
{
  // fill in cell local dofs according to twist index

  int j=0;
  // cout << "Twist0: " << TwistIndexF0 << endl;
  const int *CurrentTwistPerm = twist_permutation[TwistIndexF0].data();
  for(int i=0;i<n_dof0;i++)
  {
    aux[j] = I_KF0 + IndicesF0[CurrentTwistPerm[i]];
    // cout << IndicesF0[i] << " - " << IndicesF0[CurrentTwistPerm[i]];
    // cout << " - " << aux[j] << endl;
    j++;
  }

  // cout << "Twist1: " << TwistIndexF1 << endl;
  CurrentTwistPerm = twist_permutation[TwistIndexF1].data();
  for(int i=0;i<n_dof0;i++)
  {
    aux[j] = I_KF1 + IndicesF1[CurrentTwistPerm[i]];
    // cout << IndicesF1[i] << " - " << IndicesF1[CurrentTwistPerm[i]];
    // cout << " - " << Aux[j] << endl;
    j++;
  }

  // cout << "Twist2: " << TwistIndexF2 << endl;
  CurrentTwistPerm = twist_permutation[TwistIndexF2].data();
  for(int i=0;i<n_dof0;i++)
  {
    aux[j] = I_KF2 + IndicesF2[CurrentTwistPerm[i]];
    // cout << IndicesF2[i] << " - " << IndicesF2[CurrentTwistPerm[i]];
    // cout << " - " << Aux[j] << endl;
    j++;
  }

  // cout << "Twist3: " << TwistIndexF3 << endl;
  CurrentTwistPerm = twist_permutation[TwistIndexF3].data();
  for(int i=0;i<n_dof0;i++)
  {
    aux[j] = I_KF3 + IndicesF3[CurrentTwistPerm[i]];
    // cout << IndicesF3[i] << " - " << IndicesF3[CurrentTwistPerm[i]];
    // cout << " - " << Aux[j] << endl;
    j++;
  }

  // cout << "MapType: " << TwistIndexC << endl;
  CurrentTwistPerm = twist_permutation[TwistIndexC].data();
  for(int i=0;i<n_dof1;i++)
  {
    aux[j] = I_KC + IndicesC[CurrentTwistPerm[i]];
    // cout << IndicesC[i] << " - " << IndicesC[CurrentTwistPerm[i]];
    // cout << " - " << Aux[j] << endl;
    j++;
  }

  for(int i=0;i<n_pairs;i++)
  {
    if(CurrentPairs[i].first!=-1)
    {
      // cout << "local: " << CurrentPairs[i].first << "  ";
      // cout << CurrentPairs[i].second << endl;
      // cout << "global: " << Aux[CurrentPairs[i].first] << "  ";
      // cout << Aux[CurrentPairs[i].second] << endl;
      map_single_dof(Global, aux[CurrentPairs[i].first],
                     aux[CurrentPairs[i].second], Counter);
    }
  }

  for(int i=0;i<N_Hanging;i++)
  {
    int w = find_in_global_numbers(Global, aux[Hanging[i]]);

    // cout << i << " type: " << HangingTypes[i] << endl;
    // cout << i << " Global[w]: " << Global[w] << " " << w << endl;

    // create a hanging node only once
    if(Global[w] == HANGINGNODE) continue;

    if(Global[w] != -1)
      if(Global[w] > DirichletBound) continue;

    Global[w]=HANGINGNODE;

    auto it = hn_descriptors.find(HangingTypes[i]);
    if(it == hn_descriptors.end())
    {
      auto pair = std::pair<HNDesc, THNDesc *>(HangingTypes[i],
                                               new THNDesc(HangingTypes[i]));
      it = hn_descriptors.insert(pair).first;
    }

    std::vector<int> partner_dofs;
    partner_dofs.reserve(N_Hanging-1);
    for(int k = 0; k < N_Hanging; ++k)
    {
      if(k != i)
      {
        partner_dofs.push_back(aux[Hanging[k]]);
      }
    }
    auto hn = new THangingNode(it->second, aux.data(), coupling[i].data(),
                               partner_dofs);
    vect.push_back(hn);
    numbers.push_back(aux[Hanging[i]]);
    // cout << "HN: " << Aux[Hanging[i]] << endl;
  }

  for(int i=0;i<N_NoOpposite;i++)
  {
    int w = find_in_global_numbers(Global, aux[CurrentNoOp[i]]);
    if( Global[w] == -1 )
    {
      Counter--;
      Global[aux[CurrentNoOp[i]]]=Counter;
    }
  } // endfor i

}
