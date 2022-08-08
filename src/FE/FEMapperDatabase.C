#include "FEMapperDatabase.h"
#include "FE2DMapper1Reg.h"
#include "FE3DMapper1Reg.h"
#include "MooNMD_Io.h"
#include <limits>

const FEMapper * FEMapperDatabase::get_mapper(FEMapper_type type, bool one_reg,
                                              bool two_d)
{
  auto it = mappers.find(type);
  if(it == mappers.end())
  {
    FEMapper * mapper_pointer;
    if(one_reg)
    {
      if(two_d)
        mapper_pointer = new FEMapper1Reg2D(type);
      else
        mapper_pointer = new FEMapper1Reg3D(type);
    }
    else
    {
      mapper_pointer = new FEMapper(type);
    }
    auto it_success = mappers.insert(
        {type, std::unique_ptr<FEMapper>{mapper_pointer}});
    if(it_success.second == false)
    {
      ErrThrow("map insertion failed for fe mapper type ", type);
    }
    it = it_success.first;
  }
  return it->second.get();
}

/// @brief methods to map from two fe types to a mapper type.
/// 
/// The four methods below define which mapper type is used for two given fe
/// types. That means here one considers a joint and two neighboring cells, 
/// where each cell has a possibly different finite element type. 
///
/// Special care needs to be taken when hanging vertices are involved. In 
/// ParMooN only one hanging vertex per joint is allowed, this is leads to 
/// so-called one-regular grids and requires special mappers (suffix
/// '_1regular').
/// 
/// A few notes on the implementation: In c++ it is not possible to to switch
/// over a pair, instead the `switch_pair()` methods below construct an unsigned
/// int which is unique for each possible combination of two finite element 
/// types. In order to use this in a switch statement, one needs these methods
/// to be `constexpr`.
///@{
FEMapper_type get_mapper_type(FEDescriptor_type FE1, FEDescriptor_type FE2);
FEMapper_type get_mapper_type_1regular2d(FEDescriptor_type FE1,
                                         FEDescriptor_type FE2);
FEMapper_type get_mapper_type_1regular3d(FEDescriptor_type FE1,
                                         FEDescriptor_type FE2);
///@}

const FEMapper * FEMapperDatabase::get_fe_mapper(FEDescriptor_type FE1,
                                                 FEDescriptor_type FE2)
{
  return get_mapper(get_mapper_type(FE1, FE2), false, false);
}

const FEMapper1Reg2D * FEMapperDatabase::get_fe_mapper_one_regular2d(
  FEDescriptor_type FE1, FEDescriptor_type FE2)
{
  return static_cast<const FEMapper1Reg2D*>(
     get_mapper(get_mapper_type_1regular2d(FE1, FE2), true, true));
}

const FEMapper1Reg3D * FEMapperDatabase::get_fe_mapper_one_regular3d(
  FEDescriptor_type FE1, FEDescriptor_type FE2)
{
  return static_cast<const FEMapper1Reg3D*>(
    get_mapper(get_mapper_type_1regular3d( FE1, FE2), true, false));
}


// compute the smallest number n such that 2^n >= value
constexpr int next_exponent(unsigned int value, size_t pow = 0)
{
  return value <= 1 ? 1 
                    : (((value-1) >> pow) ? next_exponent(value, pow + 1) 
                                          : pow);
}

// store the two finite element types into one unsigned number where the first
// bits are used for first type and the last bits for the second. The function
// `next_exponent()` is used to make sure that enough bits are reserved for a
// finite element type.
// see https://stackoverflow.com/a/46236919
constexpr unsigned int switch_pair( FEDescriptor_type FE1, FEDescriptor_type FE2)
{
  static_assert(
    std::numeric_limits<unsigned int>::digits > next_exponent(n_FEDescriptors),
    "not enough bits in unsigned int"
  );
  return ((unsigned int)FE1 << next_exponent(n_FEDescriptors)) + (unsigned int)FE2;
}

FEMapper_type get_mapper_type(FEDescriptor_type FE1, FEDescriptor_type FE2)
{
  using mt = FEMapper_type;
  switch(switch_pair(FE1, FE2))
  {
    case switch_pair(FE_C_T_P00_2D, FE_C_T_P00_2D): return mt::C0_2_C0_2_2D;
    case switch_pair(FE_C_T_P00_2D, FE_C_Q_Q00_2D): return mt::C0_2_C0_2_2D;
    case switch_pair(FE_C_Q_Q00_2D, FE_C_T_P00_2D): return mt::C0_2_C0_2_2D;
    case switch_pair(FE_C_Q_Q00_2D, FE_C_Q_Q00_2D): return mt::C0_2_C0_2_2D;

    case switch_pair(FE_C_T_P0_2D, FE_C_T_P0_2D): return mt::C0_2_C0_2_2D;
    case switch_pair(FE_C_T_P0_2D, FE_C_Q_Q0_2D): return mt::C0_2_C0_2_2D;
    case switch_pair(FE_C_Q_Q0_2D, FE_C_T_P0_2D): return mt::C0_2_C0_2_2D;
    case switch_pair(FE_C_Q_Q0_2D, FE_C_Q_Q0_2D): return mt::C0_2_C0_2_2D;

    case switch_pair(FE_C_T_P1_2D, FE_C_T_P1_2D): return mt::C1_2_C1_2_2D;
    case switch_pair(FE_C_T_P1_2D, FE_C_Q_Q1_2D): return mt::C1_2_C1_2_2D;
    case switch_pair(FE_C_Q_Q1_2D, FE_C_T_P1_2D): return mt::C1_2_C1_2_2D;
    case switch_pair(FE_C_Q_Q1_2D, FE_C_Q_Q1_2D): return mt::C1_2_C1_2_2D;
    case switch_pair(FE_C_T_P1_2D, FE_C_T_UL1_2D): return mt::C1_2_C1_2_2D;
    case switch_pair(FE_C_T_UL1_2D, FE_C_T_P1_2D): return mt::C1_2_C1_2_2D;
    
    case switch_pair(FE_C_T_P2_2D, FE_C_T_P2_2D): return mt::C2_2_C2_2_2D;
    case switch_pair(FE_C_T_P2_2D, FE_C_Q_Q2_2D): return mt::C2_2_C2_2_2D;
    case switch_pair(FE_C_Q_Q2_2D, FE_C_T_P2_2D): return mt::C2_2_C2_2_2D;
    case switch_pair(FE_C_Q_Q2_2D, FE_C_Q_Q2_2D): return mt::C2_2_C2_2_2D;
    case switch_pair(FE_C_T_P2_2D, FE_C_T_UL2_2D): return mt::C2_2_C2_2_2D;
    case switch_pair(FE_C_T_UL2_2D, FE_C_T_P2_2D): return mt::C2_2_C2_2_2D;
    
    case switch_pair(FE_C_T_B2_2D, FE_C_T_P2_2D): return mt::C2_2_C2_2_2D;
    case switch_pair(FE_C_T_B2_2D, FE_C_Q_Q2_2D): return mt::C2_2_C2_2_2D;
    case switch_pair(FE_C_T_P2_2D, FE_C_T_B2_2D): return mt::C2_2_C2_2_2D;
    case switch_pair(FE_C_Q_Q2_2D, FE_C_T_B2_2D): return mt::C2_2_C2_2_2D;
    case switch_pair(FE_C_T_B2_2D, FE_C_T_B2_2D): return mt::C2_2_C2_2_2D;

    case switch_pair(FE_C_T_SV2_2D, FE_C_T_P2_2D): return mt::C2_2_C2_2_2D;
    case switch_pair(FE_C_T_SV2_2D, FE_C_Q_Q2_2D): return mt::C2_2_C2_2_2D;
    case switch_pair(FE_C_T_P2_2D, FE_C_T_SV2_2D): return mt::C2_2_C2_2_2D;
    case switch_pair(FE_C_Q_Q2_2D, FE_C_T_SV2_2D): return mt::C2_2_C2_2_2D;
    case switch_pair(FE_C_T_SV2_2D, FE_C_T_SV2_2D): return mt::C2_2_C2_2_2D;

    case switch_pair(FE_C_T_P3_2D, FE_C_T_P3_2D): return mt::C3_2_C3_2_2D;
    case switch_pair(FE_C_T_P3_2D, FE_C_Q_Q3_2D): return mt::C3_2_C3_2_2D;
    case switch_pair(FE_C_Q_Q3_2D, FE_C_T_P3_2D): return mt::C3_2_C3_2_2D;
    case switch_pair(FE_C_Q_Q3_2D, FE_C_Q_Q3_2D): return mt::C3_2_C3_2_2D;
    case switch_pair(FE_C_T_P3_2D, FE_C_T_UL3_2D): return mt::C3_2_C3_2_2D;
    case switch_pair(FE_C_T_UL3_2D, FE_C_T_P3_2D): return mt::C3_2_C3_2_2D;

    case switch_pair(FE_C_T_B3_2D, FE_C_T_P3_2D): return mt::C3_2_C3_2_2D;
    case switch_pair(FE_C_T_B3_2D, FE_C_Q_Q3_2D): return mt::C3_2_C3_2_2D;
    case switch_pair(FE_C_T_P3_2D, FE_C_T_B3_2D): return mt::C3_2_C3_2_2D;
    case switch_pair(FE_C_Q_Q3_2D, FE_C_T_B3_2D): return mt::C3_2_C3_2_2D;
    case switch_pair(FE_C_T_B3_2D, FE_C_T_B3_2D): return mt::C3_2_C3_2_2D;

    case switch_pair(FE_C_T_P4_2D, FE_C_T_P4_2D): return mt::C4_2_C4_2_2D;
    case switch_pair(FE_C_T_P4_2D, FE_C_Q_Q4_2D): return mt::C4_2_C4_2_2D;
    case switch_pair(FE_C_Q_Q4_2D, FE_C_T_P4_2D): return mt::C4_2_C4_2_2D;
    case switch_pair(FE_C_Q_Q4_2D, FE_C_Q_Q4_2D): return mt::C4_2_C4_2_2D;
    case switch_pair(FE_C_T_P4_2D, FE_C_T_UL4_2D): return mt::C4_2_C4_2_2D;
    case switch_pair(FE_C_T_UL4_2D, FE_C_T_P4_2D): return mt::C4_2_C4_2_2D;

    case switch_pair(FE_C_T_B4_2D, FE_C_T_P4_2D): return mt::C4_2_C4_2_2D;
    case switch_pair(FE_C_T_B4_2D, FE_C_Q_Q4_2D): return mt::C4_2_C4_2_2D;
    case switch_pair(FE_C_T_P4_2D, FE_C_T_B4_2D): return mt::C4_2_C4_2_2D;
    case switch_pair(FE_C_Q_Q4_2D, FE_C_T_B4_2D): return mt::C4_2_C4_2_2D;
    case switch_pair(FE_C_T_B4_2D, FE_C_T_B4_2D): return mt::C4_2_C4_2_2D;

    case switch_pair(FE_C_T_P5_2D, FE_C_T_P5_2D): return mt::C5_2_C5_2_2D;
    case switch_pair(FE_C_T_P5_2D, FE_C_Q_Q5_2D): return mt::C5_2_C5_2_2D;
    case switch_pair(FE_C_Q_Q5_2D, FE_C_T_P5_2D): return mt::C5_2_C5_2_2D;
    case switch_pair(FE_C_Q_Q5_2D, FE_C_Q_Q5_2D): return mt::C5_2_C5_2_2D;
    case switch_pair(FE_C_T_P5_2D, FE_C_T_UL5_2D): return mt::C5_2_C5_2_2D;
    case switch_pair(FE_C_T_UL5_2D, FE_C_T_P5_2D): return mt::C5_2_C5_2_2D;

    case switch_pair(FE_C_T_P6_2D, FE_C_T_P6_2D): return mt::C6_2_C6_2_2D;
    case switch_pair(FE_C_T_P6_2D, FE_C_Q_Q6_2D): return mt::C6_2_C6_2_2D;
    case switch_pair(FE_C_Q_Q6_2D, FE_C_T_P6_2D): return mt::C6_2_C6_2_2D;
    case switch_pair(FE_C_Q_Q6_2D, FE_C_Q_Q6_2D): return mt::C6_2_C6_2_2D;

    case switch_pair(FE_C_T_P7_2D, FE_C_T_P7_2D): return mt::C7_2_C7_2_2D;
    case switch_pair(FE_C_T_P7_2D, FE_C_Q_Q7_2D): return mt::C7_2_C7_2_2D;
    case switch_pair(FE_C_Q_Q7_2D, FE_C_T_P7_2D): return mt::C7_2_C7_2_2D;
    case switch_pair(FE_C_Q_Q7_2D, FE_C_Q_Q7_2D): return mt::C7_2_C7_2_2D;

    case switch_pair(FE_C_T_P8_2D, FE_C_T_P8_2D): return mt::C8_2_C8_2_2D;
    case switch_pair(FE_C_T_P8_2D, FE_C_Q_Q8_2D): return mt::C8_2_C8_2_2D;
    case switch_pair(FE_C_Q_Q8_2D, FE_C_T_P8_2D): return mt::C8_2_C8_2_2D;
    case switch_pair(FE_C_Q_Q8_2D, FE_C_Q_Q8_2D): return mt::C8_2_C8_2_2D;

    case switch_pair(FE_C_T_P9_2D, FE_C_T_P9_2D): return mt::C9_2_C9_2_2D;
    case switch_pair(FE_C_T_P9_2D, FE_C_Q_Q9_2D): return mt::C9_2_C9_2_2D;
    case switch_pair(FE_C_Q_Q9_2D, FE_C_T_P9_2D): return mt::C9_2_C9_2_2D;
    case switch_pair(FE_C_Q_Q9_2D, FE_C_Q_Q9_2D): return mt::C9_2_C9_2_2D;

    case switch_pair(FE_N_T_RT0_2D, FE_N_T_RT0_2D): return mt::N1_2_N1_2_2D;
    case switch_pair(FE_N_T_RT0_2D, FE_N_Q_RT0_2D): return mt::N1_2_N1_2_2D;
    case switch_pair(FE_N_T_RT1_2D, FE_N_T_RT1_2D): return mt::N2_2_N2_2_2D;
    case switch_pair(FE_N_T_RT2_2D, FE_N_T_RT2_2D): return mt::N3_2_N3_2_2D;
    case switch_pair(FE_N_T_RT3_2D, FE_N_T_RT3_2D): return mt::N4_2_N4_2_2D;

    case switch_pair(FE_N_T_BDM1_2D, FE_N_T_BDM1_2D): return mt::N2_2_N2_2_2D;
    case switch_pair(FE_N_T_BDM2_2D, FE_N_T_BDM2_2D): return mt::N3_2_N3_2_2D;
    case switch_pair(FE_N_T_BDM3_2D, FE_N_T_BDM3_2D): return mt::N4_2_N4_2_2D;

    case switch_pair(FE_N_Q_RT0_2D, FE_N_Q_RT0_2D): return mt::N1_2_N1_2_2D;
    case switch_pair(FE_N_Q_RT0_2D, FE_N_T_RT0_2D): return mt::N1_2_N1_2_2D;
    case switch_pair(FE_N_Q_RT1_2D, FE_N_Q_RT1_2D): return mt::N2_2_N2_2_2D;
    case switch_pair(FE_N_Q_RT2_2D, FE_N_Q_RT2_2D): return mt::N3_2_N3_2_2D;
    case switch_pair(FE_N_Q_RT3_2D, FE_N_Q_RT3_2D): return mt::N4_2_N4_2_2D;

    case switch_pair(FE_N_Q_BDM1_2D, FE_N_Q_BDM1_2D): return mt::N2_2_N2_2_2D;
    case switch_pair(FE_N_Q_BDM2_2D, FE_N_Q_BDM2_2D): return mt::N3_2_N3_2_2D;
    case switch_pair(FE_N_Q_BDM3_2D, FE_N_Q_BDM3_2D): return mt::N4_2_N4_2_2D;

    case switch_pair(FE_N_T_P1_2D, FE_N_T_P1_2D): return mt::N1_2_N1_2_2D;
    case switch_pair(FE_N_T_P1_2D, FE_N_Q_Q1_2D): return mt::N1_2_N1_2_2D;
    case switch_pair(FE_N_Q_Q1_2D, FE_N_T_P1_2D): return mt::N1_2_N1_2_2D;
    case switch_pair(FE_N_Q_Q1_2D, FE_N_Q_Q1_2D): return mt::N1_2_N1_2_2D;

    case switch_pair(FE_D_Q_P1_2D, FE_D_Q_P1_2D): return mt::C0_2_C0_2_2D;
    case switch_pair(FE_D_T_P1_2D, FE_D_T_P1_2D): return mt::C0_2_C0_2_2D;
    case switch_pair(FE_D_T_P1_2D, FE_D_Q_P1_2D): return mt::C0_2_C0_2_2D;
    case switch_pair(FE_D_Q_P1_2D, FE_D_T_P1_2D): return mt::C0_2_C0_2_2D;

    case switch_pair(FE_D_Q_P2_2D, FE_D_Q_P2_2D): return mt::C0_2_C0_2_2D;
    case switch_pair(FE_D_T_P2_2D, FE_D_T_P2_2D): return mt::C0_2_C0_2_2D;
    case switch_pair(FE_D_T_P2_2D, FE_D_Q_P2_2D): return mt::C0_2_C0_2_2D;
    case switch_pair(FE_D_Q_P2_2D, FE_D_T_P2_2D): return mt::C0_2_C0_2_2D;

    case switch_pair(FE_D_T_P2_2D, FE_D_Q_Q2_2D): return mt::C0_2_C0_2_2D;
    case switch_pair(FE_D_Q_Q2_2D, FE_D_T_P2_2D): return mt::C0_2_C0_2_2D;

    case switch_pair(FE_D_Q_P3_2D, FE_D_Q_P3_2D): return mt::C0_2_C0_2_2D;
    case switch_pair(FE_D_T_P3_2D, FE_D_T_P3_2D): return mt::C0_2_C0_2_2D;
    case switch_pair(FE_D_T_P3_2D, FE_D_Q_P3_2D): return mt::C0_2_C0_2_2D;
    case switch_pair(FE_D_Q_P3_2D, FE_D_T_P3_2D): return mt::C0_2_C0_2_2D;

    case switch_pair(FE_D_Q_P4_2D, FE_D_Q_P4_2D): return mt::C0_2_C0_2_2D;
    case switch_pair(FE_D_T_P4_2D, FE_D_T_P4_2D): return mt::C0_2_C0_2_2D;
    case switch_pair(FE_D_T_P4_2D, FE_D_Q_P4_2D): return mt::C0_2_C0_2_2D;
    case switch_pair(FE_D_Q_P4_2D, FE_D_T_P4_2D): return mt::C0_2_C0_2_2D;

    case switch_pair(FE_D_Q_P5_2D, FE_D_Q_P5_2D): return mt::C0_2_C0_2_2D;
    case switch_pair(FE_D_Q_P6_2D, FE_D_Q_P6_2D): return mt::C0_2_C0_2_2D;
    case switch_pair(FE_D_Q_P7_2D, FE_D_Q_P7_2D): return mt::C0_2_C0_2_2D;

    case switch_pair(FE_D_Q_Q1_2D, FE_D_Q_Q1_2D): return mt::C0_2_C0_2_2D;
    case switch_pair(FE_D_Q_Q2_2D, FE_D_Q_Q2_2D): return mt::C0_2_C0_2_2D;
    case switch_pair(FE_D_Q_Q3_2D, FE_D_Q_Q3_2D): return mt::C0_2_C0_2_2D;
    case switch_pair(FE_D_Q_Q4_2D, FE_D_Q_Q4_2D): return mt::C0_2_C0_2_2D;

    case switch_pair(FE_N_T_P1MOD_2D, FE_N_T_P1MOD_2D): return mt::N2_2_N2_2_2D;

    case switch_pair(FE_C_T_P1MINI_2D, FE_C_T_P1MINI_2D): return mt::C1_2_C1_2_2D;

    case switch_pair(FE_N_Q_Q2_2D, FE_N_Q_Q2_2D): return mt::N2_2_N2_2_2D;
    case switch_pair(FE_N_Q_Q2_2D, FE_N_T_P2_2D): return mt::N2_2_N2_2_2D;
    case switch_pair(FE_N_T_P2_2D, FE_N_Q_Q2_2D): return mt::N2_2_N2_2_2D;
    case switch_pair(FE_N_T_P2_2D, FE_N_T_P2_2D): return mt::N2_2_N2_2_2D;

    case switch_pair(FE_N_Q_Q3_2D, FE_N_Q_Q3_2D): return mt::N3_2_N3_2_2D;
    case switch_pair(FE_N_Q_Q3_2D, FE_N_T_P3_2D): return mt::N3_2_N3_2_2D;
    case switch_pair(FE_N_T_P3_2D, FE_N_Q_Q3_2D): return mt::N3_2_N3_2_2D;
    case switch_pair(FE_N_T_P3_2D, FE_N_T_P3_2D): return mt::N3_2_N3_2_2D;

    case switch_pair(FE_N_Q_Q4_2D, FE_N_Q_Q4_2D): return mt::N4_2_N4_2_2D;
    case switch_pair(FE_N_Q_Q4_2D, FE_N_T_P4_2D): return mt::N4_2_N4_2_2D;
    case switch_pair(FE_N_T_P4_2D, FE_N_Q_Q4_2D): return mt::N4_2_N4_2_2D;
    case switch_pair(FE_N_T_P4_2D, FE_N_T_P4_2D): return mt::N4_2_N4_2_2D;

    case switch_pair(FE_N_Q_Q5_2D, FE_N_Q_Q5_2D): return mt::N5_2_N5_2_2D;
    case switch_pair(FE_N_Q_Q5_2D, FE_N_T_P5_2D): return mt::N5_2_N5_2_2D;
    case switch_pair(FE_N_T_P5_2D, FE_N_Q_Q5_2D): return mt::N5_2_N5_2_2D;
    case switch_pair(FE_N_T_P5_2D, FE_N_T_P5_2D): return mt::N5_2_N5_2_2D;


    case switch_pair(FE_B_Q_IB2_2D, FE_B_Q_IB2_2D): return mt::C0_2_C0_2_2D;

    case switch_pair(FE_D_Q_D2_2D, FE_D_Q_D2_2D): return mt::C0_2_C0_2_2D;
    //========LOCALPROJECTION=============
    case switch_pair(FE_C_Q_UL1_2D, FE_C_Q_Q1_2D): return mt::C1_2_C1_2_2D;
    case switch_pair(FE_C_Q_UL2_2D, FE_C_Q_Q2_2D): return mt::C2_2_C2_2_2D;
    case switch_pair(FE_C_Q_UL3_2D, FE_C_Q_Q3_2D): return mt::C3_2_C3_2_2D;
    case switch_pair(FE_C_Q_UL4_2D, FE_C_Q_Q4_2D): return mt::C4_2_C4_2_2D;
    case switch_pair(FE_C_Q_UL5_2D, FE_C_Q_Q5_2D): return mt::C5_2_C5_2_2D;

    case switch_pair(FE_C_Q_Q1_2D, FE_C_Q_UL1_2D): return mt::C1_2_C1_2_2D;
    case switch_pair(FE_C_Q_Q2_2D, FE_C_Q_UL2_2D): return mt::C2_2_C2_2_2D;
    case switch_pair(FE_C_Q_Q3_2D, FE_C_Q_UL3_2D): return mt::C3_2_C3_2_2D;
    case switch_pair(FE_C_Q_Q4_2D, FE_C_Q_UL4_2D): return mt::C4_2_C4_2_2D;
    case switch_pair(FE_C_Q_Q5_2D, FE_C_Q_UL5_2D): return mt::C5_2_C5_2_2D;

    case switch_pair(FE_C_Q_UL1_2D, FE_C_Q_UL1_2D): return mt::C1_2_C1_2_2D;
    case switch_pair(FE_C_Q_UL2_2D, FE_C_Q_UL2_2D): return mt::C2_2_C2_2_2D;
    case switch_pair(FE_C_Q_UL3_2D, FE_C_Q_UL3_2D): return mt::C3_2_C3_2_2D;
    case switch_pair(FE_C_Q_UL4_2D, FE_C_Q_UL4_2D): return mt::C4_2_C4_2_2D;
    case switch_pair(FE_C_Q_UL5_2D, FE_C_Q_UL5_2D): return mt::C5_2_C5_2_2D;

    case switch_pair(FE_C_T_UL1_2D, FE_C_Q_UL1_2D): return mt::C1_2_C1_2_2D;
    case switch_pair(FE_C_T_UL2_2D, FE_C_Q_UL2_2D): return mt::C2_2_C2_2_2D;
    case switch_pair(FE_C_T_UL3_2D, FE_C_Q_UL3_2D): return mt::C3_2_C3_2_2D;
    case switch_pair(FE_C_T_UL4_2D, FE_C_Q_UL4_2D): return mt::C4_2_C4_2_2D;
    case switch_pair(FE_C_T_UL5_2D, FE_C_Q_UL5_2D): return mt::C5_2_C5_2_2D;

    case switch_pair(FE_C_Q_UL1_2D, FE_C_T_UL1_2D): return mt::C1_2_C1_2_2D;
    case switch_pair(FE_C_Q_UL2_2D, FE_C_T_UL2_2D): return mt::C2_2_C2_2_2D;
    case switch_pair(FE_C_Q_UL3_2D, FE_C_T_UL3_2D): return mt::C3_2_C3_2_2D;
    case switch_pair(FE_C_Q_UL4_2D, FE_C_T_UL4_2D): return mt::C4_2_C4_2_2D;
    case switch_pair(FE_C_Q_UL5_2D, FE_C_T_UL5_2D): return mt::C5_2_C5_2_2D;

    case switch_pair(FE_C_T_UL1_2D, FE_C_T_UL1_2D): return mt::C1_2_C1_2_2D;
    case switch_pair(FE_C_T_UL2_2D, FE_C_T_UL2_2D): return mt::C2_2_C2_2_2D;
    case switch_pair(FE_C_T_UL3_2D, FE_C_T_UL3_2D): return mt::C3_2_C3_2_2D;
    case switch_pair(FE_C_T_UL4_2D, FE_C_T_UL4_2D): return mt::C4_2_C4_2_2D;
    case switch_pair(FE_C_T_UL5_2D, FE_C_T_UL5_2D): return mt::C5_2_C5_2_2D;

    case switch_pair(FE_C_Q_UL2S_2D, FE_C_Q_UL2S_2D): return mt::C2_2_C2_2_2D;
    case switch_pair(FE_C_Q_UL3S_2D, FE_C_Q_UL3S_2D): return mt::C3_2_C3_2_2D;
    case switch_pair(FE_C_Q_UL4S_2D, FE_C_Q_UL4S_2D): return mt::C4_2_C4_2_2D;
    case switch_pair(FE_C_Q_UL5S_2D, FE_C_Q_UL5S_2D): return mt::C5_2_C5_2_2D;
    case switch_pair(FE_C_Q_UL6S_2D, FE_C_Q_UL6S_2D): return mt::C6_2_C6_2_2D;
    case switch_pair(FE_C_Q_UL7S_2D, FE_C_Q_UL7S_2D): return mt::C7_2_C7_2_2D;
    case switch_pair(FE_C_Q_UL8S_2D, FE_C_Q_UL8S_2D): return mt::C8_2_C8_2_2D;
    case switch_pair(FE_C_Q_UL9S_2D, FE_C_Q_UL9S_2D): return mt::C9_2_C9_2_2D;

    case switch_pair(FE_C_Q_UL2SE_2D, FE_C_Q_UL2SE_2D): return mt::C2_2_C2_2_2D;
    case switch_pair(FE_C_Q_UL3SE_2D, FE_C_Q_UL3SE_2D): return mt::C3_2_C3_2_2D;
    case switch_pair(FE_C_Q_UL4SE_2D, FE_C_Q_UL4SE_2D): return mt::C4_2_C4_2_2D;
    case switch_pair(FE_C_Q_UL5SE_2D, FE_C_Q_UL5SE_2D): return mt::C5_2_C5_2_2D;
    case switch_pair(FE_C_Q_UL6SE_2D, FE_C_Q_UL6SE_2D): return mt::C6_2_C6_2_2D;
    case switch_pair(FE_C_Q_UL7SE_2D, FE_C_Q_UL7SE_2D): return mt::C7_2_C7_2_2D;
    case switch_pair(FE_C_Q_UL8SE_2D, FE_C_Q_UL8SE_2D): return mt::C8_2_C8_2_2D;
    case switch_pair(FE_C_Q_UL9SE_2D, FE_C_Q_UL9SE_2D): return mt::C9_2_C9_2_2D;

    case switch_pair(FE_C_Q_M2_2D, FE_C_Q_M2_2D): return mt::C2_2_C2_2_2D;
    case switch_pair(FE_C_Q_M3_2D, FE_C_Q_M3_2D): return mt::C3_2_C3_2_2D;
    case switch_pair(FE_C_Q_M4_2D, FE_C_Q_M4_2D): return mt::C4_2_C4_2_2D;
    case switch_pair(FE_C_Q_M5_2D, FE_C_Q_M5_2D): return mt::C5_2_C5_2_2D;
    case switch_pair(FE_C_Q_M6_2D, FE_C_Q_M6_2D): return mt::C6_2_C6_2_2D;
    case switch_pair(FE_C_Q_M7_2D, FE_C_Q_M7_2D): return mt::C7_2_C7_2_2D;
    case switch_pair(FE_C_Q_M8_2D, FE_C_Q_M8_2D): return mt::C8_2_C8_2_2D;
    case switch_pair(FE_C_Q_M9_2D, FE_C_Q_M9_2D): return mt::C9_2_C9_2_2D;

    case switch_pair(FE_C_T_UL2_2D, FE_C_Q_UL2S_2D): return mt::C2_2_C2_2_2D;
    case switch_pair(FE_C_T_UL3_2D, FE_C_Q_UL3S_2D): return mt::C3_2_C3_2_2D;
    case switch_pair(FE_C_T_UL4_2D, FE_C_Q_UL4S_2D): return mt::C4_2_C4_2_2D;
    case switch_pair(FE_C_T_UL5_2D, FE_C_Q_UL5S_2D): return mt::C5_2_C5_2_2D;

    case switch_pair(FE_C_Q_UL2S_2D, FE_C_T_UL2_2D): return mt::C2_2_C2_2_2D;
    case switch_pair(FE_C_Q_UL3S_2D, FE_C_T_UL3_2D): return mt::C3_2_C3_2_2D;
    case switch_pair(FE_C_Q_UL4S_2D, FE_C_T_UL4_2D): return mt::C4_2_C4_2_2D;
    case switch_pair(FE_C_Q_UL5S_2D, FE_C_T_UL5_2D): return mt::C5_2_C5_2_2D;

    case switch_pair(FE_C_Q_EL1_2D, FE_C_Q_EL1_2D): return mt::C1_2_C1_2_2D;
    case switch_pair(FE_C_Q_UL1_2D, FE_C_Q_EL1_2D): return mt::C1_2_C1_2_2D;
    case switch_pair(FE_C_Q_EL1_2D, FE_C_Q_UL1_2D): return mt::C1_2_C1_2_2D;
    case switch_pair(FE_C_Q_EL1_2D, FE_C_Q_Q1_2D): return mt::C1_2_C1_2_2D;
    case switch_pair(FE_C_Q_Q1_2D, FE_C_Q_EL1_2D): return mt::C1_2_C1_2_2D;
// ========================================================================
// 3D
// ========================================================================

    case switch_pair(FE_C_T_P00_3D, FE_C_T_P00_3D): return mt::D_D_3D;
    case switch_pair(FE_C_T_P0_3D, FE_C_T_P0_3D): return mt::D_D_3D;
    case switch_pair(FE_C_T_P1_3D, FE_C_T_P1_3D): return mt::P1_P1_3D;
    case switch_pair(FE_C_T_P2_3D, FE_C_T_P2_3D): return mt::P2_P2_3D;
    case switch_pair(FE_C_T_P3_3D, FE_C_T_P3_3D): return mt::P3_P3_3D;

    case switch_pair(FE_C_T_B2_3D, FE_C_T_B2_3D): return mt::P2B_P2B_3D;

    case switch_pair(FE_C_H_Q00_3D, FE_C_H_Q00_3D): return mt::D_D_3D;
    case switch_pair(FE_C_H_Q0_3D, FE_C_H_Q0_3D): return mt::D_D_3D;
    case switch_pair(FE_C_H_Q1_3D, FE_C_H_Q1_3D): return mt::Q1_Q1_3D;
    case switch_pair(FE_C_H_Q2_3D, FE_C_H_Q2_3D): return mt::Q2_Q2_3D;
    case switch_pair(FE_C_H_Q3_3D, FE_C_H_Q3_3D): return mt::Q3_Q3_3D;
    case switch_pair(FE_C_H_Q4_3D, FE_C_H_Q4_3D): return mt::Q4_Q4_3D;

    case switch_pair(FE_N_H_Q1_3D, FE_N_H_Q1_3D): return mt::N1_N1_3D;
    case switch_pair(FE_N_T_P1_3D, FE_N_T_P1_3D): return mt::NP1_NP1_3D;

    case switch_pair(FE_D_T_P1_3D, FE_D_T_P1_3D): return mt::D_D_3D;
    case switch_pair(FE_D_T_P2_3D, FE_D_T_P2_3D): return mt::D_D_3D;
    case switch_pair(FE_D_T_P3_3D, FE_D_T_P3_3D): return mt::D_D_3D;

    case switch_pair(FE_D_H_P1_3D, FE_D_H_P1_3D): return mt::D_D_3D;
    case switch_pair(FE_D_H_P2_3D, FE_D_H_P2_3D): return mt::D_D_3D;
    case switch_pair(FE_D_H_P3_3D, FE_D_H_P3_3D): return mt::D_D_3D;
    
    case switch_pair(FE_D_H_Q1_3D, FE_D_H_Q1_3D): return mt::D_D_3D;
    case switch_pair(FE_D_H_Q2_3D, FE_D_H_Q2_3D): return mt::D_D_3D;

    case switch_pair(FE_B_H_IB2_3D, FE_B_H_IB2_3D): return mt::D_D_3D;

    case switch_pair(FE_N_T_P2_3D, FE_N_T_P2_3D): return mt::P1_P1_3D;
    
    case switch_pair(FE_N_H_Q2_3D, FE_N_H_Q2_3D): return mt::N2_N2_3D;
    case switch_pair(FE_N_H_Q3_3D, FE_N_H_Q3_3D): return mt::N3_N3_3D;
    case switch_pair(FE_N_H_Q4_3D, FE_N_H_Q4_3D): return mt::N4_N4_3D;
    
    //========LOCALPROJECTION==============
    case switch_pair(FE_C_H_UL1_3D, FE_C_H_Q1_3D): return mt::Q1_Q1_3D;
    case switch_pair(FE_C_H_UL2_3D, FE_C_H_Q2_3D): return mt::Q2_Q2_3D;
    case switch_pair(FE_C_H_Q1_3D, FE_C_H_UL1_3D): return mt::Q1_Q1_3D;
    case switch_pair(FE_C_H_Q2_3D, FE_C_H_UL2_3D): return mt::Q2_Q2_3D;
    case switch_pair(FE_C_H_UL1_3D, FE_C_H_UL1_3D): return mt::Q1_Q1_3D;
    case switch_pair(FE_C_H_UL2_3D, FE_C_H_UL2_3D): return mt::Q2_Q2_3D;

    case switch_pair(FE_C_H_UL3_3D, FE_C_H_Q3_3D): return mt::Q3_Q3_3D;
    case switch_pair(FE_C_H_Q3_3D, FE_C_H_UL3_3D): return mt::Q3_Q3_3D;
    case switch_pair(FE_C_H_UL3_3D, FE_C_H_UL3_3D): return mt::Q3_Q3_3D;
    //=====================================

    // ========================================================================
    // discontinuous finite element spaces
    // ========================================================================
    case switch_pair(FE_C_H_Q0_3D, FE_D_H_P1_3D): return mt::D_D_3D;
    case switch_pair(FE_D_H_P1_3D, FE_C_H_Q0_3D): return mt::D_D_3D;
    case switch_pair(FE_C_H_Q0_3D, FE_D_H_P2_3D): return mt::D_D_3D;
    case switch_pair(FE_D_H_P2_3D, FE_C_H_Q0_3D): return mt::D_D_3D;
    case switch_pair(FE_C_H_Q0_3D, FE_D_H_P3_3D): return mt::D_D_3D;
    case switch_pair(FE_D_H_P3_3D, FE_C_H_Q0_3D): return mt::D_D_3D;
    case switch_pair(FE_D_H_P1_3D, FE_D_H_P2_3D): return mt::D_D_3D;
    case switch_pair(FE_D_H_P2_3D, FE_D_H_P1_3D): return mt::D_D_3D;
    case switch_pair(FE_D_H_P1_3D, FE_D_H_P3_3D): return mt::D_D_3D;
    case switch_pair(FE_D_H_P3_3D, FE_D_H_P1_3D): return mt::D_D_3D;
    case switch_pair(FE_D_H_P2_3D, FE_D_H_P3_3D): return mt::D_D_3D;
    case switch_pair(FE_D_H_P3_3D, FE_D_H_P2_3D): return mt::D_D_3D;
    case switch_pair(FE_C_H_Q00_3D, FE_D_H_P1_3D): return mt::D_D_3D;
    case switch_pair(FE_D_H_P1_3D, FE_C_H_Q00_3D): return mt::D_D_3D;
    case switch_pair(FE_C_H_Q00_3D, FE_D_H_P2_3D): return mt::D_D_3D;
    case switch_pair(FE_D_H_P2_3D, FE_C_H_Q00_3D): return mt::D_D_3D;
    case switch_pair(FE_C_H_Q00_3D, FE_D_H_P3_3D): return mt::D_D_3D;
    case switch_pair(FE_D_H_P3_3D, FE_C_H_Q00_3D): return mt::D_D_3D;
    case switch_pair(FE_C_H_Q00_3D, FE_C_H_Q0_3D): return mt::D_D_3D;
    case switch_pair(FE_C_H_Q0_3D, FE_C_H_Q00_3D): return mt::D_D_3D;
    case switch_pair(FE_C_T_P0_3D, FE_C_T_P00_3D): return mt::D_D_3D;
    case switch_pair(FE_C_T_P00_3D, FE_C_T_P0_3D): return mt::D_D_3D;
    case switch_pair(FE_C_T_P00_3D, FE_D_T_P1_3D): return mt::D_D_3D;
    case switch_pair(FE_D_T_P1_3D, FE_C_T_P00_3D): return mt::D_D_3D;
    case switch_pair(FE_C_T_P0_3D, FE_D_T_P1_3D): return mt::D_D_3D;
    case switch_pair(FE_D_T_P1_3D, FE_C_T_P0_3D): return mt::D_D_3D;
    
    case switch_pair(FE_N_T_RT0_3D, FE_N_T_RT0_3D): return mt::NP1_NP1_3D;
    case switch_pair(FE_N_T_RT1_3D, FE_N_T_RT1_3D): return mt::NP2_NP2_3D;
    case switch_pair(FE_N_T_RT2_3D, FE_N_T_RT2_3D): return mt::P2_P2_3D;
    case switch_pair(FE_N_T_RT3_3D, FE_N_T_RT3_3D): return mt::P3_P3_3D;
    case switch_pair(FE_N_T_BDDF1_3D, FE_N_T_BDDF1_3D): return mt::NP2_NP2_3D;
    case switch_pair(FE_N_T_BDDF2_3D, FE_N_T_BDDF2_3D): return mt::P2_P2_3D;
    case switch_pair(FE_N_T_BDDF3_3D, FE_N_T_BDDF3_3D): return mt::P3_P3_3D;

    case switch_pair(FE_N_H_RT0_3D, FE_N_H_RT0_3D): return mt::N1_N1_3D;
    case switch_pair(FE_N_H_RT1_3D, FE_N_H_RT1_3D): return mt::Q1_Q1_3D;
    case switch_pair(FE_N_H_RT2_3D, FE_N_H_RT2_3D): return mt::Q2_Q2_3D;
    
    case switch_pair(FE_N_H_BDDF1_3D, FE_N_H_BDDF1_3D): return mt::N2_N2_3D;
    case switch_pair(FE_N_H_BDDF2_3D, FE_N_H_BDDF2_3D): return mt::N3_N3_3D;
    case switch_pair(FE_N_H_BDDF3_3D, FE_N_H_BDDF3_3D): return mt::N4_N4_3D;
    default:
      ErrThrow("There is no mapper type assigned to the 3D finite elements ",
               FE1, " and ", FE2);
  }
}

FEMapper_type get_mapper_type_1regular2d(FEDescriptor_type FE1,
                                         FEDescriptor_type FE2)
{
  using mt = FEMapper_type;
  switch(switch_pair(FE1, FE2))
  {
    // ========================================================================
    // ONE regular grid, same pattern on both sides
    // ========================================================================
    case switch_pair(FE_C_T_P0_2D, FE_C_T_P0_2D): return mt::C0_2_C0_2_1Reg_2D;
    case switch_pair(FE_C_T_P0_2D, FE_C_Q_Q0_2D): return mt::C0_2_C0_2_1Reg_2D;
    case switch_pair(FE_C_Q_Q0_2D, FE_C_T_P0_2D): return mt::C0_2_C0_2_1Reg_2D;
    case switch_pair(FE_C_Q_Q0_2D, FE_C_Q_Q0_2D): return mt::C0_2_C0_2_1Reg_2D;

    case switch_pair(FE_C_T_P1_2D, FE_C_T_P1_2D): return mt::C1_2_C1_2_1Reg_2D;
    case switch_pair(FE_C_T_P1_2D, FE_C_Q_Q1_2D): return mt::C1_2_C1_2_1Reg_2D;
    case switch_pair(FE_C_Q_Q1_2D, FE_C_T_P1_2D): return mt::C1_2_C1_2_1Reg_2D;
    case switch_pair(FE_C_Q_Q1_2D, FE_C_Q_Q1_2D): return mt::C1_2_C1_2_1Reg_2D;

    case switch_pair(FE_C_T_P2_2D, FE_C_T_P2_2D): return mt::C2_2_C2_2_1Reg_2D;
    case switch_pair(FE_C_T_P2_2D, FE_C_Q_Q2_2D): return mt::C2_2_C2_2_1Reg_2D;
    case switch_pair(FE_C_Q_Q2_2D, FE_C_T_P2_2D): return mt::C2_2_C2_2_1Reg_2D;
    case switch_pair(FE_C_Q_Q2_2D, FE_C_Q_Q2_2D): return mt::C2_2_C2_2_1Reg_2D;

    case switch_pair(FE_C_T_B2_2D, FE_C_T_P2_2D): return mt::C2_2_C2_2_1Reg_2D;
    case switch_pair(FE_C_T_B2_2D, FE_C_Q_Q2_2D): return mt::C2_2_C2_2_1Reg_2D;
    case switch_pair(FE_C_T_P2_2D, FE_C_T_B2_2D): return mt::C2_2_C2_2_1Reg_2D;
    case switch_pair(FE_C_Q_Q2_2D, FE_C_T_B2_2D): return mt::C2_2_C2_2_1Reg_2D;
    case switch_pair(FE_C_T_B2_2D, FE_C_T_B2_2D): return mt::C2_2_C2_2_1Reg_2D;

    case switch_pair(FE_C_T_SV2_2D, FE_C_T_P2_2D): return mt::C2_2_C2_2_1Reg_2D;
    case switch_pair(FE_C_T_SV2_2D, FE_C_Q_Q2_2D): return mt::C2_2_C2_2_1Reg_2D;
    case switch_pair(FE_C_T_P2_2D, FE_C_T_SV2_2D): return mt::C2_2_C2_2_1Reg_2D;
    case switch_pair(FE_C_Q_Q2_2D, FE_C_T_SV2_2D): return mt::C2_2_C2_2_1Reg_2D;
    case switch_pair(FE_C_T_SV2_2D, FE_C_T_SV2_2D): return mt::C2_2_C2_2_1Reg_2D;

    case switch_pair(FE_C_T_P3_2D, FE_C_T_P3_2D): return mt::C3_2_C3_2_1Reg_2D;
    case switch_pair(FE_C_T_P3_2D, FE_C_Q_Q3_2D): return mt::C3_2_C3_2_1Reg_2D;
    case switch_pair(FE_C_Q_Q3_2D, FE_C_T_P3_2D): return mt::C3_2_C3_2_1Reg_2D;
    case switch_pair(FE_C_Q_Q3_2D, FE_C_Q_Q3_2D): return mt::C3_2_C3_2_1Reg_2D;

    case switch_pair(FE_C_T_B3_2D, FE_C_T_P3_2D): return mt::C3_2_C3_2_1Reg_2D;
    case switch_pair(FE_C_T_B3_2D, FE_C_Q_Q3_2D): return mt::C3_2_C3_2_1Reg_2D;
    case switch_pair(FE_C_T_P3_2D, FE_C_T_B3_2D): return mt::C3_2_C3_2_1Reg_2D;
    case switch_pair(FE_C_Q_Q3_2D, FE_C_T_B3_2D): return mt::C3_2_C3_2_1Reg_2D;
    case switch_pair(FE_C_T_B3_2D, FE_C_T_B3_2D): return mt::C3_2_C3_2_1Reg_2D;

    // fourth order
    case switch_pair(FE_C_T_P4_2D, FE_C_T_P4_2D): return mt::C4_2_C4_2_1Reg_2D;
    case switch_pair(FE_C_T_P4_2D, FE_C_Q_Q4_2D): return mt::C4_2_C4_2_1Reg_2D;
    case switch_pair(FE_C_Q_Q4_2D, FE_C_T_P4_2D): return mt::C4_2_C4_2_1Reg_2D;
    case switch_pair(FE_C_Q_Q4_2D, FE_C_Q_Q4_2D): return mt::C4_2_C4_2_1Reg_2D;

    case switch_pair(FE_C_T_B4_2D, FE_C_T_P4_2D): return mt::C4_2_C4_2_1Reg_2D;
    case switch_pair(FE_C_T_B4_2D, FE_C_Q_Q4_2D): return mt::C4_2_C4_2_1Reg_2D;
    case switch_pair(FE_C_T_P4_2D, FE_C_T_B4_2D): return mt::C4_2_C4_2_1Reg_2D;
    case switch_pair(FE_C_Q_Q4_2D, FE_C_T_B4_2D): return mt::C4_2_C4_2_1Reg_2D;
    case switch_pair(FE_C_T_B4_2D, FE_C_T_B4_2D): return mt::C4_2_C4_2_1Reg_2D;

    // fifth order
    case switch_pair(FE_C_T_P5_2D, FE_C_T_P5_2D): return mt::C5_2_C5_2_1Reg_2D;
    case switch_pair(FE_C_T_P5_2D, FE_C_Q_Q5_2D): return mt::C5_2_C5_2_1Reg_2D;
    case switch_pair(FE_C_Q_Q5_2D, FE_C_T_P5_2D): return mt::C5_2_C5_2_1Reg_2D;
    case switch_pair(FE_C_Q_Q5_2D, FE_C_Q_Q5_2D): return mt::C5_2_C5_2_1Reg_2D;

    case switch_pair(FE_N_T_P1_2D, FE_N_T_P1_2D): return mt::N1_2_N1_2_1Reg_2D;
    case switch_pair(FE_N_T_P1_2D, FE_N_Q_Q1_2D): return mt::N1_2_N1_2_1Reg_2D;
    case switch_pair(FE_N_Q_Q1_2D, FE_N_T_P1_2D): return mt::N1_2_N1_2_1Reg_2D;
    case switch_pair(FE_N_Q_Q1_2D, FE_N_Q_Q1_2D): return mt::N1_2_N1_2_1Reg_2D;

    // second order nonconforming
    case switch_pair(FE_N_Q_Q2_2D, FE_N_Q_Q2_2D): return mt::N2_2_N2_2_1Reg_2D;
    case switch_pair(FE_N_T_P2_2D, FE_N_T_P2_2D): return mt::N2_2_N2_2_1Reg_2D;
    case switch_pair(FE_N_Q_Q2_2D, FE_N_T_P2_2D): return mt::N2_2_N2_2_1Reg_2D;
    case switch_pair(FE_N_T_P2_2D, FE_N_Q_Q2_2D): return mt::N2_2_N2_2_1Reg_2D;

    // third order nonconforming
    case switch_pair(FE_N_Q_Q3_2D, FE_N_Q_Q3_2D): return mt::N3_2_N3_2_1Reg_2D;
    case switch_pair(FE_N_T_P3_2D, FE_N_T_P3_2D): return mt::N3_2_N3_2_1Reg_2D;
    case switch_pair(FE_N_Q_Q3_2D, FE_N_T_P3_2D): return mt::N3_2_N3_2_1Reg_2D;
    case switch_pair(FE_N_T_P3_2D, FE_N_Q_Q3_2D): return mt::N3_2_N3_2_1Reg_2D;

    // fourth order nonconforming
    case switch_pair(FE_N_Q_Q4_2D, FE_N_Q_Q4_2D): return mt::N4_2_N4_2_1Reg_2D;
    case switch_pair(FE_N_T_P4_2D, FE_N_T_P4_2D): return mt::N4_2_N4_2_1Reg_2D;
    case switch_pair(FE_N_Q_Q4_2D, FE_N_T_P4_2D): return mt::N4_2_N4_2_1Reg_2D;
    case switch_pair(FE_N_T_P4_2D, FE_N_Q_Q4_2D): return mt::N4_2_N4_2_1Reg_2D;

    // fifth order nonconforming
    case switch_pair(FE_N_Q_Q5_2D, FE_N_Q_Q5_2D): return mt::N5_2_N5_2_1Reg_2D;
    case switch_pair(FE_N_T_P5_2D, FE_N_T_P5_2D): return mt::N5_2_N5_2_1Reg_2D;
    case switch_pair(FE_N_Q_Q5_2D, FE_N_T_P5_2D): return mt::N5_2_N5_2_1Reg_2D;
    case switch_pair(FE_N_T_P5_2D, FE_N_Q_Q5_2D): return mt::N5_2_N5_2_1Reg_2D;

    case switch_pair(FE_D_Q_P1_2D, FE_D_Q_P1_2D): return mt::C0_2_C0_2_1Reg_2D;
    case switch_pair(FE_D_T_P1_2D, FE_D_T_P1_2D): return mt::C0_2_C0_2_1Reg_2D;
    case switch_pair(FE_D_T_P1_2D, FE_D_Q_P1_2D): return mt::C0_2_C0_2_1Reg_2D;
    case switch_pair(FE_D_Q_P1_2D, FE_D_T_P1_2D): return mt::C0_2_C0_2_1Reg_2D;
    
    case switch_pair(FE_D_T_SV1_2D, FE_D_T_SV1_2D): return mt::C0_2_C0_2_1Reg_2D;
    case switch_pair(FE_D_T_SV1_2D, FE_D_Q_P1_2D): return mt::C0_2_C0_2_1Reg_2D;
    case switch_pair(FE_D_Q_P1_2D, FE_D_T_SV1_2D): return mt::C0_2_C0_2_1Reg_2D;
    case switch_pair(FE_D_Q_P2_2D, FE_D_Q_P2_2D): return mt::C0_2_C0_2_1Reg_2D;
    case switch_pair(FE_D_T_P2_2D, FE_D_T_P2_2D): return mt::C0_2_C0_2_1Reg_2D;
    case switch_pair(FE_D_T_P2_2D, FE_D_Q_P2_2D): return mt::C0_2_C0_2_1Reg_2D;
    case switch_pair(FE_D_Q_P2_2D, FE_D_T_P2_2D): return mt::C0_2_C0_2_1Reg_2D;

    // ========================================================================
    // ONE regular grid, different pattern
    // ========================================================================
    // coarse: second order, fine: first order
    case switch_pair(FE_N_Q_Q2_2D, FE_N_Q_Q1_2D): return mt::N2_2_N1_2_1Reg_2D;
    case switch_pair(FE_N_T_P2_2D, FE_N_T_P1_2D): return mt::N2_2_N1_2_1Reg_2D;
    case switch_pair(FE_N_Q_Q2_2D, FE_N_T_P1_2D): return mt::N2_2_N1_2_1Reg_2D;
    case switch_pair(FE_N_T_P2_2D, FE_N_Q_Q1_2D): return mt::N2_2_N1_2_1Reg_2D;

    // coarse: third order, fine: second order
    case switch_pair(FE_N_Q_Q3_2D, FE_N_Q_Q2_2D): return mt::N3_2_N2_2_1Reg_2D;
    case switch_pair(FE_N_T_P3_2D, FE_N_T_P2_2D): return mt::N3_2_N2_2_1Reg_2D;
    case switch_pair(FE_N_Q_Q3_2D, FE_N_T_P2_2D): return mt::N3_2_N2_2_1Reg_2D;
    case switch_pair(FE_N_T_P3_2D, FE_N_Q_Q2_2D): return mt::N3_2_N2_2_1Reg_2D;

    // coarse: fourth order, fine: third order
    case switch_pair(FE_N_Q_Q4_2D, FE_N_Q_Q3_2D): return mt::N4_2_N3_2_1Reg_2D;
    case switch_pair(FE_N_T_P4_2D, FE_N_T_P3_2D): return mt::N4_2_N3_2_1Reg_2D;
    case switch_pair(FE_N_Q_Q4_2D, FE_N_T_P3_2D): return mt::N4_2_N3_2_1Reg_2D;
    case switch_pair(FE_N_T_P4_2D, FE_N_Q_Q3_2D): return mt::N4_2_N3_2_1Reg_2D;

    // coarse: fifth order, fine: fourth order
    case switch_pair(FE_N_Q_Q5_2D, FE_N_Q_Q4_2D): return mt::N5_2_N4_2_1Reg_2D;
    case switch_pair(FE_N_T_P5_2D, FE_N_T_P4_2D): return mt::N5_2_N4_2_1Reg_2D;
    case switch_pair(FE_N_Q_Q5_2D, FE_N_T_P4_2D): return mt::N5_2_N4_2_1Reg_2D;
    case switch_pair(FE_N_T_P5_2D, FE_N_Q_Q4_2D): return mt::N5_2_N4_2_1Reg_2D;

    default:
      ErrThrow("There is no one-regular mapper type assigned to the 2D finite "
               "elements ", FE1, " and ", FE2);
      break;
  }
}

FEMapper_type get_mapper_type_1regular3d(FEDescriptor_type FE1,
                                         FEDescriptor_type FE2)
{
  using mt = FEMapper_type;
  switch(switch_pair(FE1, FE2))
  {
    case switch_pair(FE_C_T_P1_3D, FE_C_T_P1_3D): return mt::P1_P1_1Reg_3D;
    case switch_pair(FE_C_T_P2_3D, FE_C_T_P2_3D): return mt::P2_P2_1Reg_3D;
    case switch_pair(FE_C_T_P3_3D, FE_C_T_P3_3D): return mt::P3_P3_1Reg_3D;
    case switch_pair(FE_C_H_Q1_3D, FE_C_H_Q1_3D): return mt::Q1_Q1_1Reg_3D;
    case switch_pair(FE_C_H_Q2_3D, FE_C_H_Q2_3D): return mt::Q2_Q2_1Reg_3D;
    case switch_pair(FE_C_H_Q3_3D, FE_C_H_Q3_3D): return mt::Q3_Q3_1Reg_3D;

    case switch_pair(FE_N_T_P1_3D, FE_N_T_P1_3D): return mt::NP1_NP1_1Reg_3D;
    case switch_pair(FE_N_T_P2_3D, FE_N_T_P2_3D): return mt::NP2_NP2_1Reg_3D;
    
    case switch_pair(FE_C_H_UL1_3D, FE_C_H_UL1_3D): return mt::Q1_Q1_1Reg_3D;
    case switch_pair(FE_C_H_UL1_3D, FE_C_H_Q1_3D): return mt::Q1_Q1_1Reg_3D;
    case switch_pair(FE_C_H_Q1_3D, FE_C_H_UL1_3D): return mt::Q1_Q1_1Reg_3D;
    default:
      ErrThrow("There is no one-regular mapper type assigned to the 3D finite "
               "elements ", FE1, " and ", FE2);
      break;
  }
}
