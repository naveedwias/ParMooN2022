#ifndef __FE3DMAPPER1REG__
#define __FE3DMAPPER1REG__

#include "FEMapper.h"
#include <map>

/** find out which of the given local degress of freedom, are equivalent to the
 * same global degree of freedom for 1 regular case
 * 
 * @ruleof0
 */
class FEMapper1Reg3D : public FEMapper
{
  protected:
    /** rule how to permute the local dofs due to twist index */
    std::vector<std::vector<int>> twist_permutation;

    /** select which set of pairs will be used, due to MapType */
    const std::pair<int, int> *CurrentPairs;

    /** select which set of noopposite dof will be used, due to MapType */
    const int *CurrentNoOp;
    
    /** number of dof without an opposite */
    int N_NoOpposite;

    /** list of dof without an opposite */
    std::vector<std::vector<int>> no_opposite;

    /** number of hanging nodes */
    int N_Hanging;

    /** indices of hanging nodes */
    std::vector<int> Hanging;

    /** type of hanging nodes */
    std::vector<HNDesc> HangingTypes;

    /** numbers of the degrees of freedom in coupling */
    std::vector<std::vector<int>> coupling;

  public:

    /** @brief construct an FEMapper of the given type. 
     * 
     * Note that certain types are not supported here and need to be constructed
     * in the base class or in another derived class.
     */
    FEMapper1Reg3D(FEMapper_type t);

    /** @brief delete all members, no memory on heap is used.*/
    ~FEMapper1Reg3D() = default;
    /** @brief default copy constructor */
    FEMapper1Reg3D(const FEMapper1Reg3D&) = default;
    /** @brief default move constructor */
    FEMapper1Reg3D(FEMapper1Reg3D&&) = default;
    /** @brief default copy assignment */
    FEMapper1Reg3D & operator=(const FEMapper1Reg3D&) = default;
    /** @brief default move assignment */
    FEMapper1Reg3D& operator=(FEMapper1Reg3D&&) = default;

    /** map the given local degrees of freedom,
        coarse cell has lower number,
        0 is coarser side 0; 1,2 are on finer side 1 */
    void Map(int *Global,
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
             std::map<HNDesc, THNDesc *>& hn_descriptors) const;

};

#endif
