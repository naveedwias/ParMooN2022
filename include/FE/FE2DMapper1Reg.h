#ifndef __FE2DMAPPER1REG__
#define __FE2DMAPPER1REG__

#include "FEMapper.h"
#include <map>

/** find out which of the given local degress of freedom, are equivalent to the
 * same global degree of freedom for 1 regular case
 * 
 * @ruleof0
 */
class FEMapper1Reg2D : public FEMapper
{
  protected:
    /** number of local degrees of second element on second side */
    int N_DOF2;

    /** number of degrees of freedom at midpoint of long edge */
    int N_Mid;

    /** indices of degrees of freedom at midpoint of long edge (size: 2*N_Mid)*/
    std::vector<int> Mid;

    /** number of real degrees of freedom with no opposite */
    int N_NoOpposite;

    /** numbers of the eal degrees of freedom with no opposite */
    std::vector<int> NoOpposite;

// =======================================================================
// local hanging nodes = coupling only this node on this edge
// =======================================================================
    /** number of hanging nodes in this pattern */
    int N_Hanging;

    /** indices of hanging nodes */
    std::vector<int> Hanging;

    /** type of hanging nodes */
    std::vector<HNDesc> HangingTypes;

    /** numbers of the degrees of freedom in coupling */
    std::vector<std::vector<int>> coupling;

    /** @brief get a hanging node descriptor with a given type. If not found
     * in the map of descriptors it is created there. */
    THNDesc* get_hn_descriptor(HNDesc type,
                               std::map<HNDesc, THNDesc *>& hn_descriptors)
      const;

  public:

    /** @brief construct an FEMapper of the given type. 
     * 
     * Note that certain types are not supported here and need to be constructed
     * in the base class or in another derived class.
     */
    FEMapper1Reg2D(FEMapper_type t);

    /** @brief delete all members, no memory on heap is used.*/
    ~FEMapper1Reg2D() = default;
    /** @brief default copy constructor */
    FEMapper1Reg2D(const FEMapper1Reg2D&) = default;
    /** @brief default move constructor */
    FEMapper1Reg2D(FEMapper1Reg2D&&) = default;
    /** @brief default copy assignment */
    FEMapper1Reg2D & operator=(const FEMapper1Reg2D&) = default;
    /** @brief default move assignment */
    FEMapper1Reg2D& operator=(FEMapper1Reg2D&&) = default;

    /** map the given local degrees of freedom,
        if `coarse_fine` is true, coarse cell has lower number,
        0 is coarser side 0; 1,2 are on finer side 1 */
    void Map_1Reg(bool coarse_fine, int *Global, int I_K0, int I_K1, int I_K2,
             const int *Indices0, const int *Indices1, const int *Indices2,
             int &Counter, int LowerFirstChild,
             std::vector<THangingNode*>& vect,
             std::vector<int>& numbers,
             std::map<HNDesc, THNDesc *>& hn_descriptors) const;

};

#endif
