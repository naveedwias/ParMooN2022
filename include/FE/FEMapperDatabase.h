#include "Enumerations_fe.h"
#include "FEMapper.h"
#include <map>
#include <memory>

class FEMapper1Reg2D;
class FEMapper1Reg3D;

/**
 * @brief helper class only used during construction of finite element spaces
 * 
 * This essentially caches the objects of type FEMapper, FEMapper1Reg2D, and
 * FEMapper1Reg3D. That means it creates these objects upon request and stores
 * them.
 */
class FEMapperDatabase
{
  private:

    /** @brief store all the finite element mappers which are used */
    std::map<FEMapper_type, std::unique_ptr<FEMapper>> mappers;

    /** @brief get mapper of a given type
     * 
     * If the type is not yet stored in `FEMapperDatabase::mappers` it is 
     * created there.
     * 
     * @param type the unique type id
     * @param one_reg get the mapper for a one-regular setting (with a hanging
     *        vertex)
     * @param two_d the space dimension, which is only relevant if `one_reg` is
     *        true
     */
    const FEMapper* get_mapper(FEMapper_type type, bool one_reg, bool two_d);

  public:

    /** @brief the default constructor, does not create any mappers */
    FEMapperDatabase() = default;

    /** @brief destructor deleting all mappers */
    ~FEMapperDatabase() = default;

    /** @brief default copy constructor */
    FEMapperDatabase(const FEMapperDatabase&) = default;
    /** @brief default move constructor */
    FEMapperDatabase(FEMapperDatabase&&) = default;
    /** @brief default copy assignment */
    FEMapperDatabase & operator=(const FEMapperDatabase&) = default;
    /** @brief default move assignment */
    FEMapperDatabase& operator=(FEMapperDatabase&&) = default;

    /** @brief return FE mapper for a pair of finite elements
     * 
     * This does not work for one-regular mappers, for those see 
     * `get_fe_mapper_one_regular()`.
     */
    const FEMapper *get_fe_mapper(FEDescriptor_type FE1, FEDescriptor_type FE2);

    /** @brief return a one-regular fe mapper in 2D
     * 
     * This only works for one-regular mappers. For regular mappers, see 
     * `get_fe_mapper()`.
     */
    const FEMapper1Reg2D *get_fe_mapper_one_regular2d(FEDescriptor_type FE1,
                                                    FEDescriptor_type FE2);

    /** @brief return a one-regular fe mapper in 3D
     * 
     * This only works for one-regular mappers. For regular mappers, see 
     * `get_fe_mapper()`.
     */
    const FEMapper1Reg3D *get_fe_mapper_one_regular3d(FEDescriptor_type FE1,
                                                      FEDescriptor_type FE2);
};
