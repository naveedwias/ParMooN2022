#ifndef __FE3D__
#define __FE3D__

#include "Enumerations_fe.h"
#include "BaseFunctions.h"
#include "NodalFunctional.h"
#include "FEDescriptor.h"
#include <tuple>

/** @brief store all information for one finite element */
class FiniteElement
{
  protected:

    /** @brief type of finite element */
    FE_type type;

    /** @brief set of basis function */
    BaseFunctions base_function;

    /** @brief set of nodal functional */
    NodalFunctional nodal_functional;

    /** @brief element description */
    FEDescriptor fe_descriptor;

    /** @brief ID for reference transformation */
    ReferenceTransformation_type reference_transformation_type;

    /** @brief this constructor facilitates the implementation a bit */
    FiniteElement (FE_type id,
          std::tuple<BaseFunction_type, NodalFunctional_type, FEDescriptor_type,
                     ReferenceTransformation_type> ids);

  public:

    explicit FiniteElement ( FE_type id);

    /** @brief return type of base functions */
    BaseFunction_type GetBaseFunct_ID() const
    { return base_function.GetID(); }

    /** @brief return BaseFunctions object */
    const BaseFunctions *GetBaseFunct() const
    { return &base_function; }

    /** @brief return NodalFunctional2D_ID */
    NodalFunctional_type GetNodalFunctional_ID() const
    { return nodal_functional.GetID(); }

    /** @brief return NodalFunctional3D */
    const NodalFunctional *GetNodalFunctional() const
    { return &nodal_functional; }

    /** @brief return RefTransID */
    ReferenceTransformation_type GetRefTransID() const
    { return reference_transformation_type; }

    /** @brief return FEDesc3D_ID */
    FEDescriptor_type GetFEDesc_ID() const
    { return fe_descriptor.GetID(); }

    /** @brief return FEDesc3D */
    const FEDescriptor *GetFEDesc() const
    { return &fe_descriptor; }

    /** @brief return number of degrees of freedom */
    int GetN_DOF() const
    { return base_function.GetDimension(); }

    /** @brief check N[i](b[j]) = delta[ij] */
    void CheckNFandBF()  const;

    /** @brief return type of this finite element */
    FE_type GetID() const
    { return type; }
    
    /** @brief the space dimension this finite element is defined in*/
    int get_space_dim() const;
};

#endif
