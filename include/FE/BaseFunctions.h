#ifndef __BASEFUNCT2D__
#define __BASEFUNCT2D__

#include "Constants.h"
#include "Enumerations_fe.h"
#include "Point.h"
#include <ostream>

class TBaseCell;
class TCollection;
class TQuadFormula;

/** @brief Transforms 1D quadrature points to 2D on reference cell.
 * @details
 * Transform a given point of a 1D quadrature rule to the two dimensional
 * quadrature point on a given joint, i.e. edge.
 * @param[in] RefElement reference element, see #BFRefElements
 * @param[in] joint_nr joint number of edge onto which the quadrature point is
 * transformed to
 * @param[in] zeta 1D quadrature point
 * \return a pair of 2D coordinates xi, eta on reference cell
 */
parmoon::Point transform(BFRefElements RefElement, int joint_nr,
                                    double zeta);

/** @brief Transforms 2D quadrature points to 3D on reference cell.
 * @details
 * Transform a given point of a 2D quadrature rule to the three dimensional
 * quadrature point on a given joint, i.e. face.
 * @param[in] RefElement reference element, see #BFRefElements
 * @param[in] joint_nr joint number of face onto which the quadrature point is
 * transformed to
 * @param[in] t x coordinate of 2D quadrature point
 * @param[in] s y coordinate of 2D quadrature point
 * \return a tuple of 3D coordinates xi, eta, iota on reference cell
 */
parmoon::Point transform(BFRefElements RefElement, int joint_nr,
                         double t, double s);

constexpr int MaxN_BaseFunctions1D = 4;
constexpr int MaxN_BaseFunctions2D = 100;
constexpr int MaxN_BaseFunctions3D = 125;

/** @brief set of all basis functions on the reference element for a finite 
    element in two dimensions */
class BaseFunctions
{
  protected:
    /** @brief number of base functions = dimension of local space */
    int Dimension;

    /** @brief Id for this set of base functions */
    BaseFunction_type BaseFunct;

    /** @brief array for all functions and derivatives */
    DoubleFunct1D *Functions1D[N_MultiIndices1D];

    /** @brief array for all functions and derivatives */
    DoubleFunct2D *Functions[N_MultiIndices2D];

    /** @brief array for all functions and derivatives */
    DoubleFunct3D *Functions3D[N_MultiIndices3D];

    /** @brief reference element used for this set of base functions */
    BFRefElements RefElement;

    /** @brief polynomial degree */
    int PolynomialDegree;

    /** @brief accuracy */
    int Accuracy;

    /** @brief number of basis functions per joint where the sign
        has to be changed if needed */
    int N_BF2Change;

    /** @brief indices of basis functions with changeable sign,
        sorted by joints */
    int *** BF2Change;
 
    /** @brief Dimension of the vector basis function */
    int BaseVectDim;
    
  public:

    explicit BaseFunctions ( BaseFunction_type id);
    
    /** @brief return the dimension of local space */
    int GetDimension() const
    { return Dimension; }

    /** @brief return type */
    BaseFunction_type GetID() const
    { return BaseFunct; }

    /** @brief return the values for derivative MultiIndex at xi */
    void GetDerivatives(MultiIndex1D MultiIndex, double xi, double *values)
      const
    { Functions1D[static_cast<int>(MultiIndex)](xi, values); };

    /** @brief return the values for derivative MultiIndex at all
        quadrature points */
    void GetDerivatives(MultiIndex1D MultiIndex, const TQuadFormula *formula,
                        double **values) const;


    /** @brief return the values for derivative MultiIndex at (xi,eta) */
    void GetDerivatives(MultiIndex2D MultiIndex, double xi,
                        double eta, double *values) const
    { Functions[static_cast<int>(MultiIndex)](xi, eta, values); };

    /** @brief return the values for derivative MultiIndex at all
        quadrature points */
    void GetDerivatives(MultiIndex2D MultiIndex, const TQuadFormula *formula,
                        double **values) const;
    
    /** @brief return values on given joint for all quadrature points */
    void GetDerivatives(MultiIndex2D MultiIndex, const TQuadFormula *formula,
                        int joint, double **Values) const;

    /** @brief return derivatives on joint i */
    void GetDerivatives(MultiIndex2D MultiIndex, int N_Points,
                        const double *zeta, int i, double **Values) const;


    /** @brief return the values for derivative MultiIndex at (xi,eta) */
    void GetDerivatives(MultiIndex3D MultiIndex, double xi,
                        double eta, double zeta, double *values) const
    { Functions3D[static_cast<int>(MultiIndex)](xi, eta, zeta, values); };

    /** @brief return the values for derivative MultiIndex at all
        quadrature points */
    void GetDerivatives(MultiIndex3D MultiIndex, 
                        const TQuadFormula *formula, double **values) const;

    /** @brief return the values for derivative MultiIndex at all
        quadrature points on 'joint' */
    void GetDerivatives(MultiIndex3D MultiIndex, const TQuadFormula *formula,
                        int joint, double **values) const;

    /** @brief return reference element */
    BFRefElements GetRefElement() const
    { return RefElement; };

    /** @brief return polynomial degree */
    int GetPolynomialDegree() const
      { return PolynomialDegree; };

    /** @brief return accuracy */
    int GetAccuracy() const
      { return Accuracy; };

    /** @brief change basis functions on cell if needed */
    void ChangeBF(const TCollection *Coll, const TBaseCell *Cell,
                  double *Values) 
      const;

    /** @brief change basis functions on cell in all points if needed */
    void ChangeBF(const TCollection *Coll, const TBaseCell *Cell, int N_Points,
                  double **Values) const;
    
    /** @brief return the dimension of the vector basis function */
    int GetBaseVectDim() const
    { return BaseVectDim; }
};

#endif
