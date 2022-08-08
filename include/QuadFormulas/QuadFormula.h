#ifndef __QUAD_FORMULA__
#define __QUAD_FORMULA__

#include "Enumerations_quadrature_formula.h"
#include "Point.h"
#include <vector>

/// @brief maximum polynomial degree which can be integrated exactly.
constexpr int MAXDEGREE = 24;
/// @brief maximum number of quadrature points for all 1D quadrature formulas
constexpr int MaxN_QuadPoints_1D = 16;
/// @brief maximum number of quadrature points for all 2D quadrature formulas
constexpr int MaxN_QuadPoints_2D = 81;
/// @brief maximum number of quadrature points for all 3D quadrature formulas
constexpr int MaxN_QuadPoints_3D = 729;

/** base class for quadrature formulas */
class TQuadFormula
{
  protected:
    
    /** @brief type of this quadrature formula */
    QuadratureFormula_type type;
    
    /** @brief vector of pairs of quadrature points and its weights*/
    std::vector<std::pair<double, parmoon::Point>> quadraturePairs;
    
    /// @name redundant data, for alternative access.
    /// @todo can we get rid of this redundancy?
    //@{
    /** weights for the formula */
    std::vector<double> Weights = {};
    /** first coordinate for the quadrature formula */
    std::vector<double> Xi;
    /** second coordinate for the quadrature formula */
    std::vector<double> Eta;
    /** third coordinate for the quadrature formula */
    std::vector<double> Zeta;
    //@}

  public:
    /// @brief create a quadrature formula of a given type on the reference cell
    explicit TQuadFormula(QuadratureFormula_type id);
    
    /** @brief return number of quadrature points */
    int GetN_QuadPoints() const
    { return quadraturePairs.size(); }
    
    /** @brief return the i-th weight of the formula */
    double get_weight(unsigned int i) const
    { return quadraturePairs[i].first; }

    /** @brief return the i-th evalutation point of the formula */
    parmoon::Point get_point(unsigned int i) const
    { return quadraturePairs[i].second; }
    
    /** @brief return accuracy of this formula */
    int GetAccuracy() const;
    
    /// @brief space dimension for which this quadrature formula is made
    int get_dimension() const;
    
    /** @brief bool is true if the quadrature formula is o a line, a triangle or a
     * tetrahedron and false otherwise*/
    bool is_on_simplices() const;
    
    /** @brief this prints information (type, accuracy, quadrature points,
     * weights) of the quadrature formula*/
    void info() const;
    
    /** @brief return the type of this quadrature formula */
    QuadratureFormula_type get_type() const
    { return type; }
    
    /** @brief update one of the quadrature points and its weight.
     * 
     * This is used whenever you want to use this object on a mesh cell, i.e.,
     * not on the reference cell. You have to update all points and weights for
     * this to make sense.
     */
    void update_pair(unsigned int i, std::pair<double, parmoon::Point> pair)
    { 
      quadraturePairs[i] = pair;
      Weights[i] = pair.first;
      Xi[i] = pair.second.x;
      if(!Eta.empty())
        Eta[i] = pair.second.y;
      if(!Zeta.empty())
        Zeta[i] = pair.second.z;
    }
    
    ///@brief return individual arrays to the coordinate of the points
    //@{
    const double* get_xi() const
    { return Xi.data(); }
    const double* get_eta() const
    { return Eta.data(); }
    const double* get_zeta() const
    { return Zeta.data(); }
    //@}
    
    /** return all data of the quadrature formula */
    void GetFormulaData(int &n_points, const double* &weights,
                        const double* &xi) const;
    void GetFormulaData(int &n_points, const double* &weights, 
                        const double* &xi, const double* &eta) const;
    void GetFormulaData(int &n_points, const double* &weights, 
                        const double* &xi, const double* &eta,
                        const double* &zeta) const;
};

#endif
