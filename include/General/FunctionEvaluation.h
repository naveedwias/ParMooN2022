#ifndef __FUNCTIONEVALUATION__
#define __FUNCTIONEVALUATION__

#include <vector>

/** ************************************************************************ 
*
* @class      FunctionEvaluation
* @brief      represent an evaluation of a set of (basis) functions and their 
*             derivatives at one (quadrature) point
*
* Accessing individual functions is realized via get and set methods. Accessing
* multiple function values for a certain derivative are realized via begin/end
* methods which return iterators to the correct location in the vector stored
* by this class.
*
* @author     Ulrich Wilbrandt & Alfonso Caiazzo & Swetlana Schyschlowa
* @date       22.09.2014
*
****************************************************************************/

#include <ostream>
/** @brief all supported multi indices for derivatives in 1D, 2D, or 3D 
 * @todo what about a MultiIndex for the time derivative, e.g. Dt?
 */
enum class MultiIndex
{
  D0, D1, D2,
  D00, D10, D01, D20, D11, D02,
  D000, D100, D010, D001, D200, D110, D101, D020, D011, D002
};
/// @brief print the MultiIndex as they are written in the definition above
std::ostream& operator<<(std::ostream& out, const MultiIndex value);
/// @brief return the dimension to which a MultiIndex belongs to
unsigned int space_dimension_of_MultiIndex(MultiIndex m);

class FunctionEvaluation
{
 private:
  
  /** @brief store function values and derivatives for a set of functions
   * 
   * The size of this vector is a multiple of either 3 (1D), 6 (2D), or 10 (3D) 
   * which corresponds to the number of multi-indices up to order two available 
   * in the respective space dimension. The order of the multi-indices is 
   * described in this->get(unsigned int, unsigned int, unsigned int). For each 
   * multi-index i all evaluations (values or derivatives respectively) are 
   * stored. Its size is the number of functions and equal for all valid
   * multi-indices i. So the overall size of this vector is the product of 
   * 3 (1D), 6 (2D), or 10 (3D) and the number of functions.
   * 
   * If this has been constructed using vector valued basis functions, the size
   * of the vector this->values is d (space dimension) times larger than 
   * described above. So for each multi-index there are d components. 
   */
  std::vector<double> values;
  
  /** @brief number of evaluations per multi index and component */
  unsigned int n_funct;
  
  /** @brief number of components
   * 
   * Usually this is one, only for vector valued basis functions this is 2 or 3
   * depending on the space dimension.
   */
  unsigned int n_components;
 public:
  /// @brief constructor with a given number of functions for a certain space 
  ///        dimension
  FunctionEvaluation(unsigned int space_dimension, unsigned int n_functions,
                     unsigned int n_comp = 1);
  
  /// @brief copy constructor
  FunctionEvaluation(const FunctionEvaluation&) = default;
  
  /// @brief move constructor
  FunctionEvaluation(FunctionEvaluation&&) noexcept = default;
  
  // destructor
  ~FunctionEvaluation() noexcept = default;
  
  /// @brief copy-assignment
  FunctionEvaluation& operator=(const FunctionEvaluation&) = default;
  
  /// @brief move-assignment
  FunctionEvaluation& operator=(FunctionEvaluation&&) = default;
  
  
  /** @brief get value or derivative of the i-th function 
   * 
   * The component is needed for vector valued basis functions, in that case 
   * you can use component=1 or 2 (in 3D) as well.
   * 
   * This is the same as get(unsigned int, unsigned int, unsigned int) but the 
   * MultiIndex is mapped (hard coded) to the correct index.
   */
  double get(unsigned int i, MultiIndex, unsigned int component = 0) const;
  
  /// @brief a convenience method for less typing
  double operator()(unsigned int i, MultiIndex mi, unsigned int component = 0)
  const
  { return this->get(i, mi, component); }
  
  /** @brief get value or derivative of the i-th function 
   * 
   * The unsigned int m describes the MultiIndex. In 1D m < 3 and 
   *  - m==0 means MultiIndex::D0
   *  - m==1 means MultiIndex::D1
   *  - m==2 means MultiIndex::D2
   * In 2D m < 6 and
   *  - m==0 means MultiIndex::D00
   *  - m==1 means MultiIndex::D10
   *  - m==2 means MultiIndex::D01
   *  - m==3 means MultiIndex::D20
   *  - m==4 means MultiIndex::D11
   *  - m==5 means MultiIndex::D02
   * In 3D m < 10 and
   *  - m==0 means MultiIndex::D000
   *  - m==1 means MultiIndex::D100
   *  - m==2 means MultiIndex::D010
   *  - m==3 means MultiIndex::D001
   *  - m==4 means MultiIndex::D200
   *  - m==5 means MultiIndex::D110
   *  - m==6 means MultiIndex::D101
   *  - m==7 means MultiIndex::D020
   *  - m==8 means MultiIndex::D011
   *  - m==9 means MultiIndex::D002
   * 
   * The component is needed for vector valued basis functions, in that case 
   * you can use component=1 or 2 (in 3D) as well.
   */
  double get(unsigned int i, unsigned int m, unsigned int component = 0) const;
  /** @brief exactly like get(), but without any checks, so this should be 
   *         faster */
  double get_no_checks(unsigned int i, unsigned int m,
                       unsigned int component = 0) const;
  
  /** @brief return the beginning iterator to access the derivatives of all 
   * functions
   */
  std::vector<double>::iterator begin(MultiIndex);
  /** @brief return the beginning const_iterator to access the derivatives of
   * all functions
   */
  std::vector<double>::const_iterator begin(MultiIndex) const;
  /** @brief return the ending iterator to access the derivatives of all
   * functions
   */
  std::vector<double>::iterator end(MultiIndex);
  /** @brief return the ending const_iterator to access the derivatives of all
   * functions
   */
  std::vector<double>::const_iterator end(MultiIndex) const;
  
  /** @brief set the multi-index of the i-th function to 'a'
   * 
   * See also this->get.
   */
  void set(double a, unsigned int i, MultiIndex, unsigned int component = 0);
  
  /** @brief set the m-th multi-index of the i-th function to a
   * 
   * See also this->get.
   */
  void set(double a, unsigned int i, unsigned int m,
           unsigned int component = 0);
  
  /** @brief scale the i-th basis function by the factor a
   * 
   * This multiplies the corresponding basis function and all its derivatives
   */
  void scale(unsigned int i, double a);
  
  /** @brief add another FunctionEvaluation scaled by a to this
   * 
   *  For each derivative, do 'this = this + a*other'.
   */
  void add(const FunctionEvaluation& other, double a = 1.0);
  
  /** @brief merge all functions to just one function
   * 
   * For each MultiIndex all functions stored in this object are added into one 
   * function.
   */
  void merge();
  
  /** @brief the number of multi indices for which function evaluations are 
   *         stored in this object. 
   */
  unsigned int n_multi_indices() const
  { return this->values.size()/(this->n_funct * this->n_components); }
  
  /** @brief return the space dimension
   *
   * This uses the number of multi-indices.
   */
  unsigned int space_dimension() const;
  
  /** @brief return the number of functions
   * 
   * If the basis functions are vector valued then this function returns the 
   * product of the number of functions and the space dimension.
   */
  unsigned int n_functions() const
  { return n_funct * n_components; }
  
  /// @brief set all entries to zero (this is similar to creating a new object)
  void reset();
};



#endif //__FUNCTIONEVALUATION__
