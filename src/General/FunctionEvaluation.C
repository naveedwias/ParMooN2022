#include <FunctionEvaluation.h>
#include <MooNMD_Io.h>

std::ostream& operator<<(std::ostream& out, const MultiIndex value)
{
  const char* s = 0;
#define PROCESS_VAL(p) case(MultiIndex::p): s = #p; break;
  switch(value)
  {
    PROCESS_VAL(D0); PROCESS_VAL(D1); PROCESS_VAL(D2);
    PROCESS_VAL(D00); PROCESS_VAL(D10); PROCESS_VAL(D01);
    PROCESS_VAL(D20); PROCESS_VAL(D11); PROCESS_VAL(D02);
    PROCESS_VAL(D000); PROCESS_VAL(D100); PROCESS_VAL(D010); PROCESS_VAL(D001);
    PROCESS_VAL(D200); PROCESS_VAL(D110); PROCESS_VAL(D101); PROCESS_VAL(D020);
    PROCESS_VAL(D011); PROCESS_VAL(D002); 
    default: s = "unknown MultiIndex type"; break;
  }
#undef PROCESS_VAL
  return out << s;
}

/* ************************************************************************* */
unsigned int space_dimension_of_MultiIndex(MultiIndex m)
{
  switch(m)
  {
    case MultiIndex::D0:
    case MultiIndex::D1:
    case MultiIndex::D2:
      return 1;
    case MultiIndex::D10:
    case MultiIndex::D00:
    case MultiIndex::D01:
    case MultiIndex::D20:
    case MultiIndex::D11:
    case MultiIndex::D02:
      return 2;
    case MultiIndex::D000:
    case MultiIndex::D100:
    case MultiIndex::D010:
    case MultiIndex::D001:
    case MultiIndex::D200:
    case MultiIndex::D110:
    case MultiIndex::D101:
    case MultiIndex::D020:
    case MultiIndex::D011:
    case MultiIndex::D002:
      return 3;
    default:
      ErrThrow("unknown MultiIndex ", m);
  }
}

/* ************************************************************************* */
FunctionEvaluation::FunctionEvaluation(unsigned int space_dimension, 
                                       unsigned int n_functions, 
                                       unsigned int n_comp)
 : n_funct(n_functions), n_components(n_comp)
{
  if(this->n_funct == 0)
    ErrThrow("unable to create a FunctionEvaluation with zero functions");
  if(this->n_components == 0 || this->n_components > space_dimension)
    ErrThrow("unable to create a FunctionEvaluation with ", this->n_components,
             " components");
  // the length of each sub-vector
  unsigned int n = this->n_funct * this->n_components;
  switch(space_dimension)
  {
    case 1:
      this->values.resize(3 * n, 0.); // D0, D1, D2
      break;
    case 2:
      this->values.resize(6 * n, 0.); // D00, D10, D01, D20, D11, D02
      break;
    case 3:
      // D000, D100, D010, D001, D200, D110, D101, D020, D011, D002
      this->values.resize(10* n, 0.);
      break;
    default:
      ErrThrow("The space dimension ", space_dimension, " is unsupported");
  }
}

/* ************************************************************************* */
unsigned int MultiIndex_to_vector_index(MultiIndex m)
{
  switch(m)
  {
    case MultiIndex::D0:
    case MultiIndex::D00:
    case MultiIndex::D000:
      return 0;
    case MultiIndex::D1:
    case MultiIndex::D10:
    case MultiIndex::D100:
      return 1;
    case MultiIndex::D2:
    case MultiIndex::D01:
    case MultiIndex::D010:
      return 2;
    case MultiIndex::D20:
    case MultiIndex::D001:
      return 3;
    case MultiIndex::D11:
    case MultiIndex::D200:
      return 4;
    case MultiIndex::D02:
    case MultiIndex::D110:
      return 5;
    case MultiIndex::D101:
      return 6;
    case MultiIndex::D020:
      return 7;
    case MultiIndex::D011:
      return 8;
    case MultiIndex::D002:
      return 9;
    default:
      ErrThrow("unknown MultiIndex ", m);
  }
}


/* ************************************************************************* */
unsigned int index_in_vector(unsigned int i, MultiIndex m,
                             unsigned int component, unsigned int n_functions,
                             unsigned int n_components)
{
  if(i >= n_functions)
    ErrThrow("unable to locate entry of function ", i, ". There are only ",
             n_functions, " functions");
  if(component >= n_components)
    ErrThrow("unable to locate entry of component ", component, ". The number "
             "of components is ", n_components);
  unsigned int n = n_functions * n_components;
  return MultiIndex_to_vector_index(m)*n + n_functions * component + i;
}

/* ************************************************************************* */
double FunctionEvaluation::get(unsigned int i, MultiIndex m, 
                               unsigned int component) const
{
  if(this->space_dimension() != space_dimension_of_MultiIndex(m))
    ErrThrow("The MultiIndex ", m, " is not suitable in ",
             this->space_dimension(), "D");
  auto index = index_in_vector(i, m, component, this->n_funct,
                               this->n_components);
  return this->values[index];
}

/* ************************************************************************* */
double FunctionEvaluation::get(unsigned int i, unsigned int m, 
                               unsigned int component) const
{
  unsigned int n = this->n_funct * this->n_components;
  if(m >= this->values.size()/n) // either 3 (1D), 6 (2D), or 10 (3D)
  {
    ErrThrow("unknown index for the MultiIndex ", m);
  }
  if(i >= this->n_funct)
  {
    ErrThrow("index out of bounds in FunctionEvaluation::get ", i);
  }
  if(component >= this->n_components)
  {
    ErrThrow("component out of bounds in FunctionEvaluation::get ", i);
  }
  return this->values[m*n + i + this->n_funct * component];
}
/* ************************************************************************* */
double FunctionEvaluation::get_no_checks(unsigned int i, unsigned int m,
                                         unsigned int component) const
{
  return this->values[(m*this->n_components + component)*this->n_funct + i];
}

/* ************************************************************************* */
std::vector<double>::iterator FunctionEvaluation::begin(MultiIndex m)
{
  if(this->space_dimension() != space_dimension_of_MultiIndex(m))
    ErrThrow("The MultiIndex ", m, " is not suitable in ",
             this->space_dimension(), "D");
  unsigned int index = index_in_vector(0, m, 0, this->n_funct,
                                       this->n_components);
  return this->values.begin() + index;
}

/* ************************************************************************* */
std::vector<double>::const_iterator FunctionEvaluation::begin(MultiIndex m) 
  const
{
  if(this->space_dimension() != space_dimension_of_MultiIndex(m))
    ErrThrow("The MultiIndex ", m, " is not suitable in ",
             this->space_dimension(), "D");
  unsigned int index = index_in_vector(0, m, 0, this->n_funct,
                                       this->n_components);
  return this->values.cbegin() + index;
}

/* ************************************************************************* */
std::vector<double>::iterator FunctionEvaluation::end(MultiIndex m)
{
  if(this->space_dimension() != space_dimension_of_MultiIndex(m))
    ErrThrow("The MultiIndex ", m, " is not suitable in ",
             this->space_dimension(), "D");
  unsigned int n = this->n_funct * this->n_components;
  return this->values.begin() + (MultiIndex_to_vector_index(m)+1)*n;
}

/* ************************************************************************* */
std::vector<double>::const_iterator FunctionEvaluation::end(MultiIndex m) 
  const
{
  if(this->space_dimension() != space_dimension_of_MultiIndex(m))
    ErrThrow("The MultiIndex ", m, " is not suitable in ",
             this->space_dimension(), "D");
  unsigned int n = this->n_funct * this->n_components;
  return this->values.begin() + (MultiIndex_to_vector_index(m)+1)*n;
}

/* ************************************************************************* */
void FunctionEvaluation::set(double a, unsigned int i, MultiIndex m,
                             unsigned int component)
{
  if(this->space_dimension() != space_dimension_of_MultiIndex(m))
    ErrThrow("The MultiIndex ", m, " is not suitable in ",
             this->space_dimension(), "D");
  unsigned int index = index_in_vector(i, m, component, this->n_funct,
                                       this->n_components);
  this->values[index] = a;
}


/* ************************************************************************* */
void FunctionEvaluation::set(double a, unsigned int i, unsigned int m,
                             unsigned int component)
{
  unsigned int n = this->n_funct * this->n_components;
  if(m >= this->values.size() / n) // either 3 (1D), 6 (2D), or 10 (3D)
  {
    ErrThrow("unknown index for the MultiIndex ", m);
  }
  if(i >= this->n_funct)
  {
    ErrThrow("index out of bounds in FunctionEvaluation::get ", i);
  }
  if(component >= this->n_components)
  {
    ErrThrow("component out of bounds in FunctionEvaluation::get ", i);
  }
  this->values[m*n + i + this->n_funct * component] = a;
}

/* ************************************************************************* */
void FunctionEvaluation::scale(unsigned int i, double a)
{
  if(i >= this->n_funct)
  {
    ErrThrow("there are only ", n_funct, " function evaluations, not ", i);
  }
  // loop over all possible multi indices
  for(unsigned int m = 0, n_m = this->n_multi_indices(); m < n_m; ++m)
  {
    // loop over all components
    for(unsigned int c = 0; c < this->n_components; ++c)
    {
      this->values[(m * this->n_components  + c) * this->n_funct + i] *= a;
    }
  }
}

/* ************************************************************************* */
void FunctionEvaluation::add(const FunctionEvaluation& other, double a)
{
  if(this->n_functions() != other.n_functions())
  {
    ErrThrow("cannot add two FunctionEvaluation objects of different size");
  }
  if(this->values.size() != other.values.size())
  {
    ErrThrow("incompatible FunctionEvaluations to add");
  }
  for(unsigned int i = 0, n_i = this->values.size(); i < n_i; ++i)
  {
    this->values[i]  += a * other.values[i];
  }
}

/* ************************************************************************* */
void FunctionEvaluation::merge()
{
  // loop over all possible multi indices
  unsigned int n_m = this->n_multi_indices();
  for(unsigned int m = 0; m < n_m; ++m)
  {
    // loop over all components
    for(unsigned int c = 0; c < this->n_components; ++c)
    {
      // starting index of values in 'other'
      unsigned int index = (m * this->n_components + c) * this->n_funct;
      // set value from first function, add the rest in a loop
      this->values[m * this->n_components + c] = this->values[index];
      // loop over all functions
      for(unsigned int f = 1; f < this->n_funct; ++f)
      {
        this->values[m * this->n_components + c] += this->values[index + f];
      }
    }
  }
  this->n_funct = 1;
  this->values.resize(n_m * this->n_components); // reduce to only one function
}

/* ************************************************************************* */
unsigned int FunctionEvaluation::space_dimension() const
{
  switch(this->n_multi_indices())
  {
    case 3: return 1;
    case 6: return 2;
    case 10: return 3;
  }
  ErrThrow("cannot determine space dimension in FunctionEvaluation object");
}

/* ************************************************************************* */
void FunctionEvaluation::reset()
{
  std::fill(this->values.begin(), this->values.end(), 0.0);
}
