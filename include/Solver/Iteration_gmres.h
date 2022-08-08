#ifndef __ITERATION_GMRES__
#define __ITERATION_GMRES__

#include <IterativeMethod.h>
#include <vector>

// not nice: pollute global namespace
enum class gmres_type {left, right, flexible};

template <class LinearOperator, class Vector>
class Iteration_gmres : public IterativeMethod<LinearOperator, Vector>
{
  public:
    /** constructor */
    Iteration_gmres(std::shared_ptr<Preconditioner<Vector>> p,
                    gmres_type t = gmres_type::flexible);
    
    /** destructor */
    virtual ~Iteration_gmres() = default;
    
    /** iterate routine */
    std::pair<unsigned int, double> iterate(const LinearOperator & A, 
                                            const Vector & rhs,
                                            Vector & solution) override final;
    
  protected:
    
    gmres_type type;
    
    // temporaries, stored only to avoid reallocation:
    // these have size restart+1
    std::vector<double> s;
    std::vector<double> cs;
    std::vector<double> sn;
    std::vector<Vector> v;
    // only for flexible gmres
    std::vector<Vector> z;
    
    std::pair<unsigned int, double> left_gmres(const LinearOperator & A, 
                                               const Vector & rhs,
                                               Vector & solution);
    std::pair<unsigned int, double> right_gmres(const LinearOperator & A, 
                                                const Vector & rhs,
                                                Vector & solution);
    std::pair<unsigned int, double> flexible_gmres(const LinearOperator & A, 
                                                   const Vector & rhs,
                                                   Vector & solution);
};


#endif // __ITERATION_GMRES__
