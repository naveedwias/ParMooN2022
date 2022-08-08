#ifndef __SOLVER_H__
#define __SOLVER_H__

#include <ParameterDatabase.h>
#include <BlockFEMatrix.h>
#include <BlockVector.h>
#include <memory>

// forward declarations
template <class LinearOperator, class Vector> class IterativeMethod;
class DirectSolver;
template <class Vector> class Preconditioner;
class Multigrid;
class PETScSolver;
template <typename T> class LoopInfo;

/** @brief Solve a linear system
 *
 * This class can solve linear systems of the \f$Ax=b\f$. How exactly it is
 * solved is controlled entirely through the ParameterDatabase object given to
 * the constructor. After construction there is no way of changing the behavior
 * of the class. Typically you only need to call `solve` with a matrix, right
 * hand side and a solution vector.
 *
 * Currently this class is instantiated for `<BlockMatrix, BlockVector>` and
 * `<BlockFEMatrix, BlockVector>`.
 *
 * @todo make members direct_solver, iterative_method and preconditioner a
 * unique_ptr.
 */
template <class LinearOperator = BlockFEMatrix, class Vector = BlockVector>
class Solver
{
  public:
    /// @brief create a solver object
    ///
    /// All solver related parameters are set to default values except those
    /// which also exist in the given ParameterDatabase.
    ///
    /// If you want to solve a system only once, you can use the method `solve`
    /// which has the linear operator as well as the right hand side and
    /// solution as arguments. In case of a direct solver you can use the
    /// already computed factorization a second time through calling `solve`
    /// with only the right hand side and the solution as arguments.
    explicit Solver(const ParameterDatabase& param_db);

    /// @brief update the solver and preconditioner objects
    ///
    /// This methods prepares the members. It either creates new ones or calls
    /// an appropriate update method on the objects. For example a direct_solver
    /// object will be created (deleting the old one, if existing) in case a
    /// direct solver is used at all. A Saddle_point_preconditioner on the other
    /// hand only needs to be created once, and then update is enough. This
    /// method does this.
    ///
    /// In case of a multigrid preconditioner, the multigrid object is updated
    /// within this method - all its levels and their smoothers get updated.
    /// This behaviour is based on the assumption, that at the time this method
    /// is called all matrices on all levels have been reassembled.
    void update_matrix(const LinearOperator& matrix);

    /// @brief solve after calling `Solver::update_matrix`
    ///
    /// You need to call update matrix, to use this method. How exactly this is
    //// solved is determined by the ParameterDatabase.
    ///
    /// For direct solvers, where the factorization is stored, you can solve
    /// many times for different right hand sides.
    void solve(const Vector& rhs, Vector& solution);


    void solve_augmented(const Vector& rhs, Vector& solution);

    /// @brief solve the sytem matrix*solution = rhs
    ///
    /// This is only a wrapper for
    ///     update_matrix(matrix);
    ///     solve(rhs, solution);
    void solve(const LinearOperator& matrix, const Vector& rhs,
               Vector& solution);

    /// @brief return a constant reference to the local ParameterDatabase
    ///
    /// Note that you can not change the behavior of this class after
    /// construction. This method only lets you inspect the solver parameters
    const ParameterDatabase& get_db() const;

    /// @brief Return true if multigrid is expected as preconditioner.
    bool is_using_multigrid() const;

    /// @brief return the multigrid object, throws if not using multigrid
    std::shared_ptr<const Multigrid> get_multigrid() const;
    /// @brief return the multigrid object, throws if not using multigrid
    std::shared_ptr<Multigrid> get_multigrid();

    /// @brief If this is an iterative solver, return a const reference to
    ///        its loop info object. Otherwise throw.
    const LoopInfo<double>& get_solver_loop_info() const;

    /// @brief return a default solver parameter database
    ///
    /// Using the Solver class requires these parameters.
    static ParameterDatabase default_solver_database();

    void set_silent(bool silent);

  protected:

    /// @brief the ParameterDatabase which controls the entire solving process
    ParameterDatabase db;

    /// @brief a pointer to the linear operator
    ///
    /// This is changed in the method update_matrix. Note that this is only a
    /// pointer, so there is no way to check if the object itself changed.
    const LinearOperator* linear_operator;

    /// @brief this object is only created if needed.
    std::shared_ptr<DirectSolver> direct_solver;
    /// @brief this object is only created if needed.
    std::shared_ptr<IterativeMethod<LinearOperator, Vector>> iterative_method;
    /// @brief this object is only created if needed.
    std::shared_ptr<Preconditioner<Vector>> preconditioner;
    /// @brief this object is only created if needed.
    /// @todo the multigrid object should be part of the preconditioner object
    std::shared_ptr<Multigrid> multigrid;
    /// @brief this object is only created if needed.
    std::shared_ptr<PETScSolver> petsc_solver;
};

#endif // __SOLVER_H__
