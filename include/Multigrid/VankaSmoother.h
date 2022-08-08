/**
 * @file VankaSmoother.h
 *
 * @date 2016/05/16
 * @author Clemens Bartsch
 */

#ifndef INCLUDE_MULTIGRID_VANKASMOOTHERNEW_H_
#define INCLUDE_MULTIGRID_VANKASMOOTHERNEW_H_

#include <DofPatch.h>
#include <Smoother.h>

#include <memory>

enum class VankaType {NODAL, CELL, PATCH, CELL_JACOBI};

//forward declaration
class DenseMatrix;
class TMatrix;
class TFESpace;
class TFESpace3D;

/**
 * A class for Vanka smoother as used in the multgrid method. The
 * implementation and nomenclature are adapted to (Navier--)Stokes saddle
 * point problems - but extending it should be relatively simple.
 *
 * The Smoothers are now also available for Convection_difusion equations.
 * The only thing a matrix must fulfil in order to use a Vanka smoother is,
 * that all its spaces but the last one must be identical. This should describe
 * a saddle point problem, but might easily degenerate to a matrix where all
 * spaces are the same (e.g. a 1x1 Convection-Difusion-Reaction matrix).
 *
 * There is a variety of different Vanka smoothers available by now. The basic
 * three are nodal, cell and patch vanka. Note that cell vanka works only
 * properly for discontinuous pressure approximations in NSE. Nodal Vanka is
 * useless for CDR problems, because there the local systems degenerate to 1x1 size.
 * In SEQ, these three can be used in a "_store" version, too, where the local
 * matrices are stored. You must be aware, that this tends to eating a LOT of
 * memory, therefore it should be used for small and medium sized problems only.
 * But there it is often the fastest option.
 * Of the cell vanka a blockwise Jacobi version is available (cell_vanka_jacobi),
 * which performs the update of the global solution in a Jacobi- instead of a
 * Gauss-Seidel manner. Its main advantage ist, that it is "embarassingly parallel",
 * and therefore its MPI counterpart (which is available!) gives the exact same
 * results. In MPI, we also have the standard (Gauss Seidel) cell vanka available.
 * It did not work bad for test examples, although usually its convergence is a
 * bit worse than when using the SEQ version - in one MPI Gauss Seidel sweep
 * information transport is restricted to the process parts, while in SEQ information
 * gets transported across the entire domain.
 *
 * Note that the main issue which held us from parallelizing the nodal Vanka is
 * the issue of 'toxic systems'. In MPI, local systems stemming from a pressure
 * dof at the interface will necessarily contain halo dof. In the process local
 * matrix, halo rows are incorrect, since up to now FEMatrices in ParMooN are not
 * stored fully consistent, but only in what we call "rowwise level 1 consistency".
 * I.e. all halo rows are incorrect, which will lead to erroneous local systems
 * at the interface in the Vanka smoothers. Until this matrix storage issue is
 * solved, parallel nodal vanka is out of reach.
 *
 * Note that all Vankas can be used as preconditioners (see class Preconditioner_vanka),
 * although not all of them are already enabled there. This is but only a
 * technical issue, which you can fix in seconds. Don't expect them to work
 * well as preconditioners, but they were really useful for debugging the
 * smoothers, and they might be useful when you have to solve a big saddle point
 * problem in 3D and no grid hierarchy available for multigrid.
 *
 */
class VankaSmoother : public Smoother
{
  public:
    /** Default constructor.
     * @param type The Vanka type to use. Currently nodal, cell and patch
     * Vanke are available.
     * @param damp_factor The damping factor to use when updating the solution
     * of the global system by the solution of the local system.
     * @param store Whether to store the local systems or not. This can considerably
     * speed up the computation (especially if the grid level is traversed very
     * often) but is very memory intensive. You should use it for test cases only.
     */
    VankaSmoother(VankaType type, double damp_factor,  bool store = false);

    /// Perform one step of Vanka smoothing (solve all local systems).
    void smooth(const BlockVector& rhs, BlockVector& solution ) override;

    /// Update the local matrices. Must be called whenever the global matrix has changed.
    void update(const BlockFEMatrix&) override;

    /* ************* *
     * Special member functions. Declared but not defined, since it
     * is not yet clear whether to shallow or deep copy here.
     * ************* */
    //! Default constructor.
    VankaSmoother() = delete;

    //! Copy constructor.
    VankaSmoother( const VankaSmoother& );

    //! Move constructor.
    VankaSmoother( VankaSmoother&& );

    //! Copy assignment operator.
    VankaSmoother& operator=( const VankaSmoother& );

    //! Move assignment operator.
    VankaSmoother& operator=( VankaSmoother&& );

    ~VankaSmoother();

  private:
    /// The type of Vanka smoother (nodal, cell, patch)
    VankaType type_;

    /// The dimension of the velocity equation (usually 2 or 3).
    size_t dimension_;

    /// Damping factor to be used int the local-to-global updates.
    double damp_factor_;

    /// The corresponding global matrix.
    std::shared_ptr<TMatrix> matrix_global_;

#ifdef _MPI
    /// Vector of communicators belonging to the same matrix as matrix_global_.
    /// Necessary for MPI updates after a smoothing step.
    std::vector<const TParFECommunicator3D*> comms_;
#endif

    /// The collection of pressure dofs for the local systems.
    std::vector<DofPatch> press_dofs_local_;

    /// The collection of velocity dofs for the local systems.
    std::vector<DofPatch> velo_dofs_local_;

    //Store the matrices which were assembled in order to reuse their LU factorization.
    std::vector<DenseMatrix*> local_systems_;
    bool store_systems_;

// TODO Change these to weak pointers as soon as FESpaces are stored in a
// shared_ptr fashion.
#ifdef __2D__
    const TFESpace2D* pressure_space_;

    const TFESpace2D* velocity_space_;
#elif __3D__
    const TFESpace3D* pressure_space_;

    const TFESpace3D* velocity_space_;
#endif


    void set_up_pressure_patches(const TFESpace& pressureSpace);

    void set_up_velocity_patches(const TMatrix& pressureVelocityMatrix,
                                 const TFESpace& velocitySpace);


};


#endif /* INCLUDE_MULTIGRID_VANKASMOOTHERNEW_H_ */
