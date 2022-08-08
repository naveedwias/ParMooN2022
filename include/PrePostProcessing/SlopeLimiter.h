#ifndef __SLOPELIMITER_H__
#define __SLOPELIMITER_H__

#include "ParameterDatabase.h"
#ifdef __3D__
 #include "FEFunction3D.h"
#else
 #include "FEFunction2D.h"
#endif

/**
 * @class TSlopeLimiter
 * @author Derk Frerichs-Mihov
 * @brief A class for applying slope limiter as post-processing method to a FE
 * function.
 *
 * @details This class allows to post-process a FE function using so-called
 * slope limiters. Usually the FE function is the solution of your system which
 * shall be post-processed for instance due to spurious oscillations.
 * To be precise this class is invented to reduce unphysical
 * values in convection dominated convection-diffusion problems.
 * @note Note that this method changes the given FE function. Until now it can
 * be only used for DG-elements, i.e. completely discontinuous elements. As of
 * today, 14.02.2022, they do have the ANSATZ NUMBERS -11,-12, and so on.
 *
 * Up to now many limiters are implemented which are described in the papers
 * \"<a href="https://doi.org/10.1016/j.cam.2021.113487">On reducing spurious
 * oscillations in discontinuous Galerkin (DG) methods for steady-state
 * convection–diffusion equations</a>\" and \"<a
 * href="https://doi.org/10.1016/j.aml.2022.107969">On a technique for reducing
 * spurious oscillations in DG solutions of convection–diffusion equations</a>\"
 * by Derk Frerichs-Mihov and Volker John. It is highly recommended to have look
 * in these paper if you want to understand what this class does. All limiters
 * are explained in detail in these publications.
 *
 * The first limiter has the name "LinTriaReco" and is the implementation of the
 * book \"<a href="https://doi.org/10.1137/1.9780898717440">Discontinuous
 * Galerkin methods for solving elliptic and parabolic equations theory and
 * implementation</a>\" by Beatrice Riviere, chapter 4.3.2 "slope limiters". It
 * is the triangle case of the chapter mentioned above. However, it can be
 * applied to quadrilaterals as well. In the paper the quadrilateral case is
 * called LinQuadReco, so this is how the limiter can be accessed as well.
 *
 * The third one has the name "ConstTriaReco" and is basically the previous
 * mentioned limiter except for some changes: Instead of considering the
 * value of the FE function at the edge midpoint the integral mean along the
 * edge is examined. Second, instead of a linear reconstruction the FE
 * function on marked cells is always changed to equal its integral mean
 * over the cell. And third, for boundary cells it reflects the boundary cell
 * along the domain boundary and uses the continuation of the FE function of the
 * boundary cell to the reflected cell in the algorithm. Also this limiter can
 * be applied for quadrilateral cells which can be accessed by "ConstQuadReco".
 *
 * The next one has the name "LinQuadDeriv" and is the implementation of the
 * quadrilateral case of the same chapter of the above mentioned book. Please
 * see in the book what the constants M_lim and gamma are, that can be and have
 * to be specified by either calling the appropriate constructor or calling
 * set_m_lim() and set_gamma().
 *
 * The next one has the name "ConstQuadDeriv" and is basically the previous
 * mentioned limiter except that in the case that of a cell where the
 * function shall be limited, it is always changed to its integral mean over
 * the cell.
 *
 * Another limiter can be addressed using "ConstJump" and is based on
 *\"<a href="https://www.doi.org/10.21136/MB.2002.134171">On discontinuous
 * Galerkin methods for nonlinear convection-diffusion problems and compressible
 * flow</a>\" by Vit Doljesi, Miroslav Feistauer and Christoph Schwab.
 *
 * A modification of the previous limiter can be accessed by choosing
 * "ConstJumpMod". It is the version of the second paper mentioned in the introduction of this
 * class. The respective attributes of the class are called \link
 * characteristic_length
 * \endlink, \link characteristic_solution_scale \endlink and \link C0_CJM
 * \endlink.
 *
 * Last but not least, there are also three limiter, called ConstJumpL1Norm,
 * ConstJumpL2Norm and ConstJumpLinftyNorm. The constant \f$C_1\f$ that controls
 * how many standard deviations are added or subtracted is called \link C1_CJN
 * \endlink and can be controlled using \link set_C1_CJN \endlink
 * . To ensure that continuous solutions are not treated, no cells are
 * marked if this reference value is smaller or equal to \f$10^{-13}\f$.
 *
 * @note Note that on all cells the same limiter is used. Up to now it is
 * not possible to switch the limiter accordingly to the cells. Therefore,
 * the limiters "LinQuadDeriv" and "ConstQuadDeriv" work only on meshes that
 * consist of only quadrilaterals.
 *
 * @todo implement / check limiter in 3D
 */
class TSlopeLimiter
{
  public:

#ifdef __2D__
    using FEFunction = TFEFunction2D;
#else
    using FEFunction = TFEFunction3D;
#endif

    /**
     * @brief constructs a slope limiter using the information given in the
     * database
     * @details For this constructor to work the parameters "limiter_name",
     * "M_lim", "gamma_limiter", "alpha_ref" have to be specified in the
     * database, e.g. by calling default_slope_limiter_database() first.
     * The possible limiters are
     * "Galerkin" (no limiter at all, see also p. 12 in the paper mentioned in
     * the class description), "LinTriaReco" (see Section 3.1), "ConstTriaReco"
     * (see Section 3.1), "LinQuadReco" (see Section 3.1), "ConstQuadReco"
     * (see Section 3.1), "LinQuadDeriv" (see Section
     * 3.2), "ConstQuadDeriv" (see Section 3.2), "ConstJump" (see Section 3.3),
     * "ConstJumpMod" (see Section 3.3).
     */
    TSlopeLimiter(const ParameterDatabase& param_db);


    // Special member functions. Disable copy/move, set destructor to default.
    // Will be changed only when the underlying classes follow rule of 0/5.

    //! default copy constructor
    TSlopeLimiter(const TSlopeLimiter&) = default;

    //! Delete move constructor.
    TSlopeLimiter(TSlopeLimiter&&) = delete;

    //! Delete copy assignment operator.
    TSlopeLimiter& operator=(const TSlopeLimiter&) = delete;

    //! Delete move assignment operator.
    TSlopeLimiter& operator=(TSlopeLimiter&&) = delete;

    //! Destructor.
    ~TSlopeLimiter() = default;

    /**
     * @brief constructs all the necessary parameters needed for the class
     * TSlopeLimiter
     */
    static ParameterDatabase default_slope_limiter_database();

    /**
     * @brief gets \link limiter_name \endlink
     */
    std::string get_limiter_name() {return limiter_name;};

    /**
     * @brief returns \link m_lim \endlink
     */
    double get_m_lim() {return m_lim;};

    /**
     * @brief returns \link gamma \endlink
     */
    double get_gamma() {return gamma;};

    /**
     * @brief returns \link alpha_ref \endlink
     */
    double get_alpha_ref() {return alpha_ref;};

    /**
     * @brief returns \link characteristic_length \endlink
     */
    double get_char_length() {return characteristic_length;};

    /**
     * @brief returns \link characteristic_solution_scale \endlink
     */
    double get_char_sol_scale() {return characteristic_solution_scale;};

    /**
     * @brief returns \link C0_CJM \endlink
     */
    double get_C0_CJM() {return C0_CJM;};

    /**
     * @brief returns \link C1_CJN \endlink
     */
    double get_C1_CJN() {return C1_CJN;};

    /**
     * @brief applies the limiter to the given FE function
     * @details Here the slope limiter is actually applied. It modifies the
     * values of the FE function.
     * @note Note that this works only if the FE function is a FE function for
     * discontinuous elements, i.e. if ANSATZ FUNCTION is -11, -12, ...
     */
    void limit_function(FEFunction* const fe_function);

    /**
     * @brief return a vector of booleans indicating which cells are marked
     * @details During the application of the limiter first cells are marked
     * where the function has to be modified. This is stored in an array that
     * can be accessed by this method. The i-th entry in the vector specifies
     * whether the i-th cell of the collection that underlies the given
     * FE function is marked by the algorithm or not.
     */
    std::vector<bool> get_is_cell_to_limit(const FEFunction* const fe_function);

    /**
     * @brief return the feature vector.
     * @details During the limiting process features are computed on which the
     * marking is based. This method computes the feature vector and returns it.
     * The i-th entry in the vector corresponds to the features of the i-th cell
     * in the collection that underlies the given FE function. For details of
     * the features have a look in the implementation.
     */
    std::vector<std::vector<double>> get_features(const FEFunction*
        const fe_function);

    /**
     * @brief prints some infos about the limiter.
     */
    void info();

    /** @brief list of implemented limiter
    */
    static std::array<std::string, 12> possible_limiters;

  private:
    /** @brief controls if the limiter is applied.
     * @details This boolean controls the behaviour of the class. If this
     * boolean is set to false, limit_function(), get_features() and
     * get_is_cell_to_limit() do nothing.
     */
    bool apply_limiter;

    /** @brief name of the limiter
     * @details Either one of
     * Galerkin, LinTriaReco, ConstTriaReco, LinQuadReco, ConstTriaReco,
     * LinQuadDeriv, ConstQuadDeriv, ConstJump, ConstJumpMod.
     */
    std::string limiter_name;

    /** verifies that the given limiter name is meaningful, i.e. implemented at
     * all.
     */
    void check_limiter_name(const std::string& limiter_name);

    /** @brief is the parameter \f$M_{\mathrm{lim}}\f$ that is needed for
     * LinQuadDeriv and ConstQuadDeriv.
     * @details See the paper mentioned in the description of the class for
     * details what \f$M_{\mathrm{lim}}\f$ is.
     */
    double m_lim;

    /** @brief is the parameter \f$\gamma\f$ that is needed for LinQuadDeriv
     * and ConstQuadDeriv.
     * @details See the paper mentioned in the description
     * of the class for details what \f$\gamma\f$ is.
     */
    double gamma;

    /** @brief is the parameter \f$\alpha_{\mathrm{ref}}\f$ that is needed for
     * ConstJumpMod.
     * @details See the paper mentioned in the description
     * of the class for details what \f$\alpha_{\mathrm{ref}}\f$ is.
     */
    double alpha_ref;

    /** @brief is the parameter \f$L\f$ that is needed for ConstJumpMod.
     *
     * @details See the description of the class for information what \f$L\f$ is.
     */
    double characteristic_length;

    /** @brief is the parameter \f$u_{\mathrm{ref}}\f$ that is needed for
     * ConstJumpMod.
     *
     * @details See the description of the class for information what
     * \f$u_{\mathrm{ref}}\f$ is.
     */
    double characteristic_solution_scale;

    /** @brief is the parameter \f$C_0\f$ that is needed for ConstJumpMod.
     *
     * @details See the description of the class for information what \f$C_0\f$
     * is.
     */
    double C0_CJM;

    /** @brief is the parameter \f$C_1\f$ that is needed for ConstJumpNorm.
     *
     * @details See the description of the class for information what \f$C_1\f$
     * is.
     */
    double C1_CJN;

    /**
     * @brief feature vector of the limiter. It is only filled if
     * compute_features() or limit_fe_function() is called
     * @details This vector contains the information needed for marking the
     * cells.
     * @note Be careful, during the marking and the actual limiting process the
     * feature vector is changed again. So if you only want to have the pure
     * features of the limiter you have to compute them by compute_features().
     */
    std::vector<std::vector<double>> features;

    /**
     * @brief the vector of marked cells. It is only filled if
     * compute_cells_to_limit() or limit_fe_function() is called
     */
    std::vector<bool> is_cell_to_limit;

    /**
     * @brief sets \link m_lim \endlink
     */
    void set_m_lim(const double& value);

    /**
     * @brief sets \link gamma \endlink
     */
    void set_gamma(const double& value);

    /**
     * @brief sets \link alpha_ref \endlink
     */
    void set_alpha_ref(const double& value);

    /**
     * @brief sets \link characteristic_length \endlink
     */
    void set_char_length(const double& value);

    /**
     * @brief sets \link characteristic_solution_scale \endlink
     */
    void set_char_solution_scale(const double& value);

    /**
     * @brief sets \link C0_CJM \endlink
     */
    void set_C0_CJM(const double& value);

    /**
     * @brief sets \link C1_CJN \endlink
     */
    void set_C1_CJN(const double& value);

    /** @brief populates the feature vector of the limiter based on the given FE
     * function
     */
    void compute_features(const FEFunction* const fe_function);

    /** @brief populates the is_cell_to_limit vector of the limiter based on the
     * given FE function
     */
    void compute_cells_to_limit(const FEFunction* const fe_function);

    /** @brief limits the given FE function, i.e. changes the entries according
     * to the limiter.
     */
    void limit_fe_function(FEFunction* const fe_function);
};

#endif
