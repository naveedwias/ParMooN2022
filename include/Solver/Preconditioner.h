#ifndef __PRECONDITIONER_H__
#define __PRECONDITIONER_H__

/** @brief an abstract base class to describe preconditioners
 * 
 */
template <class Vector>
class Preconditioner
{
  public:
    /// @brief use this object as a preconditioner
    virtual void apply(const Vector & z, Vector & r) const = 0;
    
    /** @brief Method to use this preconditioner in FGMRES.
     *
     * So far i and j get ignored and simply solve(z,r) is called.
     *
     * @param i Number of current iteration since last restart in FGMRES.
     * @param j Number of current iteration in FGMRES.
     * @param z The right hand side of the preconditioning.
     * @param r The obtained vector.
     */
    void apply(unsigned int i, unsigned int j, const Vector &z, Vector &r)
      const;
    
    /** @brief update this preconditioner
     * 
     * This sometimes saves computation time, and/or reduces reallocation. In 
     * general creating a new Preconditioner should work as well. Some
     * preconditioners do not need this, so there is a default implementation 
     * here.
     */
    virtual void update() {};
    
    virtual ~Preconditioner() = default;
};

template <class Vector>
class NoPreconditioner : public Preconditioner<Vector>
{
  public:
    
    virtual void apply(const Vector & z, Vector & r) const override final
    { r = z; }
    
    /** @brief Method to use this preconditioner in FGMRES.
     *
     * So far i and j get ignored and simply solve(z,r) is called.
     *
     * @param i Number of current iteration since last restart in FGMRES.
     * @param j Number of current iteration in FGMRES.
     * @param z The right hand side of the preconditioning.
     * @param r The obtained vector.
     */
    virtual void apply(unsigned int i, unsigned int j, const Vector &z,
                       Vector &r) const final;

    virtual ~NoPreconditioner() = default;
};


// implementations
template <class Vector>
void Preconditioner<Vector>::apply(unsigned int, unsigned int, const Vector &z,
                                   Vector &r) const
{
  apply(z, r);
}

template <class Vector>
void NoPreconditioner<Vector>::apply(unsigned int, unsigned int,
                                     const Vector &z, Vector &r) const
{
  apply(z, r);
}


#endif // __PRECONDITIONER_H__
