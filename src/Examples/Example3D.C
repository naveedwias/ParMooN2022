#include <Example3D.h>


ParameterDatabase Example3D::default_example_database()
{
  Output::print<5>("creating a default Example3D parameter database");
  // we use a parmoon default database because this way these parameters are
  // available in the default NSE3D database as well.
  ParameterDatabase db("Example3D parameter database");

  db.add("example", 0,
      "Choose which example to run. \nNote that depending on the type of "
      "problem you want to solve, different values are meaningful here. See "
      "the derived classes of 'Example3D'.", -5, 200);

  /** TDatabase::ParamDB->RE_NR */
  db.add("reynolds_number", 1.,
      "Reynolds number: dimensionless number which describes how viscous the  "
      "flow is (laminar, turbulent). Re = U.L/nu, where nu is the kinematic "
      "viscosity (=mu/rho)."
      "The higher it is, the more turbulent "
      "the flow is. Reynolds number can often be in the order of "
      "magnitude of millions. Then, the classical NSE is not relevant anymore"
      "to describe the flow. One should in these cases use turbulence Models,"
      "like Smagorinsky or k-eps. Maximum value is 1000,"
      "which already corresponds to a turbulent flow. Default value is 1."
      "Note that this is also equal to 1/eps, where eps is the DIMENSIONLESS "
      "viscosity (sometimes misleadingly named as nu) "
      "which sometimes appears in the dimensionless Navier-Stokes equations"
      "instead of 1/Re.",
      0., 1000000.);

  /** TDatabase::ParamDB->PE_NR */
  db.add("diffusion_coefficient", 1.,
      "Diffusion coefficient: a factor in front of the diffusion term.",
      0., 1.);

  /** TDatabase::ParamDB->PERMEABILITY */
  db.add("permeability", (double) 1.,
      "permeability coefficient: a factor in front of the resistance term.",
      (double) 0., (double) 1000000.);
  /** TDatabase::ParamDB->VISCOSITY */
  db.add("viscosity", (double) 1.,
      "viscosity coefficient: a factor in front of the Laplacian or the resistance term.",
      (double) 0., (double) 1000000.);
  /** TDatabase::ParamDB->EFFECTIVE_VISCOSITY */
  db.add("effective_viscosity", (double) 1.,
      "effective_viscosity coefficient: a factor in front of the Laplacian term.",
      (double) 0., (double) 1000000.);


  // NEW LB 10.10.18
  db.add("inverse_permeability", (double) 0.,
      "viscosity/permeability: a factor in front of the resistance term.",
      (double) 0., (double) 1000000.);

  ///@todo boundary-related parameters could be moved soon to a dedicated class
  // Neumann BC
  db.add("n_neumann_bd",0u,
   "Number of boundary components where Neumann BC are imposed"
   " This variable is used as a preliminary flag. The number of Neumann BC"
   " must match the length of the parameter neumann_id.",0u,9u);

  db.add("neumann_id", {1u,1u,1u,1u,1u,1u,1u,1u,1u,1u},
   "Component ID of Neumann boundaries. On these edges the terms (hat p,v.n)_E"
   " will be assembled. The values -hat p- are given via the parameters "
   " neumann_values. "
   " Attention: at the moment max. 10 boundaries are allowed. This can be "
   " easily changed in the file " __FILE__,0u,999u);

  db.add("neumann_value", {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
   " Neumann values to be assigned to the respetive boundary."
   " Warning: the length must be the same as the one given in neumann_boundary_id");

  // Nitsche BC
  db.add("n_nitsche_bd",0u,
   " Number of boundary components where Nitsche BC are imposed"
   " This variable is used as a preliminary flag. The number of Nitsche BC"
   " must match the length of the parameter nitsche_id",0u,9u);

  db.add("nitsche_id", {1u,1u,1u,1u,1u,1u,1u,1u,1u,1u},
   "Component ID of boundaries where the Nitsche method is used to "
   " impose (weak) Dirichlet BC."
   " Attention: at the moment max. 10 boundaries are allowed. This can be "
   " easily changed in the file " __FILE__,0u,999u);

  db.add("nitsche_penalty", {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
   " Nitsche penalty (gamma) to be imposed on each boundary");


  db.add("symmetric_nitsche_u",-1,
   " Coefficient for the term -(u,mu dv/dn) "
   " Symmetric version (symmetric_nitsche_u = 1): assemble -(mu du/dn,v)-(u,mu dv/dn), "
   " non-sym (symmetric_nitsche_u = -1): assemble -(mu du/dn,v)+(u,mu dv/dn)");

  db.add("symmetric_nitsche_p",-1,
   " Coefficient for the term (u.n,q) "
   " Symmetric version (symmetric_nitsche_p = 1): assemble (p,v.n)+(u.n,q), "
   " non-sym (symmetric_nitsche_u = -1): assemble (p,v.n)-(u.n,q)");

  // Windkessel (Lumped-0D models) BC
  db.add("n_windkessel_bd", 0u,
         " Number of boundary components where Windkessel (RCR) BC are imposed"
         " This variable is used as a preliminary flag. The number of"
         " Windkessel BC must match the length of the parameter windkessel_id.",
         0u, 9u);

  db.add("windkessel_Rp", {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
         " Proximal resistance to be assigned to the respective Windkessel"
         " boundary. Warning: the length must be the same as the one given in"
         " windkessel_id.",
         -100., 100.);

  db.add("windkessel_Rd", {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
         " Distal resistance to be assigned to the respective Windkessel"
         " boundary. Warning: the length must be the same as the one given in"
         " windkessel_id.", -100., 100.);

  db.add("windkessel_C", {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
         " Capacitance to be assigned to the respective Windkessel"
         " boundary. Warning: the length must be the same as the one given in"
         " windkessel_id.", -100., 100.);

  db.add("windkessel_id", {1u,1u,1u,1u,1u,1u,1u,1u,1u,1u},
         "Component ID of windkessel boundaries. On these edges the terms"
         " (hat p,v.n)_E will be assembled. The values -hat p- are computed"
         " from a windkessel model. Attention: at the moment max. 10 boundaries"
         " are allowed. This can be easily changed in the file " __FILE__,
         0u, 999u);

  db.add("windkessel_distal_pressure", 0.,
         " Pressure reference to compute proximal pressure drop in Windkessel"
         " model: dP = (P - windkessel_distal_pressure) = Q * Rp.");

  db.add("windkessel_initial_distal_pressure", {0.,0., 0., 0., 0., 0., 0., 0., 0., 0.},
         " Initial distal pressure for Windkessel model");

  db.add("windkessel_damping_limiter", 0.,
         " Limiter in (0,1] to improve convergence with Windkessel BC."
         " |P_new - P_old| < windkessel_damping_limiter * P_old"
         " Zero means no limiter is used.", 0., 1.);

  db.add("windkessel_theta", 1.0,
    "Use extrapolated velocity field to compute outflows when the right-hand "
    "side is not recomputed for every nonlinear iteration.",
    0.0, 1.0);

  // NSE: directional do-nothing BC

  db.add("n_directional_do_nothing_bd", 0u,
    "Number of boundary components where directional Neumann BC are imposed. "
    "This is *only* the directional component. You must also set Neumann or "
    "Windkessel BC on the same boundary components.\n"
    "This adds g/2 (u \\cdot n)_{-} u to the RHS of the boundary condition, "
    "where the value of g is given by the directional_value array (but "
    "defaults to 1)\n"
    "Currently only implemented in 3D TNSE.",
    0u, 9u);

  db.add("directional_do_nothing_id", {1u, 1u, 1u, 1u},
    "Component ID of directional Neumann boundaries.",
    0u, 999u);

  db.add("directional_do_nothing_value", {1.0},
    "Directional component scaling.");

  db.add("directional_do_nothing_delta", 1.0e-3,
    "Sharpness of the smooth step if directional_do_nothing_type is "
    "\"smooth\". A value of zero corresponds to a sharp step.",
    0.0, 1.0);

  db.add("directional_do_nothing_D0", 0.0,
    "Scale of the time derivative term if directional_do_nothing_type is "
    "\"smooth\".",
    0.0, 10.e+3);

  db.add("directional_do_nothing_type", "sharp",
    "Directional component type.",
    { "sharp", "smooth" });

  db.add("directional_do_nothing_theta", 0.0,
    "Extrapolation parameter for the velocity vector to use in the "
    "directional do-nothing RHS. Blends between the previous step's velocity "
    "at theta = 0 and an extrapolation based on the previous two steps at "
    "theta = 1.",
    0.0, 1.0);

  db.add("darcy_permeability_jump", 1., "Some examples for Darcy equations "
         "allow jumping permeability coefficient. Typically this affects a "
         "certain part of the domain. This parameter is the factor with wich "
         " the permeability is scaled inside that subdomain.", 0., 1.e10);

   db.add("degree_polynomial", 0u, "Parameter defining the maximal degree of "
       "the polynomial used in example Polynomial.h. Specify the degree with "
       "this parameter. In principle also degrees greater 50 would be "
       "possible. Since the computations then are very slow the degree is "
       "restricted to 50.", 0u, 50u);
  return db;
}

  Example3D::Example3D(const ParameterDatabase & db)
: example_database(Example3D::default_example_database()), exact_solution(),
  boundary_conditions(), boundary_data(), problem_coefficients(nullptr)
{
  this->example_database.merge(db, false);
}

Example3D::Example3D(const std::vector<DoubleFunct3D*>& exact,
                     const std::vector<BoundCondFunct3D*>& bc,
                     const std::vector<BoundValueFunct3D*>& bd,
                     const CoeffFct3D& coeffs)
: example_database(Example3D::default_example_database()), exact_solution(exact),
  boundary_conditions(bc), boundary_data(bd), problem_coefficients(coeffs)
{ 
}
