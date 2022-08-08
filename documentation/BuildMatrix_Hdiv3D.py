import numpy as np
import sympy as sy

"""
This script calculates the coefficients of basis functions for a number of
H(div) conforming elements. In particular with the methods build_matrix_Hdiv_Q,
and build_matrix_Hdiv_T you can get the coefficients for Raviart-Thomas (RT)
and Brezzi-Douglas-Marini (BDM) elements on the reference cube [-1,1]^3 and
the reference tetrahedron whose vertices are (0,0), (1,0), and (0,1).

For each element the space consists of a number of basis functions (mainly
monomials) who are represented by expressions whose name starts with 'V'. The
nodal functionals (the degrees of freedom) are integrals over either the edges
of the reference cell or over the cell itself and are represented by methods
whose names start with NQ or NT. These integrals have an additional weighting
function which is specific to each degree of freedom. This weight functions
for the boundary degrees of freedom are the Legendre polynomials and are
represented by the methods q1, q2, q2, q4. The weights for the inner degrees of
freedom have names starting with 'qi'.

Let N_i be the degrees of freedom (functionals) and v_j the monomials for one
of the considered elements. Then this script computes (in
compute_matrix_and_inverse) the matrix with entries N_i(v_j) and inverts this
matrix. Its inverse can then be copied into ParMooN, see e.g. BF_N_Q_RT2_2D.h.
Its columns are the coefficients (with respect to the above monomials) needed
to define the corresponding basis function. The order of the monomials is
important and should be the same in ParMooN as it is here! Similarly the
degrees of freedom have to be defined in the same order in ParMooN (see e.g.
NF_N_Q_RT2_2D.h) as they are defined here.

You can call this script via
python3 BuildMatrix_Hdiv.py
It will print all the matrices which are used ParMoon.
"""

x, y, z, t, s = sy.symbols('x y z t s', real=True)

# %% Defining trivial basis functions (monomials)

V0 = sy.S.Zero
V1 = sy.S.One
# first order
Vx = x
Vy = y
Vz = z
# second order
Vxx = x**2
Vxy = x*y
Vxz = x*z
Vyy = y**2
Vyz = y*z
Vzz = z**2
# third order
Vxxx = x**3
Vxxy = x**2*y
Vxxz = x**2*z
Vxyy = x*y**2
Vxyz = x*y*z
Vxzz = x*z**2
Vyyy = y**3
Vyyz = y**2*z
Vyzz = y*z**2
Vzzz = z**3
# fourth order
Vxxxx = x**4
Vxxxy = x**3*y
Vxxxz = x**3*z
Vxxyy = x**2*y**2
Vxxyz = x**2*y*z
Vxxzz = x**2*z**2
Vxyyy = x*y**3
Vxyyz = x*y**2*z
Vxyzz = x*y*z**2
Vxzzz = x*z**3
Vyyyy = y**4
Vyyyz = y**3*z
Vyyzz = y**2*z**2
Vyzzz = y*z**3
Vzzzz = z**4

# fifth order
Vxxxxx = x**5
Vxxxxy = x**4*y
Vxxxxz = x**4*z
Vxxxyy = x**3*y**2
Vxxxyz = x**3*y*z
Vxxxzz = x**3*z**2
Vxxyyy = x**2*y**3
Vxxyyz = x**2*y**2*z
Vxxyzz = x**2*y*z**2
Vxxzzz = x**2*z**3
Vxyyyy = x*y**4
Vxyyyz = x*y**3*z
Vxyyzz = x*y**2*z**2
Vxyzzz = x*y*z**3
Vxzzzz = x*z**4
Vyyyyy = y**5
Vyyyyz = y**4*z
Vyyyzz = y**3*z**2
Vyyzzz = y**2*z**3
Vyzzzz = y*z**4
Vzzzzz = z**5
# sixth order (not all)
Vxxxyyz = x**3*y**2*z
Vxxxyzz = x**3*y*z**2
Vxxyyyz = x**2*y**3*z
Vxxyyzz = x**2*y**2*z**2
Vxxyzzz = x**2*y*z**3
Vxyyyzz = x*y**3*z**2
Vxyyzzz = x*y**2*z**3
# seventh order (not all)
Vxxxyyzz = x**3*y**2*z**2
Vxxyyyzz = x**2*y**3*z**2
Vxxyyzzz = x**2*y**2*z**3


# %% defining the test polynomials for degrees of freedom on edges
# Legendre polynomials.
def q(k, l):
    return sy.legendre(k, t) * sy.legendre(l, s)


# defining test polynomials for inner degrees of freedom
qi100 = (sy.S.One, sy.S.Zero, sy.S.Zero)
qi010 = (sy.S.Zero, sy.S.One, sy.S.Zero)
qi001 = (sy.S.Zero, sy.S.Zero, sy.S.One)

qix00 = (x, sy.S.Zero, sy.S.Zero)
qiy00 = (y, sy.S.Zero, sy.S.Zero)
qiz00 = (z, sy.S.Zero, sy.S.Zero)
qi0x0 = (sy.S.Zero, x, sy.S.Zero)
qi0y0 = (sy.S.Zero, y, sy.S.Zero)
qi0z0 = (sy.S.Zero, z, sy.S.Zero)
qi00x = (sy.S.Zero, sy.S.Zero, x)
qi00y = (sy.S.Zero, sy.S.Zero, y)
qi00z = (sy.S.Zero, sy.S.Zero, z)

qixx00 = (x**2, sy.S.Zero, sy.S.Zero)
qixy00 = (x*y, sy.S.Zero, sy.S.Zero)
qixz00 = (x*z, sy.S.Zero, sy.S.Zero)
qiyy00 = (y**2, sy.S.Zero, sy.S.Zero)
qiyz00 = (y*z, sy.S.Zero, sy.S.Zero)
qizz00 = (z**2, sy.S.Zero, sy.S.Zero)
qi0xx0 = (sy.S.Zero, x**2, sy.S.Zero)
qi0xy0 = (sy.S.Zero, x*y, sy.S.Zero)
qi0xz0 = (sy.S.Zero, x*z, sy.S.Zero)
qi0yy0 = (sy.S.Zero, y**2, sy.S.Zero)
qi0yz0 = (sy.S.Zero, y*z, sy.S.Zero)
qi0zz0 = (sy.S.Zero, z**2, sy.S.Zero)
qi00xx = (sy.S.Zero, sy.S.Zero, x**2)
qi00xy = (sy.S.Zero, sy.S.Zero, x*y)
qi00xz = (sy.S.Zero, sy.S.Zero, x*z)
qi00yy = (sy.S.Zero, sy.S.Zero, y**2)
qi00yz = (sy.S.Zero, sy.S.Zero, y*z)
qi00zz = (sy.S.Zero, sy.S.Zero, z**2)

qixxx00 = (x**3, sy.S.Zero, sy.S.Zero)
qixxy00 = (x**2*y, sy.S.Zero, sy.S.Zero)
qixxz00 = (x**2*z, sy.S.Zero, sy.S.Zero)
qixyy00 = (x*y**2, sy.S.Zero, sy.S.Zero)
qixyz00 = (x*y*z, sy.S.Zero, sy.S.Zero)
qixzz00 = (x*z**2, sy.S.Zero, sy.S.Zero)
qiyyy00 = (y**3, sy.S.Zero, sy.S.Zero)
qiyyz00 = (y**2*z, sy.S.Zero, sy.S.Zero)
qiyzz00 = (y*z**2, sy.S.Zero, sy.S.Zero)
qizzz00 = (z**3, sy.S.Zero, sy.S.Zero)
qi0xxx0 = (sy.S.Zero, x**3, sy.S.Zero)
qi0xxy0 = (sy.S.Zero, x**2*y, sy.S.Zero)
qi0xxz0 = (sy.S.Zero, x**2*z, sy.S.Zero)
qi0xyy0 = (sy.S.Zero, x*y**2, sy.S.Zero)
qi0xyz0 = (sy.S.Zero, x*y*z, sy.S.Zero)
qi0xzz0 = (sy.S.Zero, x*z**2, sy.S.Zero)
qi0yyy0 = (sy.S.Zero, y**3, sy.S.Zero)
qi0yyz0 = (sy.S.Zero, y**2*z, sy.S.Zero)
qi0yzz0 = (sy.S.Zero, y*z**2, sy.S.Zero)
qi0zzz0 = (sy.S.Zero, z**3, sy.S.Zero)
qi00xxx = (sy.S.Zero, sy.S.Zero, x**3)
qi00xxy = (sy.S.Zero, sy.S.Zero, x**2*y)
qi00xxz = (sy.S.Zero, sy.S.Zero, x**2*z)
qi00xyy = (sy.S.Zero, sy.S.Zero, x*y**2)
qi00xyz = (sy.S.Zero, sy.S.Zero, x*y*z)
qi00xzz = (sy.S.Zero, sy.S.Zero, x*z**2)
qi00yyy = (sy.S.Zero, sy.S.Zero, y**3)
qi00yyz = (sy.S.Zero, sy.S.Zero, y**2*z)
qi00yzz = (sy.S.Zero, sy.S.Zero, y*z**2)
qi00zzz = (sy.S.Zero, sy.S.Zero, z**3)

qi00xxyy = (sy.S.Zero, sy.S.Zero, x**2*y**2)
qi0xxyz0 = (sy.S.Zero, x**2*y*z, sy.S.Zero)
qi00xxyz = (sy.S.Zero, sy.S.Zero, x**2*y*z)
qi0xxzz0 = (sy.S.Zero, x**2*z**2, sy.S.Zero)
qixyyz00 = (x*y**2*z, sy.S.Zero, sy.S.Zero)
qi00xyyz = (sy.S.Zero, sy.S.Zero, x*y**2*z)
qixyzz00 = (x*y*z**2, sy.S.Zero, sy.S.Zero)
qi0xyzz0 = (sy.S.Zero, x*y*z**2, sy.S.Zero)
qiyyzz00 = (y**2*z**2, sy.S.Zero, sy.S.Zero)
qixyyzz00 = (x*y**2*z**2, sy.S.Zero, sy.S.Zero)
qi0xxyzz0 = (sy.S.Zero, x**2*y*z**2, sy.S.Zero)
qi00xxyyz = (sy.S.Zero, sy.S.Zero, x**2*y**2*z)


# %% defining the nodal functionals on hexahedra
def NH1(v, q):
    integrand = -v[2].subs({z: -1}) * q.subs({t: x, s: y})
    return sy.integrate(integrand, (x, -1, 1), (y, -1, 1))


def NH2(v, q):
    integrand = -v[1].subs({y: -1}) * q.subs({t: z, s: x})
    return sy.integrate(integrand, (x, -1, 1), (z, -1, 1))


def NH3(v, q):
    integrand = v[0].subs({x: 1}) * q.subs({t: y, s: z})
    return sy.integrate(integrand, (y, -1, 1), (z, -1, 1))


def NH4(v, q):
    integrand = v[1].subs({y: 1}) * q.subs({t: x, s: z})
    return sy.integrate(integrand, (x, -1, 1), (z, -1, 1))


def NH5(v, q):
    integrand = -v[0].subs({x: -1}) * q.subs({t: y, s: z})
    return sy.integrate(integrand, (y, -1, 1), (z, -1, 1))


def NH6(v, q):
    integrand = v[2].subs({z: 1}) * q.subs({t: x, s: y})
    return sy.integrate(integrand, (x, -1, 1), (y, -1, 1))


def NHI(v, q):
    integrand = v[0]*q[0] + v[1]*q[1] + v[2]*q[2]
    return sy.integrate(integrand, (x, -1, 1), (y, -1, 1), (z, -1, 1))


# defining the nodal functionals on tetrahedra
def NT1(v, q):
    integrand = -v[2].subs({x: t, y: s, z: 0})*q
    return 2*sy.integrate(integrand, (s, 0, 1-t), (t, 0, 1))


def NT2(v, q):
    integrand = -v[1].subs({x: t, y: 0, z: s})*q
    return 2*sy.integrate(integrand, (s, 0, 1-t), (t, 0, 1))


def NT3(v, q):
    integrand = ((v[0] + v[1] + v[2])*q).subs({x: t, y: 1-t-s, z: s})
    return 2*sy.integrate(integrand, (s, 0, 1-t), (t, 0, 1))


def NT4(v, q):
    integrand = -v[0].subs({x: 0, y: t, z: s})*q
    return 2*sy.integrate(integrand, (s, 0, 1-t), (t, 0, 1))


def NTI(v, q):
    integrand = v[0]*q[0] + v[1]*q[1] + v[2]*q[2]
    return sy.integrate(sy.integrate(integrand, (x, 0, 1-y-z)), (y, 0, 1-z),
                        (z, 0, 1))


def null(A, eps=1e-15):
    u, s, vh = np.linalg.svd(A)
    null_mask = (s <= eps)
    null_space = np.compress(null_mask, vh, axis=0)
    return np.transpose(null_space)


# %%
def compute_matrix_and_inverse(dof_list, mon_list, test_polynomials):
    # first a few checks
    dim = len(dof_list)
    assert dim == len(mon_list)
    assert dim == len(test_polynomials)
    # Compute the matrix
    A = np.matrix(np.zeros(shape=(dim, dim)))
    for i in range(dim):
        # print("i", i)
        for j in range(dim):
            val = dof_list[i](mon_list[j], test_polynomials[i])
            # print("j", j, val)
            A[i, j] = val
    np.set_printoptions(linewidth=420, threshold=20000, suppress=True,
                        precision=8)
    # print(27*A)
    # print("determinant ", np.linalg.det(A))
    # print(null(A))
    Ainv = np.linalg.inv(A)
    # set very small entries to exactly zero
    Ainv[abs(Ainv) < 1e-10] = 0
    return (A, Ainv)


# %%
def build_matrix_Hdiv_H(k):
    """
    compute Matrix of coefficients and its inverse to get the Raviart-Thomas
    (k>=0) or Brezzi-Douglas-Duran-Fortin (k<0) basis functions on the cube
    H = [-1,1]^3.

    Input: k denotes the order (k>=0 for RT, k<0 for BDDF)
    Output: Matrix A=(a_ij) with a_ij = N_i(v_j) and its inverse. Here N_i
            are the degrees of freedom (functionals) and v_j are monomial
            basis vectors. The Raviart-Thomas or Brezzi-Douglas-Duran-Fortin
            basis functions phi_j are then defined as phi_j = sum_i c_ij v_i.
            The monomial basis vectors are stored in the local variable
            mon_list.
    """

    if(k >= 3 or k <= -3):
        print("So far the order k can only be 0, 1, or 2 for ",
              "Raviart-Thomas elements and -1, -2, or -3 for ",
              "Brezzi-Douglas-Duran-Fortin elements")
        print("No matrix computed!")
        A = 0
        Ainv = 0
        return (A, Ainv)

    if(k == 0):  # 6 dofs
        dof_list = [NH1, NH2, NH3, NH4, NH5, NH6]
        mon_list = [(V1, V0, V0), (V0, V1, V0), (V0, V0, V1),
                    (Vx, V0, V0), (V0, Vy, V0), (V0, V0, Vz)]
        test_polynomials = [q(0, 0)]*6
    elif(k == 1):  # 36 dofs
        dof_list = [NH1]*4 + [NH2]*4 + [NH3]*4 + [NH4]*4 + [NH5]*4 + [NH6]*4
        dof_list += [NHI]*12
        mon_list = [(V1, V0, V0), (V0, V1, V0), (V0, V0, V1),
                    (Vx, V0, V0), (V0, Vx, V0), (V0, V0, Vx),
                    (Vy, V0, V0), (V0, Vy, V0), (V0, V0, Vy),
                    (Vz, V0, V0), (V0, Vz, V0), (V0, V0, Vz),
                    (Vxx, V0, V0), (V0, Vyy, V0), (V0, V0, Vzz),
                    (Vxy, V0, V0), (V0, Vxy, V0), (V0, V0, Vxy),
                    (Vxz, V0, V0), (V0, Vxz, V0), (V0, V0, Vxz),
                    (Vyz, V0, V0), (V0, Vyz, V0), (V0, V0, Vyz),
                    (Vxxy, V0, V0), (V0, Vxyy, V0), (V0, V0, Vxzz),
                    (Vxxz, V0, V0), (V0, Vyyz, V0), (V0, V0, Vyzz),
                    (Vxyz, V0, V0), (V0, Vxyz, V0), (V0, V0, Vxyz),
                    (Vxxyz, V0, V0), (V0, Vxyyz, V0), (V0, V0, Vxyzz)]
        test_polynomials = [q(0, 0), q(1, 0), q(0, 1), q(1, 1)]*6
        test_polynomials = [1+t+s+t*s,  1+t-s-t*s, 1-t-s+t*s, 1-t+s-t*s]*6
        test_polynomials += [qi100, qi010, qi001, qiy00, qiz00, qiyz00,
                             qi0x0, qi0z0, qi0xz0, qi00x, qi00y, qi00xy]
    elif(k == 2):  # 108 dofs
        dof_list = [NH1]*9 + [NH2]*9 + [NH3]*9 + [NH4]*9 + [NH5]*9 + [NH6]*9
        dof_list += [NHI]*54
        mon_list = [(V1, V0, V0), (V0, V1, V0), (V0, V0, V1),
                    (Vx, V0, V0), (V0, Vx, V0), (V0, V0, Vx),
                    (Vy, V0, V0), (V0, Vy, V0), (V0, V0, Vy),
                    (Vz, V0, V0), (V0, Vz, V0), (V0, V0, Vz),
                    (Vxx, V0, V0), (V0, Vxx, V0), (V0, V0, Vxx),
                    (Vxy, V0, V0), (V0, Vxy, V0), (V0, V0, Vxy),
                    (Vxz, V0, V0), (V0, Vxz, V0), (V0, V0, Vxz),
                    (Vyy, V0, V0), (V0, Vyy, V0), (V0, V0, Vyy),
                    (Vyz, V0, V0), (V0, Vyz, V0), (V0, V0, Vyz),
                    (Vzz, V0, V0), (V0, Vzz, V0), (V0, V0, Vzz),
                    (Vxxx, V0, V0), (V0, Vyyy, V0), (V0, V0, Vzzz),
                    (Vxxy, V0, V0), (V0, Vxxy, V0), (V0, V0, Vxxy),
                    (Vxxz, V0, V0), (V0, Vxxz, V0), (V0, V0, Vxxz),
                    (Vxyy, V0, V0), (V0, Vxyy, V0), (V0, V0, Vxyy),
                    (Vxyz, V0, V0), (V0, Vxyz, V0), (V0, V0, Vxyz),
                    (Vxzz, V0, V0), (V0, Vxzz, V0), (V0, V0, Vxzz),
                    (Vyyz, V0, V0), (V0, Vyyz, V0), (V0, V0, Vyyz),
                    (Vyzz, V0, V0), (V0, Vyzz, V0), (V0, V0, Vyzz),
                    (Vxxxy, V0, V0), (V0, Vxyyy, V0), (V0, V0, Vxzzz),
                    (Vxxxz, V0, V0), (V0, Vyyyz, V0), (V0, V0, Vyzzz),
                    (Vxxyy, V0, V0), (V0, Vxxyy, V0), (V0, V0, Vxxyy),
                    (Vxxyz, V0, V0), (V0, Vxxyz, V0), (V0, V0, Vxxyz),
                    (Vxxzz, V0, V0), (V0, Vxxzz, V0), (V0, V0, Vxxzz),
                    (Vxyyz, V0, V0), (V0, Vxyyz, V0), (V0, V0, Vxyyz),
                    (Vxyzz, V0, V0), (V0, Vxyzz, V0), (V0, V0, Vxyzz),
                    (Vyyzz, V0, V0), (V0, Vyyzz, V0), (V0, V0, Vyyzz),
                    (Vxxxyy, V0, V0), (V0, Vxxyyy, V0), (V0, V0, Vxxzzz),
                    (Vxxxyz, V0, V0), (V0, Vxyyyz, V0), (V0, V0, Vxyzzz),
                    (Vxxxzz, V0, V0), (V0, Vyyyzz, V0), (V0, V0, Vyyzzz),
                    (Vxxyyz, V0, V0), (V0, Vxxyyz, V0), (V0, V0, Vxxyyz),
                    (Vxxyzz, V0, V0), (V0, Vxxyzz, V0), (V0, V0, Vxxyzz),
                    (Vxyyzz, V0, V0), (V0, Vxyyzz, V0), (V0, V0, Vxyyzz),
                    (Vxxxyyz, V0, V0), (V0, Vxxyyyz, V0), (V0, V0, Vxxyzzz),
                    (Vxxxyzz, V0, V0), (V0, Vxyyyzz, V0), (V0, V0, Vxyyzzz),
                    (Vxxyyzz, V0, V0), (V0, Vxxyyzz, V0), (V0, V0, Vxxyyzz),
                    (Vxxxyyzz, V0, V0), (V0, Vxxyyyzz, V0), (V0, V0, Vxxyyzzz)]
        test_polynomials = [q(0, 0), q(1, 0), q(2, 0),
                            q(0, 1), q(1, 1), q(2, 1),
                            q(0, 2), q(1, 2), q(2, 2)]*6
        test_polynomials += [qi100, qix00, qiy00, qiz00,
                             qixy00, qixz00, qiyy00, qiyz00, qizz00,
                             qixyy00, qixyz00, qixzz00, qiyyz00, qiyzz00,
                             qixyyz00, qixyzz00, qiyyzz00, qixyyzz00,
                             qi010, qi0x0, qi0y0, qi0z0,
                             qi0xx0, qi0xy0, qi0xz0, qi0yz0, qi0zz0,
                             qi0xxy0, qi0xxz0, qi0xyz0, qi0xzz0, qi0yzz0,
                             qi0xxyz0, qi0xxzz0, qi0xyzz0, qi0xxyzz0,
                             qi001, qi00x, qi00y, qi00z,
                             qi00xx, qi00xy, qi00xz, qi00yy, qi00yz,
                             qi00xxy, qi00xxz, qi00xyy, qi00xyz, qi00yyz,
                             qi00xxyy, qi00xxyz, qi00xyyz, qi00xxyyz]
    elif(k == -1):  # 18 dofs
        dof_list = [NH1]*3 + [NH2]*3 + [NH3]*3 + [NH4]*3 + [NH5]*3 + [NH6]*3
        mon_list = [(V1, V0, V0), (V0, V1, V0), (V0, V0, V1),
                    (Vx, V0, V0), (V0, Vx, V0), (V0, V0, Vx),
                    (Vy, V0, V0), (V0, Vy, V0), (V0, V0, Vy),
                    (Vz, V0, V0), (V0, Vz, V0), (V0, V0, Vz),
                    (Vxz, -Vyz, V0), (2*Vxy, -Vyy, V0),
                    (-Vxy, V0, Vyz), (-Vxx, V0, 2*Vxz),
                    (V0, Vxy, -Vxz), (V0, 2*Vyz, Vzz)]
        test_polynomials = [sy.S.One, t, s]*6
    elif(k == -2):  # 39 dofs
        dof_list = [NH1]*6 + [NH2]*6 + [NH3]*6 + [NH4]*6 + [NH5]*6 + [NH6]*6
        dof_list += [NHI]*3
        mon_list = [(V1, V0, V0), (V0, V1, V0), (V0, V0, V1),
                    (Vx, V0, V0), (V0, Vx, V0), (V0, V0, Vx),
                    (Vy, V0, V0), (V0, Vy, V0), (V0, V0, Vy),
                    (Vz, V0, V0), (V0, Vz, V0), (V0, V0, Vz),
                    (Vxx, V0, V0), (V0, Vxx, V0), (V0, V0, Vxx),
                    (Vxy, V0, V0), (V0, Vxy, V0), (V0, V0, Vxy),
                    (Vxz, V0, V0), (V0, Vxz, V0), (V0, V0, Vxz),
                    (Vyy, V0, V0), (V0, Vyy, V0), (V0, V0, Vyy),
                    (Vyz, V0, V0), (V0, Vyz, V0), (V0, V0, Vyz),
                    (Vzz, V0, V0), (V0, Vzz, V0), (V0, V0, Vzz),
                    (Vxzz, -Vyzz, V0), (2*Vxyz, -Vyyz, V0),
                    (3*Vxyy, -Vyyy, V0),
                    (-Vxyy, V0, Vyyz), (-Vxxy, V0, 2*Vxyz),
                    (-Vxxx, V0, 3*Vxxz),
                    (V0, Vxxy, -Vxxz), (V0, 2*Vxyz, -Vxzz),
                    (V0, 3*Vyzz, -Vzzz)]
        test_polynomials = [sy.S.One, t, s, t**2, t*s, s**2]*6
        test_polynomials += [qi100, qi010, qi001]
    elif(k == -3):  # 22 dofs
        dof_list = [NH1]*10 + [NH2]*10 + [NH3]*10 + [NH4]*10 + [NH5]*10
        dof_list += [NH6]*10
        dof_list += [NHI]*12
        mon_list = [(V1, V0, V0), (V0, V1, V0), (V0, V0, V1),
                    (Vx, V0, V0), (V0, Vx, V0), (V0, V0, Vx),
                    (Vy, V0, V0), (V0, Vy, V0), (V0, V0, Vy),
                    (Vz, V0, V0), (V0, Vz, V0), (V0, V0, Vz),
                    (Vxx, V0, V0), (V0, Vxx, V0), (V0, V0, Vxx),
                    (Vxy, V0, V0), (V0, Vxy, V0), (V0, V0, Vxy),
                    (Vxz, V0, V0), (V0, Vxz, V0), (V0, V0, Vxz),
                    (Vyy, V0, V0), (V0, Vyy, V0), (V0, V0, Vyy),
                    (Vyz, V0, V0), (V0, Vyz, V0), (V0, V0, Vyz),
                    (Vzz, V0, V0), (V0, Vzz, V0), (V0, V0, Vzz),
                    (Vxxx, V0, V0), (V0, Vxxx, V0), (V0, V0, Vxxx),
                    (Vxxy, V0, V0), (V0, Vxxy, V0), (V0, V0, Vxxy),
                    (Vxxz, V0, V0), (V0, Vxxz, V0), (V0, V0, Vxxz),
                    (Vxyy, V0, V0), (V0, Vxyy, V0), (V0, V0, Vxyy),
                    (Vxyz, V0, V0), (V0, Vxyz, V0), (V0, V0, Vxyz),
                    (Vxzz, V0, V0), (V0, Vxzz, V0), (V0, V0, Vxzz),
                    (Vyyy, V0, V0), (V0, Vyyy, V0), (V0, V0, Vyyy),
                    (Vyyz, V0, V0), (V0, Vyyz, V0), (V0, V0, Vyyz),
                    (Vyzz, V0, V0), (V0, Vyzz, V0), (V0, V0, Vyzz),
                    (Vzzz, V0, V0), (V0, Vzzz, V0), (V0, V0, Vzzz),
                    (Vxzzz, -Vyzzz, V0), (2*Vxyzz, -Vyyzz, V0),
                    (3*Vxyyz, -Vyyyz, V0), (4*Vxyyy, -Vyyyy, V0),
                    (-Vxyyy, V0, Vyyyz), (-Vxxyy, V0, 2*Vxyyz),
                    (-Vxxxy, V0, 3*Vxxyz), (-Vxxxx, V0, 4*Vxxxz),
                    (V0, Vxxxy, -Vxxxz), (V0, 2*Vxxyz, -Vxxzz),
                    (V0, 3*Vxyzz, -Vxzzz), (V0, 4*Vyzzz, -Vzzzz)]
        test_polynomials = [sy.S.One, t, s, t**2, t*s, s**2,
                            t**3, t**2*s, t*s**2, s**3]*6
        test_polynomials += [qi100, qi010, qi001, qix00, qi0x0, qi00x,
                             qiy00, qi0y0, qi00y, qiz00, qi0z0, qi00z]

    a, ainv = compute_matrix_and_inverse(dof_list, mon_list, test_polynomials)
    return a, ainv, mon_list

###############################################################################
# triangles


# %%
def build_matrix_Hdiv_T(k):
    """
    compute Matrix of coefficients and its inverse to get the Raviart-Thomas
    (k>=0) or Brezzi-Douglas-Marini (k<0) basis functions on the tetrahedron
    with vertices (0,0,0), (1,0,0), (0,1,0), (0,0,1).

    Input: k denotes the order (k>=0 for RT, k<0 for BDM)
    Output: Matrix A=(a_ij) with a_ij = N_i(v_j) and its inverse. Here N_i
            are the degrees of freedom (functionals) and v_j are monomial
            basis vectors. The Raviart-Thomas or Brezzi-Douglas-Marini basis
            functions phi_j are then defined as phi_j = sum_i c_ij v_i. The
            monomial basis vectors are stored in the local variable mon_list.
    """

    if(k >= 4 or k <= -4):
        print("So far the order k can only be 0, 1, 2, or 3 for",
              "Raviart-Thomas elements and -1, -2, or -3 for",
              "Brezzi-Douglas-Marini elements")
        print("No matrix computed!")
        A = 0
        Ainv = 0
        return (A, Ainv)

    if(k == 0):  # 4 dofs
        dof_list = [NT1, NT2, NT3, NT4]
        mon_list = [(V1, V0, V0), (V0, V1, V0), (V0, V0, V1), (Vx, Vy, Vz)]
        test_polynomials = [sy.S.One]*4
    elif(k == 1):  # 15 dofs
        dof_list = [NT1]*3 + [NT2]*3 + [NT3]*3 + [NT4]*3 + [NTI]*3
        mon_list = [(V1, V0, V0), (V0, V1, V0), (V0, V0, V1),
                    (Vx, V0, V0), (V0, Vx, V0), (V0, V0, Vx),
                    (Vy, V0, V0), (V0, Vy, V0), (V0, V0, Vy),
                    (Vz, V0, V0), (V0, Vz, V0), (V0, V0, Vz),
                    (Vxx, Vxy, Vxz), (Vxy, Vyy, Vyz), (Vxz, Vyz, Vzz)]
        test_polynomials = [sy.S.One, t, s]*4 + [qi100, qi010, qi001]
    elif(k == 2):  # 36 dofs
        dof_list = [NT1]*6 + [NT2]*6 + [NT3]*6 + [NT4]*6 + [NTI]*12
        mon_list = [(V1, V0, V0), (V0, V1, V0), (V0, V0, V1),
                    (Vx, V0, V0), (V0, Vx, V0), (V0, V0, Vx),
                    (Vy, V0, V0), (V0, Vy, V0), (V0, V0, Vy),
                    (Vz, V0, V0), (V0, Vz, V0), (V0, V0, Vz),
                    (Vxx, V0, V0), (V0, Vxx, V0), (V0, V0, Vxx),
                    (Vxy, V0, V0), (V0, Vxy, V0), (V0, V0, Vxy),
                    (Vxz, V0, V0), (V0, Vxz, V0), (V0, V0, Vxz),
                    (Vyy, V0, V0), (V0, Vyy, V0), (V0, V0, Vyy),
                    (Vyz, V0, V0), (V0, Vyz, V0), (V0, V0, Vyz),
                    (Vzz, V0, V0), (V0, Vzz, V0), (V0, V0, Vzz),
                    (Vxxx, Vxxy, Vxxz), (Vxxy, Vxyy, Vxyz), (Vxxz, Vxyz, Vxzz),
                    (Vxyy, Vyyy, Vyyz), (Vxyz, Vyyz, Vyzz), (Vxzz, Vyzz, Vzzz)]
        test_polynomials = [sy.S.One, t, s, t**2, t*s, s**2]*4
        test_polynomials += [qi100, qi010, qi001, qix00, qi0x0, qi00x,
                             qiy00, qi0y0, qi00y, qiz00, qi0z0, qi00z]
    elif(k == 3):  # 70 dofs
        dof_list = [NT1]*10 + [NT2]*10 + [NT3]*10 + [NT4]*10 + [NTI]*30
        mon_list = [(V1, V0, V0), (V0, V1, V0), (V0, V0, V1),
                    (Vx, V0, V0), (V0, Vx, V0), (V0, V0, Vx),
                    (Vy, V0, V0), (V0, Vy, V0), (V0, V0, Vy),
                    (Vz, V0, V0), (V0, Vz, V0), (V0, V0, Vz),
                    (Vxx, V0, V0), (V0, Vxx, V0), (V0, V0, Vxx),
                    (Vxy, V0, V0), (V0, Vxy, V0), (V0, V0, Vxy),
                    (Vxz, V0, V0), (V0, Vxz, V0), (V0, V0, Vxz),
                    (Vyy, V0, V0), (V0, Vyy, V0), (V0, V0, Vyy),
                    (Vyz, V0, V0), (V0, Vyz, V0), (V0, V0, Vyz),
                    (Vzz, V0, V0), (V0, Vzz, V0), (V0, V0, Vzz),
                    (Vxxx, V0, V0), (V0, Vxxx, V0), (V0, V0, Vxxx),
                    (Vxxy, V0, V0), (V0, Vxxy, V0), (V0, V0, Vxxy),
                    (Vxxz, V0, V0), (V0, Vxxz, V0), (V0, V0, Vxxz),
                    (Vxyy, V0, V0), (V0, Vxyy, V0), (V0, V0, Vxyy),
                    (Vxyz, V0, V0), (V0, Vxyz, V0), (V0, V0, Vxyz),
                    (Vxzz, V0, V0), (V0, Vxzz, V0), (V0, V0, Vxzz),
                    (Vyyy, V0, V0), (V0, Vyyy, V0), (V0, V0, Vyyy),
                    (Vyyz, V0, V0), (V0, Vyyz, V0), (V0, V0, Vyyz),
                    (Vyzz, V0, V0), (V0, Vyzz, V0), (V0, V0, Vyzz),
                    (Vzzz, V0, V0), (V0, Vzzz, V0), (V0, V0, Vzzz),
                    (Vxxxx, Vxxxy, Vxxxz), (Vxxxy, Vxxyy, Vxxyz),
                    (Vxxxz, Vxxyz, Vxxzz), (Vxxyy, Vxyyy, Vxyyz),
                    (Vxxyz, Vxyyz, Vxyzz), (Vxxzz, Vxyzz, Vxzzz),
                    (Vxyyy, Vyyyy, Vyyyz), (Vxyyz, Vyyyz, Vyyzz),
                    (Vxyzz, Vyyzz, Vyzzz), (Vxzzz, Vyzzz, Vzzzz)]
        test_polynomials = [sy.S.One, t, s, t**2, t*s,
                            s**2, t**3, t**2*s, t*s**2, s**3]*4
        test_polynomials += [qi100, qi010, qi001, qix00, qi0x0, qi00x,
                             qiy00, qi0y0, qi00y, qiz00, qi0z0, qi00z,
                             qixx00, qi0xx0, qi00xx, qixy00, qi0xy0, qi00xy,
                             qixz00, qi0xz0, qi00xz, qiyy00, qi0yy0, qi00yy,
                             qiyz00, qi0yz0, qi00yz, qizz00, qi0zz0, qi00zz]
    elif(k == -1):  # 12 dofs
        dof_list = [NT1]*3 + [NT2]*3 + [NT3]*3 + [NT4]*3
        mon_list = [(V1, V0, V0), (V0, V1, V0), (V0, V0, V1),
                    (Vx, V0, V0), (V0, Vx, V0), (V0, V0, Vx),
                    (Vy, V0, V0), (V0, Vy, V0), (V0, V0, Vy),
                    (Vz, V0, V0), (V0, Vz, V0), (V0, V0, Vz)]
        test_polynomials = [sy.S.One, t, s]*4
    elif(k == -2):  # 30 dofs
        dof_list = [NT1]*6 + [NT2]*6 + [NT3]*6 + [NT4]*6 + [NTI]*6
        mon_list = [(V1, V0, V0), (V0, V1, V0), (V0, V0, V1),
                    (Vx, V0, V0), (V0, Vx, V0), (V0, V0, Vx),
                    (Vy, V0, V0), (V0, Vy, V0), (V0, V0, Vy),
                    (Vz, V0, V0), (V0, Vz, V0), (V0, V0, Vz),
                    (Vxx, V0, V0), (V0, Vxx, V0), (V0, V0, Vxx),
                    (Vxy, V0, V0), (V0, Vxy, V0), (V0, V0, Vxy),
                    (Vxz, V0, V0), (V0, Vxz, V0), (V0, V0, Vxz),
                    (Vyy, V0, V0), (V0, Vyy, V0), (V0, V0, Vyy),
                    (Vyz, V0, V0), (V0, Vyz, V0), (V0, V0, Vyz),
                    (Vzz, V0, V0), (V0, Vzz, V0), (V0, V0, Vzz)]
        test_polynomials = [sy.S.One, t, s, t**2, t*s, s**2]*4
        test_polynomials += [qi100, qi010, qi001,
                             tuple(map(lambda x, y: -x + y, qi0z0, qi00y)),
                             tuple(map(lambda x, y: x - y, qiz00, qi00x)),
                             tuple(map(lambda x, y: -x + y, qiy00, qi0x0))]
    elif(k == -3):  # 60 dofs
        dof_list = [NT1]*10 + [NT2]*10 + [NT3]*10 + [NT4]*10 + [NTI]*20
        mon_list = [(V1, V0, V0), (V0, V1, V0), (V0, V0, V1),
                    (Vx, V0, V0), (V0, Vx, V0), (V0, V0, Vx),
                    (Vy, V0, V0), (V0, Vy, V0), (V0, V0, Vy),
                    (Vz, V0, V0), (V0, Vz, V0), (V0, V0, Vz),
                    (Vxx, V0, V0), (V0, Vxx, V0), (V0, V0, Vxx),
                    (Vxy, V0, V0), (V0, Vxy, V0), (V0, V0, Vxy),
                    (Vxz, V0, V0), (V0, Vxz, V0), (V0, V0, Vxz),
                    (Vyy, V0, V0), (V0, Vyy, V0), (V0, V0, Vyy),
                    (Vyz, V0, V0), (V0, Vyz, V0), (V0, V0, Vyz),
                    (Vzz, V0, V0), (V0, Vzz, V0), (V0, V0, Vzz),
                    (Vxxx, V0, V0), (V0, Vxxx, V0), (V0, V0, Vxxx),
                    (Vxxy, V0, V0), (V0, Vxxy, V0), (V0, V0, Vxxy),
                    (Vxxz, V0, V0), (V0, Vxxz, V0), (V0, V0, Vxxz),
                    (Vxyy, V0, V0), (V0, Vxyy, V0), (V0, V0, Vxyy),
                    (Vxyz, V0, V0), (V0, Vxyz, V0), (V0, V0, Vxyz),
                    (Vxzz, V0, V0), (V0, Vxzz, V0), (V0, V0, Vxzz),
                    (Vyyy, V0, V0), (V0, Vyyy, V0), (V0, V0, Vyyy),
                    (Vyyz, V0, V0), (V0, Vyyz, V0), (V0, V0, Vyyz),
                    (Vyzz, V0, V0), (V0, Vyzz, V0), (V0, V0, Vyzz),
                    (Vzzz, V0, V0), (V0, Vzzz, V0), (V0, V0, Vzzz)]
        test_polynomials = [sy.S.One, t, s, t**2, t*s,
                            s**2, t**3, t**2*s, t*s**2, s**3]*4
        test_polynomials += [qi100, qi010, qi001, qix00, qi0x0, qi00x,
                             qiy00, qi0y0, qi00y, qiz00, qi0z0, qi00z,
                             tuple(map(lambda x, y: -x + y, qi0xz0, qi00xy)),
                             tuple(map(lambda x, y: x - y, qixz00, qi00xx)),
                             tuple(map(lambda x, y: -x + y, qixy00, qi0xx0)),
                             tuple(map(lambda x, y: -x + y, qi0yz0, qi00yy)),
                             tuple(map(lambda x, y: x - y, qiyz00, qi00xy)),
                             tuple(map(lambda x, y: -x + y, qiyy00, qi0xy0)),
                             tuple(map(lambda x, y: -x + y, qi0zz0, qi00yz)),
                             tuple(map(lambda x, y: x - y, qizz00, qi00xz))]

    a, ainv = compute_matrix_and_inverse(dof_list, mon_list, test_polynomials)
    return a, ainv, mon_list


def print_parmoon_file(ainv, mon_list, identifyer, file_name, tetra, comment,
                       factor=1.):
    n_dof = len(mon_list)
    assert n_dof == np.shape(ainv)[0]
    with open(file_name, "w") as f:
        f.write(f"// {comment}\n")
        # write the number of degrees of freedom
        n_dof_variable = f"{identifyer}_n_dof"
        coeff_matrix_name = f"{identifyer}_CM"
        f.write(f"constexpr int {n_dof_variable} = {n_dof};\n")
        # write the matrix of coefficients
        f.write("constexpr double {0}[{1}*{1}] = {{"
                .format(coeff_matrix_name, n_dof_variable))
        f.write("\n")
        for i in range(n_dof):
            for j in range(n_dof):
                if factor == 1.0:
                    f.write("{}, ".format(ainv[i, j]))
                else:
                    f.write("{}, ".format(np.round(ainv[i, j], 8)))
            f.write("\n")
        f.write("};\n")

        def function_writer(name, mon_list):
            f.write("static void {}(double xi, double eta, double zeta,\n"
                    .format(name))
            f.write("{}double * values)\n".format(' '*(13+len(name))))
            f.write("{\n")
            f.write("  // monomials, x-, y- and z-component\n")
            xi, eta, zeta = sy.symbols('xi eta zeta', real=True)
            subs = {x: xi, y: eta, z: zeta}
            x_mon = [sy.printing.cxxcode(i[0].subs(subs)) for i in mon_list]
            y_mon = [sy.printing.cxxcode(i[1].subs(subs)) for i in mon_list]
            z_mon = [sy.printing.cxxcode(i[2].subs(subs)) for i in mon_list]
            f.write("  double mon_x[{}] = {{\n  ".format(n_dof_variable))
            f.write("  {} }};\n".format(', '.join(x_mon)))
            f.write("  double mon_y[{}] = {{\n  ".format(n_dof_variable))
            f.write("  {} }};\n".format(', '.join(y_mon)))
            f.write("  double mon_z[{}] = {{\n  ".format(n_dof_variable))
            f.write("  {} }};\n".format(', '.join(z_mon)))
            f.write("  memset(values, 0.0, 3*{}*SizeOfDouble); // "
                    "3 is the space dimension\n".format(n_dof_variable))
            f.write(f"  for(int i = 0; i < {n_dof_variable}; ++i)\n  {{\n")
            f.write(f"    for(int j = 0; j < {n_dof_variable}; ++j)\n    {{\n")
            f.write(f"      double v = {coeff_matrix_name}[i+j*{n_dof_variable}];\n")
            f.write("      values[i     {}] += v*mon_x[j];\n"
                    .format(' '*len(n_dof_variable)))
            f.write(f"      values[i +   {n_dof_variable}] += v*mon_y[j];\n")
            f.write(f"      values[i + 2*{n_dof_variable}] += v*mon_z[j];\n")
            f.write("    }\n  }\n")
            if factor != 1.0:
                f.write(f"  for(int i = 0; i < 3*{n_dof_variable}; ++i)\n")
                f.write("  {\n")
                f.write(f"    values[i] /= {factor};\n")
                f.write("  }\n")
            f.write("}\n")

        # write the method for the evaluation of the basis functions
        func_name = "{}_Func".format(identifyer)
        x_name = "{}_DeriveXi".format(identifyer)
        y_name = "{}_DeriveEta".format(identifyer)
        z_name = "{}_DeriveZeta".format(identifyer)
        xx_name = "{}_DeriveXiXi".format(identifyer)
        xy_name = "{}_DeriveXiEta".format(identifyer)
        xz_name = "{}_DeriveXiZeta".format(identifyer)
        yy_name = "{}_DeriveEtaEta".format(identifyer)
        yz_name = "{}_DeriveEtaZeta".format(identifyer)
        zz_name = "{}_DeriveZetaZeta".format(identifyer)
        function_writer(func_name, mon_list)
        function_writer(x_name, [(sy.diff(i[0], x), sy.diff(i[1], x),
                                  sy.diff(i[2], x)) for i in mon_list])
        function_writer(y_name, [(sy.diff(i[0], y), sy.diff(i[1], y),
                                  sy.diff(i[2], y)) for i in mon_list])
        function_writer(z_name, [(sy.diff(i[0], z), sy.diff(i[1], z),
                                  sy.diff(i[2], z)) for i in mon_list])
        function_writer(xx_name, [(sy.diff(i[0], x, x), sy.diff(i[1], x, x),
                                   sy.diff(i[2], x, x)) for i in mon_list])
        function_writer(xy_name, [(sy.diff(i[0], x, y), sy.diff(i[1], x, y),
                                   sy.diff(i[2], x, y)) for i in mon_list])
        function_writer(xz_name, [(sy.diff(i[0], x, z), sy.diff(i[1], x, z),
                                   sy.diff(i[2], x, z)) for i in mon_list])
        function_writer(yy_name, [(sy.diff(i[0], y, y), sy.diff(i[1], y, y),
                                   sy.diff(i[2], y, y)) for i in mon_list])
        function_writer(yz_name, [(sy.diff(i[0], y, z), sy.diff(i[1], y, z),
                                   sy.diff(i[2], y, z)) for i in mon_list])
        function_writer(zz_name, [(sy.diff(i[0], z, z), sy.diff(i[1], z, z),
                                   sy.diff(i[2], z, z)) for i in mon_list])
        f.write("TBaseFunct3D * BF_{0}_Obj = \n".format(identifyer))
        reference_ele = "BFUnitTetrahedron" if tetra else "BFUnitHexahedron"
        f.write("new TBaseFunct3D({1}, BF_{0}, {2},\n"
                .format(identifyer, n_dof_variable, reference_ele))
        f.write("{}{}, {},\n".format(' '*17, func_name, x_name))
        f.write("{}{}, {},\n".format(' '*17, y_name, z_name))
        f.write("{}{}, {},\n".format(' '*17, xx_name, xy_name))
        f.write("{}{}, {},\n".format(' '*17, xz_name, yy_name))
        f.write("{}{}, {},\n".format(' '*17, yz_name, zz_name))
        # here we only check the first component to get the polynomial degree
        polynomial_degree = max([sy.degree(p[0].subs({y: x, z: x}), x)
                                 for p in mon_list])
        accuracy = 1
        f.write("{}{}, {}, 0, nullptr, 3);\n".format(' '*17, polynomial_degree,
                                                     accuracy))


if __name__ == '__main__':
    # np.set_printoptions(linewidth=180, precision=10)
    # print("Raviart-Thomas 0 on hexahedra")
    # (a, ainv, mon_list) = build_matrix_Hdiv_H(0)
    # print_parmoon_file(ainv, mon_list, "N_H_RT0_3D", "BF_N_H_RT0_3D.h",
    #                    False,
    #                    "Raviart-Thomas element of order zero on hexahedra")

    print("Raviart-Thomas 1 on hexahedra")
    (a, ainv, mon_list) = build_matrix_Hdiv_H(1)
    print_parmoon_file(64*ainv, mon_list, "N_H_RT1_3D", "BF_N_H_RT1_3D.h",
                       False,
                       "Raviart-Thomas element of order one on hexahedra",
                       factor=64.)

    # np.set_printoptions(linewidth=400, threshold=2000, suppress=True,
    #                     precision=10)
    # print("Raviart-Thomas 2 on hexahedra")
    # (a, ainv, mon_list) = build_matrix_Hdiv_H(2)
    # print_parmoon_file(102.4*ainv, mon_list, "N_H_RT2_3D", "BF_N_H_RT2_3D.h",
    #                    False,
    #                    "Raviart-Thomas element of order two on hexahedra",
    #                    factor=102.4)

    # print("Brezzi-Douglas-Duran-Fortin 1 on hexahedra")
    # (a, ainv, mon_list) = build_matrix_Hdiv_H(-1)
    # print(16*ainv)

    # print("Brezzi-Douglas-Duran-Fortin 2 on hexahedra")
    # (a, ainv, mon_list) = build_matrix_Hdiv_H(-2)
    # print(32*ainv)

    # print("Brezzi-Douglas-Duran-Fortin 3 on hexahedra")
    # (a, ainv, mon_list) = build_matrix_Hdiv_H(-3)
    # print(32*ainv)

    # Triangles:

    # print("Raviart-Thomas 0 on tetrahedra")
    # (a, ainv, mon_list) = build_matrix_Hdiv_T(0)
    # print(ainv)

    # print("Raviart-Thomas 1 on tetrahedra")
    # (a, ainv, mon_list) = build_matrix_Hdiv_T(1)
    # print(ainv)

    # print("Raviart-Thomas 2 on tetrahedra")
    # (a, ainv, mon_list) = build_matrix_Hdiv_T(2)
    # print(ainv)

    # print("Raviart-Thomas 3 on tetrahedra")
    # (a, ainv, mon_list) = build_matrix_Hdiv_T(3)
    # np.set_printoptions(linewidth=220, threshold=5000, suppress=True,
    #                     precision=8)
    # ainv = np.round(3*ainv, decimals=8)/3.
    # print(3*ainv)
    # print("inserse property")
    # print(np.sum(a*ainv-np.eye(70)))

    # print("Brezzi-Douglas-Marini 1 on tetrahedra")
    # (a, ainv, mon_list) = build_matrix_Hdiv_T(-1)
    # print(ainv)

    # print("Brezzi-Douglas-Marini 2 on tetrahedra")
    # (a, ainv, mon_list) = build_matrix_Hdiv_T(-2)
    # print(ainv)

    # print("Brezzi-Douglas-Marini 3 on tetrahedra")
    # (a, ainv, mon_list) = build_matrix_Hdiv_T(-3)
    # print(ainv)
