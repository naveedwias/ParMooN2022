import numpy as np
import sympy as sy

"""
This script calculates the coefficients of basis functions for a number of
H(div) conforming elements. In particular with the methods build_matrix_Hdiv_Q,
and build_matrix_Hdiv_T you can get the coefficients for Raviart-Thomas (RT)
and Brezzi-Douglas-Marini (BDM) elements on the reference square [-1,1]^2 and
the reference triangle whose vertices are (0,0), (1,0), and (0,1).

For each element the space consists of a number of basis functions (mainly
monomials) who are represented by expressions whose name starts with 'V'. The
nodal functionals (the degrees of freedom) are integrals over either the edges
of the reference cell or over the cell itself and are represented by methods
whose names start with NQ or NT. These integrals have an additional weighting
function which is specific to each degree of freedom. This weight function
for the boundary degrees of freedom are the Legendre polynomials and are
represented by the methods q(0), q(1), q(2), q(3). The weights for the inner
degrees of freedom have names starting with 'qi'.

Let N_i be the degrees of freedom (functionals) and v_j the monomials for one
of the considered elements. Then this script computes (in
compute_matrix_and_inverse) the matrix with entries N_i(v_j) and inverts this
matrix. Its inverse can then be copied into ParMooN, see e.g. BF_N_Q_RT2_2D.h.
Its columns are the coefficients (with respect to the above monomials) needed
to define the corresponding basis function. The order of the monomials is
important and should be the same in ParMooN as it is here! Similarly, the
degrees of freedom have to be defined in the same order in ParMooN (see e.g.
NF_N_Q_RT2_2D.h) as they are defined here.

You can call this script via
python3 BuildMatrix_Hdiv.py
It will print all the matrices which are used in ParMoon.
"""

x, y, t, s = sy.symbols('x y t s', real=True)

# %% Defining trivial basis functions (monomials)

V0 = sy.S.Zero
V1 = sy.S.One
# first order
Vx = x
Vy = y
# second order
Vxx = x**2
Vyy = y**2
Vxy = x*y
V2xy = 2*x*y
Vm2xy = -2*x*y
Vmyy = -y**2
# third order
Vxxx = x**3
Vxxy = x**2*y
Vxyy = x*y**2
Vyyy = y**3
V3xyy = 3*x*y**2
Vm3xxy = -3*x**2*y
Vmyyy = -y**3
# fourth order
Vxxxx = x*x*x*x
Vxxxy = x*x*x*y
Vxxyy = x*x*y*y
Vxyyy = x*y*y*y
Vyyyy = y*y*y*y
V4xyyy = 4*x*y*y*y
Vm4xxxy = -4*x*x*x*y
Vmyyyy = -y*y*y*y
# xy * third order
Vxxxxy = x*x*x*x*y
Vxxxyy = x*x*x*y*y
Vxxyyy = x*x*y*y*y
Vxyyyy = x*y*y*y*y
# xxyy * second order
Vxxxxyy = x*x*x*x*y*y
Vxxxyyy = x*x*x*y*y*y
Vxxyyyy = x*x*y*y*y*y
# xxxyyy * first order
Vxxxxyyy = x*x*x*x*y*y*y
Vxxxyyyy = x*x*x*y*y*y*y


# %% defining the test polynomials for degrees of freedom on edges
# Legendre polynomials.
def q(k):
    return sy.legendre(k, t)


# defining test polynomials for inner degrees of freedom
qi10 = (sy.S.One, sy.S.Zero)
qi01 = (sy.S.Zero, sy.S.One)

qix0 = (x, sy.S.Zero)
qiy0 = (y, sy.S.Zero)
qi0x = (sy.S.Zero, x)
qi0y = (sy.S.Zero, y)
qi2x0 = (2*x, sy.S.Zero)
qi02y = (sy.S.Zero, 2*y)
qiyx = (y, x)

qixx0 = (x**2, sy.S.Zero)
qixy0 = (x*y, sy.S.Zero)
qiyy0 = (y**2, sy.S.Zero)
qi0xx = (sy.S.Zero, x**2)
qi0xy = (sy.S.Zero, x*y)
qi0yy = (sy.S.Zero, y**2)

qixxx0 = (x**3, sy.S.Zero)
qixxy0 = (x**2*y, sy.S.Zero)
qixyy0 = (x*y**2, sy.S.Zero)
qiyyy0 = (y**3, sy.S.Zero)
qi0xxx = (sy.S.Zero, x**3)
qi0xxy = (sy.S.Zero, x**2*y)
qi0xyy = (sy.S.Zero, x*y**2)
qi0yyy = (sy.S.Zero, y**3)
qicb = (x-x**2-2*x*y, -y+2*x*y+y**2)  # curl of bubble function (x*y*(1-x-y))
qicbx = (x*x-x*x*x-2*x*x*y, -2*x*y+3*x*x*y+2*x*y*y)
qicby = (2*x*y-2*x*x*y-3*x*y*y, -y*y+2*x*y*y+y*y*y)

qixxxx0 = (x*x*x*x, sy.S.Zero)
qixxxy0 = (x*x*x*y, sy.S.Zero)
qixxyy0 = (x*x*y*y, sy.S.Zero)
qixyyy0 = (x*y*y*y, sy.S.Zero)
qiyyyy0 = (y*y*y*y, sy.S.Zero)
qi0xxxx = (sy.S.Zero, x*x*x*x)
qi0xxxy = (sy.S.Zero, x*x*x*y)
qi0xxyy = (sy.S.Zero, x*x*y*y)
qi0xyyy = (sy.S.Zero, x*y*y*y)
qi0yyyy = (sy.S.Zero, y*y*y*y)

qixxyyy0 = (x*x*y*y*y, sy.S.Zero)
qi0xxxyy = (sy.S.Zero, x*x*x*y*y)


# %% defining the nodal functionals on quads
def NQ1(v, q):
    integrand = -v[1].subs({y: -1}) * q.subs({t: x})
    return sy.integrate(integrand, (x, -1, 1))


def NQ2(v, q):
    integrand = v[0].subs({x: 1}) * q.subs({t: y})
    return sy.integrate(integrand, (y, -1, 1))


def NQ3(v, q):
    integrand = v[1].subs({x: -s, y: 1}) * q.subs({t: s})
    return sy.integrate(integrand, (s, -1, 1))


def NQ4(v, q):
    integrand = -v[0].subs({x: -1, y: -s}) * q.subs({t: s})
    return sy.integrate(integrand, (s, -1, 1))


def NQI(v, q):
    integrand = v[0]*q[0] + v[1]*q[1]
    return sy.integrate(integrand, (x, -1, 1), (y, -1, 1))


# defining the nodal functionals on triangles
def NT1(v, q):
    integrand = -v[1].subs({x: s, y: 0})*q.subs({t: 2*s-1})
    return sy.integrate(integrand, (s, 0, 1))


def NT2(v, q):
    integrand = ((v[0] + v[1])*q).subs({x: 1-s, y: s, t: 2*s-1})
    return sy.integrate(integrand, (s, 0, 1))


def NT3(v, q):
    integrand = -v[0].subs({x: 0, y: 1-s})*q.subs({t: 2*s-1})
    return sy.integrate(integrand, (s, 0, 1))


def NTI(v, q):
    integrand = v[0]*q[0] + v[1]*q[1]
    return sy.integrate(sy.integrate(integrand, (x, 0, 1-y)), (y, 0, 1))


# %%
def compute_matrix_and_inverse(dim, dof_list, mon_list, test_polynomials):
    # Compute the matrix
    A = np.matrix(np.zeros(shape=(dim, dim)))
    for i in range(dim):
        # print("i", i)
        for j in range(dim):
            val = dof_list[i](mon_list[j], test_polynomials[i])
            # print("j", j, val)
            A[i, j] = val
    Ainv = np.linalg.inv(A)
    # set very small entries to exactly zero
    Ainv[abs(Ainv) < 1e-10] = 0
    return (A, Ainv)


###############################################################################
# quads

def build_matrix_Hdiv_Q(k):
    """
    compute Matrix of coefficients and its inverse to get the Raviart-Thomas
    (k>=0) or Brezzi-Douglas-Marini (k<0) basis functions on the square
    Q = [-1,1]^2.

    Input: k denotes the order (k>=0 for RT, k<0 for BDM)
    Output: Matrix A=(a_ij) with a_ij = N_i(v_j) and its inverse. Here N_i
            are the degrees of freedom (functionals) and v_j are monomial
            basis vectors. The Raviart-Thomas or Brezzi-Douglas-Marini basis
            functions phi_j are then defined as phi_j = sum_i c_ij v_i. The
            monomial basis vectors are stored in the local variable Mon_List.
    """

    if(k >= 4 or k <= -4):
        print("So far the order k can only be 0, 1, 2, or 3 for ",
              "Raviart-Thomas elements and -1, -2, or -3 for ",
              "Brezzi-Douglas-Marini elements")
        print("No matrix computed!")
        A = 0
        Ainv = 0
        return (A, Ainv)

    # the dimension of the local RT space
    dim = 2*(k+1)*(k+2)
    if(k < 0):  # BDM elements
        dim = 2 + (-k+1)*(-k+2)

    if(k == 0):  # 4 dofs
        dof_list = (NQ1, NQ2, NQ3, NQ4)
        mon_list = [(V1, V0), (V0, V1), (Vx, V0), (V0, Vy)]
        test_polynomials = (q(0), q(0), q(0), q(0))
    elif(k == 1):  # 12 dofs
        dof_list = (NQ1, NQ1, NQ2, NQ2, NQ3, NQ3, NQ4, NQ4, NQI, NQI, NQI, NQI)
        mon_list = ((V1, V0), (V0, V1), (Vx, V0), (V0, Vx), (Vy, V0), (V0, Vy),
                    (Vxy, V0), (V0, Vxy), (Vxx, V0), (V0, Vyy), (Vxxy, V0),
                    (V0, Vxyy))
        test_polynomials = (q(0), q(1), q(0), q(1), q(0), q(1), q(0), q(1),
                            qi10, qi01, qiy0, qi0x)
    elif(k == 2):  # 24 dofs
        dof_list = (NQ1, NQ1, NQ1, NQ2, NQ2, NQ2, NQ3, NQ3, NQ3, NQ4, NQ4, NQ4,
                    NQI, NQI, NQI, NQI, NQI, NQI, NQI, NQI, NQI, NQI, NQI, NQI)
        mon_list = ((V1, V0), (V0, V1), (Vx, V0), (V0, Vx), (Vy, V0), (V0, Vy),
                    (Vxx, V0), (V0, Vxx), (Vxy, V0), (V0, Vxy),
                    (Vyy, V0), (V0, Vyy),
                    (Vxxx, V0), (V0, Vyyy), (Vxxy, V0), (V0, Vxxy),
                    (Vxyy, V0), (V0, Vxyy), (Vxxxy, V0), (V0, Vxyyy),
                    (Vxxyy, V0), (V0, Vxxyy), (Vxxxyy, V0), (V0, Vxxyyy))
        test_polynomials = (q(0), q(1), q(2), q(0), q(1), q(2),
                            q(0), q(1), q(2), q(0), q(1), q(2),
                            qi10, qi01, qix0, qi0x, qiy0, qi0y,
                            qi0xx, qiyy0, qixy0, qi0xy, qi0xxy, qixyy0)
    elif(k == 3):  # 40 dofs
        dof_list = (NQ1, NQ1, NQ1, NQ1, NQ2, NQ2, NQ2, NQ2, NQ3, NQ3, NQ3, NQ3,
                    NQ4, NQ4, NQ4, NQ4,
                    NQI, NQI, NQI, NQI, NQI, NQI, NQI, NQI, NQI, NQI, NQI, NQI,
                    NQI, NQI, NQI, NQI, NQI, NQI, NQI, NQI, NQI, NQI, NQI, NQI)
        mon_list = ((V1, V0), (V0, V1), (Vx, V0), (V0, Vx), (Vy, V0), (V0, Vy),
                    (Vxx, V0), (V0, Vxx), (Vxy, V0), (V0, Vxy),
                    (Vyy, V0), (V0, Vyy),
                    (Vxxx, V0), (V0, Vxxx), (Vxxy, V0), (V0, Vxxy),
                    (Vxyy, V0), (V0, Vxyy), (Vyyy, V0), (V0, Vyyy),
                    (Vxxxx, V0), (V0, Vxxxy), (Vxxxy, V0), (V0, Vxxyy),
                    (Vxxyy, V0), (V0, Vxyyy), (Vxyyy, V0), (V0, Vyyyy),
                    (Vxxxxy, V0), (V0, Vxxxyy), (Vxxxyy, V0), (V0, Vxxyyy),
                    (Vxxyyy, V0), (V0, Vxyyyy), (Vxxxxyy, V0), (V0, Vxxxyyy),
                    (Vxxxyyy, V0), (V0, Vxxyyyy), (Vxxxxyyy, V0),
                    (V0, Vxxxyyyy))
        test_polynomials = (q(0), q(1), q(2), q(3), q(0), q(1), q(2), q(3),
                            q(0), q(1), q(2), q(3), q(0), q(1), q(2), q(3),
                            qi10, qi01, qix0, qi0x, qiy0, qi0y,
                            qixx0, qi0xx, qiyy0, qi0yy, qixy0, qi0xy,
                            qiyyy0, qi0xxx, qixxy0, qi0xxy, qixyy0, qi0xyy,
                            qixxyy0, qi0xxyy, qixyyy0, qi0xxxy,
                            qixxyyy0, qi0xxxyy)
    elif(k == -1):  # 8 dofs
        dof_list = (NQ1, NQ1, NQ2, NQ2, NQ3, NQ3, NQ4, NQ4)
        mon_list = ((V1, V0), (V0, V1), (Vx, V0), (V0, Vx), (Vy, V0), (V0, Vy),
                    (Vxx, Vm2xy), (V2xy, Vmyy))
        test_polynomials = (q(0), q(1), q(0), q(1), q(0), q(1), q(0), q(1))
    elif(k == -2):  # 14 dofs
        dof_list = (NQ1, NQ1, NQ1, NQ2, NQ2, NQ2, NQ3, NQ3, NQ3, NQ4, NQ4, NQ4,
                    NQI, NQI)
        mon_list = ((V1, V0), (V0, V1), (Vx, V0), (V0, Vx), (Vy, V0), (V0, Vy),
                    (Vxx, V0), (V0, Vxx), (Vxy, V0), (V0, Vxy),
                    (Vyy, V0), (V0, Vyy), (Vxxx, Vm3xxy), (V3xyy, Vmyyy))
        test_polynomials = (q(0), q(1), q(2), q(0), q(1), q(2),
                            q(0), q(1), q(2), q(0), q(1), q(2), qi10, qi01)
    elif(k == -3):  # 22 dofs
        dof_list = (NQ1, NQ1, NQ1, NQ1, NQ2, NQ2, NQ2, NQ2, NQ3, NQ3, NQ3, NQ3,
                    NQ4, NQ4, NQ4, NQ4, NQI, NQI, NQI, NQI, NQI, NQI)
        mon_list = ((V1, V0), (V0, V1), (Vx, V0), (V0, Vx), (Vy, V0), (V0, Vy),
                    (Vxx, V0), (V0, Vxx), (Vxy, V0), (V0, Vxy),
                    (Vyy, V0), (V0, Vyy),
                    (Vxxx, V0), (V0, Vxxx), (Vxxy, V0), (V0, Vxxy),
                    (Vxyy, V0), (V0, Vxyy), (Vyyy, V0), (V0, Vyyy),
                    (Vxxxx, Vm4xxxy), (V4xyyy, Vmyyyy))
        test_polynomials = (q(0), q(1), q(2), q(3), q(0), q(1), q(2), q(3),
                            q(0), q(1), q(2), q(3), q(0), q(1), q(2), q(3),
                            qi10, qi01, qix0, qi0x, qiy0, qi0y)

    return compute_matrix_and_inverse(dim, dof_list, mon_list,
                                      test_polynomials)

###############################################################################
# triangles


def build_matrix_Hdiv_T(k):
    """
    compute Matrix of coefficients and its inverse to get the Raviart-Thomas
    (k>=0) or Brezzi-Douglas-Marini (k<0) basis functions on the triangle with
    vertices (0,0), (1,0), (0,1).

    Input: k denotes the order (k>=0 for RT, k<0 for BDM)
    Output: Matrix A=(a_ij) with a_ij = N_i(v_j) and its inverse. Here N_i
            are the degrees of freedom (functionals) and v_j are monomial
            basis vectors. The Raviart-Thomas or Brezzi-Douglas-Marini basis
            functions phi_j are then defined as phi_j = sum_i c_ij v_i. The
            monomial basis vectors are stored in the local variable Mon_List.
    """

    if(k >= 4 or k <= -4):
        print("So far the order k can only be 0, 1, 2, or 3 for",
              "Raviart-Thomas elements and -1, -2, or -3 for",
              "Brezzi-Douglas-Marini elements")
        print("No matrix computed!")
        A = 0
        Ainv = 0
        return (A, Ainv)

    # the dimension of the local RT space
    dim = (k+1)*(k+3)
    if(k < 0):  # BDM elements
        dim = (-k+1)*(-k+2)

    if(k == 0):  # 3 dofs
        dof_list = (NT1, NT2, NT3)
        mon_list = ((V1, V0), (V0, V1), (Vx, Vy))
        test_polynomials = (q(0), q(0), q(0))
    elif(k == 1):  # 8 dofs
        dof_list = (NT1, NT1, NT2, NT2, NT3, NT3, NTI, NTI)
        mon_list = ((V1, V0), (V0, V1), (Vx, V0), (V0, Vx), (Vy, V0), (V0, Vy),
                    (Vxx, Vxy), (Vxy, Vyy))
        test_polynomials = (q(0), q(1), q(0), q(1), q(0), q(1), qi10, qi01)
    elif(k == 2):  # 15 dofs
        dof_list = (NT1, NT1, NT1, NT2, NT2, NT2, NT3, NT3, NT3,
                    NTI, NTI, NTI, NTI, NTI, NTI)
        mon_list = ((V1, V0), (V0, V1), (Vx, V0), (V0, Vx), (Vy, V0), (V0, Vy),
                    (Vxx, V0), (V0, Vxx), (Vxy, V0), (V0, Vxy),
                    (Vyy, V0), (V0, Vyy),
                    (Vxxx, Vxxy), (Vxxy, Vxyy), (Vxyy, Vyyy))
        test_polynomials = (q(0), q(1), q(2), q(0), q(1), q(2),
                            q(0), q(1), q(2),
                            qi10, qi01, qix0, qi0x, qiy0, qi0y)
    elif(k == 3):  # 24 dofs
        dof_list = (NT1, NT1, NT1, NT1, NT2, NT2, NT2, NT2, NT3, NT3, NT3, NT3,
                    NTI, NTI, NTI, NTI, NTI, NTI, NTI, NTI, NTI, NTI, NTI, NTI)
        mon_list = ((V1, V0), (V0, V1), (Vx, V0), (V0, Vx), (Vy, V0), (V0, Vy),
                    (Vxx, V0), (V0, Vxx), (Vxy, V0), (V0, Vxy),
                    (Vyy, V0), (V0, Vyy),
                    (Vxxx, V0), (V0, Vxxx), (Vxxy, V0), (V0, Vxxy),
                    (Vxyy, V0), (V0, Vxyy), (Vyyy, V0), (V0, Vyyy),
                    (Vxxxx, Vxxxy), (Vxxxy, Vxxyy), (Vxxyy, Vxyyy),
                    (Vxyyy, Vyyyy))
        test_polynomials = (q(0), q(1), q(2), q(3), q(0), q(1), q(2), q(3),
                            q(0), q(1), q(2), q(3),
                            qi10, qi01, qix0, qi0x, qiy0, qi0y,
                            qixx0, qi0xx, qiyy0, qi0yy, qixy0, qi0xy)
    elif(k == -1):  # 6 dofs
        dof_list = (NT1, NT1, NT2, NT2, NT3, NT3)
        mon_list = ((V1, V0), (V0, V1), (Vx, V0), (V0, Vx), (Vy, V0), (V0, Vy))
        test_polynomials = (q(0), q(1), q(0), q(1), q(0), q(1))
    elif(k == -2):  # 12 dofs
        dof_list = (NT1, NT1, NT1, NT2, NT2, NT2, NT3, NT3, NT3, NTI, NTI, NTI)
        mon_list = ((V1, V0), (V0, V1), (Vx, V0), (V0, Vx), (Vy, V0), (V0, Vy),
                    (Vxx, V0), (V0, Vxx), (Vxy, V0), (V0, Vxy),
                    (Vyy, V0), (V0, Vyy))
        test_polynomials = (q(0), q(1), q(2), q(0), q(1), q(2),
                            q(0), q(1), q(2), qi10, qi01, qicb)
    elif(k == -3):  # 20 dofs
        dof_list = (NT1, NT1, NT1, NT1, NT2, NT2, NT2, NT2, NT3, NT3, NT3, NT3,
                    NTI, NTI, NTI, NTI, NTI, NTI, NTI, NTI)
        mon_list = ((V1, V0), (V0, V1), (Vx, V0), (V0, Vx), (Vy, V0), (V0, Vy),
                    (Vxx, V0), (V0, Vxx), (Vxy, V0), (V0, Vxy),
                    (Vyy, V0), (V0, Vyy),
                    (Vxxx, V0), (V0, Vxxx), (Vxxy, V0), (V0, Vxxy),
                    (Vxyy, V0), (V0, Vxyy), (Vyyy, V0), (V0, Vyyy))
        test_polynomials = (q(0), q(1), q(2), q(3), q(0), q(1), q(2), q(3),
                            q(0), q(1), q(2), q(3),
                            qi10, qi01, qi2x0, qi02y, qiyx, qicb, qicbx, qicby)

    return compute_matrix_and_inverse(dim, dof_list, mon_list,
                                      test_polynomials)


if __name__ == '__main__':
    np.set_printoptions(linewidth=180, precision=10)
    print("Raviart-Thomas 0 on quads")
    (a, ainv) = build_matrix_Hdiv_Q(0)
    # print(a)
    print(4*ainv)

    print("Raviart-Thomas 1 on quads")
    (a, ainv) = build_matrix_Hdiv_Q(1)
    # print(a)
    print(8*ainv)

    print("Raviart-Thomas 2 on quads")
    (a, ainv) = build_matrix_Hdiv_Q(2)
    # print(45*a)
    print(32*ainv)

    print("Raviart-Thomas 3 on quads")
    (a, ainv) = build_matrix_Hdiv_Q(3)
    np.set_printoptions(linewidth=220, threshold=2000, suppress=True,
                        precision=8)
    # print(3*a)
    print(768*ainv)
    print(a*ainv)

    print("Brezzi-Douglas-Marini 1 on quads")
    (a, ainv) = build_matrix_Hdiv_Q(-1)
    # print(3*a)
    print(16*ainv)

    print("Brezzi-Douglas-Marini 2 on quads")
    (a, ainv) = build_matrix_Hdiv_Q(-2)
    # print(15*a)
    print(8*ainv)

    print("Brezzi-Douglas-Marini 3 on quads")
    (a, ainv) = build_matrix_Hdiv_Q(-3)
    # print(45*a)
    print(32*ainv)

    # Triangles:

    print("Raviart-Thomas 0 on triangles")
    (a, ainv) = build_matrix_Hdiv_T(0)
    # print(a)
    print(ainv)

    print("Raviart-Thomas 1 on triangles")
    (a, ainv) = build_matrix_Hdiv_T(1)
    # print(a)
    print(ainv)

    print("Raviart-Thomas 2 on triangles")
    (a, ainv) = build_matrix_Hdiv_T(2)
    np.set_printoptions(linewidth=220, threshold=2000, suppress=True,
                        precision=8)
    # print(a)
    print(2*ainv)

    print("Raviart-Thomas 3 on triangles")
    (a, ainv) = build_matrix_Hdiv_T(3)
    np.set_printoptions(linewidth=220, threshold=2000, suppress=True,
                        precision=8)
    # print(a)
    print(ainv)

    print("Brezzi-Douglas-Marini 1 on triangles")
    (a, ainv) = build_matrix_Hdiv_T(-1)
    # print(a)
    print(ainv)

    print("Brezzi-Douglas-Marini 2 on triangles")
    (a, ainv) = build_matrix_Hdiv_T(-2)
    # print(a)
    print(ainv)

    print("Brezzi-Douglas-Marini 3 on triangles")
    (a, ainv) = build_matrix_Hdiv_T(-3)
    # print(a)
    print(ainv)
