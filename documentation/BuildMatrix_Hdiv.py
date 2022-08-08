import numpy as np
from scipy.integrate import quad, nquad # numerical quadrature

"""
This script calculates the coefficients of basis functions for a number of 
H(div) conforming elements. In particular with the methods build_matrix_Hdiv_Q,
and build_matrix_Hdiv_P you can get the coefficients for Raviart-Thomas (RT) and
Brezzi-Douglas-Marini (BDM) elements on the reference square [-1,1]^2 and the 
reference triangle whose vertices are (0,0), (1,0), and (0,1).

For each element the space consists of a number of basis functions (mainly
monomials) who are represented by methods whose name starts with 'V'. The nodal 
functionals (the degrees of freedom) are integrals over either the edges of the 
reference cell or over the cell itself and are represented by methods whose 
names start with NQ or NT. These integrals have an additional weighting 
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


# Defining trivial basis functions (monomials)
def V0(x,y):
  return np.zeros(np.size(x))
def V1(x,y):
  return np.ones(np.size(x))
# first order
def Vx(x,y):
  return x
def Vy(x,y):
  return y
# second order
def Vxx(x,y):
  return x*x
def Vyy(x,y):
  return y*y
def Vxy(x,y):
  return x*y
def V2xy(x,y):
  return 2*x*y
def Vm2xy(x,y):
  return -2*x*y
def Vmyy(x,y):
  return -y*y
# third order
def Vxxx(x,y):
  return x*x*x
def Vxxy(x,y):
  return x*x*y
def Vxyy(x,y):
  return x*y*y
def Vyyy(x,y):
  return y*y*y
def V3xyy(x,y):
  return 3*x*y*y
def Vm3xxy(x,y):
  return -3*x*x*y
def Vmyyy(x,y):
  return -y*y*y
# fourth order
def Vxxxx(x,y):
  return x*x*x*x
def Vxxxy(x,y):
  return x*x*x*y
def Vxxyy(x,y):
  return x*x*y*y
def Vxyyy(x,y):
  return x*y*y*y
def Vyyyy(x,y):
  return y*y*y*y
def V4xyyy(x,y):
  return 4*x*y*y*y
def Vm4xxxy(x,y):
  return -4*x*x*x*y
def Vmyyyy(x,y):
  return -y*y*y*y
# xy * third order
def Vxxxxy(x,y):
  return x*x*x*x*y
def Vxxxyy(x,y):
  return x*x*x*y*y
def Vxxyyy(x,y):
  return x*x*y*y*y
def Vxyyyy(x,y):
  return x*y*y*y*y
# xxyy * second order
def Vxxxxyy(x,y):
  return x*x*x*x*y*y
def Vxxxyyy(x,y):
  return x*x*x*y*y*y
def Vxxyyyy(x,y):
  return x*x*y*y*y*y
# xxxyyy * first order
def Vxxxxyyy(x,y):
  return x*x*x*x*y*y*y
def Vxxxyyyy(x,y):
  return x*x*x*y*y*y*y

#def V(x,y,ex,ey)
#  return x**ex * y**ey



# defining the test polynomials for degrees of freedom on edges
# Legendre polynomials. t is an array with entries in [-1,1]
def q1(t):
  return np.polynomial.legendre.legval(t, [1])
def q2(t):
  return np.polynomial.legendre.legval(t, [0, 1])
def q3(t):
  return np.polynomial.legendre.legval(t, [0, 0, 1])
def q4(t):
  return np.polynomial.legendre.legval(t, [0, 0, 0, 1])

# defining test polynomials for inner degrees of freedom
def qi10(x,y):
  return np.array([np.ones(np.size(x)), np.zeros(np.size(x))])
def qi01(x,y):
  return np.array([np.zeros(np.size(x)), np.ones(np.size(x))])

def qix0(x,y):
  return np.array([x, np.zeros(np.size(x))])
def qiy0(x,y):
  return np.array([y, np.zeros(np.size(x))])
def qi0x(x,y):
  return np.array([np.zeros(np.size(x)), x])
def qi0y(x,y):
  return np.array([np.zeros(np.size(x)), y])
def qi2x0(x,y):
  return 2*qix0(x,y)
def qi02y(x,y):
  return 2*qi0y(x,y)
def qiyx(x,y):
  return np.array([y, x])

def qixx0(x,y):
  return np.array([x*x, np.zeros(np.size(x))])
def qixy0(x,y):
  return np.array([x*y, np.zeros(np.size(x))])
def qiyy0(x,y):
  return np.array([y*y, np.zeros(np.size(x))])
def qi0xx(x,y):
  return np.array([np.zeros(np.size(x)), x*x])
def qi0xy(x,y):
  return np.array([np.zeros(np.size(x)), x*y])
def qi0yy(x,y):
  return np.array([np.zeros(np.size(x)), y*y])

def qixxx0(x,y):
  return np.array([x*x*x, np.zeros(np.size(x))])
def qixxy0(x,y):
  return np.array([x*x*y, np.zeros(np.size(x))])
def qixyy0(x,y):
  return np.array([x*y*y, np.zeros(np.size(x))])
def qiyyy0(x,y):
  return np.array([y*y*y, np.zeros(np.size(x))])
def qi0xxx(x,y):
  return np.array([np.zeros(np.size(x)), x*x*x])
def qi0xxy(x,y):
  return np.array([np.zeros(np.size(x)), x*x*y])
def qi0xyy(x,y):
  return np.array([np.zeros(np.size(x)), x*y*y])
def qi0yyy(x,y):
  return np.array([np.zeros(np.size(x)), y*y*y])
def qicb(x,y): # curl of bubble function (x*y*(1-x-y))
  return np.array([x-x*x-2*x*y, -y+2*x*y+y*y])
def qicbx(x,y):
  return np.array([x*x-x*x*x-2*x*x*y, -2*x*y+3*x*x*y+2*x*y*y])
def qicby(x,y):
  return np.array([2*x*y-2*x*x*y-3*x*y*y, -y*y+2*x*y*y+y*y*y])

def qixxxx0(x,y):
  return np.array([x*x*x*x, np.zeros(np.size(x))])
def qixxxy0(x,y):
  return np.array([x*x*x*y, np.zeros(np.size(x))])
def qixxyy0(x,y):
  return np.array([x*x*y*y, np.zeros(np.size(x))])
def qixyyy0(x,y):
  return np.array([x*y*y*y, np.zeros(np.size(x))])
def qiyyyy0(x,y):
  return np.array([y*y*y*y, np.zeros(np.size(x))])
def qi0xxxx(x,y):
  return np.array([np.zeros(np.size(x)), x*x*x*x])
def qi0xxxy(x,y):
  return np.array([np.zeros(np.size(x)), x*x*x*y])
def qi0xxyy(x,y):
  return np.array([np.zeros(np.size(x)), x*x*y*y])
def qi0xyyy(x,y):
  return np.array([np.zeros(np.size(x)), x*y*y*y])
def qi0yyyy(x,y):
  return np.array([np.zeros(np.size(x)), y*y*y*y])

def qixxyyy0(x,y):
  return np.array([x*x*y*y*y, np.zeros(np.size(x))])
def qi0xxxyy(x,y):
  return np.array([np.zeros(np.size(x)), x*x*x*y*y])


# defining the nodal functionals on quads
def NQ1(v1, v2, q):
  result = quad(lambda x :  -v2(x,-1)*q(x), -1,1)
  if(abs(result[0]) < 1e-12):
    return 0.
  return result[0]

def NQ2(v1, v2, q):
  result = quad(lambda x :  v1(1,x)*q(x), -1,1)
  if(abs(result[0]) < 1e-12):
    return 0.
  return result[0]

def NQ3(v1, v2, q):
  result = quad(lambda x :  v2(-x,1)*q(x), -1,1)
  if(abs(result[0]) < 1e-12):
    return 0.
  return result[0]

def NQ4(v1, v2, q):
  result = quad(lambda x : -v1(-1,-x)*q(x), -1,1)
  if(abs(result[0]) < 1e-12):
    return 0.
  return result[0]

def NQI(v1, v2, q):
  def integrand(x,y):
    return np.dot(np.array([v1(x,y), v2(x,y)]).reshape(-1), q(x,y).reshape(-1))
  result = nquad(integrand, [[-1,1],[-1,1]])
  if(abs(result[0]) < 1e-12):
    return 0.
  return result[0]

# defining the nodal functionals on triangles
def NT1(v1, v2, q):
  result = quad(lambda x :  -v2(x,0)*q(2*x-1), 0,1)
  if(abs(result[0]) < 1e-12):
    return 0.
  return result[0]

def NT2(v1, v2, q):
  result = quad(lambda x :  (v1(1-x,x) + v2(1-x,x))*q(2*x-1), 0,1)
  if(abs(result[0]) < 1e-12):
    return 0.
  return result[0]

def NT3(v1, v2, q):
  result = quad(lambda x :  -v1(0,1-x)*q(2*x-1), 0,1)
  if(abs(result[0]) < 1e-12):
    return 0.
  return result[0]

def NTI(v1, v2, q):
  def bounds_y():
    return [0, 1]
  def bounds_x(y):
    return [0, 1-y]
  def integrand(x,y):
    return np.dot(np.array([v1(x,y), v2(x,y)]).reshape(-1), q(x,y).reshape(-1))
  result = nquad(integrand, [bounds_x, bounds_y])
  if(abs(result[0]) < 1e-12):
    return 0.
  return result[0]



def compute_matrix_and_inverse(dim, dof_list, mon_list, test_polynomials):
  # Compute the matrix
  A = np.matrix(np.zeros(shape=(dim,dim)))
  for i in range(dim):
    for j in range(dim):
        A[i,j] = dof_list[i](mon_list[j], mon_list[j+dim], test_polynomials[i])
  Ainv = np.linalg.inv(A)
  # set very small entries to exactly zero
  Ainv[abs(Ainv)<1e-10] = 0
  return (A, Ainv)


###############################################################################
# quads

def build_matrix_Hdiv_Q(k):
  """ compute Matrix of coefficients and its inverse to get the Raviart-Thomas
      (k>=0) or Brezzi-Douglas-Marini (k<0) basis functions on the square
      Q = [-1,1]^2. 

      Input: k denotes the order (k>=0 for RT, k<0 for BDM)
      Output: Matrix A=(a_ij) with a_ij = N_i(v_j) and its inverse. Here N_i
              are the degrees of freedom (functionals) and v_j are monomial
              basis vectors. The Raviart-Thomas or Brezzi-Douglas-Marini basis
              functions phi_j are then defined as phi_j = sum_i c_ij v_i. The
              monomial basis vectors are stored in the local variable Mon_List.
  """

  if(k >= 4 or k <=-4):
    print("So far the order k can only be 0, 1, 2, or 3 for Raviart-Thomas ",
          "elements and -1, -2, or -3 for Brezzi-Douglas-Marini elements")
    print("No matrix computed!")
    A = 0
    Ainv = 0
    return (A, Ainv)

  # the dimension of the local RT space
  dim = 2*(k+1)*(k+2)
  if(k < 0): # BDM elements
    dim = 2 + (-k+1)*(-k+2)

  if(k == 0): # 4 dofs
    dof_list = (NQ1, NQ2, NQ3, NQ4)
    mon_list = (V1, V0, Vx, V0,
                V0, V1, V0, Vy)
    test_polynomials = (q1, q1, q1, q1)
  elif(k == 1): # 12 dofs
    dof_list = (NQ1, NQ1, NQ2, NQ2, NQ3, NQ3, NQ4, NQ4, NQI, NQI, NQI, NQI)
    mon_list = (V1, V0, Vx, V0, Vy, V0, Vxy, V0, Vxx, V0, Vxxy, V0,
                V0, V1, V0, Vx, V0, Vy, V0, Vxy, V0, Vyy, V0, Vxyy)
    test_polynomials = (q1, q2, q1, q2, q1, q2, q1, q2, qi10, qi01, qiy0, qi0x)
  elif(k == 2): # 24 dofs
    dof_list = (NQ1, NQ1, NQ1, NQ2, NQ2, NQ2, NQ3, NQ3, NQ3, NQ4, NQ4, NQ4,
                NQI, NQI, NQI, NQI, NQI, NQI, NQI, NQI, NQI, NQI, NQI, NQI)
    mon_list = (V1,V0, Vx,V0, Vy,V0, Vxx,V0, Vxy,V0, Vyy,V0, Vxxx,V0, Vxxy,V0, Vxyy,V0, Vxxxy,V0, Vxxyy,V0, Vxxxyy,V0,
                V0,V1, V0,Vx, V0,Vy, V0,Vxx, V0,Vxy, V0,Vyy, V0,Vyyy, V0,Vxxy, V0,Vxyy, V0,Vxyyy, V0,Vxxyy, V0,Vxxyyy)
    test_polynomials = (q1, q2, q3, q1, q2, q3, q1, q2, q3, q1, q2, q3,
                        qi10, qi01, qix0, qi0x, qiy0, qi0y, qi0xx, qiyy0, qixy0, qi0xy, qi0xxy, qixyy0)
  elif(k == 3): # 40 dofs
    dof_list = (NQ1, NQ1, NQ1, NQ1, NQ2, NQ2, NQ2, NQ2, NQ3, NQ3, NQ3, NQ3, NQ4, NQ4, NQ4, NQ4,
                NQI, NQI, NQI, NQI, NQI, NQI, NQI, NQI, NQI, NQI, NQI, NQI,
                NQI, NQI, NQI, NQI, NQI, NQI, NQI, NQI, NQI, NQI, NQI, NQI)
    mon_list = (V1,V0, Vx,V0, Vy,V0, Vxx,V0, Vxy,V0, Vyy,V0, Vxxx,V0, Vxxy,V0, Vxyy,V0, Vyyy,V0, Vxxxx,V0, Vxxxy,V0, Vxxyy,V0, Vxyyy,V0, Vxxxxy,V0, Vxxxyy,V0, Vxxyyy,V0, Vxxxxyy,V0, Vxxxyyy,V0, Vxxxxyyy,V0,
                V0,V1, V0,Vx, V0,Vy, V0,Vxx, V0,Vxy, V0,Vyy, V0,Vxxx, V0,Vxxy, V0,Vxyy, V0,Vyyy, V0,Vxxxy, V0,Vxxyy, V0,Vxyyy, V0,Vyyyy, V0,Vxxxyy, V0,Vxxyyy, V0,Vxyyyy, V0,Vxxxyyy, V0,Vxxyyyy, V0,Vxxxyyyy)
    test_polynomials = (q1, q2, q3, q4, q1, q2, q3, q4, q1, q2, q3, q4, q1, q2, q3, q4,
                        qi10, qi01, qix0, qi0x, qiy0, qi0y, qixx0, qi0xx, qiyy0, qi0yy, qixy0, qi0xy,
                        qiyyy0, qi0xxx, qixxy0, qi0xxy, qixyy0, qi0xyy,
                        qixxyy0, qi0xxyy, qixyyy0, qi0xxxy, qixxyyy0, qi0xxxyy)
  elif(k == -1): # 8 dofs
    dof_list = (NQ1, NQ1, NQ2, NQ2, NQ3, NQ3, NQ4, NQ4)
    mon_list = (V1, V0, Vx, V0, Vy, V0, Vxx,   V2xy,
                V0, V1, V0, Vx, V0, Vy, Vm2xy, Vmyy)
    test_polynomials = (q1, q2, q1, q2, q1, q2, q1, q2)
  elif(k == -2): # 14 dofs
    dof_list = (NQ1, NQ1, NQ1, NQ2, NQ2, NQ2, NQ3, NQ3, NQ3, NQ4, NQ4, NQ4, NQI, NQI)
    mon_list = (V1,V0, Vx,V0, Vy,V0, Vxx,V0, Vxy,V0, Vyy,V0, Vxxx,V3xyy,
                V0,V1, V0,Vx, V0,Vy, V0,Vxx, V0,Vxy, V0,Vyy, Vm3xxy,Vmyyy)
    test_polynomials = (q1, q2, q3, q1, q2, q3, q1, q2, q3, q1, q2, q3, qi10, qi01)
  elif(k == -3): # 22 dofs
    dof_list = (NQ1, NQ1, NQ1, NQ1, NQ2, NQ2, NQ2, NQ2, NQ3, NQ3, NQ3, NQ3, NQ4, NQ4, NQ4, NQ4,
                NQI, NQI, NQI, NQI, NQI, NQI)
    mon_list = (V1,V0, Vx,V0, Vy,V0, Vxx,V0, Vxy,V0, Vyy,V0, Vxxx,V0, Vxxy,V0, Vxyy,V0, Vyyy,V0, Vxxxx,V4xyyy, 
                V0,V1, V0,Vx, V0,Vy, V0,Vxx, V0,Vxy, V0,Vyy, V0,Vxxx, V0,Vxxy, V0,Vxyy, V0,Vyyy, Vm4xxxy,Vmyyyy)
    test_polynomials = (q1, q2, q3, q4, q1, q2, q3, q4, q1, q2, q3, q4, q1, q2, q3, q4, 
                        qi10, qi01, qix0, qi0x, qiy0, qi0y)
   
  return compute_matrix_and_inverse(dim, dof_list, mon_list, test_polynomials)



###############################################################################
# triangles


def build_matrix_Hdiv_P(k):
  """ compute Matrix of coefficients and its inverse to get the Raviart-Thomas
      (k>=0) or Brezzi-Douglas-Marini (k<0) basis functions on the triangle with
      vertices (0,0), (1,0), (0,1).

      Input: k denotes the order (k>=0 for RT, k<0 for BDM)
      Output: Matrix A=(a_ij) with a_ij = N_i(v_j) and its inverse. Here N_i
              are the degrees of freedom (functionals) and v_j are monomial
              basis vectors. The Raviart-Thomas or Brezzi-Douglas-Marini basis
              functions phi_j are then defined as phi_j = sum_i c_ij v_i. The
              monomial basis vectors are stored in the local variable Mon_List.
  """

  if(k >= 4 or k <=-4):
    print("So far the order k can only be 0, 1, 2, or 3 for Raviart-Thomas ",
          "elements and -1, -2, or -3 for Brezzi-Douglas-Marini elements")
    print("No matrix computed!")
    A = 0
    Ainv = 0
    return (A, Ainv)
  
  # the dimension of the local RT space
  dim = (k+1)*(k+3)
  if(k < 0): # BDM elements
    dim = (-k+1)*(-k+2)

  if(k == 0): # 3 dofs
    dof_list = (NT1, NT2, NT3)
    mon_list = (V1, V0, Vx,
                V0, V1, Vy)
    test_polynomials = (q1, q1, q1)
  elif(k == 1): # 8 dofs
    dof_list = (NT1, NT1, NT2, NT2, NT3, NT3, NTI, NTI)
    mon_list = (V1, V0, Vx, V0, Vy, V0, Vxx, Vxy,
                V0, V1, V0, Vx, V0, Vy, Vxy, Vyy)
    test_polynomials = (q1, q2, q1, q2, q1, q2, qi10, qi01)
  elif(k == 2): # 15 dofs
    dof_list = (NT1, NT1, NT1, NT2, NT2, NT2, NT3, NT3, NT3, NTI, NTI, NTI, NTI, NTI, NTI)
    mon_list = (V1,V0, Vx,V0, Vy,V0, Vxx,V0, Vxy,V0, Vyy,V0, Vxxx,Vxxy,Vxyy,
                V0,V1, V0,Vx, V0,Vy, V0,Vxx, V0,Vxy, V0,Vyy, Vxxy,Vxyy,Vyyy)
    test_polynomials = (q1, q2, q3, q1, q2, q3, q1, q2, q3,
                        qi10, qi01, qix0, qi0x, qiy0, qi0y)
  elif(k == 3): # 24 dofs
    dof_list = (NT1, NT1, NT1, NT1, NT2, NT2, NT2, NT2, NT3, NT3, NT3, NT3,
                NTI, NTI, NTI, NTI, NTI, NTI, NTI, NTI, NTI, NTI, NTI, NTI)
    mon_list = (V1,V0, Vx,V0, Vy,V0, Vxx,V0, Vxy,V0, Vyy,V0, Vxxx,V0, Vxxy,V0, Vxyy,V0, Vyyy,V0, Vxxxx,Vxxxy,Vxxyy,Vxyyy,
                V0,V1, V0,Vx, V0,Vy, V0,Vxx, V0,Vxy, V0,Vyy, V0,Vxxx, V0,Vxxy, V0,Vxyy, V0,Vyyy, Vxxxy,Vxxyy,Vxyyy,Vyyyy)
    test_polynomials = (q1, q2, q3, q4, q1, q2, q3, q4, q1, q2, q3, q4,
                        qi10, qi01, qix0, qi0x, qiy0, qi0y, 
                        qixx0, qi0xx, qiyy0, qi0yy, qixy0, qi0xy)
  elif(k == -1): # 6 dofs
    dof_list = (NT1, NT1, NT2, NT2, NT3, NT3)
    mon_list = (V1, V0, Vx, V0, Vy, V0,
                V0, V1, V0, Vx, V0, Vy)
    test_polynomials = (q1, q2, q1, q2, q1, q2,)
  elif(k == -2): # 12 dofs
    dof_list = (NT1, NT1, NT1, NT2, NT2, NT2, NT3, NT3, NT3, NTI, NTI, NTI)
    mon_list = (V1,V0, Vx,V0, Vy,V0, Vxx,V0, Vxy,V0, Vyy,V0,
                V0,V1, V0,Vx, V0,Vy, V0,Vxx, V0,Vxy, V0,Vyy,)
    test_polynomials = (q1, q2, q3, q1, q2, q3, q1, q2, q3, qi10, qi01, qicb)
  elif(k == -3): # 20 dofs
    dof_list = (NT1, NT1, NT1, NT1, NT2, NT2, NT2, NT2, NT3, NT3, NT3, NT3,
                NTI, NTI, NTI, NTI, NTI, NTI, NTI, NTI)
    mon_list = (V1,V0, Vx,V0, Vy,V0, Vxx,V0, Vxy,V0, Vyy,V0, Vxxx,V0, Vxxy,V0, Vxyy,V0, Vyyy,V0,
                V0,V1, V0,Vx, V0,Vy, V0,Vxx, V0,Vxy, V0,Vyy, V0,Vxxx, V0,Vxxy, V0,Vxyy, V0,Vyyy)
    test_polynomials = (q1, q2, q3, q4, q1, q2, q3, q4, q1, q2, q3, q4,
                        qi10, qi01, qi2x0, qi02y, qiyx, qicb, qicbx, qicby)
   
  return compute_matrix_and_inverse(dim, dof_list, mon_list, test_polynomials)



if __name__ == '__main__':
  np.set_printoptions(linewidth=180,precision=4)
  print("Raviart-Thomas 0 on quads")
  (a, ainv) = build_matrix_Hdiv_Q(0)
  #print(a)
  print(4*ainv)

  print("Raviart-Thomas 1 on quads")
  (a, ainv) = build_matrix_Hdiv_Q(1)
  #print(a)
  print(8*ainv)

  print("Raviart-Thomas 2 on quads")
  (a, ainv) = build_matrix_Hdiv_Q(2)
  #print(45*a)
  print(32*ainv)
  
  print("Raviart-Thomas 3 on quads")
  (a, ainv) = build_matrix_Hdiv_Q(3)
  #np.set_printoptions(linewidth=220,threshold=2000,suppress=True, precision=14)
  #print(3*a)
  print(3840*ainv)
  
  
  print("Brezzi-Douglas-Marini 1 on quads")
  (a, ainv) = build_matrix_Hdiv_Q(-1)
  #print(3*a)
  print(16*ainv)
  
  print("Brezzi-Douglas-Marini 2 on quads")
  (a, ainv) = build_matrix_Hdiv_Q(-2)
  #print(15*a)
  print(8*ainv)
  
  print("Brezzi-Douglas-Marini 3 on quads")
  (a, ainv) = build_matrix_Hdiv_Q(-3)
  #print(45*a)
  print(32*ainv)
  
  
  
  print("Raviart-Thomas 0 on triangles")
  (a, ainv) = build_matrix_Hdiv_P(0)
  #print(a)
  print(ainv)
  
  print("Raviart-Thomas 1 on triangles")
  (a, ainv) = build_matrix_Hdiv_P(1)
  #print(a)
  print(ainv)
  
  print("Raviart-Thomas 2 on triangles")
  (a, ainv) = build_matrix_Hdiv_P(2)
  #print(a)
  print(2*ainv)
  
  print("Raviart-Thomas 3 on triangles")
  (a, ainv) = build_matrix_Hdiv_P(3)
  #np.set_printoptions(linewidth=220,threshold=2000,suppress=True, precision=14)
  #print(a)
  print(ainv)
  
  
  print("Brezzi-Douglas-Marini 1 on triangles")
  (a, ainv) = build_matrix_Hdiv_P(-1)
  #print(a)
  print(ainv)
  
  print("Brezzi-Douglas-Marini 2 on triangles")
  (a, ainv) = build_matrix_Hdiv_P(-2)
  #print(a)
  print(ainv)
  
  print("Brezzi-Douglas-Marini 3 on triangles")
  (a, ainv) = build_matrix_Hdiv_P(-3)
  #print(a)
  print(ainv)
