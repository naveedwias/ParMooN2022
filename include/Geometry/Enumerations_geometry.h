#ifndef ENUMERATIONS_GEOMETRY_H
#define ENUMERATIONS_GEOMETRY_H
#include <iosfwd>

enum BoundTypes {Line = 1, Circle, Spline2D, Polygon, NonUniformSpline2D,
                 Plane = 10, Sphere, Cylinder, 
                 Wall = 100,
                 NoPRM = 4711 };

constexpr int N_ITERATORS = 9;
enum Iterators {It_EQ, It_LE, It_Finest, It_EQLevel, It_LELevel,
                It_Between, It_OCAF};

constexpr int N_MAPPER = 34;
enum Mapper {MapTriReg0,  MapTriReg1, MapTriReg2,
             MapTriBis00, MapTriBis01, MapTriBis02,
             MapTriBis10, MapTriBis11, MapTriBis12,
             MapTriBis20, MapTriBis21, MapTriBis22,
             MapTriBis010, MapTriBis011, MapTriBis012,
             MapTriBis020, MapTriBis021, MapTriBis022,
             MapTriBis100, MapTriBis101, MapTriBis102,
             MapTriBis120, MapTriBis121, MapTriBis122,
             MapTriBis200, MapTriBis201, MapTriBis202,
             MapTriBis210, MapTriBis211, MapTriBis212,
             MapQuadReg0, MapQuadReg1, MapQuadReg2, MapQuadReg3};

constexpr int N_REFDESC = 76;
// We use TriBis0 - Quad2Conf3 for generation of closures. Do not
// put any other descriptors between these.
enum Refinements {NoRef, LineReg, TriReg, QuadReg, ParallReg, RectReg,
                  TriBary, TriBis0, TriBis1, TriBis2, 
                  TriBis01, TriBis02, TriBis10, TriBis12, TriBis20, TriBis21,
                  QuadBis0, QuadBis1,
                  Quad1Conf0, Quad1Conf1, Quad1Conf2, Quad1Conf3,
                  Quad2Conf0, Quad2Conf1, Quad2Conf2, Quad2Conf3,
                  QuadToTri0, QuadToTri1, TetraReg, TetraBary, 
                  TetraReg0, TetraReg1, TetraReg2, 
                  TetraBis0, TetraBis1, TetraBis2, TetraBis3, TetraBis4, TetraBis5,
                  TetraBis01, TetraBis02, TetraBis03, TetraBis04, TetraBis05,
                  TetraBis10, TetraBis12, TetraBis13, TetraBis14, TetraBis15,
                  TetraBis20, TetraBis21, TetraBis23, TetraBis24, TetraBis25,
                  TetraBis30, TetraBis32, TetraBis34, TetraBis35,
                  TetraBis40, TetraBis41, TetraBis43, TetraBis45,
                  TetraBis51, TetraBis52, TetraBis53, TetraBis54,
                  TetraQuad0, TetraQuad1, TetraQuad2, TetraQuad3, TetraQuad4, TetraQuad5,
                  HexaReg, BrickReg};

enum RefinementMarks { NoRefinement, Refinement, DeRefinement };

enum class RefinementStrategyType
{
  // this refinement strategy marks cells if the error is above some tolerance
  // and relaxes the tolerance if not enough cells have been marked yet
  MaximalEstimatedLocalError = 0,
  // this refinement strategy marks cells if the error is above some tolerance
  // and relaxes the tolerance if not enough cells have been marked yet or the
  // accumulated error of the cells is not some prescribed portion of the global
  // error
  PortionOfEstimatedGlobalError,
  // this refinement strategy is computationally more heavy than the other
  // strategies and aims at refining as few cells as possible while sustaining
  // uniform convergence of the adaptive process
  EquilibrationStrategy
};

constexpr int N_SHAPES = 8;
enum Shapes {S_Line, Triangle, Quadrangle, Parallelogram, Rectangle,
             Tetrahedron, Hexahedron, Brick};

// enable output of refinement strategy type
std::ostream &operator<<(std::ostream &os, RefinementStrategyType type);

#endif //ENUMERATIONS_GEOMETRY_H
