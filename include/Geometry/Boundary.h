/** ************************************************************************ 
*
* @class     Boundary
* @brief     A class for handling analytical boundaries (as in PRM files)
*            In ParMooN (2D), the boundary is composed of several boundary 
*            parts, and each part can contain multiple boundary components
*            (line, circle, ....)
*            This class contains the whole boundary as a vector of
*            vectors(boundarycomponents). Each component knows its 
*            own boundary part
* @author    Alfonso Caiazzo 
* @date      30.03.16
* @History 
 ************************************************************************  */

#include <vector>
#include <string>

#ifndef __BOUNDARY__
#define __BOUNDARY__

enum boundaryType {line, circle};

///@brief a (temporary) object to improve the boundary implementation  (Alfonso)
struct BoundaryComponent
{
  boundaryType type;
  ///@brief global id of the component w.r.t to the whole boundaries
  int globalID;
  ///@brief id of the boundary part the component belongs to
  int partID;
  ///@brief id of the component w.r.t to the its part
  int localID;
  ///@brief id for adding physical (e.g. boundary conditions) properties
  int physicalID;
  ///@brief geometrical parameters describing the component (read from PRM file)
  std::vector<double> parameters;
  int nParameters;
  ///@todo handle inner interfaces
};

class Boundary
{
 public:
  Boundary();
  explicit Boundary(const std::string& prmfile);
  ~Boundary(){};

  ///@ initialize the list of bd components from a PRM file
  void initFromFile(const std::string& prmfile);

  /**
     @brief return the boundary component of point (x,y)
     @param[out] t is the local abscissa (between 0 and 1)
  */
  int isOnComponent(double x,double y, double& t);

  void info();

  ///@todo the list of part should be protected
  std::vector< std::vector< BoundaryComponent > > parts;
  
};

#endif
