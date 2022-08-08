#include <fstream>
#include <cmath>

#include "Boundary.h"
#include "MooNMD_Io.h"

Boundary::Boundary()
{
  parts.resize(0);
}

Boundary::Boundary(const std::string& PRM)
{
  parts.resize(0);
  initFromFile(PRM);
}

// read a PRM file (2D)
void Boundary::initFromFile(const std::string& PRM)
{
  std::ifstream ifile;
  ifile.open(PRM.c_str());

  if (!ifile)
  {
    Output::print(" *** Error(Boundary::initFromFile) I could not open ",PRM);
    exit(-1);
  };

  std::string theline;
  int nBoundaryParts, nBoundaryComponent,currentPart;
  int iType,nSpline,nParameters;

  int nBoundaryComponentTot = 0;
  
  getline(ifile,theline,'\n'); // NBCT
  ifile >> nBoundaryParts;
  //Output::print("n parts ",nBoundaryParts);

  getline(ifile,theline,'\n'); 
  parts.resize(nBoundaryParts);
  
  for (int i=0; i<nBoundaryParts; i++)  {
    getline(ifile,theline,'\n'); // IBCT
    ifile >> currentPart; // = i+1
    getline(ifile,theline,'\n'); 
    getline(ifile,theline,'\n'); // NCOMP
    ifile >> nBoundaryComponent;
    //Output::print("n comp ",nBoundaryComponent);
    getline(ifile,theline,'\n'); 
    parts[i].resize(nBoundaryComponent);
    
    getline(ifile,theline,'\n'); // ITYP NSPLINE NPAR
    for (int k=0; k<nBoundaryComponent; k++) {
      
      ifile >> iType;
      ifile >> nSpline;
      ifile >> nParameters;
      //Output::print("types ",iType,",",nSpline," ",nParameters);
      
      switch (iType) {
      case 1:
	parts[i][k].type = line;
	break;
      case 2:
	parts[i][k].type = circle;
	break;
      default:
	Output::print(" *** Error(Boundary::initFromFile) Boundary Component type ",
		      iType, " unknown");
	exit(-1);
	break;
      }
      parts[i][k].parameters.resize(2*nParameters);
      parts[i][k].partID = i;
      parts[i][k].localID = k;
      parts[i][k].globalID = k+nBoundaryComponentTot;

    }
    nBoundaryComponentTot += nBoundaryComponent;
    getline(ifile,theline,'\n'); 

  }
  getline(ifile,theline,'\n'); // PARAMETERS
  //Output::print(theline);
  for (int i=0; i<nBoundaryParts; i++) {
    for (unsigned int k=0; k<parts[i].size(); k++) {
      for (unsigned int p=0; p<parts[i][k].parameters.size()/2; p++) {

      ifile >> parts[i][k].parameters[2*p];
      ifile >> parts[i][k].parameters[2*p+1];
      //Output::print("param[2] ", parts[i][k].parameters[2*p], " ",
      //              parts[i][k].parameters[2*p+1]);
      }
      getline(ifile,theline,'\n');   
    }
  }
  
  Output::print<4>("  Boundary::initFromFile() finished reading ", PRM);
}

// return the Id of the boundary component of the point (x,y) and its local parameter (if any)
// return -1 if component is not found (inner node)
int Boundary::isOnComponent(double x, double y, double &t)
{
  t = -1;
  for (unsigned int i=0; i<parts.size(); i++)  {
    for (unsigned int k=0; k<parts[i].size(); k++)  {
      if (parts[i][k].type == line) {
	
	// check if (x,y) lies on the line (p1x,p1y) -- (p1x+dx,p1y+dy)
	// we do this checking the distances between the 3 points
	// P1-P2
	double p1p2 = parts[i][k].parameters[2]*parts[i][k].parameters[2] +
	  parts[i][k].parameters[3]*parts[i][k].parameters[3];
	// (x,y)-P1
	double vp1 = (x-parts[i][k].parameters[0])*(x-parts[i][k].parameters[0]) +
	  (y-parts[i][k].parameters[1])*(y-parts[i][k].parameters[1]);
	// (x,y)-P1 + (x,y)-P2
	double vp2 = (x-parts[i][k].parameters[0]-parts[i][k].parameters[2])*
	  (x-parts[i][k].parameters[0]-parts[i][k].parameters[2]) +
	  (y-parts[i][k].parameters[1]-parts[i][k].parameters[3])*
	  (y-parts[i][k].parameters[1]-parts[i][k].parameters[3]);
	//Output::print("distances: ",p1p2," ",vp1," ",vp2, " => ",std::abs(std::sqrt(vp1)+std::sqrt(vp2)-std::sqrt(p1p2)));
	if ( std::abs(std::sqrt(vp1)+std::sqrt(vp2)-std::sqrt(p1p2))<1e-10 ) {
	  t = std::sqrt(vp1/p1p2);
	  // if t=1, the vertex shall belong to the next line (with t=0)
	  if ( t<(1-1e-8)) {
	    t = t+k;
	    return i;
	    }
	  //else { cout << " --- " << t  << endl;}
	}
	
      } else if (parts[i][k].type == circle) {
	// check if (x,y) lies on the circle of center (p1x,p1y), radius p2x
	// and arc between p3x and p3y
	// step 1: check if the point is on the circle
	double vp = (x-parts[i][k].parameters[0])*(x-parts[i][k].parameters[0]) +
	  (y-parts[i][k].parameters[1])*(y-parts[i][k].parameters[1]);
	double r2 = parts[i][k].parameters[2]*parts[i][k].parameters[2];
	//Output::print("point ",x,",",y,"vp =",vp," r=",r2,"diff",std::abs(vp-r2));

	/** 
	    @warning to decide whether the point is on the circle might depend on precision
	    of the original geometry. We use now 1e-6
	*/
	if ( std::abs(vp-r2)<1e-6){
	  // step 2: check if the points is within the considered circle arc
	  double pi=acos(-1.0);
	  double angle0 = parts[i][k].parameters[4];
	  double angle1 = parts[i][k].parameters[5];
	  double theta=std::atan2(y-parts[i][k].parameters[1],x-parts[i][k].parameters[0]);
    bool full_circle = (2*pi - std::abs(angle0-angle1)) < 1.e-6;

	  if (angle0<angle1) {
	    // the boundary has "positive" sign
	    while (theta<angle0) {
	      theta= theta + 2*pi;
	    }
      if(full_circle && std::abs(theta - angle0) < 1.e-6)
      {
        t = k;
        return i;
      }
	    if (theta<angle1) {
	      // return the value of the angle renormalized between 0 and 1
	      t = k + (theta-angle0)/(angle1-angle0);
	      return i;
	    }
	  }
	  
	  else { // (angle0>angle1) 
	    // the boundary has "negative" sign 
	    while (theta<angle1) {
	      theta = theta + 2*pi;
	    }
      if(full_circle && std::abs(theta - angle1) < 1.e-6)
      {
        t = k;
        return i;
      }
	    if (theta<angle0) {
	      // return the value of the angle renormalized between 0 and 1
	      t = k + (theta-angle0)/(angle1-angle0);
	      return i;
	    }
	  }

	} // if (x,y) on the circle
      } // circle
      
    }
  }
  return -1;
}
 
void Boundary::info()
{
  Output::print(" ---- BOUNDARY INFORMATION ---- ");
  Output::print(" -- number of Parts: ",parts.size());
  for (unsigned int i=0; i<parts.size(); i++)  {
    Output::print(" --  Part ", i, ", number of components: ",parts[i].size());
  }
}
