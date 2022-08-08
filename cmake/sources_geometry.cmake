# CMakeLists.txt for subdirectory Geometry of ParMooN project. 
# Use only as subproject of ParMooN.
# 
# Change history:
# 2015/08/20 Clemens Bartsch: Rework to supply 2D and 3D library at once.
#

# Include header files.
list(APPEND PARMOON_INCLUDE_DIRS "include/Geometry")

# Source files used in 2D and 3D.
list(APPEND GEO_SOURCES "${PROJECT_SOURCE_DIR}/src/Geometry/BaseCell.C")
list(APPEND GEO_SOURCES "${PROJECT_SOURCE_DIR}/src/Geometry/BdCircle.C")
list(APPEND GEO_SOURCES "${PROJECT_SOURCE_DIR}/src/Geometry/BdLine.C")
list(APPEND GEO_SOURCES "${PROJECT_SOURCE_DIR}/src/Geometry/BdNonUniformSpline.C")
list(APPEND GEO_SOURCES "${PROJECT_SOURCE_DIR}/src/Geometry/BdPolygon.C")
list(APPEND GEO_SOURCES "${PROJECT_SOURCE_DIR}/src/Geometry/BdSpline.C")
list(APPEND GEO_SOURCES "${PROJECT_SOURCE_DIR}/src/Geometry/Boundary.C")
list(APPEND GEO_SOURCES "${PROJECT_SOURCE_DIR}/src/Geometry/BoundComp2D.C")
list(APPEND GEO_SOURCES "${PROJECT_SOURCE_DIR}/src/Geometry/BoundEdge.C")
list(APPEND GEO_SOURCES "${PROJECT_SOURCE_DIR}/src/Geometry/BoundPart.C")
list(APPEND GEO_SOURCES "${PROJECT_SOURCE_DIR}/src/Geometry/Collection.C")
list(APPEND GEO_SOURCES "${PROJECT_SOURCE_DIR}/src/Geometry/Domain.C")
list(APPEND GEO_SOURCES "${PROJECT_SOURCE_DIR}/src/Geometry/Enumerations_geometry.C")
list(APPEND GEO_SOURCES "${PROJECT_SOURCE_DIR}/src/Geometry/ErrorEstimator.C")
list(APPEND GEO_SOURCES "${PROJECT_SOURCE_DIR}/src/Geometry/GridCell.C")
list(APPEND GEO_SOURCES "${PROJECT_SOURCE_DIR}/src/Geometry/InterfaceJoint.C")
list(APPEND GEO_SOURCES "${PROJECT_SOURCE_DIR}/src/Geometry/IsoBoundEdge.C")
list(APPEND GEO_SOURCES "${PROJECT_SOURCE_DIR}/src/Geometry/IsoInterfaceJoint.C")
list(APPEND GEO_SOURCES "${PROJECT_SOURCE_DIR}/src/Geometry/IsoJointEqN.C")
list(APPEND GEO_SOURCES "${PROJECT_SOURCE_DIR}/src/Geometry/It_Between.C")
list(APPEND GEO_SOURCES "${PROJECT_SOURCE_DIR}/src/Geometry/It_EQ.C")
list(APPEND GEO_SOURCES "${PROJECT_SOURCE_DIR}/src/Geometry/It_EQLevel.C")
list(APPEND GEO_SOURCES "${PROJECT_SOURCE_DIR}/src/Geometry/It_Finest.C")
list(APPEND GEO_SOURCES "${PROJECT_SOURCE_DIR}/src/Geometry/It_LE.C")
list(APPEND GEO_SOURCES "${PROJECT_SOURCE_DIR}/src/Geometry/It_LELevel.C")
list(APPEND GEO_SOURCES "${PROJECT_SOURCE_DIR}/src/Geometry/It_OCAF.C")
list(APPEND GEO_SOURCES "${PROJECT_SOURCE_DIR}/src/Geometry/It_Search.C")
list(APPEND GEO_SOURCES "${PROJECT_SOURCE_DIR}/src/Geometry/Iterator.C")
list(APPEND GEO_SOURCES "${PROJECT_SOURCE_DIR}/src/Geometry/Joint.C")
list(APPEND GEO_SOURCES "${PROJECT_SOURCE_DIR}/src/Geometry/JointCollection.C")
list(APPEND GEO_SOURCES "${PROJECT_SOURCE_DIR}/src/Geometry/JointEqN.C")
list(APPEND GEO_SOURCES "${PROJECT_SOURCE_DIR}/src/Geometry/Line.C")
list(APPEND GEO_SOURCES "${PROJECT_SOURCE_DIR}/src/Geometry/LinesEval.C")
list(APPEND GEO_SOURCES "${PROJECT_SOURCE_DIR}/src/Geometry/MacroCell.C")
list(APPEND GEO_SOURCES "${PROJECT_SOURCE_DIR}/src/Geometry/Mapper.C")
list(APPEND GEO_SOURCES "${PROJECT_SOURCE_DIR}/src/Geometry/Mesh.C")
list(APPEND GEO_SOURCES "${PROJECT_SOURCE_DIR}/src/Geometry/PeriodicDomain.C")
list(APPEND GEO_SOURCES "${PROJECT_SOURCE_DIR}/src/Geometry/PeriodicJoint.C")
list(APPEND GEO_SOURCES "${PROJECT_SOURCE_DIR}/src/Geometry/Point.C")
list(APPEND GEO_SOURCES "${PROJECT_SOURCE_DIR}/src/Geometry/Quadrangle.C")
list(APPEND GEO_SOURCES "${PROJECT_SOURCE_DIR}/src/Geometry/ReadGeo.C")
list(APPEND GEO_SOURCES "${PROJECT_SOURCE_DIR}/src/Geometry/Rectangle.C")
list(APPEND GEO_SOURCES "${PROJECT_SOURCE_DIR}/src/Geometry/RefDesc.C")
list(APPEND GEO_SOURCES "${PROJECT_SOURCE_DIR}/src/Geometry/RefinementStrategy.C")
# list(APPEND GEO_SOURCES "${PROJECT_SOURCE_DIR}/src/Geometry/ShapeDesc.C") #CB why not?
# list(APPEND GEO_SOURCES "${PROJECT_SOURCE_DIR}/src/Geometry/SubDomainEdge3D.C") #CB why not?
list(APPEND GEO_SOURCES "${PROJECT_SOURCE_DIR}/src/Geometry/SubDomainHaloJoint.C")
list(APPEND GEO_SOURCES "${PROJECT_SOURCE_DIR}/src/Geometry/SubDomainJoint.C")
list(APPEND GEO_SOURCES "${PROJECT_SOURCE_DIR}/src/Geometry/Tests.C")
list(APPEND GEO_SOURCES "${PROJECT_SOURCE_DIR}/src/Geometry/Triangle.C")
list(APPEND GEO_SOURCES "${PROJECT_SOURCE_DIR}/src/Geometry/Vertex.C")

# Source files only used in 2D
list(APPEND GEO_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Geometry/InnerInterfaceJoint.C")
list(APPEND GEO_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Geometry/Parallelogram.C")

# Source files only used in 3D
list(APPEND GEO_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Geometry/BdCylinder.C")
list(APPEND GEO_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Geometry/BDEdge3D.C")
list(APPEND GEO_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Geometry/BdNoPRM.C")
list(APPEND GEO_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Geometry/BdPlane.C")
list(APPEND GEO_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Geometry/BdSphere.C")
list(APPEND GEO_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Geometry/BdWall.C")
list(APPEND GEO_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Geometry/BoundComp3D.C")
list(APPEND GEO_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Geometry/BoundFace.C")
list(APPEND GEO_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Geometry/Brick.C")
list(APPEND GEO_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Geometry/Edge.C")
list(APPEND GEO_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Geometry/Hexahedron.C")
list(APPEND GEO_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Geometry/InnerEdge.C")
list(APPEND GEO_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Geometry/InterfaceJoint3D.C")
list(APPEND GEO_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Geometry/IsoBoundFace.C")
list(APPEND GEO_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Geometry/IsoEdge3D.C")
list(APPEND GEO_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Geometry/IsoInterfaceJoint3D.C")
list(APPEND GEO_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Geometry/Tetrahedron.C")
list(APPEND GEO_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Geometry/Parallelogram.C")


list(APPEND PARMOON_SOURCES_2D ${GEO_SOURCES} ${GEO_SOURCES_2D})
list(APPEND PARMOON_SOURCES_3D ${GEO_SOURCES} ${GEO_SOURCES_3D})
