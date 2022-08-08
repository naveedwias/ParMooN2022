#include <fstream>

#include "Mesh.h"
#include "MooNMD_Io.h"
#include <assert.h>

#define TETRA_FACE_VERTEX { {0, 1, 2}, {0, 3, 1}, {2, 1, 3}, {0, 2, 3} }
#define HEXA_FACE_VERTEX { {0, 1, 2, 3}, {0, 1, 5, 4}, {1, 2, 6, 5}, {3, 2, 6, 7}, {0, 3, 7, 4}, {4, 5, 6, 7} }

void stripSpace(std::string &str);
void check_hexahedron(const meshHexahedron& hex,
                      const std::vector<meshNode>& vertices);

// default initialization of the mesh (dimension = 0, no elements)
Mesh::Mesh()
{
  dimension = -1;
  n_boundary_faces = 0;
  hasBothTriaAndQuads = false;
}

// initialize from a file
// note: dimension is set to 2, but it will be changed (if necessary) reading the file
Mesh::Mesh(const std::string& filename)
{
  dimension = -1;
  n_boundary_faces = 0;
  hasBothTriaAndQuads = false;

  readFromFile(filename);
}

Mesh::Mesh(const std::string& filename, const std::string& filenameBoundary)
{
  dimension = -1;
  hasBothTriaAndQuads = false;
  n_boundary_faces = 0;

  readFromFile(filename);
  setBoundary(filenameBoundary);
}

Mesh::~Mesh()
{
}

// read the mesh data from a file (.mesh)
void Mesh::readFromFile(const std::string& filename)
{
  std::ifstream ifile;
  ifile.open(filename.c_str());
  if (!ifile)
  {
    ErrThrow(" *** Error(Mesh::readFromFile) I could not open ", filename);
  };

  std::string line;

  // read header and dimension
  do
  {
    // read the whole line up to [enter], if the line is empty or starting with #
    getline(ifile,line,'\n');
    stripSpace(line); 
  }
  while (line.find("Dimension") == std::string::npos);

  line.erase(0, 9);
  if (line.length() > 0)
  {
    sscanf(line.c_str(),"%d", (int*)&dimension);
  }
  else
  {
    ifile >> dimension;
  }

  Output::root_info<4>("Mesh", "Read dimension: ", dimension);

  if (dimension < 2 || dimension > 3)
  {
    ErrThrow("I could not read the mesh dimension. Check the .mesh file");
  }

  // read nodes
  unsigned int numberOfNodes;
  do
  {
    getline(ifile,line,'\n');
    stripSpace(line);
  }
  while (line.find("Vertices") == std::string::npos);

  ifile >> numberOfNodes;
  Output::root_info<4>("Mesh", "Read n of nodes: ", numberOfNodes);

  vertex.resize(numberOfNodes);
  for (unsigned int i = 0; i < numberOfNodes; i++)
  {
    if (dimension == 2)
    {
      ifile >> vertex[i].x >> vertex[i].y >> vertex[i].reference;
    }
    else
    {
      ifile >> vertex[i].x >> vertex[i].y >> vertex[i].z >> vertex[i].reference;
    }
  }

  int current_line_number = ifile.tellg();

  // read edges
  unsigned int numberOfEdges;
  do
  {
    getline(ifile,line,'\n');
    stripSpace(line); // borra espacios en blanco
  } 
  while (line.find("Edges") == std::string::npos && !ifile.eof());

  if (line.find("Edges") != std::string::npos)
  {
    ifile >> numberOfEdges;

    if (numberOfEdges)
    {
      edge.resize(numberOfEdges);

      for (unsigned int i = 0; i < numberOfEdges; i++)
      {
        ifile >> edge[i].nodes[0] >> edge[i].nodes[1] >> edge[i].reference;
      }
    }
  }

  ifile.clear();
  ifile.seekg(current_line_number);

  // Triangles
  unsigned int numberOfTriangles = 0;
  do
  {
    getline(ifile,line,'\n');
    stripSpace(line);

    if (line.find("#") != std::string::npos)
    {
      line = "####COMMENT###";
    }
  }
  while (line.find("Triangles") == std::string::npos && !ifile.eof());

  if (line.find("Triangles") != std::string::npos)
  {
    ifile >> numberOfTriangles;

    if (numberOfTriangles)
    {
      triangle.resize(numberOfTriangles);

      for (unsigned int i = 0; i < numberOfTriangles; i++)
      {
        for (unsigned int k = 0; k < 3; k++)
        {
          ifile >> triangle[i].nodes[k];
        }
        ifile >> triangle[i].reference;
      }
    }
  }

  ifile.clear();
  ifile.seekg(current_line_number);

  // Quads
  unsigned int numberOfQuads = 0;
  do
  {
    getline(ifile,line,'\n');
    stripSpace(line); 
    if (line.find("#") != std::string::npos)
    {
      line = "####COMMENT###";
    }
  }
  while (line.find("Quadrilaterals") == std::string::npos && !ifile.eof());

  if (line.find("Quadrilaterals") != std::string::npos)
  {
    ifile >> numberOfQuads;

    if (numberOfQuads)
    {
      quad.resize(numberOfQuads);

      for (unsigned int i = 0; i < numberOfQuads; i++)
      {
        for (unsigned int k = 0; k < 4; k++)
        {
          ifile >> quad[i].nodes[k];
        }
        ifile >> quad[i].reference;
      }
    }
  }

  ifile.clear();
  ifile.seekg(current_line_number);

  if (quad.size() && triangle.size())
  {
    hasBothTriaAndQuads = true;
  }

  // tetra
  unsigned int numberOfTetra = 0;
  do
  {
    getline(ifile,line,'\n');
    stripSpace(line);

    if (line.find("#") != std::string::npos)
    {
      line = "####COMMENT###";
    }
  }
  while (line.find("Tetrahedra") == std::string::npos && !ifile.eof());

  if (line.find("Tetrahedra") != std::string::npos)
  {
    ifile >> numberOfTetra;

    if (numberOfTetra)
    {
      is_tetramesh = true;
      tetra.resize(numberOfTetra);
      for (unsigned int i = 0; i < numberOfTetra; i++)
      {
        for (unsigned int k = 0; k < 4; k++)
        {
          ifile >> tetra[i].nodes[k];
        }

        ifile >> tetra[i].reference;
      }
    }
  }

  ifile.clear();
  ifile.seekg(current_line_number);

  // hexa
  unsigned int numberOfHexa = 0;
  do
  {
    getline(ifile,line,'\n');
    stripSpace(line);

    if (line.find("#") != std::string::npos)
    {
      line = "####COMMENT###";
    }
  }
  while (line.find("Hexahedra") == std::string::npos && !ifile.eof());

  if (line.find("Hexahedra") != std::string::npos)
  {
    ifile >> numberOfHexa;

    if (numberOfHexa)
    {
      is_hexamesh = true;
      hexa.resize(numberOfHexa);

      for (unsigned int i = 0; i < numberOfHexa; i++)
      {
        for (unsigned int k = 0; k < 8; k++)
        {
          ifile >> hexa[i].nodes[k];
        }

        ifile >> hexa[i].reference;
        check_hexahedron(hexa[i], vertex);
      }
    }
  }

  ifile.clear();
  ifile.seekg(current_line_number);
  ifile.close();

  if(numberOfHexa != 0 && numberOfTetra != 0)
  {
    ErrThrow("Missing implementation for reading mesh files with hexahedra "
             "and tetrahedra.");
  }

  if(numberOfTetra + numberOfHexa == 0 && dimension == 3)
  {
    // this should be 2D and is reset here after a check.
    // All vertices should have the same z-component
    // we don't support 2D meshes embedded in 3D (yet)
    for (unsigned int i = 1; i < numberOfNodes; i++)
    {
      if(vertex[i].z != vertex[0].z)
      {
        ErrThrow("Cannot read a 3D mesh with no tetrahedra and no hexahedra.");
      }
    }

    dimension = 2;
  }

  if((numberOfTriangles) && (dimension == 2))
  {
    correct_numbering_of_vertices_if_needed<meshTriangle>(triangle);
  }

  if((numberOfQuads) && (dimension == 2))
  {
    correct_numbering_of_vertices_if_needed<meshQuad>(quad);
  }
}

void Mesh::hashTriFaces()
{
  meshTrifaceHash.clear();

  std::unordered_map<int, int> bucketCount;

  for (const meshTriangle &tri: triangle)
  {
    int h = tri.nodes[0] + tri.nodes[1] + tri.nodes[2];

    if (bucketCount.count(h) > 0)
    {
      bucketCount[h]++;
    }
    else
    {
      bucketCount[h] = 1;
    }
  }

  for (auto& p: bucketCount)
  {
    meshTrifaceHash.reserve(p.first, p.second);
  }

  for (size_t i = 0; i < triangle.size(); ++i)
  {
    meshTriangle &tri = triangle.at(i);

    meshTrifaceHash.add(i, tri.nodes[0] + tri.nodes[1] + tri.nodes[2]);
  }
}

void Mesh::hashQuadFaces()
{
  meshQuadfaceHash.clear();

  std::unordered_map<int, int> bucketCount;

  for (const meshQuad &q: quad)
  {
    int h = q.nodes[0] + q.nodes[1]
          + q.nodes[2] + q.nodes[3];

    if (bucketCount.count(h) > 0)
    {
      bucketCount[h]++;
    }
    else
    {
      bucketCount[h] = 1;
    }
  }

  for (auto& p: bucketCount)
  {
    meshQuadfaceHash.reserve(p.first, p.second);
  }

  for (size_t i = 0; i < quad.size(); ++i)
  {
    meshQuad &q = quad.at(i);

    meshQuadfaceHash.add(i, q.nodes[0] + q.nodes[1]
                          + q.nodes[2] + q.nodes[3]);
  }
}


// fill the list mapping each face to the neighbors tetra
// -1 indicates that a face is on the boundary
void Mesh::createFaceToTetrahedraMap()
{
  // create the meshTrifaceHash vector (used in findTriFace)
  if (meshTrifaceHash.size() == 0)
  {
    Output::root_info("Mesh",
      "Mesh::createFaceTetrahedraMap() is creating the hash vector");
    this->hashTriFaces();
  }

  // this vector describes the order in which the local vertices of a tetrahedra
  // (0,1,2,3) appear on the faces
  const int FaceVertex[][3] = TETRA_FACE_VERTEX;

  this->faceToTetra.resize(this->triangle.size());

  for (unsigned int i = 0; i < this->triangle.size(); ++i)
  {
    // set all values to -1
    this->faceToTetra[i].resize(2);
    this->faceToTetra[i][0] = -1;
    this->faceToTetra[i][1] = -1;
  }

  for (unsigned int i = 0; i < this->tetra.size(); ++i)
  {
    for (int face = 0; face < 4; ++face)
    {
      // find a triface with the same vertices
      int triface = findTriFace(tetra[i].nodes[FaceVertex[face][0]],
        tetra[i].nodes[FaceVertex[face][1]],
        tetra[i].nodes[FaceVertex[face][2]]);

      assert (triface != -1);

      // if the element [triface,0] is still = -1, set it = to the current tetrahedra (i)
      if(this->faceToTetra[triface][0] == -1)
      {
        this->faceToTetra[triface][0] = i;
      }
      else
      {
        // if the face has been already associated to a tetrahedra,
        // e.g. faceToTetra[triface][0]=k, then set faceToTetra[triface][1] equal to i

        assert(this->faceToTetra[triface][1] == -1);

        this->faceToTetra[triface][1] = i;
      }
    }
  }
}

// fill the list mapping each face to the neighbors hexa
// -1 indicates that a face is on the boundary
void Mesh::createFaceToHexahedraMap()
{
  // create the meshQuadfaceHash vector (used in findQuadFace)
  if (meshQuadfaceHash.size() == 0)
  {
    Output::root_info("Mesh",
      "Mesh::createFaceToHexahedraMap() is creating the hash vector");
    this->hashQuadFaces();
  }

  // this vector describes the order in which the local vertices of a hexahedra
  // (0,1,2,3,4,5,) appear on the faces
  const int FaceVertex[][4] = HEXA_FACE_VERTEX;

  this->faceToHexa.resize(this->quad.size());

  for (unsigned int i = 0; i < this->quad.size(); ++i)
  {
    // set all values to -1
    this->faceToHexa[i].resize(2);
    this->faceToHexa[i][0] = -1;
    this->faceToHexa[i][1] = -1;
  }

  for (unsigned int i = 0; i < this->hexa.size(); ++i)
  {
    for (int face = 0; face < 6; ++face)
    {
      // find a quadface with the same vertices
      int quadface = findQuadFace(hexa[i].nodes[FaceVertex[face][0]],
                                  hexa[i].nodes[FaceVertex[face][1]],
                                  hexa[i].nodes[FaceVertex[face][2]],
                                  hexa[i].nodes[FaceVertex[face][3]]);

      assert (quadface != -1);

      // if the element [quadface,0] is still = -1,
      // set it = to the current hexahedra (i)
      if (this->faceToHexa[quadface][0] == -1)
      {
        this->faceToHexa[quadface][0] = i;
      }
      else
      {
        // if the face has been already associated to a Hexahedra,
        // e.g. faceToHexa[quadface][0]=k,
        // then set faceToHexa[quadface][1] equal to i

        assert(this->faceToHexa[quadface][1] == -1);

        this->faceToHexa[quadface][1]= i;
      }
    }
  }
}

// find the index of the face with vertices a,b,c
// return -1 if face is not found
int Mesh::findTriFace(int a, int b, int c)
{
  int hash = a + b + c;

  if (this->meshTrifaceHash.size() == 0)
  {
    Output::root_info("Mesh", "Mesh::findTriFace() is creating the hash vector");
    this->hashTriFaces();
  }

  if (meshTrifaceHash.contains(hash))
  {
    for (int i: meshTrifaceHash.at(hash))
    {
      int found = 0;

      for (int j = 0; j < 3; ++j)
      {
        int vertex = triangle[i].nodes[j];

        if (a == vertex || b == vertex || c == vertex)
        {
          ++found;
        }
      }

      if (found == 3)
      {
        return i;
      }
    }
  }

  /*Output::root_warn("Mesh", "Mesh::findTriFace() ** WARNING ** no face with vertices ",
                a,",",b,",",c," found.");*/
  return -1;
}


// find the index of the face with vertices a,b,c,d
// return -1 if face is not found
int Mesh::findQuadFace(int a, int b, int c, int d)
{
  int hash = a + b + c + d;

  if (this->meshQuadfaceHash.size() == 0)
  {
    Output::root_info("Mesh", "Mesh::findQuadFace() is creating the hash vector");
    this->hashQuadFaces();
  }

  if (meshQuadfaceHash.contains(hash))
  {
    for (int i: meshQuadfaceHash.at(hash))
    {
      int found = 0;

      for (int j = 0; j < 4; ++j)
      {
        int vertex = quad[i].nodes[j];

        if (a == vertex || b == vertex || c == vertex || d == vertex)
        {
          ++found;
        }
      }

      if (found == 4)
      {
        return i;
      }
    }
  }

  /*Output::root_warn("Mesh", "Mesh::findQuadFace() ** WARNING ** no face with vertices ",
                a, ", ", b, ", ", c, ", ", d, " found.");*/

  return -1;
}

void Mesh::computeNumberOfBoundaryFaces()
{
  const std::vector<std::vector<int>>* ptr_faceToCell;

  void (Mesh::*func_createFaceToCellMap)();

  size_t nb_face, nb_face_expect;
  std::string str_face, str_faceToCell;

  if (hexa.size() == 0) // mesh with only tetrahedra
  {
    ptr_faceToCell = &faceToTetra;
    func_createFaceToCellMap = &Mesh::createFaceToTetrahedraMap;
    nb_face = triangle.size();
    nb_face_expect = 2 * tetra.size(); // add n_boundary_faces/2 when computed
    str_face = "triangles";
    str_faceToCell = "faceToTetra";
  }
  else if (tetra.size() == 0)  // mesh with only hexahedra
  {
    ptr_faceToCell = &faceToHexa;
    func_createFaceToCellMap = &Mesh::createFaceToHexahedraMap;
    nb_face = quad.size();
    nb_face_expect = 3 * hexa.size(); // add n_boundary_faces/2 when computed
    str_face = "quads";
    str_faceToCell = "faceToHexa";
  }
  else
  {
    ErrThrow(" *** Error(Mesh::computeNumberOfBoundaryFaces)");
  }

  int n_inner_faces = 0;
  this->n_boundary_faces = 0;

  if (ptr_faceToCell->size() == 0)
  {
    Output::root_info("Mesh", "Mesh::computeNumberOfBoundaryFaces() is "
      "creating the ", str_faceToCell, " map");

    (this->*func_createFaceToCellMap)();
  }

  for (unsigned int i = 0; i < ptr_faceToCell->size(); i++)
  {
    // increase the number of boundary faces if face[i] is on the boundary
    if((ptr_faceToCell->at(i)[0] < 0) || (ptr_faceToCell->at(i)[1] < 0))
    {
      this->n_boundary_faces++;
    }
    else
    {
      n_inner_faces++;
    }
  }

  nb_face_expect += n_boundary_faces / 2;

  Output::root_info("Mesh", "Mesh::computeNumberOfBoundaryFaces() - I found ",
                this->n_boundary_faces, " boundary faces and ",
                n_inner_faces, " inner faces");

  // check if faces are missing
  if (nb_face != nb_face_expect)
  {
    Output::root_warn("Mesh", "Mesh::computeNumberOfBoundaryFaces() ** WARNING: "
                  "it looks like some face is missing, ",
                  str_face, " = ", nb_face,
                  ", expected = ", nb_face_expect, " **");
  }
}

/**
   This function checks if the mesh read from file contains all the faces
   (inner+boundary). If not, the missing faces are created.

   Note: this function create a local hash vector (to handle the faces
   read from the original file). This vector is deleted at the end, and
   the meshTrifaceHash has to be created later.
   
   This function should be called only once.
 */
void Mesh::createInnerFaces()
{
  Output::root_info("Mesh::createInnerFaces()");

  unsigned int nb_vertex_per_face, nb_face_per_cell;

  if (hexa.size() == 0)  // mesh with only tetrahedra
  {
    nb_vertex_per_face = 3;
    nb_face_per_cell = 4;
  }
  else if (tetra.size() == 0)  // mesh with only hexahedra
  {
    nb_vertex_per_face = 4;
    nb_face_per_cell = 6;
  }
  else
  {
    ErrThrow(" *** Error(Mesh::createInnerFaces)");
  }

  if (is_tetramesh)
  {
    std::vector<std::vector<int>> localHash;
    int n_points = vertex.size();

    localHash.resize(nb_vertex_per_face * n_points);
    for (unsigned int i = 0; i < localHash.size(); i++)
    {
      localHash[i].resize(0);
    }

    for (size_t i = 0; i < triangle.size(); ++i)
    {
      int hash = triangle[i].nodes[0]
        + triangle[i].nodes[1]
        + triangle[i].nodes[2];

      if (localHash[hash].size() == 0)
      {
        localHash[hash].resize(2);
        localHash[hash][0] = 1;
        localHash[hash][1] = i;
      }
      else
      {
        localHash[hash][0]++;
        localHash[hash].push_back(i);
      }
    }

    // now check if there are missing faces
    int n_created_faces = 0;

    // map vertex id -> face
    const int FaceVertex[][3] = TETRA_FACE_VERTEX;

    // loop over tetra
    for (unsigned int i = 0; i < tetra.size(); i++)
    {
      // loop over faces
      for (unsigned int face = 0; face < nb_face_per_cell; face++)
      {
        int a = tetra[i].nodes[FaceVertex[face][0]];
        int b = tetra[i].nodes[FaceVertex[face][1]];
        int c = tetra[i].nodes[FaceVertex[face][2]];

        int hash = a + b + c;

        // check if face j is in the list of triangles
        if (localHash[hash].size() == 0)
        {
          // face is missing
          n_created_faces++;

          meshTriangle new_triangle;
          new_triangle.nodes[0] = a;
          new_triangle.nodes[1] = b;
          new_triangle.nodes[2] = c;
          new_triangle.reference = 0;

          triangle.push_back(new_triangle);

          // update the hash vector (to avoid creating twice the same face)
          localHash[hash].resize(2);
          localHash[hash][0] = 1;
          localHash[hash][1] = triangle.size() - 1;
        }
        else
        {
          int count = localHash[hash][0];
          int faces_found = 0;
          for(int i = 1; i <= count; ++i)
          {
            int triface = localHash[hash][i];
            unsigned int found = 0;
            for (unsigned int j = 0; j < nb_vertex_per_face; ++j)
            {
              int vertex = triangle[triface].nodes[j];
              if (a == vertex || b == vertex || c == vertex)
              {
                ++found;
              }
            }

            if(found == nb_vertex_per_face)
            {
              faces_found++;
            }
          }

          if (faces_found == 0)
          {
            // face is missing
            n_created_faces++;

            meshTriangle new_triangle;
            new_triangle.nodes[0] = a;
            new_triangle.nodes[1] = b;
            new_triangle.nodes[2] = c;
            new_triangle.reference = 0;

            triangle.push_back(new_triangle);

            // update the hash vector (to avoid creating twice the same face)
            localHash[hash][0]++;
            localHash[hash].push_back(triangle.size() - 1);
          }
        }
      } // loop over faces
    } // loop over tetra

    Output::root_info("Mesh", "Mesh::createInnerFaces() finished. Created ",
                  n_created_faces, " triangles");
  }
  else if (is_hexamesh)
  {
    std::vector<std::vector<int>> localHash;

    int n_points = vertex.size();
    localHash.resize(nb_vertex_per_face * n_points);

    for (unsigned int i = 0; i < localHash.size(); i++)
    {
      localHash[i].resize(0);
    }

    for (size_t i = 0; i < quad.size(); ++i)
    {
      int hash = quad[i].nodes[0]
        + quad[i].nodes[1]
        + quad[i].nodes[2]
        + quad[i].nodes[3];

      if (localHash[hash].size() == 0)
      {
        localHash[hash].resize(2);
        localHash[hash][0] = 1;
        localHash[hash][1] = i;
      }
      else
      {
        localHash[hash][0]++;
        localHash[hash].push_back(i);
      }
    }

    // now check if there are missing faces
    int n_created_faces = 0;

    // map vertex id -> face
    const int FaceVertex[][4] = HEXA_FACE_VERTEX;

    // loop over hexa
    for (unsigned int i = 0; i < hexa.size(); i++)
    {
      // loop over faces
      for (unsigned int face = 0; face < nb_face_per_cell; face++)
      {
        int a = hexa[i].nodes[FaceVertex[face][0]];
        int b = hexa[i].nodes[FaceVertex[face][1]];
        int c = hexa[i].nodes[FaceVertex[face][2]];
        int d = hexa[i].nodes[FaceVertex[face][3]];

        int hash = a + b + c + d;

        // check if face j is in the list of quads
        if (localHash[hash].size() == 0)
        {
          // face is missing

          n_created_faces++;

          meshQuad new_quad;
          new_quad.nodes[0] = a;
          new_quad.nodes[1] = b;
          new_quad.nodes[2] = c;
          new_quad.nodes[3] = d;
          new_quad.reference = 0;

          quad.push_back(new_quad);

          // update the hash vector (to avoid creating twice the same face)
          localHash[hash].resize(2);
          localHash[hash][0] = 1;
          localHash[hash][1] = quad.size() - 1;
        }
        else
        {
          int count = localHash[hash][0];

          int faces_found = 0;
          for(int i = 1; i <= count; ++i)
          {
            int quadface = localHash[hash][i];
            unsigned int found = 0;

            for (unsigned int j = 0; j < nb_vertex_per_face; ++j)
            {
              int vertex = quad[quadface].nodes[j];
              if (a == vertex || b == vertex || c == vertex || d == vertex)
              {
                ++found;
              }
            }

            if (found == nb_vertex_per_face)
            {
              faces_found++;
            }
          }

          if (faces_found == 0)
          {
            // face is missing

            n_created_faces++;

            meshQuad new_quad;
            new_quad.nodes[0] = a;
            new_quad.nodes[1] = b;
            new_quad.nodes[2] = c;
            new_quad.nodes[3] = d;
            new_quad.reference = 0;

            quad.push_back(new_quad);

            // update the hash vector (to avoid creating twice the same face)
            localHash[hash][0]++;
            localHash[hash].push_back(quad.size() - 1);
          }
        }
      } // loop over faces
    } // loop over hexa

    Output::root_info("Mesh", "Mesh::createInnerFaces() finished. Created ",
                  n_created_faces, " quads");
  }
}

// read a PRM file into a Boundary class
void Mesh::setBoundary(const std::string& PRM)
{
  boundary.initFromFile(PRM);
}

// write the mesh onto a .mesh file
void Mesh::writeToMesh(const std::string& filename, bool includeInteriorFaces)
{
  std::ofstream ofile;
  ofile.open(filename.c_str(),std::ios::out);
  if (!ofile)
  {
    ErrThrow(" *** Error(Mesh::writeToFile) I could not open ", filename);
  }

  Output::print("writing mesh on ", filename);

  ofile.precision(18);

  // header (.mesh)
  ofile << "MeshVersionFormatted 1" << endl;
  ofile << "Dimension " << dimension << endl << endl;
  ofile << "Vertices" << endl;
  ofile << vertex.size() << endl;

  // write nodes
  for (unsigned int i = 0; i < vertex.size(); i++)
  {
    ofile << vertex[i].x << " " << vertex[i].y << " ";

    if (dimension > 2)
    {
      ofile << vertex[i].z << " ";
    }

    ofile << vertex[i].reference << endl;
  }

  // write edges only in dimension 2
  if (dimension == 2)
  {
    ofile << endl;
    ofile << "Edges" << endl;
    ofile << edge.size() << endl;
    for (unsigned int i = 0; i < edge.size(); i++)
    {
      ofile << edge[i].nodes[0] << " " <<  edge[i].nodes[1]
        << " " << edge[i].reference << endl;
    }
  }

  if (dimension != 3 || includeInteriorFaces)
  {
    // write triangles
    ofile << endl;
    ofile << "Triangles" << endl;
    ofile <<  triangle.size() << endl;
    for (auto& tri: triangle)
    {
      for (unsigned int k = 0; k < 3; k++)
      {
        ofile << tri.nodes[k] << " ";
      }

      ofile << tri.reference << endl;
    }

    // write quadrilaterals
    ofile << endl;
    ofile << "Quadrilaterals" << endl;
    ofile <<  quad.size() << endl;
    for (auto& q: quad)
    {
      for (unsigned int k = 0; k < 4; k++) 
      {
        ofile << q.nodes[k] << " ";
      }

      ofile << q.reference << endl;
    }
  }
  else
  {
    unsigned int num_bd_triangles = 0;
    unsigned int num_bd_quads = 0;

    for (auto& tri: triangle)
    {
      if (tri.reference)
      {
        ++num_bd_triangles;
      }
    }

    for (auto& q: quad)
    {
      if (q.reference)
      {
        ++num_bd_quads;
      }
    }

    // write triangles
    ofile << endl;
    ofile << "Triangles" << endl;
    ofile <<  num_bd_triangles << endl;
    for (auto& tri: triangle)
    {
      if (!tri.reference)
      {
        continue;
      }

      for (unsigned int k = 0; k < 3; k++)
      {
        ofile << tri.nodes[k] << " ";
      }

      ofile << tri.reference << endl;
    }

    // write quadrilaterals
    ofile << endl;
    ofile << "Quadrilaterals" << endl;
    ofile <<  num_bd_quads << endl;
    for (auto& q: quad)
    {
      if (!q.reference)
      {
        continue;
      }

      for (unsigned int k = 0; k < 4; k++) 
      {
        ofile << q.nodes[k] << " ";
      }

      ofile << q.reference << endl;
    }
  }

  if (dimension == 3)
  {
    // write tetra
    ofile << endl;
    ofile << "Tetrahedra" << endl;
    ofile << tetra.size() << endl;
    for (auto& cell: tetra)
    {
      for (unsigned int k = 0; k < 4; k++)
      {
        ofile << cell.nodes[k] << " ";
      }

      ofile << cell.reference << endl;
    }

    // write hexa
    ofile << endl;
    ofile << "Hexahedra" << endl;
    ofile << hexa.size() << endl;
    for (auto& cell: hexa)
    {
      for (unsigned int k = 0; k < 8; k++)
      {
        ofile << cell.nodes[k] << " ";
      }

      ofile << cell.reference << endl;
    }
  }
}

// write the mesh onto a .(x)GEO file
// note: in 2D it needs the boundary description. If the boundary
// has not been initialized, the function returns an error
void Mesh::writeToGEO(const std::string& geoFilename)
{
  bool writeXgeo = false;
  std::ofstream geofile;
  geofile.open(geoFilename.c_str(),std::ios::out);
  if (!geofile)
  {
    ErrThrow(" *** Error(Mesh::writeToGEO) I could not open ", geoFilename);
  }

  // check errors
  if ((dimension == 2) && (boundary.parts.size() == 0))
  {
    ErrThrow(" *** Error(Mesh::writeToGEO) I need a boundary description (PRM "
             "file) in order to write a 2D GEO file");
  }

  // write output (x)GEO file
  Output::print(" writing mesh on ", geoFilename);

  // check if input file is an extended geo file (.xGEO)
  unsigned int nn = 0;
  while (geoFilename[nn] != 0)
  {
    ++nn;
  }

  if (geoFilename[nn - 4] == 'x')
  {
    Output::print("  Mesh::writeToGEO: writing mesh in xGEO format (with physical references) ");
    writeXgeo = true;
  }

  unsigned int numberOfElements = 0, maxNVertexPerElem = 0;

  if (triangle.size())
  {
    numberOfElements = triangle.size();
    maxNVertexPerElem = 3;
  }

  if (quad.size())
  {
    numberOfElements += quad.size();
    maxNVertexPerElem = 4;
  }

  if ((quad.size() == 0) && (triangle.size() == 0))
  {
    ErrThrow(" Mesh:Write2GEO Error: I cannot write a mesh without elements");
  }

  unsigned int nBoundParts = boundary.parts.size();

  geofile << "Grid 2D #generated by Mesh::writeToGEO " << endl;
  geofile << "n_elem  n_vert ignored  max_n_vert_per_elem  n_bound_parts" << endl;
  geofile <<  numberOfElements << " " <<  vertex.size() << " "
	  <<  " 0 " << maxNVertexPerElem  << " " << nBoundParts << " " << endl;
  geofile << "DCORVG" << endl;

  // vertices
  /**
     @note in 2D, The vertex belonging to the boundary are written
     according to their local parametrization. I.e., we write
     compID+t 0
     where compID is the ID of the boundary component of the vertex, and t is its 
     local parametrization.
     Inner vertices are written simply using their coordinates.
   */
  double localParam = -1;
  int partID;
  for (unsigned int i = 0; i < vertex.size(); i++)
  {
    // check if the vertex is on the boundary
    partID = boundary.isOnComponent(vertex[i].x,vertex[i].y,localParam);

    if (partID >= 0)
    {
      Output::print("Vertex ", vertex[i].x, " ", vertex[i].y,
        " is on part ", partID + 1,
        "; t = ", localParam, " --> ", localParam, " 0.0");

      geofile << localParam << " 0.000000000 " << endl;

      // note: a boundary vertex takes the reference of its boundary component+1
      vertex[i].reference = partID + 1;
    }
    else
    {
      geofile << vertex[i].x << " " << vertex[i].y << endl;

      // inner nodes are marked with reference 0
      vertex[i].reference = 0;
    }
  }

  geofile << "KVERT" << endl;

  // elements
  for (unsigned int i = 0; i < triangle.size(); i++)
  {
    for (unsigned int k = 0; k < 3; k++)
    {
      geofile << triangle[i].nodes[k] << " ";
    }

    // note: for mixed meshes, we weite a 0 as triangles 4th elements
    if (maxNVertexPerElem == 4)
    {
      geofile << 0 << " ";
    }

    if (writeXgeo)
    {
      geofile << triangle[i].reference;
    }

    geofile << endl;
  }

  for (unsigned int i = 0; i < quad.size(); i++)
  {
    for (unsigned int k = 0; k < 4; k++)
    {
      geofile << quad[i].nodes[k] << " ";
    }

    if (writeXgeo)
    {
      geofile << quad[i].reference;
    }

    geofile << endl;
  }

  geofile << "KNPR" << endl;

  for (unsigned int i = 0; i < vertex.size(); i++)
  {
    geofile << vertex[i].reference << " ";

    if (((i + 1) % 10) == 0)
    {
      geofile << endl;
    }
  }

  geofile << endl;
  geofile << "KMM";
  geofile << endl;

  ///@todo what does this last row mean?
  geofile << "1  2  3  4" << endl;
  geofile.close();
}

void Mesh::fixBoundaryElements()
{
  if (is_tetramesh)
  {
    fixBoundaryTetrahedra();
  }
  else
  {
    Output::root_warn("Mesh", "Boundary element correction is implemented for "
      "tetrahedral meshes only.");
  }
}

void Mesh::fixBoundaryTetrahedra()
{
  // identify and subdivide tetrahedra with more than one boundary face

  if (faceToTetra.size() == 0)
  {
    createFaceToTetrahedraMap();
  }

  Output::root_info("Mesh",
    "Identifying and fixing multi-boundary tetrahedra.");

  // initial counts
  size_t num_tetra = tetra.size();
  size_t num_triangles = triangle.size();
  size_t num_vert = vertex.size();

  int num_cells_split = 0;

  // note that the upper bound tetra.size() (rather than num_tetra) is sound,
  // assuming the mesh is topologically nondegenerate (in the sense that each
  // edge of a boundary face is adjacent to exactly two boundary faces):
  //
  // we say a tetrahedron is "bad", specifically "k-bad", if it has k > 1
  // boundary faces.
  //
  // - centroid split removes a bad tetrahedron with no side effects on other
  //   tetrahedra.
  // - an edge split can create new bad tetrahedra in three scenarios:
  //   - splitting a 2-bad tetrahedron along an edge of only one of its
  //     boundary faces - in this case only one of the two resulting tetrahedra
  //     will be bad, so if the edge is the back edge of another bad
  //     tetrahedron, the number of bad tetrahedra still decreases.
  //   - splitting a 3-bad tetrahedron along an edge of only one of its
  //     boundary faces. this may happen in particularly unfortunate
  //     geometries. in this case the 3-bad tetrahedron is replaced by two
  //     2-bad tetrahedra, increasing the total number of bad tetrahedra.
  //     but the number of 3-bad tetrahedra decreases, so this can happen only
  //     a finite number of times.
  //   - splitting any bad tetrahedron along an edge shared by two of its
  //     boundary faces. this case is not relevant because such an edge cannot
  //     be the back edge of a 2-bad tetrahedron without a topological
  //     degeneracy.
  //
  // consequently bad tetrahedra may appear past the end of the original list,
  // but the below loop cannot create infinite additional tetrahedra.
  //
  // it may create very many in the case of especially unhappy geometries where
  // edges are shared by extremely many tetrahedra, but if n_k is the number of
  // k-bad original tetrahedra and m is the maximum number of tetrahedra
  // adjacent to any edge, then the final number of tetrahedra is
  // <= num_tetra + 2 n_4 + m * (n_2 + 2 n_3)
  // < (2 + 3 m) num_tetra.

  for (size_t c = 0; c < tetra.size(); c++)
  {
    auto& cell = tetra[c];

    // first, count boundary faces.

    int numBoundaryFaces = 0;
    int f0 = -1;
    int f1 = -1;

    for (int i = 0; i < 4; i++)
    {
      // no need to worry about face vertex order here,
      // so use [i, i+1, i+2] mod 4
      int f = findTriFace(cell.nodes[i],
        cell.nodes[(i + 1) & 3],
        cell.nodes[(i + 2) & 3]);

      if (faceToTetra[f][0] < 0 || faceToTetra[f][1] < 0)
      {
        // this is a boundary face

        if (numBoundaryFaces == 0)
        {
          f0 = i;
        }
        else
        {
          f1 = i;
        }

        ++numBoundaryFaces;
      }
    }

    if (numBoundaryFaces < 2)
    {
      // nothing to see here, move along
      continue;
    }
    else if (numBoundaryFaces == 2)
    {
      // find the back edge by finding the vertex
      // opposite each boundary face

      int a = (f0 + 3) & 3;
      int b = (f1 + 3) & 3;

      Output::root_info<4>("Mesh", "Cell ", c, " has 2 boundary faces! "
        "Splitting the back edge ",
        cell.nodes[a], ", ", cell.nodes[b], ".");

      num_cells_split += splitTetrahedraAtEdge(cell.nodes[a], cell.nodes[b]);
    }
    else
    {
      // more than two boundary faces!

      Output::root_info<4>("Mesh", "Cell ", c, " has ", numBoundaryFaces,
      " boundary faces! Splitting it up at the centroid.");

      splitTetrahedronCentroid((int)c);
      ++num_cells_split;
    }
  }

  if (num_cells_split > 0)
  {
    Output::root_info("Mesh", "Split ", num_cells_split, " cells. The mesh now has ",
      tetra.size(), " cells, ",
      triangle.size(), " triangles, ",
      vertex.size(), " vertices."
      "\nPreviously: ",
      num_tetra, " cells, ",
      num_triangles, " triangles, ",
      num_vert, " vertices.");
  }
}

void Mesh::buildEdgeSums()
{
  edgeSumToTriangle.clear();

  std::unordered_map<int, int> bucketCount;

  for (size_t f = 0; f < triangle.size(); f++)
  {
    auto& face = triangle[f];
    for (int i = 0; i < 3; i++)
    {
      int s = face.nodes[i] + face.nodes[(i + 1) % 3];

      if (bucketCount.count(s) > 0)
      {
        bucketCount[s]++;
      }
      else
      {
        bucketCount[s] = 1;
      }
    }
  }

  for (auto& p: bucketCount)
  {
    edgeSumToTriangle.reserve(p.first, p.second);
  }

  for (size_t f = 0; f < triangle.size(); f++)
  {
    auto& face = triangle[f];

    for (int i = 0; i < 3; i++)
    {
      // we will only encounter each edge sum once
      // per nondegenerate triangle:
      // WLOG a < b < c. then b + c > a + c > a + b.

      int s = face.nodes[i] + face.nodes[(i + 1) % 3];

      edgeSumToTriangle.add(f, s);
    }
  }
}

void Mesh::removeEdgeSums(int f)
{
  auto& tri = triangle[f];

  for (int i = 0; i < 3; i++)
  {
    int s = tri.nodes[i] + tri.nodes[(i + 1) % 3];

    edgeSumToTriangle.remove(f, s);
  }
}

void Mesh::updateEdgeSums(int f)
{
  auto& tri = triangle[f];

  for (int i = 0; i < 3; i++)
  {
    int s = tri.nodes[i] + tri.nodes[(i + 1) % 3];

    edgeSumToTriangle.add(f, s);
  }
}

int Mesh::splitTetrahedraAtEdge(int a, int b)
{
  if (edgeSumToTriangle.size() == 0)
  {
    // initialize the edge mapping vector
    buildEdgeSums();
  }

  // create edge midpoint vertex

  int new_node_index = 1 + (int)vertex.size();
  vertex.resize(new_node_index);

  auto& va = vertex[a - 1];
  auto& vb = vertex[b - 1];

  meshNode &v = vertex[new_node_index - 1];

  v.x = 0.5 * (va.x + vb.x);
  v.y = 0.5 * (va.y + vb.y);
  v.z = 0.5 * (va.z + vb.z);
  v.reference = std::max(va.reference, vb.reference);

  // cells affected by splitting each triangle
  std::set<int> cellsToSplit;

  // copy! if we simply iterated over edgeSumToTriangle here,
  // we might run into trouble due to changes while splitting
  const std::vector<int> affectedTriangles = edgeSumToTriangle.at(a + b);

  // split all triangles sharing this edge
  for (int f: affectedTriangles)
  {
    auto& face = triangle[f];

    for (int i = 0; i < 3; i++)
    {
      int j = (i + 1) % 3;

      if ((face.nodes[i] == a && face.nodes[j] == b)
          || (face.nodes[i] == b && face.nodes[j] == a))
      {
        splitTriangleAtEdge(f, i, new_node_index, cellsToSplit);
        break;
      }
    }
  }

  // now split the relevant cells
  for (int c: cellsToSplit)
  {
    splitTetrahedronAtEdge(c, a, b, new_node_index);
  }

  return cellsToSplit.size();
}

void Mesh::splitTriangleAtEdge(int f, int i, int new_node_index, std::set<int> &cellsToSplit)
{
  std::array<int, 3> original_nodes = triangle[f].nodes;
  int ref = triangle[f].reference;

  int a = original_nodes[i];
  int b = original_nodes[(i + 1) % 3];
  int c = original_nodes[(i + 2) % 3];

  // remove old edge sum data
  removeEdgeSums(f);

  // add a new face
  int new_face_index = (int)triangle.size();
  triangle.resize(new_face_index + 1);
  faceToTetra.emplace(faceToTetra.end(), 2, -1);

  auto& tri_0 = triangle[f];
  auto& tri_1 = triangle[new_face_index];

  // orientation-preserving node assignment

  tri_0.nodes[0] = a;
  tri_0.nodes[1] = new_node_index;
  tri_0.nodes[2] = c;
  tri_0.reference = ref;

  tri_1.nodes[0] = b;
  tri_1.nodes[1] = c;
  tri_1.nodes[2] = new_node_index;
  tri_1.reference = ref;

  // rehash old face
  unhashTriangle(f, a + b + c);
  hashNewTriangle(f);

  // hash new face
  hashNewTriangle(new_face_index);

  // write new edge sum data
  updateEdgeSums(f);
  updateEdgeSums(new_face_index);

  // queue adjacent cells for splitting

  int c0 = faceToTetra[f][0];
  int c1 = faceToTetra[f][1];

  if (c0 >= 0)
  {
    cellsToSplit.insert(c0);
  }
  else
  {
    ++n_boundary_faces;
  }

  if (c1 >= 0)
  {
    cellsToSplit.insert(c1);
  }
  else
  {
    ++n_boundary_faces;
  }

  // reset face to cell mapping since we'll be refilling these
  // when splitting the tetrahedra
  faceToTetra[f][0] = -1;
  faceToTetra[f][1] = -1;
  faceToTetra[new_face_index][0] = -1;
  faceToTetra[new_face_index][1] = -1;
}

void Mesh::splitTetrahedronAtEdge(int c, int a, int b, int new_node_index)
{
  auto& cell = tetra[c];
  std::array<int, 4> original_nodes = cell.nodes;
  int ref = cell.reference;

  const int FaceVertex[][3] = TETRA_FACE_VERTEX;

  // allocate new cell

  int new_cell_index = (int)tetra.size();
  tetra.resize(new_cell_index + 1);

  auto& old_cell = tetra[c];
  auto& new_cell = tetra[new_cell_index];

  // set vertices

  // vertices of the edge opposite to ab
  int other_a = -1;
  int other_b = -1;

  for (int i = 0; i < 4; i++)
  {
    int v = original_nodes[i];

    if (v != a && v != b)
    {
      if (other_a < 0)
      {
        other_a = v;
      }
      else
      {
        other_b = v;
      }
    }

    // this preserves orientation
    old_cell.nodes[i] = v == b ? new_node_index : v;
    new_cell.nodes[i] = v == a ? new_node_index : v;
  }

  old_cell.reference = ref;
  new_cell.reference = ref;

  // create interior face (split vertex to opposite edge)

  int new_face_index = (int)triangle.size();

  triangle.resize(new_face_index + 1);
  faceToTetra.emplace(faceToTetra.end(), 2, -1);

  auto& tri = triangle[new_face_index];

  // same order the opposite vertices originally appeared in
  tri.nodes[0] = other_a;
  tri.nodes[1] = other_b;
  tri.nodes[2] = new_node_index;

  // interior faces get reference 0
  tri.reference = 0;

  // update lookup structures
  hashNewTriangle(new_face_index);
  updateEdgeSums(new_face_index);

  // set adjacency for the new cell at index c
  for (int i = 0; i < 4; i++)
  {
    int fa = old_cell.nodes[FaceVertex[i][0]];
    int fb = old_cell.nodes[FaceVertex[i][1]];
    int fc = old_cell.nodes[FaceVertex[i][2]];

    int f = findTriFace(fa, fb, fc);

    if (f < 0)
    {
      ErrThrow("Face doesn't exist.");
    }

    if (fa != new_node_index && fb != new_node_index && fc != new_node_index)
    {
      // as one of the original cell's exterior faces, this is already
      // marked adjacent to the original cell (index c), same index as
      // the first new cell
      continue;
    }
    else
    {
      // this is one of the new split faces or the interior face

      if (faceToTetra[f][0] < 0)
      {
        faceToTetra[f][0] = c;
      }
      else
      {
        faceToTetra[f][1] = c;
      }
    }
  }

  // set adjacency for the new cell at index new_cell_index
  for (int i = 0; i < 4; i++)
  {
    int fa = new_cell.nodes[FaceVertex[i][0]];
    int fb = new_cell.nodes[FaceVertex[i][1]];
    int fc = new_cell.nodes[FaceVertex[i][2]];

    int f = findTriFace(fa, fb, fc);

    if (f < 0)
    {
      ErrThrow("Face doesn't exist.");
    }

    if (fa != new_node_index && fb != new_node_index && fc != new_node_index)
    {
      // as one of the original cell's exterior faces, this was previously
      // marked adjacent to the original cell (index c)

      if (faceToTetra[f][0] == c)
      {
        faceToTetra[f][0] = new_cell_index;
      }
      else
      {
        faceToTetra[f][1] = new_cell_index;
      }
    }
    else
    {
      // this is one of the new split faces or the interior face

      if (faceToTetra[f][0] < 0)
      {
        faceToTetra[f][0] = new_cell_index;
      }
      else
      {
        faceToTetra[f][1] = new_cell_index;
      }
    }
  }
}

void Mesh::unhashTriangle(int f, int old_hash)
{
  meshTrifaceHash.remove(f, old_hash);
}

void Mesh::hashNewTriangle(int f)
{
  auto& tri = triangle[f];

  meshTrifaceHash.add(f, tri.nodes[0] + tri.nodes[1] + tri.nodes[2]);
}

void Mesh::splitTetrahedronCentroid(int c)
{
  if (edgeSumToTriangle.size() == 0)
  {
    buildEdgeSums();
  }

  auto& cell = tetra[c];
  int cell_ref = cell.reference;
  std::array<int, 4> original_nodes = cell.nodes;

  const int FaceVertex[][3] = TETRA_FACE_VERTEX;

  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
  int vertex_ref = vertex[original_nodes[0] - 1].reference;

  for (int i = 0; i < 4; i++)
  {
    // node indices are 1-based
    auto& node = vertex[original_nodes[i] - 1];
    x += node.x;
    y += node.y;
    z += node.z;
    vertex_ref = std::min(vertex_ref, node.reference);
  }

  x *= 0.25;
  y *= 0.25;
  z *= 0.25;

  // node indices are 1-based
  int new_node_index = 1 + (int)vertex.size();
  vertex.resize(new_node_index);

  int base_tetra_index = (int)tetra.size();
  tetra.resize(base_tetra_index + 3);

  meshNode &v = vertex[new_node_index - 1];

  v.x = x;
  v.y = y;
  v.z = z;
  v.reference = vertex_ref;

  for (int i = 0; i < 4; i++)
  {
    int ii = original_nodes[FaceVertex[i][0]];
    int jj = original_nodes[FaceVertex[i][1]];
    int kk = original_nodes[FaceVertex[i][2]];

    int new_cell_index = i > 0 ? base_tetra_index + i - 1 : (int)c;

    meshTetrahedron &subcell = tetra[new_cell_index];

    subcell.nodes[0] = ii;
    subcell.nodes[1] = jj;
    subcell.nodes[2] = kk;
    subcell.nodes[3] = new_node_index;
    subcell.reference = cell_ref;

    // update faceToTetra map for original exterior face
    // (not necessary for i == 0 since new_cell_index == c)
    if (i > 0)
    {
      int f = findTriFace(ii, jj, kk);

      if (faceToTetra[f][0] == (int)c)
      {
        faceToTetra[f][0] = new_cell_index;
      }
      else
      {
        faceToTetra[f][1] = new_cell_index;
      }
    }

    // update face references for new interior faces
    for (int r = 1; r < 4; r++)
    {
      int ar = subcell.nodes[FaceVertex[r][0]];
      int br = subcell.nodes[FaceVertex[r][1]];
      int cr = subcell.nodes[FaceVertex[r][2]];

      int face_r = findTriFace(ar, br, cr);

      if (face_r < 0)
      {
        // new face! create and insert

        face_r = (int)triangle.size();

        triangle.emplace(triangle.end());

        auto& tri = triangle[face_r];

        tri.nodes[0] = ar;
        tri.nodes[1] = br;
        tri.nodes[2] = cr;

        // interior triangles are ref zero
        tri.reference = 0;

        // add to map
        hashNewTriangle(face_r);
        updateEdgeSums(face_r);

        // we're the first neighbour
        faceToTetra.emplace(faceToTetra.end(), 2, -1);
        faceToTetra[face_r][0] = new_cell_index;
      }
      else
      {
        // this one already exists, so we're the second neighbour
        faceToTetra[face_r][1] = new_cell_index;
      }
    }
  }
}

void Mesh::info()
{
  Output::print(" ---- MESH INFORMATION ---- ");
  Output::print(" -- Dimension: ",dimension);
  Output::print(" --  N. Nodes: ",vertex.size());
  Output::print(" --  N. Edges: ",edge.size());
  Output::print(" --  N. Triangles: ",triangle.size());
  Output::print(" --  N. Quandrilaterals: ",quad.size());
  Output::print(" --  N. Tetrahedra: ",tetra.size());
  Output::print(" --  N. Hexahedra: ",hexa.size());
}

template <class Element2D>
void Mesh::correct_numbering_of_vertices_if_needed(std::vector<Element2D>& elements_2d)
{
  bool numbering_was_corrected = false;

  for (unsigned int i = 0; i < elements_2d.size(); i++)
  {
    if (first_three_vertices_are_numbered_clockwisely<Element2D>(elements_2d[i]))
    {
      reverse_clockwise_numbering<Element2D>(elements_2d[i]);
      numbering_was_corrected = true;
    }
  }

  if(numbering_was_corrected)
  {
    Output::root_info("Mesh", "Numbering of some of element vertices was modified"
                  " to get counter-clockwise numbering required by ParMooN.");
  }
}

template <class Element2D>
bool Mesh::first_three_vertices_are_numbered_clockwisely(const Element2D& element_2d)
{
  double vec1[2], vec2[2];

  int verties_indices[3] = { element_2d.nodes[0] - 1,
                             element_2d.nodes[1] - 1,
                             element_2d.nodes[2] - 1 };

  vec1[0] = vertex[verties_indices[1]].x - vertex[verties_indices[0]].x;
  vec1[1] = vertex[verties_indices[1]].y - vertex[verties_indices[0]].y;

  vec2[0] = vertex[verties_indices[2]].x - vertex[verties_indices[1]].x;
  vec2[1] = vertex[verties_indices[2]].y - vertex[verties_indices[1]].y;

  return get_z_component_of_cross_product(vec1, vec2) < 0.0;
}

template <class Element2D>
void Mesh::reverse_clockwise_numbering(Element2D& element_2d)
{
  int i = element_2d.nodes[0];
  element_2d.nodes[0] = element_2d.nodes[2];
  element_2d.nodes[2] = i;
}

void check_hexahedron(const meshHexahedron& hex,
                      const std::vector<meshNode>& vertices)
{
  constexpr int n_faces = 6;
  constexpr int n_vert_per_face = 4;

  std::array<std::array<int, n_vert_per_face>, n_faces> face_vertex = 
    { HEXA_FACE_VERTEX };

  for(int f = 0; f < n_faces; ++f)
  {
    std::array<int, n_vert_per_face> face_vertices =
    {
       hex.nodes[face_vertex[f][0]]-1, hex.nodes[face_vertex[f][1]]-1,
       hex.nodes[face_vertex[f][2]]-1, hex.nodes[face_vertex[f][3]]-1
    };

    // check if the four vertices really are on one plane
    double x1 = vertices[face_vertices[1]].x - vertices[face_vertices[0]].x;
    double y1 = vertices[face_vertices[1]].y - vertices[face_vertices[0]].y;
    double z1 = vertices[face_vertices[1]].z - vertices[face_vertices[0]].z;
    double x2 = vertices[face_vertices[2]].x - vertices[face_vertices[0]].x;
    double y2 = vertices[face_vertices[2]].y - vertices[face_vertices[0]].y;
    double z2 = vertices[face_vertices[2]].z - vertices[face_vertices[0]].z;
    double x3 = vertices[face_vertices[3]].x - vertices[face_vertices[0]].x;
    double y3 = vertices[face_vertices[3]].y - vertices[face_vertices[0]].y;
    double z3 = vertices[face_vertices[3]].z - vertices[face_vertices[0]].z;

    // compute determinant of matrix (Rule of Sarrus)
    // ( x1  x2  x3 )
    // ( y1  y2  y3 )
    // ( z1  z2  z3 )
    double det = x1 * y2 * z3 + x2 * y3 * z1 + x3 * y1 * z2
      - z1 * y2 * x3 - z2 * y3 * x1 - z3 * y1 * x2;

    if (std::abs(det) > 1.e-5)
    {
      Output::print("The vertices ", face_vertices[0], ", ", face_vertices[1], 
                    ", ", face_vertices[2], ", and ", face_vertices[3], 
                    " are not on a plane. Their coordinates are:\n(",
                    vertices[face_vertices[0]].x, ",",
                    vertices[face_vertices[0]].y, ",",
                    vertices[face_vertices[0]].z, "),\n(",
                    vertices[face_vertices[1]].x, ",",
                    vertices[face_vertices[1]].y, ",",
                    vertices[face_vertices[1]].z, "),\n(",
                    vertices[face_vertices[2]].x, ",",
                    vertices[face_vertices[2]].y, ",",
                    vertices[face_vertices[2]].z, "),   and\n(",
                    vertices[face_vertices[3]].x, ",",
                    vertices[face_vertices[3]].y, ",",
                    vertices[face_vertices[3]].z, ")");

      Output::print("the computed determinant is ", det,
                    " which should have been (close to) zero.");

      ErrThrow("There is a hexahedron which has a face whose vertices are not on a plane.");
    }
  }
}

void stripSpace(std::string &str)
{
  for (size_t i = 0; i < str.length(); i++)
  {
    if (str[i] == ' ') 
    {
      str.erase(i, 1);
      i--;
    }
  }
}

bool Mesh::IndexLookup::contains(int hash) const
{
  return map.count(hash) > 0;
}

const std::vector<int>& Mesh::IndexLookup::at(int hash) const
{
  if (map.count(hash) > 0)
  {
    // .at rather than [] since we're returning const
    return map.at(hash);
  }

  ErrThrow("Bad hash");
}

void Mesh::IndexLookup::reserve(int hash, int capacity)
{
  if (map.count(hash) > 0)
  {
    map[hash].reserve(capacity);
  }
  else
  {
    map.emplace(hash, 0);
    map[hash].reserve(capacity);
  }
}

void Mesh::IndexLookup::add(int value, int hash)
{
  if (map.count(hash) > 0)
  {
    map[hash].push_back(value);
  }
  else
  {
    map.emplace(hash, 1);
    map[hash][0] = value;
  }
}

void Mesh::IndexLookup::remove(int value, int hash)
{
  if (map.count(hash) > 0)
  {
    auto& vec = map[hash];

    for (size_t i = 0; i < vec.size(); i++)
    {
      if (vec[i] == value)
      {
        vec[i] = vec.back();
        vec.erase(vec.end() - 1);
        return;
      }
    }
  }
}

void Mesh::IndexLookup::clear()
{
  map.clear();
}