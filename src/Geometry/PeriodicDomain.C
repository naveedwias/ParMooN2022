#include <PeriodicDomain.h>

#include <BaseCell.h>
#include <BoundComp.h>
#include <JointEqN.h>
#include <MooNMD_Io.h>
#include <PeriodicJoint.h>
#include <Vertex.h>

#ifdef __2D__
  #include <BoundEdge.h>
#endif
#ifdef __3D__
  #include <BoundFace.h>
#endif
#ifdef _MPI
#include <MeshPartition.h>
#endif

#include <algorithm>
#include <cmath>
#include <numeric>


/* ************************************************************************** */
void normalized_vector(std::vector<double>& v, double tolerance)
{
  if(v.empty())
  {
    return;
  }

  double v_norm = std::sqrt(std::inner_product(v.begin(),
                                               v.end(),
                                               v.begin(), 0.));
  if(v_norm<tolerance)
  {
    ErrThrow("normalized_vector(): unable to normalize null vector");
  }
  for(auto& v_i: v)
  {
    v_i /= v_norm;
  }
}


/* ************************************************************************** */
bool are_vectors_collinear(const std::vector<double>& a,
                           const std::vector<double>& b,
                           double tol)
{
  if(a.size() != b.size())
  {
    return false;
  }
  else
  {
    // normalize vectors:
    double a_norm = std::sqrt(std::inner_product(a.begin(), a.end(),
                                                 a.begin(), 0.));
    double b_norm = std::sqrt(std::inner_product(b.begin(), b.end(),
                                                 b.begin(), 0.));
    double min_norm = std::min(a_norm, b_norm);

    if(min_norm<tol)
    {
      ErrThrow("are_vectors_collinear(): norm of vector is smaller than ",
               "tolerance:", tol);
    }

    double ab_norm = std::abs(std::inner_product(a.begin(), a.end(),
                                                 b.begin(), 0.));
    if(std::abs(ab_norm/(a_norm*b_norm) -1)>tol)
    {
      return false;
    }
    else
    {
      return true;
    }
  }
}


/* ************************************************************************** */
bool are_vectors_equal(const std::vector<double>& a,
                       const std::vector<double>& b,
                       double tol)
{
  if(a.size() != b.size())
  {
    return false;
  }
  else
  {
    for(size_t i=0 ; i<a.size() ; i++)
    {
      double abs_max = std::max(std::abs(a[i]), std::abs(b[i]));
      if(abs_max<=tol)
      {
        if(std::abs(a[i] - b[i])>tol)
        {
          return false;
        }
      }
      else
      {
        if((std::abs(a[i]-b[i])/abs_max)>tol)
        {
          return false;
        }
      }
    }
    return true;
  }
}


/* ************************************************************************** */
bool have_common_element(std::vector<int>& V1, std::vector<int>& V2)
{
  auto it1 = V1.begin();
  auto it2 = V2.begin();

  while(it1 != V1.end() && it2 != V2.end())
  {
    if(*it1 < *it2)
    {
      ++it1;
    }
    else if(*it2 < *it1)
    {
      ++it2;
    }
    else
    {
      return true;
    }
  }
  return false;
}


/* ************************************************************************** */
void sort_and_select(std::vector<int>& V)
{
  std::sort(V.begin(), V.end());
  auto new_end = std::unique(V.begin(), V.end());
  V.erase(new_end, V.end());
}


/* ************************************************************************** */
std::vector<double> set_periodicDirection(const ParameterDatabase& domain_db)
{
#ifdef __2D__
  int dim = 2;
#else
  int dim = 3;
#endif

  std::vector<double> periodic_directions;
  if(domain_db.contains("periodic_translations"))
  {
    std::vector<double> directions = domain_db["periodic_translations"];
    if(directions.size()==1 && directions[0]==0.)
    {
      Output::print("set_PeriodicJoint(): no periodic_translations defined, the"
                    " direction will be supposed normal to the boundary faces");
    }
    else if(directions.size()%dim != 0)
    {
      ErrThrow("set_PeriodicJoint(): Some coordinates in parameter: "
               "periodic_translations are missing.");
    }
    else
    {
      periodic_directions = directions;
    }
  }
  return periodic_directions;
}


/* ************************************************************************** */
#ifdef __2D__
bool find_PeriodicJoint(BoundCondFunct2D* BoundaryCondition,
#else
bool find_PeriodicJoint(BoundCondFunct3D* BoundaryCondition,
#endif
                        TDomain *Domain,
                        std::vector<std::pair<int, int>>& joint_list,
                        std::vector<parmoon::Point>& barycenter)
{
  unsigned int nb_periodic = 0;
  const TCollection* coll = Domain->GetCollection(It_Finest, 0);
  const int n_cells = coll->GetN_Cells();

  for(int i=0 ; i<n_cells ; i++)
  {
    auto cell = coll->GetCell(i);
#ifdef __2D__
    int n_joints = cell->GetN_Edges();
#else
    int n_joints = cell->GetN_Faces();
#endif
    for(int j=0 ; j<n_joints ; j++)
    {
      auto joint = cell->GetJoint(j);

#ifdef __2D__
      if(joint->GetType() == BoundaryEdge)
      {
        auto boundjoint = static_cast<TBoundEdge*>(joint);
#else
      if(joint->GetType() == BoundaryFace)
      {
        auto boundjoint = static_cast<TBoundFace*>(joint);
#endif
        const TBoundComp* BoundComp = boundjoint->GetBoundComp();
        BoundCond bd_cond;
        // compute joint center:
#ifdef __2D__
        int comp = BoundComp->GetID();
        double t0, t1, Xe, Ye;
        boundjoint->GetParameters(t0, t1);
        double T = (t0+t1)/2.;
        boundjoint->GetXYofT(T, Xe, Ye);
        BoundaryCondition(comp, T, bd_cond);
        parmoon::Point p_center(Xe, Ye);
#else
        int comp                    = BoundComp->get_physical_id();
        const int *TmpFV, *TmpLen;
        int MaxLen;
        cell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);
        double t0 = 1.0/TmpLen[j];
        double xf = 0; double yf = 0; double zf = 0;
        double X, Y, Z;
        for(int iv=0 ; iv<TmpLen[j] ; iv++)
        {
          cell->GetVertex(TmpFV[j*MaxLen+iv])->GetCoords(X, Y, Z);
          xf += t0*X;
          yf += t0*Y;
          zf += t0*Z;
        }
        BoundaryCondition(comp, xf, yf, zf, bd_cond);
        parmoon::Point p_center(xf, yf, zf);
#endif
        if(bd_cond == PERIODIC)
        {
          joint_list.push_back(std::make_pair(i, j));
          barycenter.push_back(p_center);
          nb_periodic++;
        }
      }
    }
  }
  if(nb_periodic != 0)
  {
    Output::print("Found ", nb_periodic, " periodic joints");
#ifdef _MPI
    Domain->SetPeriodicDomain();
#endif
    return true;
  }
  else
  {
    return false;
  }
}


/* ************************************************************************** */
#ifdef __2D__
void set_PeriodicJoint(BoundCondFunct2D* BoundaryCondition,
#else
void set_PeriodicJoint(BoundCondFunct3D* BoundaryCondition,
#endif
                       TDomain *Domain)
{
  std::vector<std::pair<int, int>> joint_list;
  std::vector<parmoon::Point> barycenter;

  if(!find_PeriodicJoint(BoundaryCondition, Domain, joint_list, barycenter))
  {
    return;
  }
  else
  {
    std::vector<double> periodic_directions = set_periodicDirection(
                                                        Domain->get_database());
    int nb_pair = joint_list.size()/2;
    auto coll = Domain->GetCollection(It_Finest, 0);

    Output::print("Setting periodic joints");

#ifdef _MPI
    // used VertexNumbers to constructed PeriodicVertexConnect
    int N_Cells = coll->GetN_Cells();
    // this assumes that all cells have the same geometry:
    TBaseCell* first_cell = coll->GetCell(0);
    int N_VertInCell = first_cell->GetN_Vertices();
    int N_AllLocVert = N_VertInCell * N_Cells;
    std::vector<int> VertexNumbers(N_AllLocVert);
    std::vector<const TVertex*> Vertices = get_sorted_vertices(coll);
    int N_TotalVertices = get_number_vertex(Vertices, coll, VertexNumbers);
    std::vector<std::vector<int>>* PeriodicVertexConnect = Domain->
                                                get_PeriodicVertexConnect_ptr();
    PeriodicVertexConnect->resize(N_TotalVertices);
#endif

    for(int idx_1 = 0; idx_1 < nb_pair; idx_1++)
    {
      bool is_found = false;
      int idx_c1 = joint_list[idx_1].first;
      auto cell_1 = coll->GetCell(idx_c1);
      int idx_j1 = joint_list[idx_1].second;
      auto joint_1 = cell_1->GetJoint(idx_j1);
      parmoon::Point b_1 = barycenter[idx_1];

#ifdef __2D__
      const unsigned int dim = 2;
      auto boundjoint = static_cast<TBoundEdge*>(joint_1);
      std::vector<double> n_1(2, 0.);
      boundjoint->get_normal(n_1[0], n_1[1]);
      std::vector<double> n_2(2, 0.);
#else
      const unsigned int dim = 3;
      auto boundjoint = static_cast<TBoundFace*>(joint_1);
      std::vector<double> n_1(3, 0.);
      boundjoint->get_normal_vector(b_1.x, b_1.y, b_1.z,
                                    n_1[0], n_1[1], n_1[2]);
      std::vector<double> n_2(3, 0.);
#endif

      int idx_c2, idx_j2;
      parmoon::Point v_12(dim);
      TBaseCell* cell_2;
      TJoint* joint_2;

      // find the joint idx_2 periodic with idx_1
      for(unsigned int idx_2 = idx_1+1; idx_2 < joint_list.size(); idx_2++)
      {
        idx_c2 = joint_list[idx_2].first;
        cell_2 = coll->GetCell(idx_c2);
        idx_j2 = joint_list[idx_2].second;
        joint_2 = cell_2->GetJoint(idx_j2);
        parmoon::Point b_2 = barycenter[idx_2];
#ifdef __2D__
        auto boundjoint = static_cast<TBoundEdge*>(joint_2);
        boundjoint->get_normal(n_2[0], n_2[1]);
#else
        auto boundjoint = static_cast<TBoundFace*>(joint_2);
        boundjoint->get_normal_vector(b_2.x, b_2.y, b_2.z,
                                      n_2[0], n_2[1], n_2[2]);
#endif

        // test if joint idx_2 is the translated from idx_1 in
        // the periodic_direction or (else) normal direction
        if(are_vectors_collinear(n_1, n_2))
        {
          std::vector<double> b_12 = static_cast<parmoon::Point>(b_2-b_1);

          if(!periodic_directions.empty())
          {
            for(unsigned int i_d=0 ; i_d<periodic_directions.size()/dim ; i_d++)
            {
              std::vector<double> dir(dim, 0.);
              std::vector<double> dir_opp(dim, 0.);

              for(unsigned int v_i = 0; v_i < dim; v_i++)
              {
                dir[v_i] = periodic_directions[i_d*dim+v_i];
                dir_opp[v_i] = -periodic_directions[i_d*dim+v_i];
              }
              if(are_vectors_equal(dir, b_12)||are_vectors_equal(dir_opp, b_12))
              {
                is_found = true;
                v_12 = b_2 - b_1;

                joint_list.erase(joint_list.begin() + idx_2);
                barycenter.erase(barycenter.begin() + idx_2);
                break;
              }
            }
            if(is_found)
            {
              break;
            }
          }
          else
          {
            normalized_vector(b_12);
            if(are_vectors_collinear(n_1, b_12))
            {
              is_found = true;
              v_12 = b_2 - b_1;

              joint_list.erase(joint_list.begin() + idx_2);
              barycenter.erase(barycenter.begin() + idx_2);
              break;
            }
          }
        }
      } // for(unsigned int idx_2...

      if(is_found)
      {
        is_found = false;
        // delete old joints
        delete joint_1;
        delete joint_2;
        // make new joint
        auto joint = new TPeriodicJoint(cell_1, cell_2);
        // set joint
        cell_1->SetJoint(idx_j1, joint);
        cell_2->SetJoint(idx_j2, joint);

#ifdef __3D__
        // this assumes that all cells have the same geometry:
        const int *TmpFV, *TmpLen;
        int MaxLen;
        cell_1->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);
        parmoon::Point P0_1(dim);
        cell_1->GetVertex(TmpFV[idx_j1*MaxLen])->GetCoords(P0_1.x,
                                                           P0_1.y,
                                                           P0_1.z);

        // find opposite vertex P0_2 to local vertex zero P0_1 of face idx_j1
        parmoon::Point P0_2 = P0_1 + v_12;
        int maptype = -1;
        for(int i=0 ; i<TmpLen[idx_j2] ; i++)
        {
          parmoon::Point P_2(dim);
          cell_2->GetVertex(TmpFV[idx_j2*MaxLen+i])->GetCoords(P_2.x,
                                                               P_2.y,
                                                               P_2.z);
          if(P_2.is_equal(P0_2))
          {
            is_found = true;
            maptype = i;

#ifdef _MPI
            // fill PeriodicVertexConnect
            // check the orientation of the faces idx_j1 and idx_j2
            // (if the normal vectors n_1 and n_2 have the same direction)
            int coeff_orientation = std::inner_product(n_1.begin(), n_1.end(),
                                    n_2.begin(), 0.) < 0. ? -1 : 1;

            // each vertex of face idx_j1 is connected to the corresponding
            // vertex of face idx_j2
            for(int j=0 ; j<MaxLen ; j++)
            {
              // index of the corresponding vertex according to face orientation
              int k = (i + (MaxLen + coeff_orientation*j)) % MaxLen;

              // Get number of vertices in VertexNumbers
              auto Vertex_nb_1 = idx_c1*N_VertInCell + TmpFV[idx_j1*MaxLen + j];
              auto Vertex_nb_2 = idx_c2*N_VertInCell + TmpFV[idx_j2*MaxLen + k];
              // get vertices global index
              int idx_P1 = VertexNumbers[Vertex_nb_1];
              int idx_P2 = VertexNumbers[Vertex_nb_2];

              int idx_min = std::min(idx_P1, idx_P2);
              int idx_max = std::max(idx_P1, idx_P2);

              PeriodicVertexConnect->at(idx_min).push_back(idx_min);
              PeriodicVertexConnect->at(idx_min).push_back(idx_max);
            }
#endif // #ifdef _MPI

            break;
          } //  if(P_2.is_equal(P0_2))
        } // for(int i=0 ; i<TmpLen[idx_j2] ; i++)
        if(!is_found)
        {
          ErrThrow("set_PeriodicJoint(): check if the tolerance in is_equal() "
                   "is not too small");
        }
        else
        {
          joint->SetMapType(maptype);
        }
#endif // #ifdef __3D__
      }/* if(is_found) */
      else
      {
        ErrThrow("set_PeriodicJoint(): mesh on periodic face seems not to be "
                 "identical or the periodic_directions parameter is "
                 "mismatching the geometry");
      }
    } // for(int idx_1 = 0; idx_1 < nb_pair; idx_1++)

#ifdef _MPI
    // format PeriodicVertexConnect
    // 1. removing empty vectors
    auto new_end = std::remove_if(PeriodicVertexConnect->begin(),
                                  PeriodicVertexConnect->end(),
                    [](const std::vector<int>& periodicVertices)
                    { return periodicVertices.size()<1; });
    PeriodicVertexConnect->erase(new_end, PeriodicVertexConnect->end());

    // 2. sort elements in each vector and keep only one of each
    for(size_t i=0 ; i<PeriodicVertexConnect->size() ; i++)
    {
      sort_and_select(PeriodicVertexConnect->at(i));
    }

    // 3. propagate vertex connectiviy (transitivity)
    for(auto it_i=PeriodicVertexConnect->begin() ;
        it_i!=PeriodicVertexConnect->end() ; it_i++)
    {
      auto it_j = it_i+1;
      while(it_j!=PeriodicVertexConnect->end())
      {
        if(have_common_element(*it_i, *it_j))
        {
          (*it_i).insert((*it_i).end(), (*it_j).begin(), (*it_j).end());
          sort_and_select(*it_i);
          it_j = PeriodicVertexConnect->erase(it_j);
        }
        else
        {
          it_j++;
        }
      }
    }

    // 4. set Vertices as periodic and PeriodicVertIndex
    TVertex::NbPeriodicVert = PeriodicVertexConnect->size();
    for(int i=0 ; i<TVertex::NbPeriodicVert ; i++)
    {
      for(auto j:PeriodicVertexConnect->at(i))
      {
        for(int i_c=0 ; i_c<N_Cells ; i_c++)
        {
          for(int j_v=0 ; j_v<N_VertInCell ; j_v++)
          {
            if(VertexNumbers[i_c*N_VertInCell + j_v] == j)
            {
              TVertex* current_vertex = coll->GetCell(i_c)->GetVertex(j_v);
              current_vertex->SetPeriodicVertIndex(i);
              break;
            }
          }
        }
      }
    }
#endif // #ifdef _MPI
  }
}
