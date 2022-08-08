/** ************************************************************************ 
*
* @class     TGridCell
* @brief     represent geometric information of the cell
* @author    Volker Behns  09.07.97
* @History 
 ************************************************************************  */

#ifndef __GRIDCELL__
#define __GRIDCELL__

#include <BaseCell.h>

/**  @brief represent geometric information of the cell */
class TGridCell : public TBaseCell
{
  protected:
    /**  @brief field of pointer to children */
    TBaseCell **Children;
    /**  @brief pointer to parent cell */
    TBaseCell *Parent;

    /**  @brief field of all vertices */
    TVertex  **Vertices;

    /**  @brief grid level on with this cell was generated */
    int RefLevel;

    virtual TBaseCell *GetChild(int C_i) override;
    virtual TBaseCell *GetParent() override;
    
    friend class TDomain;

    /**  @brief used by IsLineCutingCell() to test/apply the constraints (i.e.
     * for x in [range_m, range_M] then range_m <= A*x + b <= range_M)
     * and restrict the interval [range_m, range_M] if necessary */
    bool test_constraint(double  A,
                         double  B, 
                         double& range_m,
                         double& range_M) const;


  public:
    // Constructor
    TGridCell(const TRefDesc *refdesc, int reflevel);

    // Destructor
    ~TGridCell();

    // Methods
    /**  @brief set the pointer to vertex with number i */
    virtual int SetVertex(int Vert_i, TVertex *Vert) override;
    /**  @brief return the pointer to vertex with number i */
    virtual TVertex *GetVertex(int Vert_i) override;
    virtual const TVertex *GetVertex(int Vert_i) const override;
    /**  @brief return field of pointers to all vertices */
    const TVertex * const *GetVertices() const
    { return Vertices; }


    /** @brief computes the centroid of 2D cell
     * @param[in,out] x x-coordinate of centroid
     * @param[in,out] y y-coordinate of centroid
     * @return centroid of 2D cell
     * @note This is only correct in the case that the cell is a polygon, i.e
     * for elements with curved edges the function returns a wrong value.
     */
    virtual void ComputeCentroid(double& x, double& y) const override;
    
    /**  @brief Implemented only for polygonal cells in 2D
     * and tetrahedral cells in 3D */
    virtual parmoon::Point getCentroid() const override;

    /** @brief computes the reflection of a given point along a given edge
     *
     * @param[in] edge_j local edge number of the edge that is part of the
     * straight where the point is reflected at
     * @param[in] x_orig x coordinate of point that shall be reflected
     * @param[in] y_orig y coordinate of point that shall be reflected
     * @param[in] x_reflected x coordinate of the reflected point
     * @param[in] y_reflected y coordinate of the reflected point
     * @return the reflected point
     * @note Note that this method assumes the edge to be a part of a straight
     * line, i.e. for curved boundaries this method does not compute the correct
     * reflection.
     */
    virtual void ComputeReflectedPoint(int edge_j, double x_orig, double y_orig,
        double& x_reflected, double& y_reflected) const override;

    /**  @brief return number of children */
    virtual int GetN_Children() const override;
    /**  @brief return number of parents */
    virtual int GetN_Parents() const override;

    /**  @brief return pointer to child cell with number C\_i */
    virtual const TBaseCell *GetChild(int C_i) const override;
    /**  @brief return pointer to neighbor cell across joint with number
     * joint\_nr
     * */
    virtual const TBaseCell *GetNeighbor(int joint_nr) const override;
    /**  @brief return pointer to parent cell */
    virtual const TBaseCell *GetParent() const override;
    /**  @brief set parent */
    virtual int SetParent(TBaseCell *parent) override;
    /**  @brief return local number of child Me */
    virtual int GetChildNumber(TBaseCell *Me) const override;

    /**  @brief put out postscript data to a file */
    virtual void PS(std::ofstream &dat, double scale, double StartX,
                    double StartY, int cell_index = -1) const override;
    
    /**  @brief derefine the cell */
    virtual int Derefine() override;
    /**  @brief refine or derefine the cell according to cell's clipboard */
    virtual int RefDeref() override;
    /**  @brief set marks in neighbour cells in order to maintain 1-regularity */
    virtual int Gen1RegMarks() override;
    /**  @brief generate conforming closures */
    virtual int MakeConfClosure() override;
    /**  @brief refine a cell */
    virtual int Refine(int RefLevel) override;

    /**  @brief generate a 2-regular grid */
    virtual int Gen1RegGrid() override;
    /**  @brief set refinement for the neighbour of your parent on joint LocJointNum */
    virtual int Ref1Reg(int LocJointNum, TBaseCell *&RefCell) override;
    /**  @brief check whether the surroundings of cell is 1-regular */
    virtual int Check1Reg() override;
    /**  @brief set RefDesc to no refinement */
    virtual int SetNoRefinement() override;
    /**  @brief set RefDesc to regular refinement */
    virtual int SetRegRefine() override;
    /**  @brief set RefDesc to adaptive refinement */
    virtual int Set1Refine(int i) override;      
    /**  @brief check whether a cell should be refined */
    virtual int IsToRefine() const override;
    /**  @brief check whether exist some children */
    virtual int ExistChildren() const override
    { return Children == nullptr ? false : true; }

    /**  @brief check if the line define by the position P and the direction is
     * intersecting this cell.
     * If so, lmin and lmax are the intersecting points */
    virtual bool IsLineCutingCell(int                  direction,
#ifdef __2D__
                                  std::array<double,2> P,
#else // __3D__
                                  std::array<double,3> P,
#endif
                                  double&              lmin,
                                  double&              lmax) const override;

#ifdef __2D__
    /**  @brief return coordinates of mid point P\_j on edge J\_i */
    virtual int LineMidXY(int J_i, int P_j, double &X, double &Y) override;
    /**  @brief return parameters on boundary of subedge SJ\_j on edge J\_i */
    virtual int LineMidT(int J_i, int SJ_j, double &T_0, double &T_1) override;
#endif

    /** @brief computes the midpoint of an edge of a 2D cell
     * @param[in] joint_i local joint number of the edge
     * @param[out] x x-coordinate of midpoint
     * @param[out] y y-coordinate of midpoint
     * @return the midpoint of the edge
     * @note This is only correct in the case that the cell is a polygon, i.e
     * for elements with curved edges the function returns a wrong value.
     */
    virtual parmoon::Point ComputeMidOfJoint(int joint_i) const override;

    /**  @brief return whether a point is inside a cell */
    virtual bool PointInCell(parmoon::Point p) const override;
    
    /**  @brief get diameter of a cell */
    virtual double GetDiameter() const override
    { return RefDesc->GetShapeDesc()->GetDiameter(Vertices); }

    /**  @brief get shortest edge of a cell */
    virtual double GetShortestEdge() const override
    { return RefDesc->GetShapeDesc()->GetShortestEdge(Vertices); }

    /**  @brief return the length of the cell defined with the reference map */
    virtual double GetLengthWithReferenceMap() const override
    { return RefDesc->GetShapeDesc()->GetLengthWithReferenceMap(Vertices); }

     /**  @brief get measure of a cell */
    virtual double GetMeasure() const override
    { return RefDesc->GetShapeDesc()->GetMeasure(Vertices); }

    /**  @brief get geometry level */
    virtual int GetGeoLevel() override;

    /**  @brief return subgrid ID */
    virtual int GetSubGridID() const override;

    /**  @brief compute number of edges at the boundary */
    virtual int get_n_boundary_joints() const;
    
    virtual void check() const override;
};

#endif
