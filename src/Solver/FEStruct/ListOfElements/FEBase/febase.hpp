/** @file febase.hpp */
#ifndef SRC_SOLVER_FESTRUCT_LISTOFELEMENTS_FEBASE_FEBASE_HPP
#define SRC_SOLVER_FESTRUCT_LISTOFELEMENTS_FEBASE_FEBASE_HPP

#include <vector>

#include "../../../../Algorithms/Math/math.hpp"
#include "../../../../Core/core.hpp"
#include "felocalinfos.hpp"

template <typename T>
using LambdaOnPoint = std::function<T (Point *)>;

#define LAMBDA_EVALDOUBLE(X) \
    [] (Point * p) -> real_t {VOID_USE(p); return X; }

#define LAMBDA_FUNCTIONDOUBLE(x) \
    [] (Point * p) -> real_t {VOID_USE(p); x }

#define LAMBDA_EVALPOINT(x, y, z) \
    [] (Point * p) -> Point {VOID_USE(p); return Point(x, y, z); }

#define LAMBDA_FUNCTIONPOINT(x) \
    [] (Point * p) -> Point {VOID_USE(p); x }

class FEBase
{
public:
    FEBase ();
    virtual ~FEBase ();

    virtual void  LocalCompute (Point * atpoint, FELocalInfos * obj) = 0;
    virtual Point TransformRefToEle (Point * p)                      = 0;
    virtual Point TransformEleToRef (Point * p)                      = 0;

    void                    CastForCell (Cell * target);
    Cell *                  GetTarget ();
    VTK_CELL_TYPE           GetCellType ();
    ul_t                    GetNumberOfPoints ();
    ul_t                    GetNumberOfEdges ();
    real_t                  GetMesureOfRef ();
    ul_t                    GetOrderOfElement ();
    ul_t                    GetDimension ();
    PHYSICAL_CELL_TYPE      GetTypePhysicalBase ();
    Point *                 GetPoint (ul_t id);
    Edge *                  GetEdge (ul_t id);
    Point *                 GetNormalToEdge (ul_t id);
    Point *                 GetTangentToEdge (ul_t id);
    LambdaOnPoint<real_t> & GetPhi (ul_t id);
    LambdaOnPoint<Point> &  GetGradPhi (ul_t id);

protected:
    VTK_CELL_TYPE                      m_type_vtk;
    real_t                             m_volref;
    ul_t                               m_order;
    ul_t                               m_dim;
    PHYSICAL_CELL_TYPE                 m_type_physical;
    std::vector<Point *>               m_coor;
    std::vector<Edge *>                m_edges;
    std::vector<Point *>               m_edges_normals;
    std::vector<Point *>               m_edges_tangents;
    std::vector<LambdaOnPoint<real_t>> m_phi;
    std::vector<LambdaOnPoint<Point>>  m_grad_phi;
    Cell *                             m_target;
    FE_CLASS_TYPE                      m_cat;

private:
    FEBase (const FEBase &) = delete;
    FEBase & operator= (const FEBase &) = delete;
};

#endif /* SRC_SOLVER_FESTRUCT_LISTOFELEMENTS_FEBASE_FEBASE_HPP */
