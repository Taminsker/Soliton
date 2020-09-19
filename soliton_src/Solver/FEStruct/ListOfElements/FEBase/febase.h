/** @file febase.h */
#ifndef FEBASE_H
#define FEBASE_H

#include <vector>
#include <Algorithms/Math/math.h>
#include <Core/core.h>

#include <Solver/FEStruct/ListOfElements/FEBase/felocalinfos.h>

template <typename T>
using LambdaOnPoint =  std::function <T(Point*)>;

#define LAMBDA_EVALDOUBLE(X)\
    [] (Point* p) -> double {(void)p;; return X;}

#define LAMBDA_FUNCTIONDOUBLE(x)\
    [] (Point* p) -> double {(void)p; x}

#define LAMBDA_EVALPOINT(x, y, z)\
    [] (Point* p) -> Point {(void)p; return Point(x, y, z);}

#define LAMBDA_FUNCTIONPOINT(x)\
    [] (Point* p) -> Point {(void)p; x}

enum class PHYSICAL_CELL_TYPE
{
    EMPTY,
    VERTEX,
    LINE,
    TRIANGLE,
    QUADRANGLE,
    TETRAHEDRON,
    FIRST = EMPTY,
    LAST = TETRAHEDRON
};

enum class FE_CLASS_TYPE
{
    EMPTY = 0,
    LAGRANGE = 1,
    HERMITE = 2,
    FIRST = EMPTY,
    LAST = HERMITE
};

class FEBase
{
public:
    FEBase ();
    virtual ~FEBase ();

    virtual void    CastForCell (Cell* target) = 0;
    virtual void    LocalCompute (Point* atpoint, FELocalInfos* obj) = 0;
    virtual Point   TransformRefToEle (Point *p) = 0;
    virtual Point   TransformEleToRef (Point* p) = 0;
    Cell*           GetTarget ();

//    VTK_CELL_TYPE   type;
//    std::size_t npts;
//    std::size_t nedges;
//    double volRef;
//    std::size_t order;
//    std::size_t dim;
//    PHYSICAL_CELL_TYPE typePhysical;
//    std::vector <Point*> coor;
//    std::vector <Edge*> edges;
//    std::vector <Point*> edges_normals;
//    std::vector <Point*> edges_tangents;
//    std::vector <LambdaOnPoint <double>> Phi;
//    std::vector <LambdaOnPoint <Point>> gradPhi;

    VTK_CELL_TYPE           GetCellType ();
    std::size_t             GetNumberOfPoints ();
    std::size_t             GetNumberOfEdges ();
    double                  GetMesureOfRef ();
    std::size_t             GetOrderOfElement ();
    std::size_t             GetDimension ();
    PHYSICAL_CELL_TYPE      GetTypePhysicalBase ();
    Point*                  GetPoint (std::size_t id);
    Edge*                   GetEdge (std::size_t id);
    Point*                  GetNormalToEdge (std::size_t id);
    Point*                  GetTangentToEdge (std::size_t id);
    LambdaOnPoint<double>   GetPhi (std::size_t id);
    double                  EvalPhi(std::size_t id, Point* atpoint);
    LambdaOnPoint<Point>    GetGradPhi (std::size_t id);
    Point                   EvalGradPhi(std::size_t id, Point* atpoint);

protected:
    VTK_CELL_TYPE                           m_type_vtk;
//    std::size_t                             m_npts;
//    std::size_t                             m_nedges;
    double                                  m_volref;
    std::size_t                             m_order;
    std::size_t                             m_dim;
    PHYSICAL_CELL_TYPE                      m_type_physical;
    std::vector <Point*>                    m_coor;
    std::vector <Edge*>                     m_edges;
    std::vector <Point*>                    m_edges_normals;
    std::vector <Point*>                    m_edges_tangents;
    std::vector <LambdaOnPoint <double>>    m_phi;
    std::vector <LambdaOnPoint <Point>>     m_grad_phi;
    Cell*                                   m_target;
    FE_CLASS_TYPE                           m_cat;

private:
    FEBase (const FEBase&) = delete;
    FEBase& operator= (const FEBase&) = delete;
};


#endif // FEBASE_H
