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
    [] (Point* p) -> double {VOID_USE(p); return X;}

#define LAMBDA_FUNCTIONDOUBLE(x)\
    [] (Point* p) -> double {VOID_USE(p); x}

#define LAMBDA_EVALPOINT(x, y, z)\
    [] (Point* p) -> Point {VOID_USE(p); return Point(x, y, z);}

#define LAMBDA_FUNCTIONPOINT(x)\
    [] (Point* p) -> Point {VOID_USE(p); x}


class FEBase
{
public:
    FEBase ();
    virtual ~FEBase ();

    virtual void            LocalCompute (Point* atpoint, FELocalInfos* obj) = 0;
    virtual Point           TransformRefToEle (Point *p) = 0;
    virtual Point           TransformEleToRef (Point* p) = 0;

    void                    CastForCell (Cell* target);
    Cell*                   GetTarget ();
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
    LambdaOnPoint<double>&  GetPhi (std::size_t id);
    LambdaOnPoint<Point>&   GetGradPhi (std::size_t id);

protected:
    VTK_CELL_TYPE                           m_type_vtk;
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
