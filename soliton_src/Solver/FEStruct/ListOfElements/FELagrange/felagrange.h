/** @file felagrange.h */

#ifndef FELAGRANGE_H
#define FELAGRANGE_H

#include <Solver/FEStruct/ListOfElements/PhysicalElements/physicalelements.h>

enum class LAGRANGE_ORDER : int
{
    ORDER_P0 = 0,
    ORDER_P1 = 1,
    ORDER_P2 = 2,
    ORDER_P3 = 3,
    ORDER_P4 = 4,
    FIRST   = ORDER_P0,
    LAST    = ORDER_P4
};

/** @class FELagrange
 * \brief Defines a class for stuff.
 * \tparam name Type to work with.
 * \tparam npts Type to work with.
 * \tparam nddl Type to work with.
 */
template <PHYSICAL_CELL_TYPE t_name, int t_npts, int t_nddl>
class FELagrange : public FEPhysicalElement<t_name, t_npts>
{
public:
    /** @brief Generic constructor for FELagrange<_name, _npts, _nddl>
     */
    FELagrange ();

    /** @brief Destructor for FELagrange<_name, _npts, _nddl>
     */
    ~FELagrange () override = default;

    /** @brief Generic update FELagrange<_name, _npts, _nddl>
     */
    void CastForCell (Cell* target) override;

    /** @brief Generic local computation of the transposed inverse Jacobian matrix and of the determinant of the transformation from the reference element to the physical element.
     */
    void LocalCompute (Point* atpoint, FELocalInfos* obj) override;

    /** @brief Transformation between the reference element and the physical element FELagrange<_name, _npts, _nddl>
     */
    Point TransformRefToEle (Point *p) override;

    /** @brief Transformation between the physical element and the reference element FELagrange<_name, _npts, _nddl>
     */
    Point TransformEleToRef (Point* p) override;

protected:
    void Build ();

private:
    FELagrange (const FELagrange&) = delete;
    FELagrange& operator= (const FELagrange&) = delete;
};


template <PHYSICAL_CELL_TYPE t_name, int t_npts, int t_nddl>
FELagrange<t_name, t_npts, t_nddl>::FELagrange () :
    FEPhysicalElement<t_name, t_npts> ()
{
    Build ();
}

template <PHYSICAL_CELL_TYPE t_name, int t_npts, int t_nddl>
Point FELagrange<t_name, t_npts, t_nddl>::TransformRefToEle (Point* p)
{
    Point out;

    if (this->m_target->GetPoints()->size() != this->m_coor.size())
    {
        ERROR << "Wrong Cell in 'TransformRefToEle'" << ENDLINE;
        return out;
    }
    for (std::size_t i = 0; i < this->m_coor.size(); ++i)
        out = out + this->m_phi[i](p) * *this->m_target->GetPoints ()->at (i);
    return out;
}

template <PHYSICAL_CELL_TYPE t_name, int t_npts, int t_nddl>
Point FELagrange<t_name, t_npts, t_nddl>::TransformEleToRef (Point *p)
{
    (void)p;

    // We try to solve F(t) = T(t) - p = 0
    // with jac_{F}(t) = Jac_{T}(t)

    // Newton method
    // X_{n+1} = X_{n} - JacInvT(X_{n}) * F(X_N)

    double eps = 1e-10;
    double err = 100.;
    int maxiter = 100;
    int it = 0;

    Point x_n = *this->m_coor[0];
    Point x_np1 = x_n + Point(1., 1., 1.);
    Matrix3x3 Jac, JacInvT;
    Point F;

    for (it = 0; it < maxiter; ++it)
    {
        // Jacobian matrix
        Jac.setZero ();
        for (std::size_t i = 0; i < this->m_coor.size(); ++i)
        {
            Jac += OuterProduct (*this->m_target->GetPoints ()->at (i), this->m_grad_phi[i](&x_n));
        }

        // Jacobian inverse transpose
        JacInvT = Jac.completeOrthogonalDecomposition ().pseudoInverse ().transpose ();

        F = this->TransformRefToEle (&x_n) - *p;

        x_np1 = x_n - JacInvT * F;

        err = (x_n - x_np1).EuclidianNorm ();
        if (err < eps)
            break;

        x_n = x_np1;
    }

    return x_np1;
}

template <PHYSICAL_CELL_TYPE t_name, int t_npts, int t_nddl>
void FELagrange<t_name, t_npts, t_nddl>::LocalCompute (Point* atpoint, FELocalInfos* obj)
{
    // Init
    obj->JacInvT.setZero ();
    obj->detJac = 0.0;

    if (this->m_target->GetPoints()->size() != this->m_coor.size())
    {
        ERROR << "internal error" << ENDLINE;
        return;
    }

    // Jacobian matrix
    obj->Jac.setZero ();
    for (std::size_t i = 0; i < this->m_coor.size(); ++i)
    {
        obj->Jac += OuterProduct (*this->m_target->GetPoints ()->at (i), this->m_grad_phi [i] (atpoint));
    }

    // Jacobian inverse transpose
    obj->JacInvT = obj->Jac.completeOrthogonalDecomposition ().pseudoInverse ().transpose ();

    // Det of jacobian matrix
    Eigen::JacobiSVD<Matrix3x3> svd(obj->Jac);
    obj->detJac = 1.;
    for (int i = 0; i < 3; ++i)
    {
        double value = svd.singularValues ().coeff (i);
        if (std::abs (value) > 0)
            obj->detJac *= value;
    }
    obj->detJac = std::abs (obj->detJac);

    // Normal vectors on edges
    obj->normalEdge_ref = {0., 0., 0.};
    obj->tangentEdge_ref = {0., 0., 0.};

    for (std::size_t idedge = 0; idedge < this->m_edges.size(); idedge++)
    {
        Edge* e = this->m_edges [idedge];
        double prod_phi = 1.;

        for (Point* p : *e->GetPoints ())
            prod_phi *= this->m_phi [static_cast<std::size_t>(p->GetGlobalIndex ())] (atpoint);

        obj->normalEdge_ref  += std::abs(prod_phi) * *this->m_edges_normals [idedge];
        obj->tangentEdge_ref += std::abs(prod_phi) * *this->m_edges_tangents [idedge];
    }

    obj->Jac_cmp = obj->Jac;
    obj->JacInvT_cmp = obj->JacInvT;

    for (int i = 0; i < 3; i++)
    {
        if (obj->Jac_cmp.col (i) == Eigen::Vector3d::Zero ())
            obj->Jac_cmp.coeffRef (i, i) = 1.;

        if (obj->JacInvT_cmp.col (i) == Eigen::Vector3d::Zero ())
            obj->JacInvT_cmp.coeffRef (i, i) = 1.;
    }

    obj->normalEdge = obj->JacInvT_cmp * obj->normalEdge_ref;
    obj->normalEdge.Normalize ();
    obj->normalEdge_ref.Normalize ();

    obj->tangentEdge = obj->Jac_cmp * obj->tangentEdge_ref;
    obj->tangentEdge.Normalize ();
    obj->tangentEdge_ref.Normalize ();

    obj->normalCell = CrossProduct (obj->normalEdge, obj->tangentEdge);
    obj->normalCell.Normalize ();
    obj->normalCell_ref = CrossProduct (obj->normalEdge_ref, obj->tangentEdge_ref);
    obj->normalCell_ref.Normalize ();

    return;
}

template <PHYSICAL_CELL_TYPE t_name, int t_npts, int t_nddl>
void FELagrange<t_name, t_npts, t_nddl>::CastForCell (Cell* target)
{
    this->m_target = target;

    return;
}

/**
 * \brief Emp0N0DDL Lagrange Element
 * \class Emp0N0DDL
 * \implements FELagrange
 * \implements Emp0N
 */
typedef FELagrange<PHYSICAL_CELL_TYPE::EMPTY, 0, 0>        Emp0N0DDL;

/**
 * \brief Ver1N1DDL Lagrange Element
 * \class Ver1N1DDL
 * \implements FELagrange
 * \implements Ver1N
 */
typedef FELagrange<PHYSICAL_CELL_TYPE::VERTEX, 1, 1>       Ver1N1DDL;

/**
 * \brief Lin2N2DDL Lagrange Element
 * \class Lin2N2DDL
 * \implements FELagrange
 * \implements Lin2N
 */
typedef FELagrange<PHYSICAL_CELL_TYPE::LINE, 2, 2>         Lin2N2DDL;

/**
 * \brief Lin3N3DDL Lagrange Element
 * \class Lin3N3DDL
 * \implements FELagrange
 * \implements Lin3N
 */
typedef FELagrange<PHYSICAL_CELL_TYPE::LINE, 3, 3>         Lin3N3DDL;

/**
 * \brief Tri3N3DDL Lagrange Element
 * \class Tri3N3DDL
 * \implements FELagrange
 * \implements Tri3N
 */
typedef FELagrange<PHYSICAL_CELL_TYPE::TRIANGLE, 3, 3>     Tri3N3DDL;

/**
 * \brief Tri6N6DDL Lagrange Element
 * \class Tri6N6DDL
 * \implements FELagrange
 * \implements Tri6N
 */
typedef FELagrange<PHYSICAL_CELL_TYPE::TRIANGLE, 6, 6>     Tri6N6DDL;

/**
 * \brief Quad4N4DDL Lagrange Element
 * \class Quad4N4DDL
 * \implements FELagrange
 * \implements Quad4N
 */
typedef FELagrange<PHYSICAL_CELL_TYPE::QUADRANGLE, 4, 4>   Quad4N4DDL;

/**
 * \brief Quad8N8DDL Lagrange Element
 * \class Quad8N8DDL
 * \implements FELagrange
 * \implements Quad8N
 */
typedef FELagrange<PHYSICAL_CELL_TYPE::QUADRANGLE, 8, 8>   Quad8N8DDL;

/**
 * \brief Quad9N9DDL Lagrange Element
 * \class Quad9N9DDL
 * \implements FELagrange
 * \implements Quad9N
 */
typedef FELagrange<PHYSICAL_CELL_TYPE::QUADRANGLE, 9, 9>   Quad9N9DDL;

/**
 * \brief Tet4N4DDL Lagrange Element
 * \class Tet4N4DDL
 * \implements FELagrange
 * \implements Tet4N
 */
typedef FELagrange<PHYSICAL_CELL_TYPE::TETRAHEDRON, 4, 4>  Tet4N4DDL;
//}

#endif // FELAGRANGE_H
