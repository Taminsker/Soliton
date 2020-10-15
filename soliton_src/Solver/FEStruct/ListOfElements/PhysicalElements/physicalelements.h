/** @file physicalelements.h **/
#ifndef PHYSICALELEMENTS_H
#define PHYSICALELEMENTS_H

#include <Solver/FEStruct/ListOfElements/FEBase/febase.h>

/** @brief Template class for Finite Element : physical part.
 * \tparam _name name of type @see PHYSICAL_CELL_TYPE.
 * \tparam _npts number of points on the physical element.
 *  \see from FEBase.
 */
template <PHYSICAL_CELL_TYPE t_name, int t_npts>
class FEPhysicalElement : public FEBase
{
public:
    FEPhysicalElement ();
    virtual ~FEPhysicalElement () = default;
private:
    void Build ();
    FEPhysicalElement (const FEPhysicalElement&) = delete;
    FEPhysicalElement& operator= (const FEPhysicalElement&) = delete;

    int count_points;
    int count_edges;

    void AddPointDefinition (double x = 0, double y = 0., double z = 0.);
    void AddEdgeElement (std::vector <std::size_t> vec, GMSH_CELL_TYPE type, Point normal, Point tangent);
};

/** @brief Generic constructor for FEPhysicalElement (default empty element).
 * \tparam _name name of type @see PHYSICAL_CELL_TYPE.
 * \tparam _npts number of points on the physical element.
 * \see from FEBase.
 */
template <PHYSICAL_CELL_TYPE t_name, int t_npts>
FEPhysicalElement <t_name, t_npts>::FEPhysicalElement () :
    FEBase (),
    count_points (0),
    count_edges (0)
{
    Build ();
}

template <PHYSICAL_CELL_TYPE t_name, int t_npts>
void FEPhysicalElement <t_name, t_npts>::AddPointDefinition (double x, double y, double z)
{
    Point* p = new Point(x, y, z);
    p->SetGlobalIndex (count_points);
    m_coor.push_back (p);
    count_points++;

    return;
}

template <PHYSICAL_CELL_TYPE t_name, int t_npts>
void FEPhysicalElement <t_name, t_npts>::AddEdgeElement (std::vector <std::size_t> vec, GMSH_CELL_TYPE gmshtype, Point normal, Point tangent)
{
    Edge* e = new Edge();
    e->SetType (gmshtype);
    e->SetGlobalIndex (count_edges);

    for (std::size_t id : vec)
        e->AddPoint (m_coor [id]);

    m_edges.push_back (e);

    m_edges_normals.push_back (new Point(normal));
    m_edges_tangents.push_back (new Point(tangent));

    count_edges++;
    return;
}

/**
 * \brief Emp0N Physical Element
 * \class Emp0N
 * \implements FEPhysicalElement
 */
typedef FEPhysicalElement<PHYSICAL_CELL_TYPE::EMPTY, 0>        Emp0N;

/**
 * \brief Ver1N Physical Element
 * \class Ver1N
 * \implements FEPhysicalElement
 */
typedef FEPhysicalElement<PHYSICAL_CELL_TYPE::VERTEX, 1>       Ver1N;

/**
 * \brief Lin2N Physical Element
 * \class Lin2N
 * \implements FEPhysicalElement
 */
typedef FEPhysicalElement<PHYSICAL_CELL_TYPE::LINE, 2>         Lin2N;

/**
 * \brief Lin3N Physical Element
 * \class Lin3N
 * \implements FEPhysicalElement
 */
typedef FEPhysicalElement<PHYSICAL_CELL_TYPE::LINE, 3>         Lin3N;

/**
 * \brief Tri3N Physical Element
 * \class Tri3N
 * \implements FEPhysicalElement
 */
typedef FEPhysicalElement<PHYSICAL_CELL_TYPE::TRIANGLE, 3>     Tri3N;

/**
 * \brief Tri6N Physical Element
 * \class Tri6N
 * \implements FEPhysicalElement
 */
typedef FEPhysicalElement<PHYSICAL_CELL_TYPE::TRIANGLE, 6>     Tri6N;

/**
 * \brief Quad4N Physical Element
 * \class Quad4N
 * \implements FEPhysicalElement
 */
typedef FEPhysicalElement<PHYSICAL_CELL_TYPE::QUADRANGLE, 4>   Quad4N;

/**
 * \brief Quad8N Physical Element
 * \class Quad8N
 * \implements FEPhysicalElement
 */
typedef FEPhysicalElement<PHYSICAL_CELL_TYPE::QUADRANGLE, 8>   Quad8N;

/**
 * \brief Quad9N Physical Element
 * \class Quad9N
 * \implements FEPhysicalElement
 */
typedef FEPhysicalElement<PHYSICAL_CELL_TYPE::QUADRANGLE, 9>   Quad9N;

/**
 * \brief Tet4N Physical Element
 * \class Tet4N
 * \implements FEPhysicalElement
 */
typedef FEPhysicalElement<PHYSICAL_CELL_TYPE::TETRAHEDRON, 4>  Tet4N;

#endif // PHYSICALELEMENTS_H
