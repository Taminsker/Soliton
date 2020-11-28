#include "quadstore.h"

QuadStore::QuadStore() :
    m_quad_empty (nullptr),
    m_quad_vertex (nullptr),
    m_quad_line ({}),
    m_quad_triangle ({}),
    m_quad_quadrangle ({})
{
    // 0D
    m_quad_empty = new QuadObject (PHYSICAL_CELL_TYPE::EMPTY, 0, 0);
    m_quad_empty->w = {};
    m_quad_empty->pts = {};

    // 1D
    m_quad_vertex = new QuadObject (PHYSICAL_CELL_TYPE::VERTEX, 1, 1);
    m_quad_vertex->w = {1};
    m_quad_vertex->pts = {{}};

    QuadObject* line_1 = new QuadObject (PHYSICAL_CELL_TYPE::LINE, 1, 1);
    line_1->w = {2};
    line_1->pts = {{}};
    m_quad_line.push_back (line_1);

    QuadObject* line_2 = new QuadObject (PHYSICAL_CELL_TYPE::LINE, 2, 3);
    line_2->w = {1., 1.};
    line_2->pts = {{-0.577350269189625}, {+0.577350269189625}};
    m_quad_line.push_back (line_2);


    QuadObject* line_3 = new QuadObject (PHYSICAL_CELL_TYPE::LINE, 3, 5);
    line_3->w = {0.555555555555556, 0.888888888888889, 0.555555555555556};
    line_3->pts = {{-0.774596669241483}, {}, {+0.774596669241483}};
    m_quad_line.push_back (line_3);

    QuadObject* line_4 = new QuadObject (PHYSICAL_CELL_TYPE::LINE, 4, 7);
    line_4->w = {0.347854845137454,
                 0.652145154862545,
                 0.652145154862545,
                 0.347854845137454};
    line_4->pts = {{-0.861136311594052},
                   {-0.339981043584856},
                   {+0.339981043584856},
                   {+0.861136311594052}};
    m_quad_line.push_back (line_4);

    QuadObject* line_5 = new QuadObject (PHYSICAL_CELL_TYPE::LINE, 5, 9);
    line_5->w = {0.236926885056189,
                 0.478628670499365,
                 0.568888889888889,
                 0.478628670499365,
                 0.236926885056189};
    line_5->pts = {{-0.906179845938664},
                   {-0.538469310105683},
                   {0.0000000000000000},
                   {+0.538469310105683},
                   {+0.906179845938664}};
    m_quad_line.push_back (line_5);


    // 2D
    QuadObject* triangle_1 = new QuadObject (PHYSICAL_CELL_TYPE::TRIANGLE, 1, 1);
    triangle_1->w = {0.5};
    triangle_1->pts = {{+0.3333333333333333, +0.3333333333333333}};
    m_quad_triangle.push_back (triangle_1);

    QuadObject* triangle_3 = new QuadObject (PHYSICAL_CELL_TYPE::TRIANGLE, 3, 2);
    triangle_3->w = {+0.1666666666666667,
                     +0.1666666666666667,
                     +0.1666666666666667};
    triangle_3->pts = {{+0.6666666666666667, +0.1666666666666667},
                       {+0.1666666666666667, +0.6666666666666667},
                       {+0.1666666666666667, +0.1666666666666667}};
    m_quad_triangle.push_back (triangle_3);

    QuadObject* triangle_4 = new QuadObject (PHYSICAL_CELL_TYPE::TRIANGLE, 4, 3);
    triangle_4->w = {-0.2812500000000000,
                     +0.2604166666666667,
                     +0.2604166666666667,
                     +0.2604166666666667,
                     +0.2604166666666667};
    triangle_4->pts = {{+0.3333333333333333, +0.3333333333333333},
                       {+0.2000000000000000, +0.2000000000000000},
                       {+0.2000000000000000, +0.6000000000000000},
                       {+0.6000000000000000, +0.2000000000000000}};
    m_quad_triangle.push_back (triangle_4);

    QuadObject* triangle_6 = new QuadObject (PHYSICAL_CELL_TYPE::TRIANGLE, 6, 4);
    triangle_6->w = {+0.1116907948390055,
                     +0.1116907948390055,
                     +0.1116907948390055,
                     +0.0549758718276610,
                     +0.0549758718276610,
                     +0.0549758718276610};
    triangle_6->pts = {{+0.108103018168070, +0.445948490915965},
                       {+0.445948490915965, +0.108103018168070},
                       {+0.445948490915965, +0.445948490915965},
                       {+0.816847572980459, +0.091576213509771},
                       {+0.091576213509771, +0.816847572980459},
                       {+0.091576213509771, +0.091576213509771}};
    m_quad_triangle.push_back (triangle_6);

    QuadObject* triangle_12 = new QuadObject (PHYSICAL_CELL_TYPE::TRIANGLE, 12, 6);
    triangle_12->w = {+0.0254224531851035,
                      +0.0254224531851035,
                      +0.0254224531851035,
                      +0.0583931378631895,
                      +0.0583931378631895,
                      +0.0583931378631895,
                      +0.0414255378091870,
                      +0.0414255378091870,
                      +0.0414255378091870,
                      +0.0414255378091870,
                      +0.0414255378091870,
                      +0.0414255378091870};
    triangle_12->pts = {{+0.873821971016996, +0.063089014491502},
                        {+0.063089014491502, +0.873821971016996},
                        {+0.063089014491502, +0.063089014491502},
                        {+0.501426509658179, +0.249286745170910},
                        {+0.249286745170910, +0.501426509658179},
                        {+0.249286745170910, +0.249286745170910},
                        {+0.636502499121399, +0.310352451033785},
                        {+0.636502499121399, +0.053145049844816},
                        {+0.310352451033785, +0.636502499121399},
                        {+0.310352451033785, +0.053145049844816},
                        {+0.053145049844816, +0.310352451033785},
                        {+0.053145049844816, +0.636502499121399}};
    m_quad_triangle.push_back (triangle_12);

    //    QuadObject* triangle_16 = new QuadObject (TYPE_TRIANGLE, 16, 8);
    //    triangle_16->w = {0.5};
    //    triangle_16->pts = {{+0.3333333333333333, +0.3333333333333333}};
    //    m_quad_triangle.push_back (triangle_16);


    m_quad_quadrangle.push_back (QuadrangleFromLine (line_1));
    m_quad_quadrangle.push_back (QuadrangleFromLine (line_2));
    m_quad_quadrangle.push_back (QuadrangleFromLine (line_3));
    m_quad_quadrangle.push_back (QuadrangleFromLine (line_4));
    m_quad_quadrangle.push_back (QuadrangleFromLine (line_5));
}

QuadStore::~QuadStore ()
{
    delete m_quad_empty;
    delete m_quad_vertex;

    for (QuadObject * obj : m_quad_line)
        delete obj;
    m_quad_line.clear ();

    for (QuadObject * obj : m_quad_triangle)
        delete obj;
    m_quad_triangle.clear ();

    for (QuadObject * obj : m_quad_quadrangle)
        delete obj;
    m_quad_quadrangle.clear ();
}

QuadStore::QuadObject* QuadStore::Get (FEBase *object)
{
    ul_t order_attempt = object->GetOrderOfElement () + 1;

    switch (object->GetTypePhysicalBase ())
    {
    default:
        return m_quad_empty;
    case PHYSICAL_CELL_TYPE::VERTEX:
        return m_quad_vertex;
    case PHYSICAL_CELL_TYPE::LINE:
        for (ul_t id = 0; id < m_quad_line.size (); ++id)
            if (m_quad_line [id]->order > order_attempt)
                return m_quad_line [id];
        break;
    case PHYSICAL_CELL_TYPE::TRIANGLE:
        for (ul_t id = 0; id < m_quad_triangle.size (); ++id)
            if (m_quad_triangle [id]->order > order_attempt)
                return m_quad_triangle [id];
        break;
    case PHYSICAL_CELL_TYPE::QUADRANGLE:
        for (ul_t id = 0; id < m_quad_quadrangle.size (); ++id)
            if (m_quad_quadrangle [id]->order > order_attempt)
                return m_quad_quadrangle [id];
        break;
    }

    return m_quad_empty;
}

QuadStore::QuadObject::QuadObject (PHYSICAL_CELL_TYPE _type, ul_t _npts, ul_t _order) :
    type (_type),
    npts (_npts),
    order (_order),
    w ({}),
    pts ({})
{}

QuadStore::QuadObject::~QuadObject ()
{
    type = PHYSICAL_CELL_TYPE::EMPTY;
    npts = 0;
    order = 0;
    w.clear ();
    pts.clear ();
}

void QuadStore::QuadObject::Print()
{
    INFOS << "Quad : type " << to_string(type) << " npts " << npts << " order " << order << ENDLINE;

    for (ul_t i = 0; i < npts; ++i)
        INFOS << "\t P : " << pts.at (i) << " --- W : " << w.at (i) << ENDLINE;

    return;
}


QuadStore::QuadObject* QuadrangleFromLine(QuadStore::QuadObject* line)
{
    ul_t npts = line->npts;
    npts *= npts;

    QuadStore::QuadObject* quadrangle = new QuadStore::QuadObject (PHYSICAL_CELL_TYPE::QUADRANGLE, npts, line->order);

    for (ul_t i = 0; i < line->npts; ++i)
        for (ul_t j = 0; j < line->npts; ++j)
        {
            real_t w = line->w.at (i) * line->w.at (j);
            Point p = {line->pts.at (i).x, line->pts.at (j).x};

            quadrangle->w.push_back (w);
            quadrangle->pts.push_back (p);
        }

    return quadrangle;
}
