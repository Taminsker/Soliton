#ifdef STAND_ALONE
#define BOOST_TEST_MODULE SolitonTests_Quadratures
#endif

#include <boost/test/unit_test.hpp>
#include <Solver/solver.h>
#include <Algorithms/Math/math.h>
#include <Core/core.h>

real_t fun1d (Point p)
{
    return p.x;
}

real_t fun2d (Point p)
{
    return p.x + p.y;
}

Point grad2d_1 (Point p)
{
    return p;
}

Point grad2d_2 (Point)
{
    return Point(1, 1, 1);
}
BOOST_AUTO_TEST_SUITE(Quadratures)

BOOST_AUTO_TEST_CASE(quadrature_1d)
{
    BOOST_TEST_MESSAGE("-- Tests quadrature 1D");
    {
        Point p0(0, 0, 0);
        Point p1(4.5, 0, 0);
        Cell cell;
        cell.AddPoint (&p0);
        cell.AddPoint (&p1);
        cell.SetType (GMSH_CELL_TYPE::GMSH_2_NODE_LINE);

        FEStore festore;
        FEBase* obj = festore.GetElementFor (&cell);
        obj->CastForCell (&cell);

        QuadStore quadstore;

        QuadStore::QuadObject* quadobj = quadstore.Get (obj);
        ul_t npts = quadobj->npts;

        real_t value = 0.;
        FELocalInfos loc;

        for (ul_t k = 0; k < npts; ++k)
        {
            Point pk = quadobj->pts [k];
            real_t wk = quadobj->w [k];

            Point intp = obj->TransformRefToEle (&pk);
            obj->LocalCompute (&pk, &loc);

            //      INFOS << loc.detJac << ENDLINE;

            value += loc.detJac * wk * fun1d(intp);
        }

        real_t realvalue = 0.5 * (p1.x * p1.x - p0.x * p0.x);

        //    std::cout << realvalue << std::endl;
        //    std::cout << value << std::endl;
        BOOST_CHECK_EQUAL (std::abs(value - realvalue) < 1e-3, true);
    }

    //  {
    //    Point p0(0, 0, 0);
    //    Point p1(3, 0, 0);
    //    Point p2 (1, 0, 0);
    //    Cell cell;
    //    cell.AddPoint (&p0);
    //    cell.AddPoint (&p1);
    //    cell.AddPoint (&p2);
    //    cell.SetType (GMSH_CELL_TYPE::GMSH_3_NODE_QUADRATIC_LINE);

    //    FEStore festore;
    //    FEBase* obj = festore.GetElementFor (&cell);
    //    obj->CastForCell (&cell);

    //    QuadStore quadstore;
    //    QuadStore::QuadObject* quadobj = quadstore.Get (obj);
    //    ul_t npts = quadobj->npts;
    ////    std::cout << "quadid : " << quadobj->npts << std::endl;

    //    real_t value = 0.;
    //    real_t detjac = 0;
    //    Matrix3x3 JacInvT;

    //    for (ul_t k = 0; k < npts; ++k)
    //    {
    //      Point pk = quadobj->pts [k];
    //      real_t wk = quadobj->w [k];

    //      Point intp = obj->TransformRefToEle (&pk);
    //      obj->LocalCompute (&pk, &JacInvT, &detjac);

    ////      std::cout << intp << detjac << std::endl;
    //      value += detjac * wk * fun1d(intp);
    //    }

    //    real_t realvalue = 0.5 * (p1.x * p1.x - p0.x * p0.x);
    ////    real_t realvalue = p1.x - p0.x;


    //    BOOST_CHECK_EQUAL (std::abs(value - realvalue) < 1e-3, true);
    //  }
}



BOOST_AUTO_TEST_CASE(quadrature_2d)
{
    BOOST_TEST_MESSAGE("-- Tests quadrature 2D");

    Point p0(0, 0, 0);
    Point p1(5, 0, 0);
    Point p2(0, 3, 0);
    Cell cell;
    cell.AddPoint (&p0);
    cell.AddPoint (&p1);
    cell.AddPoint (&p2);
    cell.SetType (GMSH_CELL_TYPE::GMSH_3_NODE_TRIANGLE);

    FEStore festore;
    FEBase* obj = festore.GetElementFor (&cell);
    obj->CastForCell (&cell);

    QuadStore quadstore;
    QuadStore::QuadObject* quadobj = quadstore.Get (obj);
    ul_t npts = quadobj->npts;

    //  INFOS << "quadobj : type " << quadobj->type << ENDLINE;
    //  INFOS << "quadobj : npts " << quadobj->npts << ENDLINE;
    //  INFOS << "quadobj : order " << quadobj->order << ENDLINE;
    //  INFOS << "quadobj : w " << ENDLINE;
    //  for (ul_t i = 0; i < quadobj->npts; ++i)
    //    std::cout << "$ " << quadobj->w [i] << ENDLINE;
    //  INFOS << "quadobj : pts " << ENDLINE;
    //  for (ul_t i = 0; i < quadobj->npts; ++i)
    //    std::cout << "$ " << quadobj->pts [i] << ENDLINE;

    std::vector<Triplet_eig> triplet;
    triplet.push_back ({0, 0, 1});
    triplet.push_back ({0, 0, 1});
    triplet.push_back ({0, 0, 0.5});
    Eigen::SparseMatrix<real_t> m(3, 3);
    m.setFromTriplets(triplet.begin(), triplet.end());

    real_t value1 = 0.;
    real_t value2 = 0.;
    FELocalInfos loc;
    std::function<Point(Point*)> gradphi_0 = obj->GetGradPhi (0);

    for (ul_t k = 0; k < npts; ++k)
    {
        Point pk = quadobj->pts [k];
        real_t wk = quadobj->w [k];
        Point intp = obj->TransformRefToEle (&pk);
        obj->LocalCompute (&pk, &loc);

        value1 += loc.detJac * wk * fun2d(intp);

        value2 += loc.detJac * wk * ((loc.JacInvT * gradphi_0(&pk)) | (loc.JacInvT * gradphi_0(&pk)));
    }

    real_t realvalue = 20.;

    //  INFOS << "value1 = " << value1 << ENDLINE;
    //  INFOS << "value2 = " << value2 << ENDLINE;

    BOOST_CHECK_EQUAL (std::abs(value1 - realvalue) < 1e-10, true);
    //  BOOST_CHECK_EQUAL (std::abs(value2 - realvalue) < 1e-10, true);




}

BOOST_AUTO_TEST_SUITE_END()
