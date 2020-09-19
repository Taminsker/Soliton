#ifdef STAND_ALONE
#define BOOST_TEST_MODULE SolitonTests_LagrangeElements
#endif

#include <boost/test/unit_test.hpp>
#include <Solver/solver.h>
#include <Core/core.h>

BOOST_AUTO_TEST_SUITE(Emp0N0DDL_felag)

BOOST_AUTO_TEST_CASE(Emp0N0DDL_constructor)
{
    BOOST_TEST_MESSAGE("-- Tests for Emp0N0DDL");
    FELagrange<PHYSICAL_CELL_TYPE::EMPTY, 0, 0> obj;

    BOOST_CHECK_EQUAL(obj.GetNumberOfPoints (), 0);
    BOOST_CHECK_EQUAL(obj.GetCellType (), VTK_CELL_TYPE::VTK_EMPTY_CELL);

}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(Ver1N1DDL_felag)

BOOST_AUTO_TEST_CASE(Ver1N1DDL_constructor)
{
    BOOST_TEST_MESSAGE("-- Tests for Ver1N1DDL");
    FELagrange<PHYSICAL_CELL_TYPE::VERTEX, 1, 1> obj;

    Point p0(1, 1, 1);
    Cell cell;
    cell.AddPoint (&p0);
    cell.SetType (GMSH_CELL_TYPE::GMSH_1_NODE_POINT);
    obj.CastForCell (&cell);

    BOOST_CHECK_EQUAL(obj.GetNumberOfPoints (), 1);
    BOOST_CHECK_EQUAL(obj.GetCellType (), VTK_CELL_TYPE::VTK_VERTEX);

    Point k0(0, 0, 0);
    Point p00 = obj.TransformRefToEle (&k0);
    BOOST_CHECK_EQUAL(p0, p00);

    FELocalInfos loc;
    obj.LocalCompute (&k0, &loc);
    BOOST_CHECK_EQUAL(std::isinf (loc.detJac), false);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(Lin2N2DDL_felag)

BOOST_AUTO_TEST_CASE(Lin2N2DDL_constructor)
{
    BOOST_TEST_MESSAGE("-- Tests for Lin2N2DDLL");
    FELagrange<PHYSICAL_CELL_TYPE::LINE, 2, 2> obj;

    Point p0(-0.25, 0, 0);
    Point p1(0.25, 0, 0);
    Cell cell;
    cell.AddPoint (&p0);
    cell.AddPoint (&p1);
    cell.SetType (GMSH_CELL_TYPE::GMSH_2_NODE_LINE);
    obj.CastForCell (&cell);

    BOOST_CHECK_EQUAL(obj.GetNumberOfPoints (), 2);
    BOOST_CHECK_EQUAL(obj.GetCellType (), VTK_CELL_TYPE::VTK_LINE);

    Point k0(-1, 0, 0);
    Point k1(1, 0, 0);
    Point p00 = obj.TransformRefToEle (&k0);
    BOOST_CHECK_EQUAL(p0, p00);

    Point p11 = obj.TransformRefToEle (&k1);
    BOOST_CHECK_EQUAL(p1, p11);

    FELocalInfos loc;
    obj.LocalCompute (&k0, &loc);
    double vol = (p1 - p0).EuclidianNorm ();
    BOOST_CHECK_EQUAL(std::isinf (loc.detJac), false);
    BOOST_CHECK_EQUAL (loc.detJac, 0.5 * vol);

    Point kpt1(0.5, 0, 0);
    Point pt1 = obj.TransformRefToEle (&kpt1);
    Point kpt11 = obj.TransformEleToRef (&pt1);

    BOOST_CHECK_EQUAL(kpt1, kpt11);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(Lin3N3DDL_felag)

BOOST_AUTO_TEST_CASE(Lin3N3DDL_constructor)
{
    BOOST_TEST_MESSAGE("-- Tests for Lin3N3DDLL");
    FELagrange<PHYSICAL_CELL_TYPE::LINE, 3, 3> obj;

    Point p0(0, 0, 0);
    Point p1(4.5, 0, 0);
    Point p2(2.3, 0, 0);
    Cell cell;
    cell.AddPoint (&p0);
    cell.AddPoint (&p1);
    cell.AddPoint (&p2);
    cell.SetType (GMSH_CELL_TYPE::GMSH_3_NODE_QUADRATIC_LINE);
    obj.CastForCell (&cell);

    BOOST_CHECK_EQUAL(obj.GetNumberOfPoints (), 3);
    BOOST_CHECK_EQUAL(obj.GetCellType (), VTK_CELL_TYPE::VTK_QUADRATIC_EDGE);

    Point k0(-1, 0, 0);
    Point p00 = obj.TransformRefToEle (&k0);
    BOOST_CHECK_EQUAL(p0, p00);

    Point k1(1, 0, 0);
    Point p11 = obj.TransformRefToEle (&k1);
    BOOST_CHECK_EQUAL(p1, p11);

    Point k2(0, 0, 0);
    Point p22 = obj.TransformRefToEle (&k2);
    BOOST_CHECK_EQUAL(p2, p22);

    FELocalInfos loc;
    obj.LocalCompute (&k0, &loc);
//        double vol = (p1 - p0).EuclidianNorm ();
    BOOST_CHECK_EQUAL(std::isinf (loc.detJac), false);
    //    BOOST_CHECK_EQUAL (detJac, 0.5 * vol);

    Point kpt1(0.5, 0, 0);
    Point pt1 = obj.TransformRefToEle (&kpt1);
    Point kpt11 = obj.TransformEleToRef (&pt1);

    BOOST_CHECK_EQUAL(kpt1, kpt11);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(Tri3N3DDL_felag)

BOOST_AUTO_TEST_CASE(Tri3N3DDL_constructor)
{
    BOOST_TEST_MESSAGE("-- Tests for Tri3N3DDL");
    Tri3N3DDL obj;

    Point p0(1., 1., 0);
    Point p1(2.5, 1., 0);
    Point p2(1.5, 1.5, 0);
    Cell cell;
    cell.AddPoint (&p0);
    cell.AddPoint (&p1);
    cell.AddPoint (&p2);
    cell.SetType (GMSH_CELL_TYPE::GMSH_3_NODE_TRIANGLE);

    obj.CastForCell (&cell);

    {
        BOOST_CHECK_EQUAL(obj.GetNumberOfPoints (), 3);
        BOOST_CHECK_EQUAL(obj.GetCellType (), VTK_CELL_TYPE::VTK_TRIANGLE);

        Point k0(0, 0, 0);
        Point p00 = obj.TransformRefToEle (&k0);
        BOOST_CHECK_EQUAL(p0, p00);

        Point k1(1, 0, 0);
        Point p11 = obj.TransformRefToEle (&k1);
        BOOST_CHECK_EQUAL(p1, p11);

        Point k2(0, 1, 0);
        Point p22 = obj.TransformRefToEle (&k2);
        BOOST_CHECK_EQUAL(p2, p22);

        FELocalInfos loc;
        obj.LocalCompute (&k0, &loc);
        BOOST_CHECK_EQUAL(std::isinf (loc.detJac), false);

        double vol = 0.5 * CrossProduct (p1-p0, p2-p0).EuclidianNorm ();

//        std::cout << loc.detJac << std::endl;
//        std::cout << obj.GetDimension () * vol << std::endl;

        BOOST_CHECK_EQUAL (std::abs(loc.detJac - static_cast<double>(obj.GetDimension ()) * vol) < 1e-7, true);

        Point kpt1(0.5, 0.25, 0.);
        Point pt1 = obj.TransformRefToEle (&kpt1);
        Point kpt11 = obj.TransformEleToRef (&pt1);

        BOOST_CHECK_EQUAL(kpt1, kpt11);
    }
    {
        BOOST_TEST_CHECKPOINT ("Elementary matrix");
        FEStore festore;
        QuadStore quadstore;
        FEBase* febase = festore.GetElementFor (&cell);
        febase->CastForCell (&cell);
        QuadStore::QuadObject* quadobject = quadstore.Get (febase);

        // compute
        Matrix3x3 mat;
        mat.setZero ();
        FELocalInfos loc;

        for (std::size_t i = 0; i < 3; ++i)
        {
//            Point* p_i = cell.GetPoints ()->at (i);
            std::function<Point(Point*)> grad_phi_i = febase->GetGradPhi (i);

            for (std::size_t j = 0; j < 3; ++j)
            {
//                Point* p_j = cell.GetPoints ()->at (j);
                std::function<Point(Point*)> grad_phi_j = febase->GetGradPhi (j);

                double coeff = 0.;

                for (std::size_t k = 0; k < quadobject->npts; ++k)
                {
                    Point pk = quadobject->pts [k];
                    double wk = quadobject->w [k];

                    febase->LocalCompute (&pk, &loc);
//                    std::cout << "JacInvT : [" << i << ", " << j << "]" << std::endl << JacInvT << std::endl;

                    coeff += loc.detJac * wk * ((loc.JacInvT * grad_phi_i(&pk)) | (loc.JacInvT * grad_phi_j(&pk)));

                }

                mat.coeffRef(int(i),int(j)) = coeff;
            }
        }

        // real
        Matrix3x3 matreal;
        matreal.setZero();
        Eigen::Matrix<double, 2, 2> middle;
        middle.col(0) << (p1 - p0).x, (p1 - p0).y ;
        middle.col(1) << (p2 - p0).x, (p2 - p0).y;

        double vol = 0.5 * CrossProduct (p1-p0, p2-p0).EuclidianNorm ();

        auto temp = middle.inverse ().transpose();
        middle = temp;

        for (std::size_t i = 0; i < 3; ++i)
        {
            Eigen::Vector2d grad_phi_i;
            grad_phi_i << febase->GetGradPhi (i) (&p0).x, febase->GetGradPhi (i) (&p0).y;
            grad_phi_i = middle * grad_phi_i;

            for (std::size_t j = 0; j < 3; ++j)
            {
                Eigen::Vector2d grad_phi_j;
                grad_phi_j << febase->GetGradPhi (j) (&p0).x, febase->GetGradPhi (j) (&p0).y;
                grad_phi_j = middle * grad_phi_j;

                matreal.coeffRef(int(i),int(j)) = vol * grad_phi_i.adjoint () * grad_phi_j;
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(Tri6N6DDL_felag)

BOOST_AUTO_TEST_CASE(Tri6N6DDL_constructor)
{
    BOOST_TEST_MESSAGE("-- Tests for Tri6N6DDL");
    FELagrange<PHYSICAL_CELL_TYPE::TRIANGLE, 6, 6> obj;

    Point p0(0, 0, 0);
    Point p1(1, 0, 0);
    Point p2(0, 1, 0);
    Point p3(0.5, 0, 0);
    Point p4(0.5, 0.5, 0);
    Point p5(0, 0.5, 0);
    Cell cell;
    cell.AddPoint (&p0);
    cell.AddPoint (&p1);
    cell.AddPoint (&p2);
    cell.AddPoint (&p3);
    cell.AddPoint (&p4);
    cell.AddPoint (&p5);
    cell.SetType (GMSH_CELL_TYPE::GMSH_6_NODE_QUADRATIC_TRIANGLE);
    obj.CastForCell (&cell);

    BOOST_CHECK_EQUAL(obj.GetNumberOfPoints (), 6);
    BOOST_CHECK_EQUAL(obj.GetCellType (), VTK_CELL_TYPE::VTK_QUADRATIC_TRIANGLE);

    Point k0(0, 0, 0);
    Point p00 = obj.TransformRefToEle (&k0);
    BOOST_CHECK_EQUAL(p0, p00);

    Point k1(1, 0, 0);
    Point p11 = obj.TransformRefToEle (&k1);
    BOOST_CHECK_EQUAL(p1, p11);

    Point k2(0, 1, 0);
    Point p22 = obj.TransformRefToEle (&k2);
    BOOST_CHECK_EQUAL(p2, p22);

    Point k3(0.5, 0, 0);
    Point p33 = obj.TransformRefToEle (&k3);
    BOOST_CHECK_EQUAL(p3, p33);

    Point k4(0.5, 0.5, 0);
    Point p44 = obj.TransformRefToEle (&k4);
    BOOST_CHECK_EQUAL(p4, p44);

    Point k5(0, 0.5, 0);
    Point p55 = obj.TransformRefToEle (&k5);
    BOOST_CHECK_EQUAL(p5, p55);

    FELocalInfos loc;
    obj.LocalCompute (&k0, &loc);
    BOOST_CHECK_EQUAL(std::isinf (loc.detJac), false);

    double vol = 0.5 * CrossProduct (p1-p0, p2-p0).EuclidianNorm ();
    BOOST_CHECK_EQUAL (loc.detJac, static_cast<double>(obj.GetDimension ()) * vol);

    Point kpt1(0.5, 0.25, 0.);
    Point pt1 = obj.TransformRefToEle (&kpt1);
    Point kpt11 = obj.TransformEleToRef (&pt1);

    BOOST_CHECK_EQUAL(kpt1, kpt11);
}

BOOST_AUTO_TEST_SUITE_END()
