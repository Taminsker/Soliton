#ifdef STAND_ALONE
#define BOOST_TEST_MODULE SolitonTests_LagrangeElements
#endif

#include <boost/test/unit_test.hpp>

#include "Core/core.hpp"
#include "Solver/solver.hpp"

#define BOOST_CHECK_EQUAL_MESSAGE(L, R, M) \
    BOOST_TEST_INFO (M);                   \
    BOOST_CHECK_EQUAL (L, R)

#define BOOST_WARN_EQUAL_MESSAGE(L, R, M) \
    BOOST_TEST_INFO (M);                  \
    BOOST_WARN_EQUAL (L, R)

BOOST_AUTO_TEST_SUITE (Emp0N0DDL_felag)

BOOST_AUTO_TEST_CASE (Emp0N0DDL_constructor)
{
    BOOST_TEST_MESSAGE ("-- Tests for Emp0N0DDL");
    FELagrange<PHYSICAL_CELL_TYPE::EMPTY, 0, 0> obj;

    BOOST_CHECK_EQUAL (obj.GetNumberOfPoints (), 0);
    BOOST_CHECK_EQUAL (int (obj.GetCellType ()), int (VTK_CELL_TYPE::VTK_EMPTY_CELL));
}

BOOST_AUTO_TEST_SUITE_END ()

BOOST_AUTO_TEST_SUITE (Ver1N1DDL_felag)

BOOST_AUTO_TEST_CASE (Ver1N1DDL_constructor)
{
    BOOST_TEST_MESSAGE ("-- Tests for Ver1N1DDL");
    FELagrange<PHYSICAL_CELL_TYPE::VERTEX, 1, 1> obj;

    Point p0 (1, 1, 1);
    Cell  cell;
    cell.AddPoint (&p0);
    cell.SetType (GMSH_CELL_TYPE::GMSH_1_NODE_POINT);
    obj.CastForCell (&cell);

    BOOST_CHECK_EQUAL (obj.GetNumberOfPoints (), 1);
    BOOST_CHECK_EQUAL (int (obj.GetCellType ()), int (VTK_CELL_TYPE::VTK_VERTEX));

    Point k0 (0, 0, 0);
    Point p00 = obj.TransformRefToEle (&k0);
    BOOST_CHECK_EQUAL (p0, p00);

    FELocalInfos loc;
    obj.LocalCompute (&k0, &loc);
    BOOST_CHECK_EQUAL (std::isinf (loc.detJac), false);
}

BOOST_AUTO_TEST_SUITE_END ()

BOOST_AUTO_TEST_SUITE (Lin2N2DDL_felag)

BOOST_AUTO_TEST_CASE (Lin2N2DDL_constructor)
{
    BOOST_TEST_MESSAGE ("-- Tests for Lin2N2DDLL");
    FELagrange<PHYSICAL_CELL_TYPE::LINE, 2, 2> obj;

    Point p0 (-0.25, 0, 0);
    Point p1 (0.25, 0, 0);
    Cell  cell;
    cell.AddPoint (&p0);
    cell.AddPoint (&p1);
    cell.SetType (GMSH_CELL_TYPE::GMSH_2_NODE_LINE);
    obj.CastForCell (&cell);

    BOOST_CHECK_EQUAL (obj.GetNumberOfPoints (), 2);
    BOOST_CHECK_EQUAL (int (obj.GetCellType ()), int (VTK_CELL_TYPE::VTK_LINE));

    Point k0 (-1, 0, 0);
    Point k1 (1, 0, 0);
    Point p00 = obj.TransformRefToEle (&k0);
    BOOST_CHECK_EQUAL (p0, p00);

    Point p11 = obj.TransformRefToEle (&k1);
    BOOST_CHECK_EQUAL (p1, p11);

    FELocalInfos loc;
    obj.LocalCompute (&k0, &loc);
    real_t vol = (p1 - p0).EuclidianNorm ();
    BOOST_CHECK_EQUAL (std::isinf (loc.detJac), false);
    BOOST_CHECK_EQUAL (loc.detJac, 0.5 * vol);

    Point kpt1 (0.5, 0, 0);
    Point pt1   = obj.TransformRefToEle (&kpt1);
    Point kpt11 = obj.TransformEleToRef (&pt1);

    BOOST_CHECK_EQUAL (kpt1, kpt11);
}

BOOST_AUTO_TEST_SUITE_END ()

BOOST_AUTO_TEST_SUITE (Lin3N3DDL_felag)

BOOST_AUTO_TEST_CASE (Lin3N3DDL_constructor)
{
    BOOST_TEST_MESSAGE ("-- Tests for Lin3N3DDLL");
    FELagrange<PHYSICAL_CELL_TYPE::LINE, 3, 3> obj;

    Point p0 (0, 0, 0);
    Point p1 (4.5, 0, 0);
    Point p2 (2.3, 0, 0);
    Cell  cell;
    cell.AddPoint (&p0);
    cell.AddPoint (&p1);
    cell.AddPoint (&p2);
    cell.SetType (GMSH_CELL_TYPE::GMSH_3_NODE_QUADRATIC_LINE);
    obj.CastForCell (&cell);

    BOOST_CHECK_EQUAL (obj.GetNumberOfPoints (), 3);
    BOOST_CHECK_EQUAL (int (obj.GetCellType ()), int (VTK_CELL_TYPE::VTK_QUADRATIC_EDGE));

    Point k0 (-1, 0, 0);
    Point p00 = obj.TransformRefToEle (&k0);
    BOOST_CHECK_EQUAL (p0, p00);

    Point k1 (1, 0, 0);
    Point p11 = obj.TransformRefToEle (&k1);
    BOOST_CHECK_EQUAL (p1, p11);

    Point k2 (0, 0, 0);
    Point p22 = obj.TransformRefToEle (&k2);
    BOOST_CHECK_EQUAL (p2, p22);

    FELocalInfos loc;
    obj.LocalCompute (&k0, &loc);
    //    real_t vol = (p1 - p0).EuclidianNorm ();
    BOOST_CHECK_EQUAL (std::isinf (loc.detJac), false);
    //  BOOST_CHECK_EQUAL (detJac, 0.5 * vol);

    Point kpt1 (0.5, 0, 0);
    Point pt1   = obj.TransformRefToEle (&kpt1);
    Point kpt11 = obj.TransformEleToRef (&pt1);

    BOOST_CHECK_EQUAL (kpt1, kpt11);
}

BOOST_AUTO_TEST_SUITE_END ()

BOOST_AUTO_TEST_SUITE (Tri3N3DDL_felag)

BOOST_AUTO_TEST_CASE (Tri3N3DDL_constructor)
{
    BOOST_TEST_MESSAGE ("-- Tests for Tri3N3DDL");
    Tri3N3DDL obj;

    Point p0 (1., 1., 0);
    Point p1 (2.5, 1., 0);
    Point p2 (1.5, 1.5, 0);
    Cell  cell;
    cell.AddPoint (&p0);
    cell.AddPoint (&p1);
    cell.AddPoint (&p2);
    cell.SetType (GMSH_CELL_TYPE::GMSH_3_NODE_TRIANGLE);

    obj.CastForCell (&cell);

    {
        BOOST_CHECK_EQUAL (obj.GetNumberOfPoints (), 3);
        BOOST_CHECK_EQUAL (int (obj.GetCellType ()), int (VTK_CELL_TYPE::VTK_TRIANGLE));

        Point k0 (0, 0, 0);
        Point p00 = obj.TransformRefToEle (&k0);
        BOOST_CHECK_EQUAL (p0, p00);

        Point k1 (1, 0, 0);
        Point p11 = obj.TransformRefToEle (&k1);
        BOOST_CHECK_EQUAL (p1, p11);

        Point k2 (0, 1, 0);
        Point p22 = obj.TransformRefToEle (&k2);
        BOOST_CHECK_EQUAL (p2, p22);

        FELocalInfos loc;
        obj.LocalCompute (&k0, &loc);
        BOOST_CHECK_EQUAL (std::isinf (loc.detJac), false);

        real_t vol = 0.5 * CrossProduct (p1 - p0, p2 - p0).EuclidianNorm ();

        //    std::cout << loc.detJac << std::endl;
        //    std::cout << obj.GetDimension () * vol << std::endl;

        BOOST_CHECK_EQUAL (std::abs (loc.detJac - static_cast<real_t> (obj.GetDimension ()) * vol) < 1e-7, true);

        Point kpt1 (0.5, 0.25, 0.);
        Point pt1   = obj.TransformRefToEle (&kpt1);
        Point kpt11 = obj.TransformEleToRef (&pt1);

        BOOST_CHECK_EQUAL (kpt1, kpt11);
    }
    {
        BOOST_TEST_CHECKPOINT ("Elementary matrix");
        FEStore   festore;
        QuadStore quadstore;
        FEBase *  febase = festore.GetElementFor (&cell);
        febase->CastForCell (&cell);
        QuadStore::QuadObject * quadobject = quadstore.Get (febase);

        // compute
        Matrix3x3 mat;
        mat.setZero ();
        FELocalInfos loc;

        for (ul_t i = 0; i < 3; ++i)
        {
            //      Point* p_i = cell.GetPoints ()->at (i);
            std::function<Point (Point *)> grad_phi_i = febase->GetGradPhi (i);

            for (ul_t j = 0; j < 3; ++j)
            {
                //        Point* p_j = cell.GetPoints ()->at (j);
                std::function<Point (Point *)> grad_phi_j = febase->GetGradPhi (j);

                real_t coeff = 0.;

                for (ul_t k = 0; k < quadobject->npts; ++k)
                {
                    Point  pk = quadobject->pts [k];
                    real_t wk = quadobject->w [k];

                    febase->LocalCompute (&pk, &loc);
                    //          std::cout << "JacInvT : [" << i << ", " << j << "]" << std::endl << JacInvT << std::endl;

                    coeff += loc.detJac * wk * ((loc.JacInvT * grad_phi_i (&pk)) | (loc.JacInvT * grad_phi_j (&pk)));
                }

                mat.coeffRef (int (i), int (j)) = coeff;
            }
        }

        // real
        Matrix3x3 matreal;
        matreal.setZero ();
        Eigen::Matrix<real_t, 2, 2> middle;
        middle.col (0) << (p1 - p0).x, (p1 - p0).y;
        middle.col (1) << (p2 - p0).x, (p2 - p0).y;

        real_t vol = 0.5 * CrossProduct (p1 - p0, p2 - p0).EuclidianNorm ();

        auto temp = middle.inverse ().transpose ();
        middle    = temp;

        for (ul_t i = 0; i < 3; ++i)
        {
            Eigen::Vector2d grad_phi_i;
            grad_phi_i << febase->GetGradPhi (i) (&p0).x, febase->GetGradPhi (i) (&p0).y;
            grad_phi_i = middle * grad_phi_i;

            for (ul_t j = 0; j < 3; ++j)
            {
                Eigen::Vector2d grad_phi_j;
                grad_phi_j << febase->GetGradPhi (j) (&p0).x, febase->GetGradPhi (j) (&p0).y;
                grad_phi_j = middle * grad_phi_j;

                matreal.coeffRef (int (i), int (j)) = vol * grad_phi_i.adjoint () * grad_phi_j;
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END ()

BOOST_AUTO_TEST_SUITE (Tri6N6DDL_felag)

BOOST_AUTO_TEST_CASE (Tri6N6DDL_constructor)
{
    BOOST_TEST_MESSAGE ("-- Tests for Tri6N6DDL");
    FELagrange<PHYSICAL_CELL_TYPE::TRIANGLE, 6, 6> obj;

    Point p0 (0, 0, 0);
    Point p1 (1, 0, 0);
    Point p2 (0, 1, 0);
    Point p3 (0.5, 0, 0);
    Point p4 (0.5, 0.5, 0);
    Point p5 (0, 0.5, 0);
    Cell  cell;
    cell.AddPoint (&p0);
    cell.AddPoint (&p1);
    cell.AddPoint (&p2);
    cell.AddPoint (&p3);
    cell.AddPoint (&p4);
    cell.AddPoint (&p5);
    cell.SetType (GMSH_CELL_TYPE::GMSH_6_NODE_QUADRATIC_TRIANGLE);
    obj.CastForCell (&cell);

    BOOST_CHECK_EQUAL (obj.GetNumberOfPoints (), 6);
    BOOST_CHECK_EQUAL (int (obj.GetCellType ()), int (VTK_CELL_TYPE::VTK_QUADRATIC_TRIANGLE));

    Point k0 (0, 0, 0);
    Point p00 = obj.TransformRefToEle (&k0);
    BOOST_CHECK_EQUAL (p0, p00);

    Point k1 (1, 0, 0);
    Point p11 = obj.TransformRefToEle (&k1);
    BOOST_CHECK_EQUAL (p1, p11);

    Point k2 (0, 1, 0);
    Point p22 = obj.TransformRefToEle (&k2);
    BOOST_CHECK_EQUAL (p2, p22);

    Point k3 (0.5, 0, 0);
    Point p33 = obj.TransformRefToEle (&k3);
    BOOST_CHECK_EQUAL (p3, p33);

    Point k4 (0.5, 0.5, 0);
    Point p44 = obj.TransformRefToEle (&k4);
    BOOST_CHECK_EQUAL (p4, p44);

    Point k5 (0, 0.5, 0);
    Point p55 = obj.TransformRefToEle (&k5);
    BOOST_CHECK_EQUAL (p5, p55);

    FELocalInfos loc;
    obj.LocalCompute (&k0, &loc);
    BOOST_CHECK_EQUAL (std::isinf (loc.detJac), false);

    real_t vol = 0.5 * CrossProduct (p1 - p0, p2 - p0).EuclidianNorm ();
    BOOST_CHECK_EQUAL (loc.detJac, static_cast<real_t> (obj.GetDimension ()) * vol);

    Point kpt1 (0.5, 0.25, 0.);
    Point pt1   = obj.TransformRefToEle (&kpt1);
    Point kpt11 = obj.TransformEleToRef (&pt1);

    BOOST_CHECK_EQUAL (kpt1, kpt11);
}

BOOST_AUTO_TEST_SUITE_END ()

BOOST_AUTO_TEST_SUITE (Quad4N4DDL_felag)

BOOST_AUTO_TEST_CASE (Quad4N4DDL_constructor)
{
    BOOST_TEST_MESSAGE ("-- Tests for Quad4N4DDL");

    constexpr PHYSICAL_CELL_TYPE physicalType = PHYSICAL_CELL_TYPE::QUADRANGLE;
    constexpr GMSH_CELL_TYPE     gmshType     = GMSH_CELL_TYPE::GMSH_4_NODE_QUADRANGLE;
    //  VTK_CELL_TYPE vtkType = Convert (gmshType);
    constexpr int npts = 4;
    constexpr int ddl  = 4;

    FELagrange<physicalType, npts, ddl> obj;

    Point p0 (0, 0, 0);
    Point p1 (1, 0, 0);
    Point p2 (1, 1, 0);
    Point p3 (0, 1, 0);

    real_t vol = 0.5 * CrossProduct (p1 - p0, p2 - p0).EuclidianNorm ();
    vol += 0.5 * CrossProduct (p1 - p3, p2 - p3).EuclidianNorm ();

    BOOST_WARN_EQUAL_MESSAGE (std::abs (vol) < EPSILON, false, "WARNING volume of element is null");

    std::vector<Point *> listPoints = {&p0, &p1, &p2, &p3};
    BOOST_CHECK_EQUAL_MESSAGE (listPoints.size (), npts, "test");

    Cell cell;
    cell.AddPoints (listPoints)->SetType (gmshType);
    obj.CastForCell (&cell);

    BOOST_CHECK_EQUAL_MESSAGE (obj.GetNumberOfPoints (), listPoints.size (), "Misconfigured number of points");
    VTK_CELL_TYPE vtkType = Convert<GMSH_CELL_TYPE, VTK_CELL_TYPE> (gmshType);
    BOOST_CHECK_EQUAL (int (obj.GetCellType ()), int (vtkType));

    FELocalInfos loc;
    for (ul_t idPoint = 0; idPoint < obj.GetNumberOfPoints (); ++idPoint)
    {
        Point pt = obj.TransformRefToEle (obj.GetPoint (idPoint));
        BOOST_CHECK_EQUAL_MESSAGE (pt, *listPoints.at (idPoint), "The point " + std::to_string (idPoint) + " is misconfigured (TransformRefToEle)");
        obj.LocalCompute (obj.GetPoint (idPoint), &loc);
        BOOST_CHECK_EQUAL_MESSAGE (std::isinf (loc.detJac), false, "The point " + std::to_string (idPoint) + " is misconfigured (detJac == INF)");
    }

    //  obj.LocalCompute (obj.GetPoint (0), &loc);
    //  BOOST_CHECK_EQUAL (loc.detJac, static_cast<real_t>(obj.GetDimension ()) * vol);
}

BOOST_AUTO_TEST_SUITE_END ()

BOOST_AUTO_TEST_SUITE (Quad8N8DDL_felag)

BOOST_AUTO_TEST_CASE (Quad8N8DDL_constructor)
{
    BOOST_TEST_MESSAGE ("-- Tests for Quad8N8DDL");

    constexpr PHYSICAL_CELL_TYPE physicalType = PHYSICAL_CELL_TYPE::QUADRANGLE;
    constexpr GMSH_CELL_TYPE     gmshType     = GMSH_CELL_TYPE::GMSH_8_NODE_QUADRATIC_QUADRANGLE;
    //  VTK_CELL_TYPE vtkType = Convert (gmshType);
    constexpr int npts = 8;
    constexpr int ddl  = 8;

    FELagrange<physicalType, npts, ddl> obj;

    Point p0 (0, 0, 0);
    Point p1 (1, 0, 0);
    Point p2 (1, 1, 0);
    Point p3 (0, 1, 0);
    Point p4 (0.5, 0, 0);
    Point p5 (1, 0.5, 0);
    Point p6 (0.5, 1, 0);
    Point p7 (0, 0.5, 0);

    real_t vol = 0.5 * CrossProduct (p1 - p0, p2 - p0).EuclidianNorm ();
    vol += 0.5 * CrossProduct (p1 - p3, p2 - p3).EuclidianNorm ();

    BOOST_WARN_EQUAL_MESSAGE (std::abs (vol) < EPSILON, false, "WARNING volume of element is null");

    std::vector<Point *> listPoints = {&p0, &p1, &p2, &p3, &p4, &p5, &p6, &p7};
    BOOST_CHECK_EQUAL_MESSAGE (listPoints.size (), npts, "test");

    Cell cell;
    cell.AddPoints (listPoints)->SetType (gmshType);
    obj.CastForCell (&cell);

    BOOST_CHECK_EQUAL_MESSAGE (obj.GetNumberOfPoints (), listPoints.size (), "Misconfigured number of points");
    VTK_CELL_TYPE vtkType = Convert<GMSH_CELL_TYPE, VTK_CELL_TYPE> (gmshType);
    BOOST_CHECK_EQUAL (int (obj.GetCellType ()), int (vtkType));

    FELocalInfos loc;
    for (ul_t idPoint = 0; idPoint < obj.GetNumberOfPoints (); ++idPoint)
    {
        Point pt = obj.TransformRefToEle (obj.GetPoint (idPoint));
        BOOST_CHECK_EQUAL_MESSAGE (pt, *listPoints.at (idPoint), "The point " + std::to_string (idPoint) + " is misconfigured (TransformRefToEle)");
        obj.LocalCompute (obj.GetPoint (idPoint), &loc);
        BOOST_CHECK_EQUAL_MESSAGE (std::isinf (loc.detJac), false, "The point " + std::to_string (idPoint) + " is misconfigured (detJac == INF)");
    }

    //  obj.LocalCompute (obj.GetPoint (0), &loc);
    //  BOOST_CHECK_EQUAL (loc.detJac, static_cast<real_t>(obj.GetDimension ()) * vol);
}

BOOST_AUTO_TEST_SUITE_END ()

BOOST_AUTO_TEST_SUITE (Quad9N9DDL_felag)

BOOST_AUTO_TEST_CASE (Quad9N9DDL_constructor)
{
    BOOST_TEST_MESSAGE ("-- Tests for Quad9N9DDL");

    constexpr PHYSICAL_CELL_TYPE physicalType = PHYSICAL_CELL_TYPE::QUADRANGLE;
    constexpr GMSH_CELL_TYPE     gmshType     = GMSH_CELL_TYPE::GMSH_9_NODE_QUADRATIC_QUADRANGLE;
    //  VTK_CELL_TYPE vtkType = Convert (gmshType);
    constexpr int npts = 9;
    constexpr int ddl  = 9;

    FELagrange<physicalType, npts, ddl> obj;

    Point p0 (0, 0, 0);
    Point p1 (1, 0, 0);
    Point p2 (1, 1, 0);
    Point p3 (0, 1, 0);
    Point p4 (0.5, 0, 0);
    Point p5 (1, 0.5, 0);
    Point p6 (0.5, 1, 0);
    Point p7 (0, 0.5, 0);
    Point p8 (0.5, 0.5, 0);

    real_t vol = 0.5 * CrossProduct (p1 - p0, p2 - p0).EuclidianNorm ();
    vol += 0.5 * CrossProduct (p1 - p3, p2 - p3).EuclidianNorm ();

    BOOST_WARN_EQUAL_MESSAGE (std::abs (vol) < EPSILON, false, "WARNING volume of element is null");

    std::vector<Point *> listPoints = {&p0, &p1, &p2, &p3, &p4, &p5, &p6, &p7, &p8};
    BOOST_CHECK_EQUAL_MESSAGE (listPoints.size (), npts, "test");

    Cell cell;
    cell.AddPoints (listPoints)->SetType (gmshType);
    obj.CastForCell (&cell);

    BOOST_CHECK_EQUAL_MESSAGE (obj.GetNumberOfPoints (), listPoints.size (), "Misconfigured number of points");
    VTK_CELL_TYPE vtkType = Convert<GMSH_CELL_TYPE, VTK_CELL_TYPE> (gmshType);
    BOOST_CHECK_EQUAL (int (obj.GetCellType ()), int (vtkType));

    FELocalInfos loc;
    for (ul_t idPoint = 0; idPoint < obj.GetNumberOfPoints (); ++idPoint)
    {
        Point pt = obj.TransformRefToEle (obj.GetPoint (idPoint));
        BOOST_CHECK_EQUAL_MESSAGE (pt, *listPoints.at (idPoint), "The point " + std::to_string (idPoint) + " is misconfigured (TransformRefToEle)");
        obj.LocalCompute (obj.GetPoint (idPoint), &loc);
        BOOST_CHECK_EQUAL_MESSAGE (std::isinf (loc.detJac), false, "The point " + std::to_string (idPoint) + " is misconfigured (detJac == INF)");
    }

    //  obj.LocalCompute (obj.GetPoint (0), &loc);
    //  BOOST_CHECK_EQUAL (loc.detJac, static_cast<real_t>(obj.GetDimension ()) * vol);
}

BOOST_AUTO_TEST_SUITE_END ()
