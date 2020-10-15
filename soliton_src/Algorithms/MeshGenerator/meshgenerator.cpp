#include "meshgenerator.h"
#include "../AlgoMesh/algomesh.h"
#include "../Hash4Edges/hash4edges.h"

#include <Core/core.h>
#include <IO/io.h>
#include <Solver/solver.h>

#include <fstream>

std::string GenerateWithGMSH (InputDatStruct* input)
{
    HEADERFUN("GenerateWithGMSH");
    BEGIN << "Generate a mesh with GMSH in msh file." << ENDLINE;

    int err = 0;
    double count = 0;
    std::string out = "temp_msh_file.msh";

    // Before
    count = (input->grid_x_p - input->grid_x_m) * (input->grid_y_p - input->grid_y_m) / (input->hsize * input->hsize);
    count *= 1.2;


    if (count >= 1e6)
    {
        STATUS << "GMSH can take a long time ... be patient." << ENDLINE;
        INFOS << COLOR_BLUE << "Prevision : points " << count << "\t cells : " << count * 2. << ENDLINE;
    }


    // GeoFile Generation

    std::ofstream geofile (".temp_geo_file.geo");

    //    geofile << "Include \"" << input->filename << "\";" << std::endl;

    geofile << "xm = " << input->grid_x_m << ";" << std::endl;
    geofile << "xp = " << input->grid_x_p << ";" << std::endl;
    geofile << "ym = " << input->grid_y_m << ";" << std::endl;
    geofile << "yp = " << input->grid_y_p << ";" << std::endl;
    geofile << "hsize = " << input->hsize << ";" << std::endl;

    geofile << "// point" << std::endl;
    geofile << "Point(1) = {xm, ym, 0, hsize};" << std::endl;
    geofile << "Point(2) = {xp, ym, 0, hsize};" << std::endl;
    geofile << "Point(3) = {xp, yp, 0, hsize};" << std::endl;
    geofile << "Point(4) = {xm, yp, 0, hsize};" << std::endl;

    geofile << "// build" << std::endl;
    geofile << "Line(1) = {4, 3};" << std::endl;
    geofile << "Line(2) = {2, 3};" << std::endl;
    geofile << "Line(3) = {4, 1};" << std::endl;
    geofile << "Line(4) = {1, 2};" << std::endl;
    geofile << "Curve Loop(1) = {1, -2, -4, -3};" << std::endl;
    geofile << "Plane Surface(1) = {1};" << std::endl;

//    geofile << "// physical" << std::endl;
//    geofile << "Physical Curve(\"wall\") = {1, 4};" << std::endl;
//    geofile << "Physical Curve(\"inlet\") = {2};" << std::endl;
//    geofile << "Physical Curve(\"outlet\") = {3};" << std::endl;
//    geofile << "Physical Surface(\"domain\") = {1};" << std::endl;

//    geofile << "Physical Curve(" << int(PHYS::WALL) << ") = {1, 4};" << std::endl;
//    geofile << "Physical Curve(" << int(PHYS::INLET) << ") = {2};" << std::endl;
//    geofile << "Physical Curve(" << int(PHYS::OUTLET) << ") = {3};" << std::endl;
    geofile << "Physical Surface(" << int(PHYS::DOMAIN) << ") = {1};" << std::endl;

    geofile << "// mesh generation" << std::endl;
    geofile << "//Mesh.SaveAll=0;" << std::endl;

    geofile << "Mesh.ElementOrder = " << input->ele_order << ";" << std::endl;
    if (input->ele_type == 1)
        geofile << "Mesh.RecombineAll = 1 ;" << std::endl;

    geofile << "Mesh.MshFileVersion = 2.2;" << std::endl;
    geofile << "Mesh 2;" << std::endl;
    geofile << "Save \"" << out << "\";" << std::endl;

    geofile.close ();

#ifdef DEBUG
    err = system("gmsh .temp_geo_file.geo");
#else
    err = system("gmsh .temp_geo_file.geo >> .log_gmsh.txt");
    err = system("rm .log_gmsh.txt");
#endif

    if (err != 0)
    {
        ERROR << "mesh is generated with GMSH : failed. Please try 'gmsh gmsh .temp_geo_file.geo' or look at the geo file. " << BLINKRETURN << ENDLINE;
        return "";
    } else
    {
        STATUS << "mesh is generated with GMSH : success." << ENDLINE;
    }

    err = system("rm .temp_geo_file.geo");

    ENDFUN;
    return out;
}


void ObjectGenerator (InputDatStruct* data, Sto4Sol* store)
{
    HEADERFUN("ObjectGenerator");

    if (store->mesh->GetNumberOfPoints () == 0)
    {
        ERROR << "the mesh is completely empty. Generate and Parse a Msh file before set the objects !" << BLINKRETURN << ENDLINE;
        return;
    }

    std::size_t count = data->objects.size ();

    for (std::size_t num = 0; num < count; ++num)
    {
        ObjectDatStruct obj = data->objects.at (num);
        Mesh* object = new Mesh();
        object->SetName (data->objects.at (num).basename);

        if (obj.filename_msh.empty ())
        {
            if (obj.algo_gen == 1)
                AlgoGen1 (&obj, object);
            else if (obj.algo_gen == 2)
                AlgoGen2 (&obj, object);
            else if (obj.algo_gen == 3)
                AlgoGen3 (&obj, object);
            else
                ERROR << "the object " << num << " is set with algo " << obj.algo_gen << ") is not supported yet... " << BLINK << "SKIP" << ENDLINE;
        }
        else
        {
            ParseMSH (object, obj.filename_msh, true);

            if (obj.rinp)
                MoveObject (object, obj.radius, {obj.x_center, obj.y_center, obj.z_center});
        }


        BuildEdgesWithHashMap(object);
        //            ComputeNormalsOnEdges (object);
        //            ComputeNormalsOnCells (object);
        //            ComputeNormalsOnPoints (object);


        store->listobjects.push_back (object);

    }

    return;
}

void AlgoGen1 (ObjectDatStruct* data, Mesh* object)
{
    HEADERFUN("AlgoGen1");
    BEGIN << "Generate the object with AlgoGen1 (circle radius) : " << COLOR_BLUE << object->GetName () << ENDLINE;

    if (data->radius <= 1e-6)
    {
        ERROR << "the object " << object->GetName () << " is set with algo " << data->algo_gen << "(radius generator) but radius is " << data->radius << ". Correct it please !" << BLINK << "SKIP" << ENDLINE;
    }

    //    if (    data->xp <= obj.xc || data->xm >= obj.xc ||
    //            data->yp <= obj.yc || data->ym >= obj.yc)
    //    {
    //        ERROR << "the object " << num << " is set with center [" << obj.xc << ", " << obj.yc << "] but it's not in the mesh window [" << data->xm << ", " << data->xp << "]x[" << data->ym << ", " << data->yp  << "]. Correct it please !" << BLINK << "SKIP" << ENDLINE;
    //    }

    double radmin = 2. * M_PI / double (data->nbpts);
    int countpt = 0;

    for (double rr = 0.; rr < 2. * M_PI; rr = rr + radmin)
    {
        Point* p = new Point();
        p->x = data->radius * std::cos (rr) + data->x_center;
        p->y = data->radius * std::sin (rr) + data->y_center;
        p->z = 0.0;
        p->SetGlobalIndex (countpt);

        object->AddPoint (p);
        countpt++;
    }

    countpt = object->GetNumberOfPoints ();

    for (int i = 0; i < countpt; ++i)
    {
        Cell* c = new Cell();
        c->SetType (GMSH_CELL_TYPE::GMSH_2_NODE_LINE);
        c->SetGlobalIndex (i);
        c->AddPoint (object->GetPoint (i));

        if (i == countpt - 1)
            c->AddPoint (object->GetPoint (0));
        else
            c->AddPoint (object->GetPoint (i+1));

        object->AddCell (c);
    }

    ENDFUN;
    return;
}

void AlgoGen2 (ObjectDatStruct* data, Mesh* object)
{
    HEADERFUN("AlgoGen2");
    BEGIN << "Generate the object with AlgoGen2 (puzzle piece) : " << COLOR_BLUE << object->GetName () << ENDLINE;

    if (data->radius <= 1e-6)
    {
        ERROR << "the object " << object->GetName () << " is set with algo " << data->algo_gen << "(radius generator) but radius is " << data->radius << ". Correct it please !" << BLINK << "SKIP" << ENDLINE;
    }

    //    if (    data->xp <= obj.xc || data->xm >= obj.xc ||
    //            data->yp <= obj.yc || data->ym >= obj.yc)
    //    {
    //        ERROR << "the object " << num << " is set with center [" << obj.xc << ", " << obj.yc << "] but it's not in the mesh window [" << data->xm << ", " << data->xp << "]x[" << data->ym << ", " << data->yp  << "]. Correct it please !" << BLINK << "SKIP" << ENDLINE;
    //    }

    double radmin = 2. * M_PI / double (data->nbpts);
    int countpt = 0;

    std::vector<Point*> Normals;

    //    for (double rr = 0.; rr < 2. * M_PI; rr += radmin)
    for (double rr = 2. * M_PI; rr > 0.; rr -= radmin)
    {
        Point* p = new Point();
        p->x = 0.6 * std::cos (rr) - 0.3 * std::cos (3. * rr);
        p->y = 0.7 * std::sin (rr) - 0.07 * std::sin (3. * rr) + 0.2 * std::sin (7. * rr);
        p->z = 0.0;

        *p = data->radius * *p + Point (data->x_center, data->y_center, 0.);
        p->SetGlobalIndex (countpt);

        object->AddPoint (p);
        countpt++;
    }

    countpt = object->GetNumberOfPoints ();

    int countcell = 0;
    for (int i = 0; i < countpt; i++)
    {
        Cell* c = new Cell();
        c->SetType (GMSH_CELL_TYPE::GMSH_2_NODE_LINE);
        c->SetGlobalIndex (countcell);
        c->AddPoint (object->GetPoint (i));

        if (i == countpt - 1)
            c->AddPoint (object->GetPoint (0));
        else
            c->AddPoint (object->GetPoint (i+1));

        object->AddCell (c);
        countcell++;
    }

    ENDFUN;
    return;
}

void AlgoGen3 (ObjectDatStruct* data, Mesh* object)
{
    HEADERFUN("AlgoGen2");
    BEGIN << "Generate the object with AlgoGen2 (puzzle piece) : " << COLOR_BLUE << object->GetName () << ENDLINE;

    if (data->radius <= 1e-6)
    {
        ERROR << "the object " << object->GetName () << " is set with algo " << data->algo_gen << "(radius generator) but radius is " << data->radius << ". Correct it please !" << BLINK << "SKIP" << ENDLINE;
    }

    //    if (    data->xp <= obj.xc || data->xm >= obj.xc ||
    //            data->yp <= obj.yc || data->ym >= obj.yc)
    //    {
    //        ERROR << "the object " << num << " is set with center [" << obj.xc << ", " << obj.yc << "] but it's not in the mesh window [" << data->xm << ", " << data->xp << "]x[" << data->ym << ", " << data->yp  << "]. Correct it please !" << BLINK << "SKIP" << ENDLINE;
    //    }

    double radmin = 2. * M_PI / double (data->nbpts);
    int countpt = 0;

    std::vector<Point*> Normals;

    for (double rr = 0.; rr < 2. * M_PI; rr = rr + radmin)
    {
        Point* p = new Point();
        p->x = 0.6 * std::cos (rr) - 0.3 * std::cos (3. * rr);
        p->y = 0.7 * std::sin (rr) - 0.07 * std::sin (3. * rr) + 0.2 * std::sin (7. * rr);
        p->z = 0.0;

        *p = data->radius * *p + Point (data->x_center, data->y_center, 0.);
        p->SetGlobalIndex (countpt);

        object->AddPoint (p);
        countpt++;
    }

    countpt = object->GetNumberOfPoints ();

    for (int i = 0; i < countpt; ++i)
    {
        Cell* c = new Cell();
        c->SetType (GMSH_CELL_TYPE::GMSH_2_NODE_LINE);
        c->SetGlobalIndex (i);
        c->AddPoint (object->GetPoint (i));

        if (i == countpt - 1)
            c->AddPoint (object->GetPoint (0));
        else
            c->AddPoint (object->GetPoint (i+1));

        object->AddCell (c);
    }

    ENDFUN;
    return;
}

