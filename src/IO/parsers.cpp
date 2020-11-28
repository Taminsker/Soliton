#include "parsers.h"

#include <Algorithms/algorithms.h>
#include <Core/core.h>
#include <Solver/solver.h>

#define ERRORPARAMS \
{\
    ERROR << "unknown parameter " << field << " = " << value << ENDLINE;\
    errs++;\
    count--;\
    continue;\
    }

void ParseInputDatFile (InputDatStruct* out, std::string filename)
{
    BEGIN << "Parse the input dat file." << ENDLINE;

    std::ifstream infile (filename);
    std::string line, field, value;
    std::vector <std::string> spl;
    unsigned int errs = 0;
    unsigned int count = 0;

    if (!infile.is_open ())
    {
        ERROR << "the file " << filename << " can not be open." << BLINKRETURN << ENDLINE;
        return;
    } else
    {
        STATUS << "the file " << COLOR_BLUE << filename << COLOR_DEFAULT << " is open." << ENDLINE;
    }

    //  out->filename_msh = filename;

    while (std::getline(infile, line))
    {
        if (line.front () != '@' && line.size () != 0)
        {
            RemoveBlankSpace (&line, &spl);

            if (spl.size () == 0)
                continue;

            if (spl.size () >= 3)
            {
                count++;

                field = spl.at(0);
                value = spl.at (2);

                if (field == "file_msh")
                    out->filename_msh = value;
                else if (field == "grid_x_m")
                    out->grid_x_m = stod(value);
                else if (field == "grid_x_p")
                    out->grid_x_p = stod(value);
                else if (field == "grid_y_m")
                    out->grid_y_m = stod(value);
                else if (field == "grid_y_p")
                    out->grid_y_p = stod(value);
                else if (field == "hsize")
                    out->hsize = stod(value);
                else if (field == "ele_type")
                    out->ele_type = stoi(value);
                else if (field == "ele_order")
                    out->ele_order = stoi(value);
                else if (field == "damping")
                    out->damping = CastToBool (&value);
                else if (field == "zeta_0")
                    out->zeta_0 = stod(value);
                else if (field == "beta_0")
                    out->beta_0 = stod(value);
                else if (field == "g")
                    out->g_variable = stod(value);
                else if (field == "colContItem")
                    out->colConcItem = CastToBool (&value);
                else if (field == "dt")
                    out->dt = stod(value);
                else if (field == "object_policy")
                    out->objectsAreFixed = CastToBool (&value);
                else if (field == "coeff_penalization")
                    out->penal = stod(value);
                else if (field == "power_penalization")
                    out->powpenalty = stoi(value);
                else if (field == "synthetize")
                    out->synthetize = CastToBool (&value);
                else
                    ERRORPARAMS
            }

            field = spl.at(0);

            if (field == "$BEGIN_OBJECT")
            {
                ObjectDatStruct obs;

                while (std::getline(infile, line))
                {
                    if (line.empty ())
                        continue;

                    RemoveBlankSpace (&line, &spl);

                    if (spl.size () == 0)
                        continue;

                    if (spl[0] == "$END_OBJECT")
                        break;

                    if (infile.eof ())
                    {
                        ERROR << "not enabled to find $END_OBJECT..." << BLINKRETURN << ENDLINE;
                        return;
                    }

                    if (line.front () != '@' && line.size () != 0)
                    {
                        RemoveBlankSpace (&line, &spl);

                        if (spl.size () >= 3)
                        {
                            count++;

                            field = spl.at(0);
                            value = spl.at (2);

                            if (field == "file_msh")
                                obs.filename_msh = value;
                            else if (field == "algo_gen")
                                obs.algo_gen = static_cast<unsigned int>(stod(value));
                            else if (field == "nbpts")
                                obs.nbpts = stoi(value);
                            else if (field == "rinp")
                                obs.rinp = CastToBool (&value);
                            else if (field == "x_center")
                                obs.x_center = stod(value);
                            else if (field == "y_center")
                                obs.y_center = stod(value);
                            else if (field == "z_center")
                                obs.z_center = stod(value);
                            else if (field == "basename")
                                obs.basename = value;
                            else if (field == "radius")
                                obs.radius = stod(value);
                            else
                                ERRORPARAMS
                        }
                    }
                }

                out->objects.push_back (obs);
            }
        }
    }

    STATUS << "EOF : " << COLOR_BLUE << count << " success" << COLOR_DEFAULT << " and " << COLOR_RED << errs << " errors." << ENDLINE;

    ENDFUN;

#ifdef DEBUG
    Print (out);
#endif
    return;
}

void ParseMSH (Mesh* mesh, std::string filename, bool keep_original)
{
    BEGIN << "Parse the msh file." << ENDLINE;

    if (mesh == nullptr)
    {
        ERROR << "the store object is not correct initialised." << BLINKRETURN << ENDLINE;
        return;
    }

    std::ifstream infile (filename);
    std::string line;
    std::string delimeter = "\"";
    std::vector <std::string> sv;
    int numpoints = 0;
    int numcells = 0;
    int err = 0;
    int count = 0;
    int depth_tag = 0;

    // OPEN
    if (!infile.is_open ())
    {
        ERROR << "the file " << filename << " can not be open." << BLINKRETURN << ENDLINE;
        return;
    } else
    {
        STATUS << "the file " << filename << " is " << COLOR_BLUE << "open." << ENDLINE;
    }

    // MESH FORMAT

    std::getline(infile, line);
    while (line != "$MeshFormat")
    {
        if (infile.eof ())
        {
            ERROR << "The $MeshFormat field does not seem to be present" << BLINKRETURN << ENDLINE;
            return;
        }

        std::getline(infile, line);
    }

#ifdef DEBUG
    STATUS << "$MeshFormat detected." << ENDLINE;
#endif
    std::getline(infile, line);
    RemoveBlankSpace (&line, &sv);

    if (sv.at (0) != "2.2")
    {
        ERROR << "the msh file format must be 2.2 format" << BLINKRETURN << ENDLINE;
        return;
    }

    // NODES
    while (line != "$Nodes")
    {
        if (infile.eof ())
        {
            ERROR << "the $Nodes field does not seem to be present" << BLINKRETURN << ENDLINE;
            return;
        }

        std::getline(infile, line);
    }

#ifdef DEBUG
    STATUS << "$Nodes detected." << ENDLINE;
#endif

    std::getline(infile, line);
    RemoveBlankSpace (&line, &sv);

    numpoints = stoi(sv.at (0));

    count = 0;
    while (count < 10 * numpoints)
    {
        if (infile.eof ())
        {
            ERROR << "the $EndNodes field does not seem to be present" << BLINKRETURN << ENDLINE;
            return;
        }

        std::getline(infile, line);
        RemoveBlankSpace (&line, &sv);

        if (line == "$EndNodes")
            break;

        if (sv.size () < 4)
        {
            ERROR << "one point doesn't respect the rule : id x y z, please check your msh file." << BLINKRETURN << ENDLINE;
            return;
        }

        auto P = new Point ();
        P->x = stod (sv.at (1));
        P->y = stod (sv.at (2));
        P->z = stod (sv.at (3));
        P->SetGlobalIndex (count);

        mesh->AddPoint (P);

        count++;
    }

    if (count != numpoints)
    {
        ERROR << "the number of points does not coincide." << BLINKRETURN << ENDLINE;
        return;
    }

#ifdef DEBUG
    STATUS << "$EndNodes detected." << ENDLINE;
#endif

    INFOS << "nodes : \t" << count << ENDLINE;

    // ELEMENTS
    while (line != "$Elements")
    {
        if (infile.eof ())
        {
            ERROR << "the $Elements field does not seem to be present" << BLINKRETURN << ENDLINE;
            return;
        }
        std::getline(infile, line);
    }

#ifdef DEBUG
    STATUS << "$Elements detected." << ENDLINE;
#endif

    std::getline(infile, line);
    RemoveBlankSpace (&line, &sv);

    numcells = stoi(sv.at (0));

    count = 0;
    while (count < 10 * numcells)
    {
        if (infile.eof ())
        {
            ERROR << "the $EndElements field does not seem to be present" << BLINKRETURN << ENDLINE;
            return;
        }

        std::getline(infile, line);
        RemoveBlankSpace (&line, &sv);

        if (line == "$EndElements")
            break;

        if (sv.size () < 3)
        {
            ERROR << "One cell doesn't respect the rule : id type ntag <tag> <points>, please check your msh file." << BLINKRETURN << ENDLINE;
            return;
        }

        auto C = new Cell ();
        C->SetType (GMSH_CELL_TYPE (stoi (sv.at (1))));
        C->SetGlobalIndex (count);

        depth_tag = 3 + stoi (sv.at (2));

        for (ul_t i = ul_t(depth_tag); i < sv.size (); ++i)
            C->AddPoint (mesh->GetPoint (stoi(sv.at (i)) - 1));

        mesh->AddCell (C);

        count++;
    }

    if (count != numcells)
    {
        ERROR << "The number of cells does not coincide." << BLINKRETURN << ENDLINE;
        return;
    }

#ifdef DEBUG
    STATUS << "$EndElements detected." << ENDLINE;
#endif

    INFOS << "elements : \t" << count << ENDLINE;
    // END READ

    if (!keep_original)
    {
        err = system( ("rm " + filename).c_str ());

        if (err != 0)
            ERROR << "can not delete the msh file..." << ENDLINE;
    }

    ENDFUN;

    //  DeleteLines (mesh);
    Build_NtoN (mesh);
    BuildEdgesWithHashMap (mesh);

    return;
}


void Print (InputDatStruct* input)
{
    INFOS << "dt    = " << input->dt << ENDLINE;
    INFOS << "zeta_0  = " << input->zeta_0<< ENDLINE;
    INFOS << "beta_0  = " << input->beta_0<< ENDLINE;

    INFOS << "file   = " << input->filename_msh << ENDLINE;
    INFOS << "xm    = " << input->grid_x_m << ENDLINE;
    INFOS << "xp    = " << input->grid_x_p << ENDLINE;
    INFOS << "ym    = " << input->grid_y_m << ENDLINE;
    INFOS << "yp    = " << input->grid_y_p << ENDLINE;
    INFOS << "hsize  = " << input->hsize << ENDLINE;
    INFOS << "eltype  = " << input->ele_type << ENDLINE;
    INFOS << "elorder = " << input->ele_order << ENDLINE;

    INFOS << "obj   = " << input->objectsAreFixed << ENDLINE;

    INFOS << "penalization_coeff = " << input->penal << ENDLINE;
    INFOS << "penalization_power = " << input->powpenalty << ENDLINE;
    INFOS << "colConcItem     = " << input->colConcItem << ENDLINE;
    INFOS << "synthetize     = " << input->synthetize << ENDLINE;


    for (ul_t i = 0; i < input->objects.size (); ++i)
    {
        ObjectDatStruct obj = input->objects.at (i);
        INFOS << SEPARATOR << ENDLINE;
        INFOS << "object" << ENDLINE;
        INFOS << "\txc    = " << obj.x_center << ENDLINE;
        INFOS << "\tyc    = " << obj.y_center << ENDLINE;
        INFOS << "\tr    = " << obj.radius << ENDLINE;
        INFOS << "\talgo   = " << obj.algo_gen << ENDLINE;
        INFOS << "\trinp    = " << obj.rinp << ENDLINE;
        INFOS << "\tnbpts  = " << obj.nbpts << ENDLINE;
        INFOS << "\tfile_msh = " << obj.filename_msh << ENDLINE;
    }

    ENDFUN;

    return;
}

void split (std::string* s, std::string* delimiter, std::vector<std::string>* out)
{
    ul_t pos_start = 0, pos_end, delim_len = delimiter->length();
    std::string token;
    out->clear ();

    while ((pos_end = s->find (*delimiter, pos_start)) != std::string::npos) {
        token = s->substr (pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        out->push_back (token);
    }

    out->push_back (s->substr (pos_start));
    return;
}

void RemoveBlankSpace (std::string* s, std::vector<std::string>* out)
{
    std::istringstream iss(*s);
    *out = {std::istream_iterator<std::string>{iss},      std::istream_iterator<std::string>{}};

    return;
}

bool CastToBool(std::string *s)
{
    return (*s == "true" ? true : false);
}
