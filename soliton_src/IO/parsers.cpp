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
    HEADERFUN("ParseInputDatFile");
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
        STATUS << "the file " << filename << " is open." << ENDLINE;
    }

    out->filename = filename;

    while (std::getline(infile, line))
    {
        if (line.front () != '#' && line.size () != 0)
        {
            RemoveBlankSpace (&line, &spl);

            if (spl.size () == 0)
                continue;

            if (spl.size () >= 3)
            {
                count++;

                field = spl.at(0);
                value = spl.at (2);

                if (field == "xm")
                    out->xm = stod(value);
                else if (field == "xp")
                    out->xp = stod(value);
                else if (field == "ym")
                    out->ym = stod(value);
                else if (field == "yp")
                    out->yp = stod(value);
                else if (field == "hsize")
                    out->hsize = stod(value);
                else if (field == "dt")
                    out->dt = stod(value);
                else if (field == "Pk")
                    out->Pk = stoi(value);
                else if (field == "mouv")
                {
                    out->mouv = false;
                    if (value == "true")
                        out->mouv = true;
                }
                else if (field == "damping")
                {
                    out->damping = false;
                    if (value == "true")
                        out->damping = true;
                }
                else if (field == "gen")
                {
                    out->gen = false;
                    if (value == "true")
                        out->gen= true;
                }
                else if (field == "eq_type")
                {
                    if (value == "poisson" || value == "wave")
                        out->eq_type = value;
                    else
                        ERRORPARAMS
                }
                else if (field == "zeta_0")
                    out->zeta_0 = stod(value);
                else if (field == "beta_0")
                    out->beta_0 = stod(value);
                else if (field == "file_msh")
                    out->filename_msh = value;
                else
                    ERRORPARAMS
            }

            field = spl.at(0);

            if (field == "$BEGIN_OBJECT")
            {
                ObjectDatStruct obs;

                while (std::getline(infile, line))
                {
                    if (line == "$END_OBJECT")
                        break;

                    if (infile.eof ())
                    {
                        ERROR << "not enabled to find $END_OBJECT..." << BLINKRETURN << ENDLINE;
                        return;
                    }

                    if (line.front () != '#' && line.size () != 0)
                    {
                        RemoveBlankSpace (&line, &spl);

                        if (spl.size () >= 3)
                        {
                            count++;

                            field = spl.at(0);
                            value = spl.at (2);

                            if (field == "type")
                            {
                                if (value != "object")
                                {
                                    ERROR << "only object type is now supported" << BLINKRETURN << ENDLINE;
                                    return;
                                }
                            }
                            else if (field == "xc")
                                obs.xc = stod(value);
                            else if (field == "yc")
                                obs.yc = stod(value);
                            else if (field == "r")
                                obs.radius = stod(value);
                            else if (field == "nbpts")
                                obs.nbpts = stoi(value);
                            else if (field == "gen")
                            {
                                obs.gen_with_algo = false;
                                if (value == "true")
                                    obs.gen_with_algo = true;
                            }
                            else if (field == "enableMover")
                            {
                                obs.enableMover = false;
                                if (value == "true")
                                    obs.enableMover = true;
                            }
                            else if (field == "algo_gen")
                                obs.algo_gen = static_cast<unsigned int>(stod(value));
                            else if (field == "file_msh")
                                obs.filename_msh = value;
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
    HEADERFUN("ParseMSH");
    BEGIN << "Parse the msh file." << ENDLINE;

    if (mesh == nullptr)
    {
        ERROR << "the store object is not correct initialised." << BLINKRETURN << ENDLINE;
        return;
    }

    std::ifstream infile (filename);
    std::string line;
    std::string delimeter = "\"";
    std::vector <std::string> sv, tagphysical;
    std::vector<int> tag;
    int numtag = 0;
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

    // PhysicalNames
    while (line != "$PhysicalNames")
    {
        if (infile.eof ())
        {
            ERROR << "the $PhysicalNames field does not seem to be present" << BLINKRETURN << ENDLINE;
            return;
        }

        std::getline(infile, line);
    }

#ifdef DEBUG
    STATUS << "$PhysicalNames detected." << ENDLINE;
#endif

    std::getline(infile, line);
    RemoveBlankSpace (&line, &sv);

    numtag = stoi(sv.at (0));

    count = 0;
    while (count < 10 * numtag)
    {
        if (infile.eof ())
        {
            ERROR << "the $EndPhysicalNames field does not seem to be present" << BLINKRETURN << ENDLINE;
            return;
        }

        std::getline(infile, line);
        RemoveBlankSpace (&line, &sv);

        if (line == "$EndPhysicalNames")
            break;

        if (sv.size () < 3)
        {
            ERROR << "one physical name doesn't respect the rule : onid id str, please check your msh file." << BLINKRETURN << ENDLINE;
            return;
        }

        split (&sv.at (2), &delimeter, &sv);

        tagphysical.push_back (sv.at (1));
        count++;
    }

    if (count != numtag)
    {
        ERROR << "the number of physical names does not coincide." << BLINKRETURN << ENDLINE;
        return;
    }

//    mesh->SetTagPhysical (tagphysical);

#ifdef DEBUG
    STATUS << "$EndPhysicalNames detected." << ENDLINE;
#endif

    INFOS << "physical names : \t" << count << ENDLINE;

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
    tag.resize (std::size_t (numcells));

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

        tag.at (std::size_t (count)) = stoi (sv.at (3));

        depth_tag = 3 + stoi (sv.at (2));

        for (std::size_t i = std::size_t(depth_tag); i < sv.size (); ++i)
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

    mesh->GetCellsData ()->GetIntArrays ()->Add (NAME_TAG_PHYSICAL, tag);

    ENDFUN;

    //    DeleteLines (mesh);
    Build_NtoN (mesh);
    BuildEdgesWithHashMap (mesh);

    return;
}


void Print (InputDatStruct* input)
{
    HEADERFUN("Print");

    INFOS << "dt       = " << input->dt << ENDLINE;
    INFOS << "mouv     = " << input->mouv << ENDLINE;
    INFOS << "eq_type  = " << input->eq_type << ENDLINE;
    INFOS << "zeta_0   = " << input->zeta_0<< ENDLINE;
    INFOS << "beta_0   = " << input->beta_0<< ENDLINE;

    INFOS << "xm       = " << input->xm << ENDLINE;
    INFOS << "xp       = " << input->xp << ENDLINE;
    INFOS << "ym       = " << input->ym << ENDLINE;
    INFOS << "yp       = " << input->yp << ENDLINE;

    INFOS << "hsize    = " << input->hsize << ENDLINE;

    for (std::size_t i = 0; i < input->objects.size (); ++i)
    {
        auto obj = input->objects.at (i);
        INFOS << SEPARATOR << ENDLINE;
        INFOS << "\txc       = " << obj.xc << ENDLINE;
        INFOS << "\tyc       = " << obj.yc << ENDLINE;
        INFOS << "\tr        = " << obj.radius << ENDLINE;
        INFOS << "\tgen      = " << obj.gen_with_algo << ENDLINE;
        INFOS << "\talgo     = " << obj.algo_gen << ENDLINE;
        INFOS << "\tnbpts    = " << obj.nbpts << ENDLINE;
        INFOS << "\tfile_msh = " << obj.filename_msh << ENDLINE;
    }
    ENDFUN;

    return;
}

void split (std::string* s, std::string* delimiter, std::vector<std::string>* out)
{
    std::size_t pos_start = 0, pos_end, delim_len = delimiter->length();
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
    *out = {std::istream_iterator<std::string>{iss},           std::istream_iterator<std::string>{}};

    return;
}
