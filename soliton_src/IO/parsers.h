#ifndef PARSERS_H
#define PARSERS_H

#include <fstream>
#include <string>
#include <sstream>
#include <iterator>
#include <vector>

#include <Core/Defs4Soliton/defs4soliton.h>

class Mesh;
class Sto4Sol;

struct ObjectDatStruct
{
    double xc = 0.0;
    double yc = 0.0;
    bool enableMover = false;
    bool gen_with_algo = false;
    unsigned int algo_gen = 0;
    std::string filename_msh = "";
    double radius  = 1.;
    int nbpts = 500;
};

struct InputDatStruct
{
    std::string filename = "";
    double xm = 0.;
    double xp = 0.;
    double ym = 0.;
    double yp = 0.;
    double hsize = 0.;
    double dt = 0.;
    bool mouv = false;
    std::string eq_type = "poisson";
    double zeta_0 = 1.0;
    double beta_0 = 0.0;

    std::vector <ObjectDatStruct> objects;
};

SOLITON_RETURN ParseInputDatFile (InputDatStruct* out, std::string filename);
SOLITON_RETURN ParseMSH (Mesh* mesh, std::string filename, bool keep_original = false);

SOLITON_RETURN Print (InputDatStruct* input);

SOLITON_RETURN split (std::string* s, std::string* delimiter, std::vector<std::string>* out);

SOLITON_RETURN RemoveBlankSpace (std::string* s, std::vector<std::string>* out);

#endif // PARSERS_H
