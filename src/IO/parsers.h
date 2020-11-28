#ifndef PARSERS_H
#define PARSERS_H

#include <Enums/enums.h>
#include <ProgDef/proddef.h>

#include <fstream>
#include <functional>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

class Mesh;
class Point;
class MeshStorage;

struct ObjectDatStruct
{
    std::string filename_msh = "";
    std::string basename     = "obj_";
    ul_t        algo_gen     = 0;
    int         nbpts        = 500;
    bool        rinp         = false;
    real_t      x_center     = 0.0;
    real_t      y_center     = 0.0;
    real_t      z_center     = 0.0;
    real_t      radius       = 1.;
};

struct InputDatStruct
{
    // Mesh part
    std::string filename_msh = "";
    real_t      grid_x_m     = 0.;
    real_t      grid_x_p     = 0.;
    real_t      grid_y_m     = 0.;
    real_t      grid_y_p     = 0.;
    real_t      hsize        = 0.;
    int         ele_type     = 0;
    int         ele_order    = 1;

    // Objects
    std::vector<ObjectDatStruct> objects;

    // Solver part
    real_t dt              = 0.;
    bool   objectsAreFixed = true;

    // Addons definitions
    bool   damping     = false;
    real_t zeta_0      = 1.0;
    real_t beta_0      = 0.0;
    real_t penal       = 1E5;
    int    powpenalty  = 1;
    bool   colConcItem = false;
    bool   synthetize  = true;
    real_t g_variable  = 9.81;
};

void
ParseInputDatFile (InputDatStruct * out, std::string filename);

void
ParseMSH (Mesh * mesh, std::string filename, bool keep_original = false);

void
Print (InputDatStruct * input);

void
split (std::string * s, std::string * delimiter, std::vector<std::string> * out);

void
RemoveBlankSpace (std::string * s, std::vector<std::string> * out);

bool
CastToBool (std::string * s);
#endif  // PARSERS_H
