#ifndef PARSERS_H
#define PARSERS_H

#include <fstream>
#include <string>
#include <sstream>
#include <iterator>
#include <vector>
#include <functional>

#include <ProgDef/proddef.h>
#include <Enums/enums.h>
#include <Solver/UpgradableSolver/ItemSolver/enum_item_solver_type.h>
#include <Solver/UpgradableSolver/UpgradableSolverBase/scheme_temporal_solver.h>

class Mesh;
class Point;
class Sto4Sol;

struct ObjectDatStruct
{
    std::string     filename_msh = "";
    std::string     basename = "obj_";
    unsigned int    algo_gen = 0;
    int             nbpts = 500;

    bool            rinp = false;
    double          x_center = 0.0;
    double          y_center = 0.0;
    double          z_center = 0.0;
    double          radius  = 1.;
};

struct ItemSolverDatStruct
{
    ITEM_T           type = ITEM_T::EMPTY;
    PHYS     tagToApply = PHYS::NONE;
    std::string      fun = "null_fun";
    int              dert = 0;
    int              id_obj = -1;
    bool             varTime = false;
};

struct InputDatStruct
{
    // Mesh part
    std::string     filename_msh = "";
    double          grid_x_m = 0.;
    double          grid_x_p = 0.;
    double          grid_y_m = 0.;
    double          grid_y_p = 0.;
    double          hsize = 0.;
    int             ele_type = 0;
    int             ele_order = 1;

    // Objects
    std::vector <ObjectDatStruct> objects;

    // Solver part
    double              dt = 0.;
    bool                objectsAreFixed = true;
    SCH_T               scheme = SCH_T::NO_TIME;
    std::vector<INTER>  listTagsSolver = {};

    std::vector <ItemSolverDatStruct> items;

    // Addons definitions
    bool            damping = false;
    double          zeta_0 = 1.0;
    double          beta_0 = 0.0;
    double          penal = 1E5;
    bool            colConcItem = false;

    // fun part
    std::function<double(Point, double)> fun_analytic;
    std::function<double(Point, double)> fun_neumann;
    std::function<double(Point, double)> fun_dirichlet;
    std::function<double(Point, double)> fun_secmember;


};

void ParseInputDatFile (InputDatStruct* out, std::string filename);
void ParseMSH (Mesh* mesh, std::string filename, bool keep_original = false);

void Print (InputDatStruct* input);

void split (std::string* s, std::string* delimiter, std::vector<std::string>* out);

void RemoveBlankSpace (std::string* s, std::vector<std::string>* out);

bool CastToBool(std::string *s);
#endif // PARSERS_H
