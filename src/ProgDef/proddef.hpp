#ifndef PRODDEF_H
#define PRODDEF_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cstddef>

/**
 * PROG DEF
 **/

#define SOLITON_INLINE       inline
#define SOLITON_DEVICE_FUNC  SOLITON_INLINE __attribute__ ((flatten)) __attribute__ ((always_inline))
#define SOLITON_RESERVE_NAME MACSOLRES_

// Math def
typedef unsigned long ul_t;
typedef double        real_t;

namespace unstd
{
template <typename T>
using matrix = std::vector<std::vector<T>>;
}

typedef Eigen::SparseMatrix<real_t, Eigen::RowMajor> SparseMatrix_eig;
typedef Eigen::Matrix<real_t, Eigen::Dynamic, 1>     PlainVector_eig;
typedef Eigen::Matrix<real_t, 3, 3>                  Matrix3x3_eig;
typedef Eigen::Triplet<real_t>                       Triplet_eig;

/**
 * NAMES DEFINITIONS
 **/

#define NAME_TAG_DAMPING_AREA        std::string ("damping_area")
#define NAME_TAG_PHYSICAL            std::string ("physical")
#define NAME_TAG_INTERSECTION        std::string ("_intersection")
#define NAME_DISPLACEMENT_VECTOR     std::string ("_d")
#define NAME_LEVELSET                std::string ("_LS")
#define NAME_CONTRIBUTIONS_ON_CELLS  std::string ("_cellsContrib")
#define NAME_CONTRIBUTIONS_ON_EDGES  std::string ("_edgesContrib")
#define NAME_CONTRIBUTIONS_ON_POINTS std::string ("_pointsContrib")
#define NAME_NORMAL_ON_CELLS         std::string ("normal_on_cells")
#define NAME_NORMAL_ON_POINTS        std::string ("normal_on_points")
#define NAME_NORMAL_ON_EDGES         std::string ("normal_on_edges")
#define SOLITONCONTAINER_ORDER_MAX   5
#define SOLITONQUEUE_MAX             6 * SOLITONCONTAINER_ORDER_MAX

/**
 * PROGRAM DEFINITIONS
 **/

#define STR(x)                  #x
#define XSTR(x)                 STR (x)
#define CONCAT_STR(x, y, z)     x##y##z
#define CONCAT_STR_2(x, y)      x##y
#define MAKE_SOLITON_RESERVE(X) CONCAT_STR_2 (MACSOLRES_, X)
#define VOID_USE(X)             (void)X
#define SPACE                   std::left << std::setw (12)
#define NEG_ZERO                -0.000000E+00
#define EPSILON                 1E-10
#define MAX_VALUE               1E25
#define MATH_PI                 3.14159265358979323846264338328
#define MATH_PI_2               1.57079632679489661923132169164
#define COLOR_BLACK             "\033[0;30m"
#define COLOR_RED               "\033[0;31m"
#define COLOR_GREEN             "\033[0;32m"
#define COLOR_YELLOW            "\033[0;33m"
#define COLOR_BLUE              "\033[0;34m"
#define COLOR_MAGENTA           "\033[0;35m"
#define COLOR_WHITE             "\033[0;97m"
#define COLOR_DEFAULT           "\033[0;0m"
#define BACK_MAGENTA            "\033[0;7m"
#define UNDERLINE               "\033[0;4m"
#define BLINK                   "\033[0;5m"
#define REVERSE                 "\033[1;7m"
#define ENDLINE                 COLOR_DEFAULT << std::endl
#define FLUSHLINE               COLOR_DEFAULT << std::flush
#define SEPARATOR               "--------------------------"
#define TREE_BRANCH             std::cout << "\u251C\u2500\u2500 "
#define BEGIN                   std::cout << REVERSE << "--> "
#define ENDFUN                  std::cout << ENDLINE
#define COUT                    std::cout

#ifdef DEBUG
#define STATUS std::cout << "(" << __FUNCTION__ << ")\t" \
                         << " STATUS : "
#define INFOS std::cout << "(" << __FUNCTION__ << ")\t" \
                        << "* "
#define ERROR std::cerr << COLOR_RED << "(" << __FUNCTION__ << ")\t" \
                        << " ERROR : "
#define WARNING std::cout << COLOR_YELLOW << "(" << __FUNCTION__ << ")\t" \
                          << " WARNING : "
#define BLINKRETURN BLINK << " RETURN "
#else
#define STATUS      std::cout << "STATUS : "
#define INFOS       std::cout << " * "
#define ERROR       std::cerr << COLOR_RED << "ERROR : "
#define WARNING     std::cout << COLOR_YELLOW << "WARNING : "
#define BLINKRETURN ""
#endif

#endif  // PRODDEF_H
