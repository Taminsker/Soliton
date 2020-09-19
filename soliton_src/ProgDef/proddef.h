#ifndef PRODDEF_H
#define PRODDEF_H


/**
  * PROG DEF
  **/

#define SOLITON_DEVICE_FUNC __attribute__((flatten)) __attribute__((always_inline))
#define SOLITON_INLINE inline

#define SOLITON_RESERVE_NAME MACSOLRES_

/**
  * NAMES DEFINITIONS
  **/

#define NAME_TAG_DAMPING_AREA           std::string("damping_area")
#define NAME_TAG_PHYSICAL               std::string("physical")
#define NAME_TAG_SURROGATE              std::string("_surrogate")
#define NAME_TAG_INTERSECTION           std::string("_intersection")
#define NAME_DISPLACEMENT_VECTOR        std::string("_d")
#define NAME_LEVELSET                   std::string("_LS")

#define NAME_NORMAL_ON_CELLS            std::string("normal_on_cells")
#define NAME_NORMAL_ON_POINTS           std::string("normal_on_points")
#define NAME_NORMAL_ON_EDGES            std::string("normal_on_edges")


#define UPSOLVER_MAX_TIME_DER           2
#define SOLQUEUE_MAX                    6 * UPSOLVER_MAX_TIME_DER


/**
  * PROGRAM DEFINITIONS
  **/

#define STR(x) #x
#define XSTR(x) STR(x)
#define CONCAT_STR(x, y, z) x##y##z
#define CONCAT_STR_2(x, y) x##y
#define MAKE_SOLITON_RESERVE(X) CONCAT_STR_2(MACSOLRES_, X)

#define SPACE std::left << std::setw(12)

#define EPSILON 1e-10
#define MATH_PI 3.14159265358979323846264338328
#define MATH_PI_2 1.57079632679489661923132169164

#define COLOR_BLACK     "\033[1;30m"
#define COLOR_RED       "\033[1;31m"
#define COLOR_GREEN     "\033[1;32m"
#define COLOR_YELLOW    "\033[1;33m"
#define COLOR_BLUE      "\033[1;34m"
#define COLOR_MAGENTA   "\033[1;35m"
#define COLOR_WHITE     "\033[1;97m"
#define COLOR_DEFAULT   "\033[0m"

#define BACK_MAGENTA    "\033[1;7m"
#define UNDERLINE       "\033[1;4m"
#define BLINK           "\033[1;5m"
#define REVERSE         "\033[1;7m"

#define ENDLINE         COLOR_DEFAULT << std::endl
#define FLUSHLINE       COLOR_DEFAULT << std::flush
#define SEPARATOR       "--------------------------"
#define BEGIN           std::cout << REVERSE << "--> "
#define ENDFUN          std::cout << std::endl

#ifdef DEBUG
    #define HEADERFUN(x) \
    std::string __NAMEFUN__ = x

    #define STATUS          std::cout                   << "(" << __NAMEFUN__ << ")\t" << " STATUS : "
    #define INFOS           std::cout                   << "(" << __NAMEFUN__ << ")\t" << "*  "
    #define ERROR           std::cerr << COLOR_RED      << "(" << __NAMEFUN__ << ")\t" << " ERROR : "
    #define WARNING         std::cout << COLOR_YELLOW   << "(" << __NAMEFUN__ << ")\t" << " WARNING : "
    #define BLINKRETURN     BLINK << " RETURN "
#else
    #define HEADERFUN(x)    (void)#x
    #define STATUS          std::cout                   << "STATUS : "
    #define INFOS           std::cout                   << " * "
    #define ERROR           std::cerr << COLOR_RED      << "ERROR : "
    #define WARNING         std::cout << COLOR_YELLOW   << "WARNING : "
    #define BLINKRETURN     ""
#endif

#define CLASS_NAME(X)   \
private:                \
    typedef X ClassName
#endif // PRODDEF_H
