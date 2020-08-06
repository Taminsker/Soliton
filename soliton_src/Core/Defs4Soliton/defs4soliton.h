#ifndef DEFS4SOLITON_H
#define DEFS4SOLITON_H

#include <iostream>
#include <string>

#define NAME_SURROGATETAG       std::string("_SurrogateTag")
#define NAME_INTERSECTIONTAG    std::string("_IntersectionTag")
#define NAME_DISPLACEMENTVECTOR std::string("_d")
#define NAME_LEVELSET           std::string("_LS")

#define NAME_NORMALONCELLS      std::string("NormalOnCells")
#define NAME_NORMALONPOINTS     std::string("NormalOnPoints")
#define NAME_NORMALONEDGES      std::string("NormalOnEdges")


typedef int SOLITON_RETURN;
#define SOLITON_SUCCESS     EXIT_SUCCESS
#define SOLITON_FAILURE     EXIT_FAILURE
#define USE_SOLITON_RETURN(x) \
    if (x != SOLITON_SUCCESS)\
        return SOLITON_FAILURE


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
#define FLUSHLINE       COLOR_DEFAULT << std::flush;
#define SEPARATOR       "--------------------------"
#define BEGIN           std::cout << REVERSE << "--> "
#define ENDFUN          std::cout << std::endl

#ifdef DEBUG
    #define HEADERFUN(x) \
    std::string __NAMEFUN__ = x

    #define STATUS          std::cout                   << "(" << __NAMEFUN__ << ")\t" << " STATUS : "
    #define INFOS           std::cout                   << "(" << __NAMEFUN__ << ")\t" << "*  "
    #define ERROR           std::cerr << COLOR_RED      << "(" << __NAMEFUN__ << ")\t" << " ERROR : "
    #define BLINKRETURN     BLINK << " RETURN "
#else
    #define HEADERFUN(x)    (void)#x
    #define STATUS          std::cout                   << "STATUS : "
    #define INFOS           std::cout                   << " * "
    #define ERROR           std::cerr << COLOR_RED      << "ERROR : "
    #define BLINKRETURN     ""
#endif


#endif // DEFS4SOLITON_H
