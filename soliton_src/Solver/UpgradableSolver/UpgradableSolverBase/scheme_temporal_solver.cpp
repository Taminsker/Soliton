//#include "scheme_temporal_solver.h"

//std::string ToString (SCHEME_TEMPORAL_SOLVER type)
//{
//    switch (type)
//    {
//    case SCHEME_TEMPORAL_SOLVER::NO_TIME:
//        return "NO_TIME";
//    case SCHEME_TEMPORAL_SOLVER::EULER_IMPLICIT:
//        return "EULER_IMPLICIT";
//    case SCHEME_TEMPORAL_SOLVER::EULER_EXPLICIT:
//        return "EULER_EXPLICIT";
//    case SCHEME_TEMPORAL_SOLVER::TVD_RK_2:
//        return "TVD_RK_2";
//    case SCHEME_TEMPORAL_SOLVER::TVD_RK_4:
//        return "TVD_RK_4";
//    case SCHEME_TEMPORAL_SOLVER::NEWMARK_SCHEME:
//        return "NEWMARK_SCHEME";
//    }

//    return "NO_TIME";
//}

//std::ostream& operator<< (std::ostream& out, SCHEME_TEMPORAL_SOLVER type)
//{
//    return out << ToString (type) << std::flush;
//}

//SCHEME_TEMPORAL_SOLVER FromStringScheme (std::string s)
//{
//    if (s == "NO_TIME")
//        return SCHEME_TEMPORAL_SOLVER::NO_TIME;
//    else if (s == "EULER_IMPLICIT")
//        return SCHEME_TEMPORAL_SOLVER::EULER_IMPLICIT;
//    else if (s == "EULER_EXPLICIT")
//        return SCHEME_TEMPORAL_SOLVER::EULER_EXPLICIT;
//    else if (s == "TVD_RK_2")
//        return SCHEME_TEMPORAL_SOLVER::TVD_RK_2;
//    else if (s == "TVD_RK_4")
//        return SCHEME_TEMPORAL_SOLVER::TVD_RK_4;
//    else if (s == "NEWMARK_SCHEME")
//        return SCHEME_TEMPORAL_SOLVER::NEWMARK_SCHEME;
//    else
//        return SCHEME_TEMPORAL_SOLVER::NO_TIME;
//}
