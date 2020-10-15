//#include "enum_item_solver_type.h"

//std::string ToString (ITEM_SOLVER_TYPE type)
//{
//    switch (type)
//    {
//    default:
//        return "EMPTY";
//    case ITEM_SOLVER_TYPE::PHI_PHI:
//        return "PHI_PHI";
//    case ITEM_SOLVER_TYPE::DAMPING_PHI_PHI:
//        return "DAMPING_PHI_PHI";
//    case ITEM_SOLVER_TYPE::GRAD_GRAD:
//        return "GRAD_GRAD";
//    case ITEM_SOLVER_TYPE::SECOND_MEMBER:
//        return "SECOND_MEMBER";
//    case ITEM_SOLVER_TYPE::DIRICHLET_BOUNDARY:
//        return "DIRICHLET_BOUNDARY";
//    case ITEM_SOLVER_TYPE::NEUMANN_BOUNDARY:
//        return "NEUMANN_BOUNDARY";
//    case ITEM_SOLVER_TYPE::DIRICHLET_BOUNDARY_SBM:
//        return "DIRICHLET_BOUNDARY_SBM";
//    case ITEM_SOLVER_TYPE::NEUMANN_BOUNDARY_SBM:
//        return "NEUMANN_BOUNDARY_SBM";
//    }
//}


//std::ostream& operator<< (std::ostream& out, ITEM_SOLVER_TYPE type)
//{
//    return out << ToString (type) << std::flush;
//}


//ITEM_SOLVER_TYPE FromStringItemSolver (std::string tag)
//{
//    if (tag == "PHI_PHI")
//        return ITEM_SOLVER_TYPE::PHI_PHI;
//     else if (tag == "DAMPING_PHI_PHI")
//        return ITEM_SOLVER_TYPE::DAMPING_PHI_PHI;
//    else if (tag == "GRAD_GRAD")
//       return ITEM_SOLVER_TYPE::GRAD_GRAD;
//    else if (tag == "SECOND_MEMBER")
//       return ITEM_SOLVER_TYPE::SECOND_MEMBER;
//    else if (tag == "DIRICHLET_BOUNDARY")
//       return ITEM_SOLVER_TYPE::DIRICHLET_BOUNDARY;
//    else if (tag == "NEUMANN_BOUNDARY")
//       return ITEM_SOLVER_TYPE::NEUMANN_BOUNDARY;
//    else if (tag == "DIRICHLET_BOUNDARY_SBM")
//       return ITEM_SOLVER_TYPE::DIRICHLET_BOUNDARY_SBM;
//    else if (tag == "NEUMANN_BOUNDARY_SBM")
//       return ITEM_SOLVER_TYPE::NEUMANN_BOUNDARY_SBM;
//    else
//        return ITEM_SOLVER_TYPE::EMPTY;
//}
