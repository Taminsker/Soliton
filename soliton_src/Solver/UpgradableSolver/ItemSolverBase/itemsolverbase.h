//#ifndef ITEMSOLVERBASE_H
//#define ITEMSOLVERBASE_H

//#include <Solver/UpgradableSolver/ItemSolver/enum_item_solver_type.h>
//#include <Core/core.h>
//#include <Solver/DFStruct/dfstruct.h>
//#include <Solver/FEStruct/festruct.h>
//#include <Solver/QuadStruct/quadstruct.h>

//#define NULL_FUNC_NAME soliton_null_fun
//double NULL_FUNC_NAME (Point a, double b);

//#define NEED_TO_BE_RECOMPUTE true

//class ItemSolverBase
//{
//public:
//    ItemSolverBase ();
//    ItemSolverBase (const ItemSolverBase& tocopy);
//    virtual ~ItemSolverBase ();

//    operator SparseMatrix* ();
//    operator PlainVector* ();

//    SparseMatrix*           CastSparseMatrix ();
//    PlainVector*            CastPlainVector ();

//    std::vector<int>*      GetApplicableCellsVector();
//    std::vector<int>*      GetApplicableEdgesVector();


//    virtual void            Compute (double time, bool view = true, bool force = false) = 0;
//    virtual ItemSolverBase& SetEmbMesh (Mesh* mesh) = 0;
//    virtual ItemSolverBase& SetTagToApply (PHYS tag) = 0;
//    virtual ItemSolverBase& SetFunction (std::function<double(Point, double)> fun) = 0;
//    virtual ItemSolverBase& SetDerT (std::size_t dert) = 0;
//    virtual ItemSolverBase& SetVariableOverTime (bool variableOverTime) = 0;
////    virtual ItemSolverBase& Set
//    void                    ConfigureFromSolver (Mesh **mesh, double *coeffpen, bool *collect);
//    std::size_t             GetTemporalDerivationOrder ();
//    ITEM_SOLVER_TYPE        GetTypeOfItem ();
//    bool                    ItsVariableOverTime ();
//    PHYS            GetTagToApply();

//    friend std::ostream& operator<< (std::ostream& out, const ItemSolverBase& item);

//protected:
//    /*!!! Set from solver */
//    Mesh                                    **m_mesh;
//    double                                  *m_coeff_pen;
//    bool                                    *m_collect;

//    /*!!! Internal objects */
//    DFStore                                 *m_dfstore;
//    FEStore                                 *m_festore;
//    QuadStore                               *m_quadstore;
//    std::vector<Triplet>                    m_triplets;
//    std::vector<Duo>                        m_duos;
//    std::vector<int>                        m_applicableCells;
//    std::vector<int>                        m_applicableEdges;

//    /*!!! Item define objects */
//    Mesh                                    *m_emb;
//public:
//    std::function<double(Point, double)>    m_fun;
//protected:
//    PHYS                            m_tagToApply;
//    std::size_t                             m_dert;
//    ITEM_SOLVER_TYPE                        m_type;
//    bool                                    m_variableOverTime;

//    /*!!! Outer */
//    SparseMatrix                            m_mat;
//    PlainVector                             m_second_member;

//};

//std::ostream& operator<< (std::ostream& out, const ItemSolverBase& item);

//#endif // ITEMSOLVERBASE_H
