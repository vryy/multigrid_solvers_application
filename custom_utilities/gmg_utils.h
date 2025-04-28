/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
 */


/* *********************************************************
 *
 *   Last Modified by:    $Author: hbui $
 *   Date:                $Date: 14/7/2018 $
 *   Revision:            $Revision: 1.1 $
 *
 * ***********************************************************/


#if !defined(KRATOS_MULTIGRID_SOLVERS_APP_GMG_UTILS )
#define  KRATOS_MULTIGRID_SOLVERS_APP_GMG_UTILS


/* System includes */
#include <algorithm>
//#include <random>


/* External includes */


/* Project includes */
#include "includes/define.h"
#include "utilities/timer.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "custom_utilities/mg_level.h"
#include "custom_utilities/matrix_based_mg_level.h"
#include "custom_utilities/multi_index.h"


namespace Kratos
{

/**@name Kratos Globals */
/*@{ */


/*@} */
/**@name Type Definitions */
/*@{ */


/*@} */


/**@name  Enum's */
/*@{ */


/*@} */
/**@name  Functions */
/*@{ */



/*@} */
/**@name Kratos Classes */
/*@{ */

/**
 * This class defines utility for aggregation of node clusters to be used in deflated solvers.
 */
template<class TSparseSpaceType, class TDenseSpaceType, class TModelPartType>
class GMGUtils
{
public:

    /**@name Type Definitions */
    /*@{ */

    KRATOS_CLASS_POINTER_DEFINITION(GMGUtils);

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TSparseSpaceType::MatrixPointerType SparseMatrixPointerType;

    typedef typename TSparseSpaceType::SizeType SizeType;

    typedef typename TSparseSpaceType::IndexType IndexType;

    typedef typename TModelPartType::DataType ValueType;

    typedef boost::numeric::ublas::unbounded_array<IndexType> IndexContainerType;

    typedef boost::numeric::ublas::unbounded_array<int> IntegerContainerType;

    typedef boost::numeric::ublas::unbounded_array<ValueType> ValueContainerType;

    typedef boost::numeric::ublas::coordinate_matrix<ValueType> LocalSparseMatrixType;

//    typedef boost::shared_ptr<LocalSparseMatrixType> LocalSparseMatrixPointerType;

//    typedef std::vector<IndexType> IndexVectorType;
//    typedef boost::numeric::ublas::vector<IndexType> IndexVectorType;
    typedef boost::numeric::ublas::vector<int> IndexVectorType; // this must be int since this is the only one support by kratos python

//    typedef boost::shared_ptr<IndexVectorType> IndexVectorPointerType;

    typedef MatrixBasedMGLevel<TSparseSpaceType, TDenseSpaceType, TModelPartType> MatrixBasedLevelType;

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TModelPartType> LinearSolverType;

    typedef BuilderAndSolver<TSparseSpaceType, TDenseSpaceType, LinearSolverType, TModelPartType> BuilderAndSolverType;

    typedef typename BuilderAndSolverType::TSchemeType SchemeType;

    /*@} */
    /**@name Life Cycle
     */
    /*@{ */

    /** Constructor.
     */


    /** Destructor.
     */

    /*@} */
    /**@name Operators
     */

    /*@} */
    /**@name Operations */
    /*@{ */

    // static void Initialize(typename SchemeType::Pointer pScheme, ModelPart::Pointer pModelPart)
    // {
    //     if (!pScheme->SchemeIsInitialized())
    //         pScheme->Initialize(pModelPart);
    //     if (!pScheme->ElementsAreInitialized())
    //         pScheme->InitializeElements(pModelPart);
    //     if (!pScheme->ConditionsAreInitialized())
    //         pScheme->InitializeConditions(pModelPart);
    // }


    // static void InitializeSolutionStep()
    // {

    // }


    /// Compute the coarse matrix on the multigrid level of GMG. The model_part is required to compute the coarse matrix, which is sparse.
    static void ComputeCoarseMatrix(typename BuilderAndSolverType::Pointer pBuilderAndSolver,
        typename SchemeType::Pointer pScheme,
        typename TModelPartType::Pointer pModelPart,
        MatrixBasedLevelType& rLevel)
    {
        VectorType dummy;
        SparseMatrixPointerType pA = rLevel.GetCoarseMatrix();
        TSparseSpaceType::Resize(dummy, TSparseSpaceType::Size1(*pA));
        pBuilderAndSolver->Build(pScheme, *pModelPart, *pA, dummy);
        pBuilderAndSolver->ApplyDirichletConditions(pScheme, *pModelPart, *pA, dummy, dummy);
    }

    /// Compute the coarse matrix on the multigrid level of GMG. The model_part is required to compute the coarse matrix, which is sparse.
    static void ComputeCoarseMatrix(typename BuilderAndSolverType::Pointer pBuilderAndSolver,
        typename SchemeType::Pointer pScheme,
        typename TModelPartType::Pointer pModelPart,
        typename MatrixBasedLevelType::Pointer pLevel)
    {
        ComputeCoarseMatrix(pBuilderAndSolver, pScheme, pModelPart, *pLevel);
    }


    /*@} */
    /**@name Access */
    /*@{ */


    /*@} */
    /**@name Inquiry */
    /*@{ */


    /*@} */
    /**@name Friends */
    /*@{ */


    /*@} */

private:
    /**@name Static Member Variables */
    /*@{ */


    /*@} */
    /**@name Member Variables */
    /*@{ */


    /*@} */
    /**@name Private Operators*/
    /*@{ */


    /*@} */
    /**@name Private Operations*/
    /*@{ */


    /*@} */
    /**@name Private  Access */
    /*@{ */


    /*@} */
    /**@name Private Inquiry */
    /*@{ */


    /*@} */
    /**@name Un accessible methods */
    /*@{ */


    /*@} */

}; /* Class ClassName */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

} /* namespace Kratos.*/

#endif /* KRATOS_MULTIGRID_SOLVERS_APP_GMG_UTILS  defined */
