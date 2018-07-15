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

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 2015 Mar 30 $
//   Revision:            $Revision: 1.1 $
//
//


#if !defined(KRATOS_MULTIGRID_SOLVERS_APP_RUGE_STUEBEN_SOLVER_FACTORY_H_INCLUDED )
#define  KRATOS_MULTIGRID_SOLVERS_APP_RUGE_STUEBEN_SOLVER_FACTORY_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <sstream>
#include <cstddef>


// External includes
//#include "external_includes/pyamg/relaxation.h"
//#include "external_includes/pyamg/ruge_stuben.h"


// Project includes
#include "includes/define.h"
#include "linear_solvers/linear_solver.h"
#include "custom_utilities/matrix_based_mg_level.h"
#include "custom_utilities/matrix_based_mg_projector.h"
#include "custom_utilities/amg_utils.h"
#include "custom_utilities/parameter_list.h"
#include "custom_utilities/multilevel_solver_factory.h"
#include "custom_linear_solvers/multilevel_solver.h"

#define DEBUG_MULTILEVEL_SOLVER_FACTORY

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

template<class TSparseSpaceType, class TDenseSpaceType>
class RugeStuebenSolverFactory : public MultilevelSolverFactory<TSparseSpaceType, TDenseSpaceType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of RugeStuebenSolverFactory
    KRATOS_CLASS_POINTER_DEFINITION(RugeStuebenSolverFactory);

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType> LinearSolverType;

    typedef MultilevelSolverFactory<TSparseSpaceType, TDenseSpaceType> BaseType;

    typedef typename LinearSolverType::Pointer LinearSolverPointerType;

    typedef MultilevelSolver<TSparseSpaceType, TDenseSpaceType> MultilevelSolverType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::MatrixPointerType SparseMatrixPointerType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    typedef typename TDenseSpaceType::VectorType DenseVectorType;

    typedef MatrixBasedMGLevel<TSparseSpaceType, TDenseSpaceType> LevelType;

    typedef typename LevelType::Pointer LevelPointerType;

    typedef std::vector<LevelPointerType> LevelContainerType;

    typedef typename LevelContainerType::iterator LevelIteratorType;

    typedef MatrixBasedMGProjector<TSparseSpaceType> MGProjectorType;

    typedef Kratos::ParameterList<std::string> ParameterListType;

    typedef AMGUtils<TSparseSpaceType> AMGUtilsType;

    typedef typename AMGUtilsType::IndexVectorType IndexVectorType;

    typedef typename AMGUtilsType::ValueContainerType ValueContainerType;

    typedef typename AMGUtilsType::IndexContainerType IndexContainerType;

    typedef boost::shared_ptr<IndexVectorType> IndexVectorPointerType;

    typedef boost::numeric::ublas::unbounded_array<VectorType> VectorContainerType;

    typedef typename TSparseSpaceType::SizeType SizeType;

    typedef typename TSparseSpaceType::IndexType IndexType;

    typedef int  IntegerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    RugeStuebenSolverFactory(ParameterListType& amg_parameter_list) : BaseType(amg_parameter_list)
    {
    }


    /// Copy constructor.
    RugeStuebenSolverFactory(const RugeStuebenSolverFactory& Other) : BaseType(Other)
    {
    }

    /// Destructor.
    virtual ~RugeStuebenSolverFactory() {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    RugeStuebenSolverFactory& operator=(const RugeStuebenSolverFactory& Other)
    {
        BaseType::operator=(Other);
        return *this;
    }


    ///@}
    ///@name Operations
    ///@{

    virtual void GenerateMultilevelSolver(MultilevelSolverType& solver, SparseMatrixType& rA) const
    {
        #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
        std::cout << "##################################" << std::endl;
        std::cout << "### Executing RugeStuebenSolverFactory::" << __FUNCTION__ << std::endl;
        #endif

        // general parameters
        ParameterListType amg_parameter_list = BaseType::mamg_parameter_list;

        IndexType max_levels = amg_parameter_list.get<int>("max_levels", 10);
        IndexType max_coarse = amg_parameter_list.get<int>("max_coarse", 500);

        // strength
        ParameterListType strength_param = amg_parameter_list.sublist("strength");
        std::string strength_name = strength_param.get("name", "classical");

        // splitting method
        std::string CF = amg_parameter_list.get("CF", "RS");

        #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
        std::cout << "ruge_stuben parameters:" << std::endl;
        std::cout << "   "; KRATOS_WATCH(max_levels)
        std::cout << "   "; KRATOS_WATCH(max_coarse)
        std::cout << "   "; KRATOS_WATCH(strength_name)
        std::cout << "   "; KRATOS_WATCH(CF)
        #endif

        SizeType last_size = TSparseSpaceType::Size1(rA);

        // set up level 0
        solver.AddLevel(typename LevelType::Pointer(new LevelType(0)));
        LevelType& first_level = dynamic_cast<LevelType&>(solver.GetLastLevel());
        SparseMatrixPointerType pA = first_level.GetCoarseMatrix();
        TSparseSpaceType::Resize(*pA, last_size, last_size);
        TSparseSpaceType::Copy(rA, *pA);

        // create other levels until the maximum coarse size is reached
        while(solver.GetNumberOfLevels() < max_levels && last_size > max_coarse)
        {
            LevelType& current_level = dynamic_cast<LevelType&>(solver.GetLastLevel());

            #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
            std::cout << "---------------" << std::endl;
            std::cout << "> generating level " << current_level.LevelDepth() << std::endl;
            #endif

            SparseMatrixPointerType A = current_level.GetCoarseMatrix();

            #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
            std::cout << "   retrieved initial matrix for level " << current_level.LevelDepth() << std::endl;
            std::cout << "    "; KRATOS_WATCH(A->size1())
            std::cout << "    "; KRATOS_WATCH(A->size2())
            std::cout << "    "; KRATOS_WATCH(last_size)
            #endif

            // compute strength of connection
            SparseMatrixType C(last_size, last_size);
            if(strength_name == std::string("classical"))
            {
                #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
                std::cout << "   calculating classical strength of connection";
                #endif

                double theta = strength_param.get("theta", 0.25);
                AMGUtilsType::ClassicalStrengthOfConnection(C, *A, theta);
            }
            else if(strength_name == std::string("symmetric"))
            {
                #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
                std::cout << "   calculating symmetric strength of connection";
                #endif

                double theta = strength_param.get("theta", 1.0);
                AMGUtilsType::SymmetricStrengthOfConnection(C, *A, theta);
            }
            else if(strength_name == std::string("None"))
            {
                #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
                std::cout << "   use coarse matrix as strength of connection" << std::endl;
                #endif
                TSparseSpaceType::Copy(*A, C);
            }
            else
                KRATOS_THROW_ERROR(std::logic_error, "strength_name is undefined or not supported:", strength_name)

            #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
            std::cout << " completed" << std::endl;
            #endif

            // compute splitting
            IndexVectorType splitting(last_size, 0);
            if(CF == std::string("RS"))
            {
                #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
                std::cout << "   computing RS splitting";
                #endif
                AMGUtilsType::RS(splitting, C);
            }
            else if(CF == std::string("PMIS"))
            {
                #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
                std::cout << "   computing PMIS splitting";
                #endif
                AMGUtilsType::PMIS(splitting, C);
            }
            else if(CF == std::string("PMISc"))
            {
                #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
                std::cout << "   computing PMISc splitting";
                #endif
                std::string& coloring_method = strength_param.get("coloring_method", "JP");
                AMGUtilsType::PMISc(splitting, C, coloring_method);
            }
            else if(CF == std::string("CLJP"))
            {
                #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
                std::cout << "   computing CLJP splitting";
                #endif
                AMGUtilsType::CLJP(splitting, C);
            }
            else if(CF == std::string("CLJPc"))
            {
                #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
                std::cout << "   computing CLJPc splitting" << std::endl;
                #endif
                AMGUtilsType::CLJPc(splitting, C);
            }

            #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
            std::cout << " completed" << std::endl;
            #endif

            // generate prolongation operator
            #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
            std::cout << "   computing prolongation operator by direct interpolation";
            #endif

            typename MGProjectorType::Pointer pPrologator
                = typename MGProjectorType::Pointer(new MGProjectorType());
            SparseMatrixPointerType P = pPrologator->GetOperator();
            AMGUtilsType::DirectInterpolation(*P, *A, C, splitting); // P will be resized here
            current_level.SetProlongationOperator(pPrologator);

            #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
            std::cout << " completed" << std::endl;
            #endif

            // generate restriction operator
            #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
            std::cout << "   computing restriction operator by transposing the prolongation operator";
            #endif

            SizeType AfterCoarsenSize  = TSparseSpaceType::Size2(*P);
            SizeType BeforeCoarsenSize = TSparseSpaceType::Size1(*P);
            typename MGProjectorType::Pointer pRestrictor
                = typename MGProjectorType::Pointer(new MGProjectorType());
            SparseMatrixPointerType R = pRestrictor->GetOperator();
            TSparseSpaceType::Resize(*R, AfterCoarsenSize, BeforeCoarsenSize);

            #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
            std::cout << "{ BeforeCoarsenSize: " << BeforeCoarsenSize;
            std::cout << ", AfterCoarsenSize: " << AfterCoarsenSize << " }";
            #endif

            AMGUtilsType::Transpose(*P, *R);

            current_level.SetRestrictionOperator(pRestrictor);

            #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
            std::cout << " completed" << std::endl;
            #endif

            #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
            std::cout << "  generate level " << current_level.LevelDepth() << " completed" << std::endl;
            #endif

//            current_level->PrintData(std::cout);

            // generate coarsen matrix for next level
            #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
            std::cout << "  compute coarse matrix";
            #endif

            solver.AddLevel(typename LevelType::Pointer(new LevelType(current_level.LevelDepth()+1)));
            LevelType& last_level = dynamic_cast<LevelType&>(solver.GetLastLevel());
            SparseMatrixType tmp(BeforeCoarsenSize, AfterCoarsenSize);
            AMGUtilsType::Mult(*A, *P, tmp);
            SparseMatrixPointerType Ac = last_level.GetCoarseMatrix();
            TSparseSpaceType::Resize(*Ac, AfterCoarsenSize, AfterCoarsenSize);
            AMGUtilsType::Mult(*R, tmp, *Ac);

            #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
            std::cout << " for level " << last_level.LevelDepth() << " completed" << std::endl;
            #endif

            last_size = AfterCoarsenSize;

            bool terminate = false;
            if (last_size <= max_coarse)
            {
                terminate = true;
                #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
                std::cout << " level " << last_level.LevelDepth() << " coarse size = " << last_size
                          << " <= max_coarse = " << max_coarse << ". The process terminated." << std::endl;
                #endif
            }

            if (last_size > max_coarse && solver.GetNumberOfLevels() >= max_levels)
            {
                terminate = true;
                #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
                std::cout << " level " << last_level.LevelDepth() << " coarse size = " << last_size
                          << " > max_coarse = " << max_coarse << ". But the maximum number of levels is reached."
                          << " The process terminated but the coarsest size is " << last_size << "." << std::endl;
                #endif
            }

            if (terminate)
            {
                typename MGProjectorType::Pointer pPrologator
                    = typename MGProjectorType::Pointer(new MGProjectorType());
                last_level.SetProlongationOperator(pPrologator);

                typename MGProjectorType::Pointer pRestrictor
                    = typename MGProjectorType::Pointer(new MGProjectorType());

                last_level.SetRestrictionOperator(pPrologator);
            }
        }

        #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
        std::cout << "Executing RugeStuebenSolverFactory::" << __FUNCTION__ << " completed" << std::endl;
        std::cout << "#######################" << std::endl;
        #endif
    }

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; // Class RugeStuebenSolverFactory

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


}  // namespace Kratos.

#undef DEBUG_MULTILEVEL_SOLVER_FACTORY

#endif // KRATOS_MULTIGRID_SOLVERS_APP_SMOOTHED_AGGREGATION_SOLVER_FACTORY_H_INCLUDED  defined


