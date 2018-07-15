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


#if !defined(KRATOS_MULTIGRID_SOLVERS_APP_SMOOTHED_AGGREGATION_SOLVER_FACTORY_H_INCLUDED )
#define  KRATOS_MULTIGRID_SOLVERS_APP_SMOOTHED_AGGREGATION_SOLVER_FACTORY_H_INCLUDED



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
#include "custom_utilities/mg_level.h"
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
class SmoothedAggregationSolverFactory : public MultilevelSolverFactory<TSparseSpaceType, TDenseSpaceType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SmoothedAggregationSolverFactory
    KRATOS_CLASS_POINTER_DEFINITION(SmoothedAggregationSolverFactory);

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType> LinearSolverType;

    typedef MultilevelSolverFactory<TSparseSpaceType, TDenseSpaceType> BaseType;

    typedef typename LinearSolverType::Pointer LinearSolverPointerType;

    typedef MultilevelSolver<TSparseSpaceType, TDenseSpaceType> MultilevelSolverType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::MatrixPointerType SparseMatrixPointerType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    typedef typename TDenseSpaceType::VectorType DenseVectorType;

    typedef MGLevel<TSparseSpaceType, TDenseSpaceType> LevelType;

    typedef typename ĹevelType::Pointer LevelPointerType;

    typedef std::vector<LevelPointerType> LevelContainerType;

    typedef typename LevelContainerType::iterator LevelIteratorType;

    typedef Kratos::ParameterList<std::string> ParameterListType;

    typedef AMGUtils<TSparseSpaceType> AMGUtilsType;

    typedef typename AMGUtilsType::IndexVectorType IndexVectorType;

    typedef typename AMGUtilsType::ValueContainerType ValueContainerType;

    typedef typename AMGUtilsType::IndexContainerType IndexContainerType;

    typedef boost::shared_ptr<IndexVectorType> IndexVectorPointerType;

    typedef boost::numeric::ublas::unbounded_array<VectorType> VectorContainerType;

    typedef std::size_t  SizeType;

    typedef unsigned int  IndexType;

    typedef int  IntegerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SmoothedAggregationSolverFactory(ParameterListType& amg_parameter_list) : BaseType(amg_parameter_list)
    {
    }


    /// Copy constructor.
    SmoothedAggregationSolverFactory(const SmoothedAggregationSolverFactory& Other) : BaseType(Other)
    {
    }

    /// Destructor.
    virtual ~SmoothedAggregationSolverFactory() {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    SmoothedAggregationSolverFactory& operator=(const SmoothedAggregationSolverFactory& Other)
    {
        BaseType::operator=(Other);
        return *this;
    }


    ///@}
    ///@name Operations
    ///@{

    virtual void GenerateMultilevelSolver(MultilevelSolverType& solver, SparseMatrixType& rA)
    {
        #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
        std::cout << "#######################" << std::endl;
        std::cout << "At generate_smoothed_aggregation" << std::endl;
        #endif

        SizeType m = TSparseSpaceType::Size1(rA);

        // general parameters
        ParameterListType& amg_parameter_list = BaseType::mamg_parameter_list;

        IndexType max_levels = amg_parameter_list.get("max_levels", 10);
        IndexType max_coarse = amg_parameter_list.get("max_coarse", 500);
        std::string& symmetry = amg_parameter_list.get("symmetry", "hermitian");

        // right near-nullspace candidates
        ParameterListType& B_list = amg_parameter_list.sublist("B");
        IndexType B_length = B_list.get("length", 0);
        bool empty_right_nullspace = false;
        if(B_length == 0)
        {
            B_length += 1;
            empty_right_nullspace = true;
        }
        VectorContainerType B(B_length);
        if(empty_right_nullspace)
        {
            VectorType b(m, 1.00);
            B[0] = b;
        }
        else
        {
            for(IndexType i = 0; i < B_length; i++)
            {
                std::stringstream tmp;
                tmp << (i+1);
                VectorType& b = B_list.template get<VectorType>(tmp.str());
                B[i] = b;
            }
        }

        // left near-nullspace candidates
        ParameterListType& BH_list = amg_parameter_list.sublist("BH");
        IndexType BH_length = BH_list.get("length", 0);
        bool empty_left_nullspace = false;
        if(BH_length == 0)
        {
            BH_length += 1;
            empty_left_nullspace = true;
        }
        VectorContainerType BH(BH_length);
        if(empty_left_nullspace)
        {
            VectorType bh(m, 1.00);
            BH[0] = bh;
        }
        else
        {
            for(IndexType i = 0; i < BH_length; i++)
            {
                std::stringstream tmp;
                tmp << (i+1);
                VectorType& bh = BH_list.template get<VectorType>(tmp.str());
                BH[i] = bh;
            }
        }

        // strength of connection
        ParameterListType& strength_params = amg_parameter_list.sublist("strength");

        // aggregate
        ParameterListType& aggregate_params = amg_parameter_list.sublist("aggregate");

        // smooth
        ParameterListType& smooth_params = amg_parameter_list.sublist("smooth");

        // Bimprove
        ParameterListType& Bimprove_params = amg_parameter_list.sublist("Bimprove");

        std::string strength_name = strength_params.get("name", "classical");
        std::string aggregate_name = aggregate_params.get("name", ""); // TODO
        std::string smooth_name = smooth_params.get("name", ""); // TODO
        std::string Bimprove_name = smooth_params.get("name", ""); // TODO

        std::cout << "smoothed_aggregation parameters:" << std::endl;
        KRATOS_WATCH(max_levels);
        KRATOS_WATCH(max_coarse);
        KRATOS_WATCH(strength_name);
        KRATOS_WATCH(aggregate_name);
        KRATOS_WATCH(smooth_name);
        KRATOS_WATCH(Bimprove_name);
        std::cout << "end of smoothed_aggregation parameters" << std::endl;

        SizeType last_size = TSparseSpaceType::Size1(rA);

//        KRATOS_WATCH("at generate_ruge_stuben");
//        KRATOS_WATCH(&rA);

        // set up level 0
        solver.CreateLevel();
        LevelType& first_level = solver.GetLastLevel();
        SparseMatrixPointerType pA = first_level.GetCoarsenMatrix();
        TSparseSpaceType::Resize(*pA, last_size, last_size);
        TSparseSpaceType::Copy(rA, *pA);

        int cnt = 0;
        while(solver.GetNumberOfLevels() < max_levels && last_size > max_coarse)
        {
            #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
            std::cout << "#######" << std::endl;
            std::cout << "...generating level " << cnt << std::endl;
            #endif

            LevelType& current_level = solver.GetLastLevel();

            SparseMatrixPointerType A = current_level.GetCoarsenMatrix();

            #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
            std::cout << "...retrieved initial matrix for level " << cnt << std::endl;
            KRATOS_WATCH(A->size1())
            KRATOS_WATCH(A->size2())
            #endif

            // compute strength of connection
            SparseMatrixType C(last_size, last_size);
            if(strength_name.compare("classical") == 0)
            {
                #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
                std::cout << "...calculating classical strength of connection";
                #endif

                double theta = strength_params.get("theta", 0.25);
                AMGUtilsType::ClassicalStrengthOfConnection(C, *A, theta);
//                AMGUtilsType::AddToDiagonal(C, 1.0); //do we need this?

                #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
                std::cout << " completed" << std::endl;
                #endif
            }
            else if(strength_name.compare("symmetric") == 0)
            {
                #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
                std::cout << "...calculating symmetric strength of connection";
                #endif

                double theta = strength_params.get("theta", 1.0);
                AMGUtilsType::SymmetricStrengthOfConnection(C, *A, theta);
//                AMGUtilsType::AddToDiagonal(C, 1.0); //do we need this?

                #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
                std::cout << " completed" << std::endl;
                #endif
            }
            else if(strength_name.compare("distance") == 0)
            {
                #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
                std::cout << "...calculating distance strength of connection";
                #endif

                AMGUtilsType::DistanceStrengthOfConnection(C, *A);

                #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
                std::cout << " completed" << std::endl;
                #endif
            }
            else if(strength_name == std::string("None"))
            {
                #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
                std::cout << "...use coarse matrix as strength of connection" << std::endl;
                #endif
                TSparseSpaceType::Copy(*A, C);
            }
            else
                KRATOS_THROW_ERROR(std::logic_error, "strength_name is undefined or not supported:", strength_name)

            //TODO


            ++cnt;
        }

        #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
        std::cout << "generate_smoothed_aggregation completed" << std::endl;
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


