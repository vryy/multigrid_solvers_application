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
//   Date:                $Date: 16/7/2018 $
//   Revision:            $Revision: 1.1 $
//
//


#if !defined(KRATOS_MULTIGRID_SOLVERS_APP_GMG_STRUCTURED_SOLVER_FACTORY_H_INCLUDED )
#define  KRATOS_MULTIGRID_SOLVERS_APP_GMG_STRUCTURED_SOLVER_FACTORY_H_INCLUDED



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
#include "custom_utilities/parameter_list.h"
#include "custom_linear_solvers/multilevel_solver.h"
#include "custom_utilities/structured_mg_prolongator.h"
#include "custom_utilities/structured_matrix_based_mg_prolongator.h"
#include "custom_utilities/structured_mg_restrictor.h"
#include "custom_utilities/structured_matrix_based_mg_restrictor.h"
#include "custom_utilities/multilevel_solver_factory.h"
#include "custom_utilities/gmg_utils.h"

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

/**
 ** Solver factory for geometric multigrid (GMG). With GMG, it's important to not delete the row/column of Dirichlet BC.
 */
template<class TSparseSpaceType, class TDenseSpaceType, std::size_t TDim>
class GMGStructuredSolverFactory : public MultilevelSolverFactory<TSparseSpaceType, TDenseSpaceType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of GMGStructuredSolverFactory
    KRATOS_CLASS_POINTER_DEFINITION(GMGStructuredSolverFactory);

    typedef MultilevelSolverFactory<TSparseSpaceType, TDenseSpaceType> BaseType;

    typedef typename BaseType::MultilevelSolverType MultilevelSolverType;

    typedef typename BaseType::SparseMatrixType SparseMatrixType;

    typedef typename BaseType::SparseMatrixPointerType SparseMatrixPointerType;

    typedef typename BaseType::ParameterListType ParameterListType;

    typedef MatrixBasedMGLevel<TSparseSpaceType, TDenseSpaceType> LevelType;

    typedef typename LevelType::Pointer LevelPointerType;

    typedef GMGUtils<TSparseSpaceType, TDenseSpaceType> GMGUtilsType;

    typedef typename GMGUtilsType::BuilderAndSolverType BuilderAndSolverType;

    typedef typename GMGUtilsType::SchemeType SchemeType;

    typedef typename BaseType::SizeType SizeType;

    typedef typename BaseType::IndexType IndexType;

    typedef typename BaseType::IntegerType IntegerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    GMGStructuredSolverFactory(ParameterListType& gmg_parameter_list)
    : BaseType(gmg_parameter_list)
    {
        IndexType nlevels = gmg_parameter_list.get("num_levels", 1);
        mpModelParts.resize(nlevels);
        mpBuilderAndSolvers.resize(nlevels);
        mpSchemes.resize(nlevels);
    }

    /// Copy constructor.
    GMGStructuredSolverFactory(const GMGStructuredSolverFactory& rOther)
    : BaseType(rOther)
    , mpModelParts(rOther.mpModelParts)
    , mpBuilderAndSolvers(rOther.mpBuilderAndSolvers)
    , mpSchemes(rOther.mpSchemes)
    {}

    /// Destructor.
    virtual ~GMGStructuredSolverFactory()
    {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    GMGStructuredSolverFactory& operator=(const GMGStructuredSolverFactory& rOther)
    {
        BaseType::operator=(rOther);
        mpModelParts = rOther.mpModelParts;
        mpBuilderAndSolvers = rOther.mpBuilderAndSolvers;
        mpSchemes = rOther.mpSchemes;
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    typename BuilderAndSolverType::Pointer GetBuilderAndSolver(const std::size_t& lvl)
    {
        return mpBuilderAndSolvers[lvl];
    }

    typename BuilderAndSolverType::Pointer GetBuilderAndSolver(const std::size_t& lvl) const
    {
        return mpBuilderAndSolvers[lvl];
    }

    void SetBuilderAndSolver(const std::size_t& lvl, typename BuilderAndSolverType::Pointer pBuilderAndSolver)
    {
        if (lvl < mpBuilderAndSolvers.size())
            mpBuilderAndSolvers[lvl] = pBuilderAndSolver;
        else
            KRATOS_THROW_ERROR(std::logic_error, "Error setting builder_and_solver for non-existing level", lvl)
    }

    typename SchemeType::Pointer GetScheme(const std::size_t& lvl)
    {
        return mpSchemes[lvl];
    }

    typename SchemeType::Pointer GetScheme(const std::size_t& lvl) const
    {
        return mpSchemes[lvl];
    }

    void SetScheme(const std::size_t& lvl, typename SchemeType::Pointer pScheme)
    {
        if (lvl < mpSchemes.size())
            mpSchemes[lvl] = pScheme;
        else
            KRATOS_THROW_ERROR(std::logic_error, "Error setting scheme for non-existing level", lvl)
    }

    ModelPart::Pointer GetModelPart(const std::size_t& lvl)
    {
        return mpModelParts[lvl];
    }

    ModelPart::Pointer GetModelPart(const std::size_t& lvl) const
    {
        return mpModelParts[lvl];
    }

    void SetModelPart(const std::size_t& lvl, ModelPart::Pointer pModelPart)
    {
        if (lvl < mpModelParts.size())
            mpModelParts[lvl] = pModelPart;
        else
            KRATOS_THROW_ERROR(std::logic_error, "Error setting model_part for non-existing level", lvl)
    }

    virtual void InitializeMultilevelSolver(MultilevelSolverType& solver) const
    {
        ParameterListType gmg_parameter_list = BaseType::mmg_parameter_list;

        IndexType nlevels = gmg_parameter_list.get("num_levels", 1);

        for (std::size_t lvl = 0; lvl < nlevels; ++lvl)
        {
            solver.AddLevel(typename LevelType::Pointer(new LevelType(lvl)));
        }
    }

    virtual void GenerateMultilevelSolver(MultilevelSolverType& solver, SparseMatrixType& rA) const
    {
        typedef StructuredMGProlongator<TSparseSpaceType, TDim> MGProlongatorType;
        typedef StructuredMGRestrictor<TSparseSpaceType, TDim> MGRestrictorType;
        // typedef StructuredMatrixBasedMGProlongator<TSparseSpaceType, TDim> MGProlongatorType;
        // typedef StructuredMatrixBasedMGRestrictor<TSparseSpaceType, TDim> MGRestrictorType;

        #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
        std::cout << "##################################" << std::endl;
        std::cout << "### Executing GMGStructuredSolverFactory::" << __FUNCTION__ << std::endl;
        #endif

        // general parameters
        ParameterListType gmg_parameter_list = BaseType::mmg_parameter_list;

        IndexType nlevels = gmg_parameter_list.get("num_levels", 1);
        IndexType block_size = gmg_parameter_list.get("block_size", 1);
        IndexType coarse_div_1 = gmg_parameter_list.get("coarse_div_1", 10);
        IndexType coarse_div_2 = gmg_parameter_list.get("coarse_div_2", 10);
        IndexType coarse_div_3 = gmg_parameter_list.get("coarse_div_3", 10);

        #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
        std::cout << "gmg parameters:" << std::endl;
        std::cout << "   "; KRATOS_WATCH(nlevels)
        std::cout << "   "; KRATOS_WATCH(block_size)
        std::cout << "   "; KRATOS_WATCH(coarse_div_1)
        std::cout << "   "; KRATOS_WATCH(coarse_div_2)
        std::cout << "   "; KRATOS_WATCH(coarse_div_3)
        #endif

        KRATOS_WATCH(solver.GetNumberOfLevels())
        for (std::size_t lvl = 0; lvl < nlevels; ++lvl)
        {
            LevelType& current_level = dynamic_cast<LevelType&>(solver.GetLevel(lvl));

            #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
            std::cout << "---------------" << std::endl;
            std::cout << "> generating level " << current_level.LevelDepth() << std::endl;
            #endif

            // generate coarsen matrix for next level
            #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
            std::cout << "   compute coarse matrix";
            #endif

            if (lvl > 0)
            {
                // set up level lvl. Level 0 coarse matrix is assumed to be rA
                // assemble the coarse matrix
                GMGUtilsType::ComputeCoarseMatrix(this->GetBuilderAndSolver(lvl),
                    this->GetScheme(lvl), this->GetModelPart(lvl), current_level);
            }

            #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
            std::cout << " for level " << current_level.LevelDepth() << " completed" << std::endl;
            #endif

            // generate prolongation operator
            #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
            std::cout << "   generating prolongation operator";
            #endif

            typename MGProlongatorType::Pointer pPrologator;
            if (lvl < nlevels-1)
                pPrologator = typename MGProlongatorType::Pointer(new MGProlongatorType(this->GetModelPart(lvl+1), this->GetModelPart(lvl)));
            else
                pPrologator = typename MGProlongatorType::Pointer(new MGProlongatorType());
            pPrologator->SetBlockSize(block_size);
            pPrologator->SetDivision(0, coarse_div_1 << (nlevels-1-lvl));
            pPrologator->SetDivision(1, coarse_div_2 << (nlevels-1-lvl));
            if (TDim == 3)
                pPrologator->SetDivision(2, coarse_div_3 << (nlevels-1-lvl));
            pPrologator->Initialize();

            current_level.SetProlongationOperator(pPrologator);

            #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
            std::cout << " completed" << std::endl;
            #endif

            // generate prolongation operator
            #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
            std::cout << "   generating restriction operator";
            #endif

            typename MGRestrictorType::Pointer pRestrictor;
            if (lvl < nlevels-1)
                pRestrictor = typename MGRestrictorType::Pointer(new MGRestrictorType(this->GetModelPart(lvl+1), this->GetModelPart(lvl)));
            else
                pRestrictor = typename MGRestrictorType::Pointer(new MGRestrictorType());
            pRestrictor->SetBlockSize(block_size);
            pRestrictor->SetDivision(0, coarse_div_1 << (nlevels-2-lvl));
            pRestrictor->SetDivision(1, coarse_div_2 << (nlevels-2-lvl));
            if (TDim == 3)
                pRestrictor->SetDivision(2, coarse_div_3 << (nlevels-2-lvl));
            pRestrictor->Initialize();

            current_level.SetRestrictionOperator(pRestrictor);

            #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
            std::cout << " completed" << std::endl;
            #endif
        }

        #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
        std::cout << "Executing GMGStructuredSolverFactory::" << __FUNCTION__ << " completed" << std::endl;
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

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Geometric multigrid solver factory";
        return  buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        BaseType::PrintData(rOStream);
    }


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

    std::vector<typename BuilderAndSolverType::Pointer> mpBuilderAndSolvers;
    std::vector<typename SchemeType::Pointer> mpSchemes;
    std::vector<ModelPart::Pointer> mpModelParts;

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

}; // Class GMGStructuredSolverFactory

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TSparseSpaceType, class TDenseSpaceType, std::size_t TDim>
inline std::istream& operator >> (std::istream& IStream, GMGStructuredSolverFactory<TSparseSpaceType, TDenseSpaceType, TDim>& rThis)
{
    return IStream;
}

/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType, std::size_t TDim>
inline std::ostream& operator << (std::ostream& rOStream, const GMGStructuredSolverFactory<TSparseSpaceType, TDenseSpaceType, TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.

#undef DEBUG_MULTILEVEL_SOLVER_FACTORY

#endif // KRATOS_MULTIGRID_SOLVERS_APP_GMG_STRUCTURED_SOLVER_FACTORY_H_INCLUDED  defined
