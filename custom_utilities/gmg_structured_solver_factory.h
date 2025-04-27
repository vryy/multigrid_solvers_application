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
#include "includes/model_part.h"
#include "linear_solvers/linear_solver.h"
#include "custom_utilities/parameter_list.h"
#include "custom_linear_solvers/multilevel_solver.h"
#include "custom_utilities/structured_mesh_mg_projector.h"
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

    typedef StructuredMeshMGProjector<TSparseSpaceType, TDim> StructuredMeshMGProjectorType;

    typedef typename LevelType::Pointer LevelPointerType;

    typedef GMGUtils<TSparseSpaceType, TDenseSpaceType, ModelPart> GMGUtilsType;

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
        mpProlongationType.resize(nlevels);
        mpRestrictionType.resize(nlevels);
        mpTransferType.resize(nlevels);
    }

    /// Copy constructor.
    GMGStructuredSolverFactory(const GMGStructuredSolverFactory& rOther)
    : BaseType(rOther)
    , mpModelParts(rOther.mpModelParts)
    {}

    /// Destructor.
    ~GMGStructuredSolverFactory() override
    {}

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    GMGStructuredSolverFactory& operator=(const GMGStructuredSolverFactory& rOther)
    {
        BaseType::operator=(rOther);
        mpModelParts = rOther.mpModelParts;
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    ModelPart::Pointer GetModelPart(const std::size_t& lvl)
    {
        return mpModelParts[lvl];
    }

    ModelPart::Pointer GetModelPart(const std::size_t& lvl) const
    {
        return mpModelParts[lvl];
    }

    /// Set the model_part at specific level
    void SetModelPart(const std::size_t& lvl, ModelPart::Pointer pModelPart)
    {
        if (lvl < mpModelParts.size())
            mpModelParts[lvl] = pModelPart;
        else
            KRATOS_ERROR << "Error setting model_part for non-existing level " << lvl;
    }

    /// Set the restriction operator for all levels
    void SetRestrictionOperator(typename StructuredMeshMGProjectorType::Pointer pProjector)
    {
        for (unsigned int lvl = 0; lvl < mpRestrictionType.size(); ++lvl)
            mpRestrictionType[lvl] = pProjector;
    }

    /// Set the transfer operator for all levels
    void SetTransferOperator(typename StructuredMeshMGProjectorType::Pointer pProjector)
    {
        for (unsigned int lvl = 0; lvl < mpTransferType.size(); ++lvl)
            mpTransferType[lvl] = pProjector;
    }

    /// Set the prolongation operator for all levels
    void SetProlongationOperator(typename StructuredMeshMGProjectorType::Pointer pProjector)
    {
        for (unsigned int lvl = 0; lvl < mpProlongationType.size(); ++lvl)
            mpProlongationType[lvl] = pProjector;
    }

    void InitializeMultilevelSolver(MultilevelSolverType& solver) const override
    {
        // create the levels
        ParameterListType gmg_parameter_list = BaseType::mmg_parameter_list;

        IndexType nlevels = gmg_parameter_list.get("num_levels", 1);

        for (std::size_t lvl = 0; lvl < nlevels; ++lvl)
        {
            solver.AddLevel(typename LevelType::Pointer(new LevelType(lvl)));
        }

        // initialize the levels
        #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
        std::cout << "##################################" << std::endl;
        std::cout << "### Executing GMGStructuredSolverFactory::" << __FUNCTION__ << std::endl;
        #endif

        // general parameters
        IndexType block_size = gmg_parameter_list.get("block_size", 1);
        IndexType coarse_div[TDim];
        for (unsigned int dim = 0; dim < TDim; ++dim)
        {
            std::stringstream ss;
            ss << "coarse_div_" << (dim+1);
            coarse_div[dim] = gmg_parameter_list.get(ss.str(), 10);
        }

        #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
        std::cout << "gmg parameters:" << std::endl;
        std::cout << "   "; KRATOS_WATCH(nlevels)
        std::cout << "   "; KRATOS_WATCH(block_size)
        for (unsigned int dim = 0; dim < TDim; ++dim)
            std::cout << "   coarse_div[" << dim << "]: " << coarse_div[dim] << std::endl;
        #endif

        KRATOS_WATCH(solver.GetNumberOfLevels())
        for (std::size_t lvl = 0; lvl < nlevels; ++lvl)
        {
            LevelType& current_level = dynamic_cast<LevelType&>(solver.GetLevel(lvl));

            #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
            std::cout << "---------------" << std::endl;
            std::cout << "> generating level " << current_level.LevelDepth() << std::endl;
            #endif

            // generate prolongation operator
            #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
            std::cout << "   generating prolongation operator";
            #endif

            typename StructuredMeshMGProjectorType::Pointer pProlongator;
            if (lvl < nlevels-1)
            {
                if (mpProlongationType[lvl] == NULL)
                {
                    std::cout << " (WARNING: the prolongation operator type for level " << lvl << " is not set. User is required to set it separately)";
                    pProlongator = typename StructuredMeshMGProjectorType::Pointer(new StructuredMeshMGProjectorType());
                }
                else
                    pProlongator = boost::dynamic_pointer_cast<StructuredMeshMGProjectorType>(
                        mpProlongationType[lvl]->Create(this->GetModelPart(lvl+1), this->GetModelPart(lvl), block_size));
            }
            else
                pProlongator = typename StructuredMeshMGProjectorType::Pointer(new StructuredMeshMGProjectorType());
            for (unsigned int dim = 0; dim < TDim; ++dim)
            {
                pProlongator->SetFineMeshSize(dim, coarse_div[dim] << (nlevels-1-lvl));
                pProlongator->SetCoarseMeshSize(dim, coarse_div[dim] << (nlevels-2-lvl));
            }
            pProlongator->Initialize();
            current_level.SetProlongationOperator(pProlongator);

            #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
            std::cout << " completed" << std::endl;
            #endif

            // generate restriction operator
            #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
            std::cout << "   generating restriction operator";
            #endif

            typename StructuredMeshMGProjectorType::Pointer pRestrictor;
            if (lvl < nlevels-1)
            {
                if (mpRestrictionType[lvl] == NULL)
                {
                    std::cout << " (WARNING: the restriction operator type for level " << lvl << " is not set. User is required to set it separately)";
                    pRestrictor = typename StructuredMeshMGProjectorType::Pointer(new StructuredMeshMGProjectorType());
                }
                else
                    pRestrictor = boost::dynamic_pointer_cast<StructuredMeshMGProjectorType>(
                        mpRestrictionType[lvl]->Create(this->GetModelPart(lvl+1), this->GetModelPart(lvl), block_size));
            }
            else
                pRestrictor = typename StructuredMeshMGProjectorType::Pointer(new StructuredMeshMGProjectorType());
            for (unsigned int dim = 0; dim < TDim; ++dim)
            {
                pRestrictor->SetFineMeshSize(dim, coarse_div[dim] << (nlevels-1-lvl));
                pRestrictor->SetCoarseMeshSize(dim, coarse_div[dim] << (nlevels-2-lvl));
            }
            pRestrictor->Initialize();
            current_level.SetRestrictionOperator(pRestrictor);

            #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
            std::cout << " completed" << std::endl;
            #endif

            // generate transfer operator
            #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
            std::cout << "   generating transfer operator";
            #endif

            typename StructuredMeshMGProjectorType::Pointer pTransferOperator;
            if (lvl < nlevels-1)
            {
                if (mpTransferType[lvl] == NULL)
                {
                    std::cout << " (WARNING: the transfer operator type for level " << lvl << " is not set. Will use restriction operator for transfer)";
                    pTransferOperator = pRestrictor;
                }
                else
                    pTransferOperator = boost::dynamic_pointer_cast<StructuredMeshMGProjectorType>(
                        mpTransferType[lvl]->Create(this->GetModelPart(lvl+1), this->GetModelPart(lvl), block_size));
            }
            else
                pTransferOperator = typename StructuredMeshMGProjectorType::Pointer(new StructuredMeshMGProjectorType());
            for (unsigned int dim = 0; dim < TDim; ++dim)
            {
                pTransferOperator->SetFineMeshSize(dim, coarse_div[dim] << (nlevels-1-lvl));
                pTransferOperator->SetCoarseMeshSize(dim, coarse_div[dim] << (nlevels-2-lvl));
            }
            pTransferOperator->Initialize();
            current_level.SetTransferOperator(pTransferOperator);

            #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
            std::cout << " completed" << std::endl;
            #endif
        }

        #ifdef DEBUG_MULTILEVEL_SOLVER_FACTORY
        std::cout << "Executing GMGStructuredSolverFactory::" << __FUNCTION__ << " completed" << std::endl;
        std::cout << "#######################" << std::endl;
        #endif
    }

    void GenerateMultilevelSolver(MultilevelSolverType& solver, SparseMatrixType& rA) const override
    {}

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
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "Geometric multigrid solver factory";
        return  buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
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

    std::vector<ModelPart::Pointer> mpModelParts;
    std::vector<typename StructuredMeshMGProjectorType::Pointer> mpProlongationType;
    std::vector<typename StructuredMeshMGProjectorType::Pointer> mpRestrictionType;
    std::vector<typename StructuredMeshMGProjectorType::Pointer> mpTransferType;

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


///@}

}  // namespace Kratos.

#undef DEBUG_MULTILEVEL_SOLVER_FACTORY

#endif // KRATOS_MULTIGRID_SOLVERS_APP_GMG_STRUCTURED_SOLVER_FACTORY_H_INCLUDED  defined
