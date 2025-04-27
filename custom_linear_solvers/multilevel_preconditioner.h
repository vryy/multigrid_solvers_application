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
//   Date:                $Date: 2013-14-1 9:34:00 $
//   Revision:            $Revision: 1.1 $
//
//

#if !defined(KRATOS_BUI_MULTILEVEL_PRECONDITIONER_H_INCLUDED)
#define KRATOS_BUI_MULTILEVEL_PRECONDITIONER_H_INCLUDED




// System includes



// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include <boost/numeric/ublas/vector.hpp>
#include "utilities/openmp_utils.h"
#include "custom_linear_solvers/multilevel_solver.h"



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

///@name  Preconditioners
///@{

/// MultilevelPreconditioner class.
/**   */
template<class TSparseSpaceType, class TDenseSpaceType>
class MultilevelPreconditioner : public Preconditioner<TSparseSpaceType, TDenseSpaceType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Block2PhasePreconditioner
    KRATOS_CLASS_POINTER_DEFINITION (MultilevelPreconditioner);

    typedef Preconditioner<TSparseSpaceType, TDenseSpaceType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::MatrixPointerType SparseMatrixPointerType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    typedef MultilevelSolver<TSparseSpaceType, TDenseSpaceType> MultilevelSolverType;

    typedef typename MultilevelSolverType::Pointer MultilevelSolverPointerType;

    typedef typename TSparseSpaceType::SizeType SizeType;

    typedef typename TSparseSpaceType::IndexType IndexType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MultilevelPreconditioner(MultilevelSolverPointerType ml_solver)
    {
        mml_solver = ml_solver;
        std::cout << "invoking multiphase_application/multilevel_preconditioner" << std::endl;
    }


    /// Copy constructor.
    MultilevelPreconditioner(const MultilevelPreconditioner& Other)
    {
        mml_solver = Other.mml_solver;
    }

    /// Destructor.
    ~MultilevelPreconditioner() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    MultilevelPreconditioner& operator=(const MultilevelPreconditioner& Other)
    {
        mml_solver = Other.mml_solver;
        return *this;
    }

    void SetMultilevelSolver(MultilevelSolverPointerType ml_solver)
    {
        mml_solver = ml_solver;
    }

    ///@}
    ///@name Operations
    ///@{

    void Initialize(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        mml_solver->Initialize(rA, rX, rB);

        mml_solver->SetMaxIterationsNumber(1);

        if(mml_solver->mpFactory != NULL)
        {

            mml_solver->mLevels.clear();
            mml_solver->mResidualNorm = 0.00;
            mml_solver->mIterationsNumber = 0;
            mml_solver->mBNorm = 0.00;

            mml_solver->mpFactory->GenerateMultilevelSolver(*mml_solver, rA);

            mml_solver->SetUpSmoothers();

            std::cout << "multilevel structure is generated" << std::endl;

        }
    }

    bool AdditionalPhysicalDataIsNeeded() override
    {
        return false;
    }

    void Mult(SparseMatrixType& rA, VectorType& rX, VectorType& rY) override
    {
        VectorType z = rX;
        TSparseSpaceType::Mult(rA, z, rY);
        ApplyLeft(rY);
    }

    /** calculate preconditioned_X = A^{-1} * X;
        @param rX  Unknows of preconditioner suystem
    */
    VectorType& ApplyLeft(VectorType& rX) override
    {
        SizeType size = TSparseSpaceType::Size(rX);

        VectorType pX(size, 0.00);

        SolveOneStep(pX, rX);

//        KRATOS_WATCH(pX);

        TSparseSpaceType::Copy(pX, rX);

        return rX;
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

    /// Return information about this object.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "MultilevelPreconditioner with ";
        buffer << mml_solver->Info();
        return buffer.str();
    }

    /// Print information about this object.
    void  PrintInfo(std::ostream& OStream) const override
    {
        OStream << Info();
    }

    void PrintData(std::ostream& OStream) const override
    {
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

    MultilevelSolverPointerType mml_solver;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void SolveOneStep(VectorType &rX, VectorType& rB)
    {
        if(mml_solver->GetNumberOfLevels() == 1)
        {
            mml_solver->GetLevel(0).Inverse(mml_solver->mpCoarseSolver, rX, rB);
        }
        else
        {
            mml_solver->RecursiveSolve(0, rX, rB, mml_solver->mCycle);
        }
    }

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Unaccessible methods
    ///@{


    ///@}

}; // Class MultilevelSolver

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}
}  // namespace Kratos.

#endif // KRATOS_BUI_MULTILEVEL_PRECONDITIONER_H_INCLUDED  defined
