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
//   Date:                $Date: 27 Sep 2018 $
//   Revision:            $Revision: 1.1 $
//
//


#if !defined(KRATOS_VANKA_SMOOTHER_H_INCLUDED )
#define  KRATOS_VANKA_SMOOTHER_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <sstream>
#include <cstddef>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "linear_solvers/reorderer.h"
#include "linear_solvers/linear_solver.h"
#include "utilities/openmp_utils.h"
#include "utilities/progress.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "structural_application/custom_utilities/sd_math_utils.h"

// #define ENABLE_PROFILING

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

/// Class for Vanka smoother.
/**
 * TODO add some information or reference
*/
template<class TSparseSpaceType, class TDenseSpaceType,
         class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class VankaSmoother : public LinearSolver<TSparseSpaceType, TDenseSpaceType, ModelPart, TReordererType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of LinearSolver
    KRATOS_CLASS_POINTER_DEFINITION(VankaSmoother);

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, ModelPart, TReordererType> BaseType;

    typedef typename BaseType::SparseMatrixType SparseMatrixType;

    typedef typename BaseType::VectorType VectorType;

    typedef typename BaseType::DenseMatrixType DenseMatrixType;

    typedef typename BaseType::DenseVectorType DenseVectorType;

    typedef typename BaseType::SizeType SizeType;
    typedef typename BaseType::IndexType IndexType;

    typedef ModelPart::ElementsContainerType ElementsArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    VankaSmoother(ModelPart::Pointer pModelPart)
    : mpModelPart(pModelPart), mLambda(0.5)
    {
    }

    VankaSmoother(ModelPart::Pointer pModelPart, const double& lambda)
    : mpModelPart(pModelPart), mLambda(lambda)
    {
    }

    /// Copy constructor.
    VankaSmoother(const VankaSmoother& Other)
    : BaseType(Other), mpModelPart(Other.mpModelPart), mLambda(Other.mLambda)
    {
    }

    /// Destructor.
    ~VankaSmoother() override {}

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    VankaSmoother& operator=(const VankaSmoother& Other)
    {
        BaseType::operator=(Other);
        mpModelPart = Other.mpModelPart;
        mLambda = Other.mLambda;
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /** This function is designed to be called as few times as possible. It creates the data structures
     * that only depend on the connectivity of the matrix (and not on its coefficients)
     * so that the memory can be allocated once and expensive operations can be done only when strictly
     * needed
    @param rA. System matrix
    @param rX. Solution vector. it's also the initial guess for iterative linear solvers.
    @param rB. Right hand side vector.
    */
    void Initialize(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        BaseType::Initialize(rA, rX, rB);
    }

    /** Normal solve method.
    Solves the linear system Ax=b and puts the result on SystemVector& rX.
    rVectorx is also th initial guess for iterative methods.
    @param rA. System matrix
    @param rX. Solution vector. it's also the initial
    guess for iterative linear solvers.
     @param rB. Right hand side vector.
    */
    bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        if(this->IsNotConsistent(rA, rX, rB))
            return false;

        ProcessInfo& CurrentProcessInfo = mpModelPart->GetProcessInfo();
        Element::EquationIdVectorType EquationId;
        DenseMatrixType LocalMatrix, InvLocalMatrix;
        VectorType rLocal, xLocal;
        VectorType r(rB.size()), aux(rB.size());

        #ifdef ENABLE_PROFILING
        Kratos::progress_display show_progress(mpModelPart->Elements().size());
        double cur_time;
        double extract_local_time = 0.0;
        double update_local_time = 0.0;
        double local_solve_time = 0.0;
        #endif

        // compute the residual
        noalias(r) = rB;
        TSparseSpaceType::Mult(rA, rX, aux);
        TSparseSpaceType::UnaliasedAdd(r, -1.0, aux);

        for (typename ElementsArrayType::ptr_iterator it = mpModelPart->Elements().ptr_begin(); it != mpModelPart->Elements().ptr_end(); ++it)
        {
            // extract the rows corresponding to the element
            (*it)->EquationIdVector(EquationId, CurrentProcessInfo);

            #ifdef ENABLE_PROFILING
            cur_time = OpenMPUtils::GetCurrentTime();
            #endif

            // extract the local residual
            if (rLocal.size() != EquationId.size())
                rLocal.resize(EquationId.size(), false);
            for (std::size_t i = 0; i < EquationId.size(); ++i)
                rLocal(i) = r(EquationId[i]);

            // extract the local matrix
            if (LocalMatrix.size1() != EquationId.size())
            {
                LocalMatrix.resize(EquationId.size(), EquationId.size(), false);
                InvLocalMatrix.resize(EquationId.size(), EquationId.size(), false);
                xLocal.resize(EquationId.size(), false);
            }

            for (std::size_t i = 0; i < EquationId.size(); ++i)
                for (std::size_t j = 0; j < EquationId.size(); ++j)
                    LocalMatrix(i, j) = rA(EquationId[i], EquationId[j]);

            #ifdef ENABLE_PROFILING
            extract_local_time += (OpenMPUtils::GetCurrentTime() - cur_time);
            cur_time = OpenMPUtils::GetCurrentTime();
            #endif

            // solve the local system
            SD_MathUtils<double>::InvertMatrix(LocalMatrix, InvLocalMatrix);
            noalias(xLocal) = prod(InvLocalMatrix, rLocal);

            #ifdef ENABLE_PROFILING
            local_solve_time += (OpenMPUtils::GetCurrentTime() - cur_time);
            cur_time = OpenMPUtils::GetCurrentTime();
            #endif

            // update the global solution
            for (std::size_t i = 0; i < EquationId.size(); ++i)
                rX(EquationId[i]) += mLambda*xLocal(i);

            // update the residual
            noalias(rLocal) = prod(LocalMatrix, xLocal);
            for (std::size_t i = 0; i < EquationId.size(); ++i)
                r(EquationId[i]) -= mLambda*rLocal(i);

            #ifdef ENABLE_PROFILING
            update_local_time += (OpenMPUtils::GetCurrentTime() - cur_time);
            ++show_progress;
            #endif
        }

        #ifdef ENABLE_PROFILING
        KRATOS_WATCH(extract_local_time)
        KRATOS_WATCH(local_solve_time)
        KRATOS_WATCH(update_local_time)
        #endif

        return true;
    }

    /** Multi solve method for solving a set of linear systems with same coefficient matrix.
    Solves the linear system Ax=b and puts the result on SystemVector& rX.
    rVectorx is also th initial guess for iterative methods.
    @param rA. System matrix
    @param rX. Solution vector. it's also the initial
    guess for iterative linear solvers.
     @param rB. Right hand side vector.
    */
    bool Solve(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB) override
    {
        return false;
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
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "Vanka smoother, model_part " << mpModelPart->Name() << ", lambda = " << mLambda;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
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

    double mLambda; // stabilization parameter
    ModelPart::Pointer mpModelPart;

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

}; // Class VankaSmoother

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

}  // namespace Kratos.

#undef ENABLE_PROFILING

#endif // KRATOS_VANKA_SMOOTHER_H_INCLUDED  defined
