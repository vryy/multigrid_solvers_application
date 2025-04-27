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
//   Date:                $Date: 2013 Jan 8 17:59:00 $
//   Revision:            $Revision: 1.1 $
//
//


#if !defined(KRATOS_GAUSS_SEIDEL_ITERATIVE_SOLVER_H_INCLUDED )
#define  KRATOS_GAUSS_SEIDEL_ITERATIVE_SOLVER_H_INCLUDED


// System includes

// External includes
#include "external_includes/pyamg/relaxation.h"

// Project includes
#include "includes/define.h"
#include "linear_solvers/reorderer.h"
#include "linear_solvers/linear_solver.h"
#include "includes/model_part.h"

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

/// Class for SOR Gauss Seidel iterative solver.
/** This class define the general interface for the linear solvers in Kratos.
    There is three template parameter:
    - TSparseSpaceType which specify type
      of the unknowns, coefficients, sparse matrix, vector of
  unknowns, right hand side vector and their respective operators.
    - TDenseMatrixType which specify type of the
      matrices used as temporary matrices or multi solve unknowns and
  right hand sides and their operators.
    - TReordererType which specify type of the Orderer that performs the reordering of matrix to optimize the solution.
*/
template<class TSparseSpaceType, class TDenseSpaceType,
         class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class GaussSeidelIterativeSolver : public LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of LinearSolver
    KRATOS_CLASS_POINTER_DEFINITION(GaussSeidelIterativeSolver);

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    typedef typename TDenseSpaceType::VectorType DenseVectorType;

    typedef std::size_t  SizeType;
    typedef unsigned int  IndexType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    GaussSeidelIterativeSolver() : mMaxIterationsNumber(1), mSweepMode("symmetric")
    {
    }

    GaussSeidelIterativeSolver(const IndexType NewMaxIterationsNumber)
        : mMaxIterationsNumber(NewMaxIterationsNumber), mSweepMode("symmetric")
    {
    }

    GaussSeidelIterativeSolver(const IndexType NewMaxIterationsNumber, const std::string SweepMode)
        : mMaxIterationsNumber(NewMaxIterationsNumber), mSweepMode(SweepMode)
    {
    }

    /// Copy constructor.
    GaussSeidelIterativeSolver(const GaussSeidelIterativeSolver& Other) : BaseType(Other),
        mMaxIterationsNumber(Other.mMaxIterationsNumber),
        mSweepMode(Other.mSweepMode)
    {
    }

    /// Destructor.
    ~GaussSeidelIterativeSolver() override {}

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    GaussSeidelIterativeSolver& operator=(const GaussSeidelIterativeSolver& Other)
    {
        BaseType::operator=(Other);
        mMaxIterationsNumber = Other.mMaxIterationsNumber;
        mSweepMode = Other.mSweepMode;
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

        if(mSweepMode.compare("forward") == 0 || mSweepMode.compare("backward") == 0)
        {
            for(IndexType i = 0; i < mMaxIterationsNumber; i++)
                SolveOneStep(rA, rX, rB, mSweepMode);
        }
        else if (mSweepMode.compare("symmetric") == 0)
        {
            for(IndexType i = 0; i < mMaxIterationsNumber; i++)
            {
                SolveOneStep(rA, rX, rB, "forward");
                SolveOneStep(rA, rX, rB, "backward");
            }
        }
        else
        {
            KRATOS_THROW_ERROR(std::logic_error, "valid sweep directions are 'forward', 'backward', and 'symmetric'", "")
        }

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
        KRATOS_ERROR << "This solver does not support multisolve";
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
        buffer << "Gauss Seidel iterative solver (" << mSweepMode << " mode)";
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

    IndexType mMaxIterationsNumber;
    std::string mSweepMode;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void SolveOneStep(SparseMatrixType& rA, VectorType& rX, VectorType& rB, const std::string& Sweep)
    {
        IndexType row_start, row_stop, row_step;

        if(Sweep.compare("forward") == 0)
        {
            row_start = 0;
            row_stop = TSparseSpaceType::Size(rX);
            row_step = 1;
        }
        else if (Sweep.compare("backward") == 0)
        {
            row_start = TSparseSpaceType::Size(rX) - 1;
            row_stop = -1;
            row_step = -1;
        }

        gauss_seidel(rA.index1_data(), rA.index2_data(), rA.value_data(), rX, rB, row_start, row_stop, row_step);
    }

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

}; // Class GaussSeidelIterativeSolver

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

}  // namespace Kratos.

#endif // KRATOS_GAUSS_SEIDEL_ITERATIVE_SOLVER_H_INCLUDED  defined
