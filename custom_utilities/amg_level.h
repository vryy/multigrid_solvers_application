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
//   Date:                $Date: 2013 Jan 9 10:53:00 $
//   Revision:            $Revision: 1.1 $
//
//


#if !defined(KRATOS_AMG_LEVEL_H_INCLUDED )
#define  KRATOS_AMG_LEVEL_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <sstream>
#include <cstddef>


// External includes


// Project includes
#include "includes/define.h"
#include "linear_solvers/linear_solver.h"
#include "linear_solvers/reorderer.h"
#include "includes/ublas_interface.h"

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
class AMGLevel
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of AMGLevel
    KRATOS_CLASS_POINTER_DEFINITION(AMGLevel);

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::MatrixPointerType SparseMatrixPointerType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    typedef typename TDenseSpaceType::VectorType DenseVectorType;

    typedef Reorderer<TSparseSpaceType, TDenseSpaceType> ReordererType;

    typedef typename LinearSolver<TSparseSpaceType, TDenseSpaceType, ReordererType>::Pointer LinearSolverPointerType;

    typedef std::size_t  SizeType;
    
    typedef unsigned int  IndexType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    AMGLevel() : mLevelDepth(0)
    {
        Initialize();
    }
    
    AMGLevel(LinearSolverPointerType pPreSmootherSolver, LinearSolverPointerType pPostSmootherSolver) : mLevelDepth(0)
    {
        mpPreSmootherSolver = pPreSmootherSolver;
        mpPostSmootherSolver = pPostSmootherSolver;
        Initialize();
    }
    
    /// Copy constructor. Implement copy constructor is important in order to pass the data to the container (i.e. std::vector)
    AMGLevel(const AMGLevel& Other) : mpPreSmootherSolver(Other.mpPreSmootherSolver),
                                      mpPostSmootherSolver(Other.mpPostSmootherSolver),
                                      mLevelDepth(Other.mLevelDepth),
                                      mpR(Other.mpR),
                                      mpA(Other.mpA),
                                      mpP(Other.mpP)
    {
    }
    
    /// Destructor.
    virtual ~AMGLevel() {}
    
    
    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator. It's also important like the Copy constructor
    AMGLevel& operator= (const AMGLevel& Other)
    {
        mpPreSmootherSolver = Other.mpPreSmootherSolver;
        mpPostSmootherSolver = Other.mpPostSmootherSolver;
        mLevelDepth = Other.mLevelDepth;
        mpR = Other.mpR;
        mpA = Other.mpA;
        mpP = Other.mpP;
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    void ApplyPreSmoother(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
    {
        mpPreSmootherSolver->Solve(rA, rX, rB);
    }

    void ApplyPostSmoother(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
    {
        mpPostSmootherSolver->Solve(rA, rX, rB);
    }
    
    void ApplyRestriction(VectorType& rX, VectorType& rY)
    {
        if(mpR == NULL)
        {
            std::stringstream ss;
            ss << "The restriction operator has not been set for " << Info();
            KRATOS_THROW_ERROR(std::logic_error, ss.str(), "");
        }
        TSparseSpaceType::Mult(*mpR, rX, rY);
    }
    
    void ApplyProlongation(VectorType& rX, VectorType& rY)
    {
        if(mpP == NULL)
        {
            std::stringstream ss;
            ss << "The prolongation operator has not been set for " << Info();
            KRATOS_THROW_ERROR(std::logic_error, ss.str(), "");
        }
        TSparseSpaceType::Mult(*mpP, rX, rY);
    }
    
    ///@}
    ///@name Access
    ///@{

    void SetPreSmoother(LinearSolverPointerType pPreSmootherSolver)
    {
        mpPreSmootherSolver = pPreSmootherSolver;
    }

    void SetPostSmoother(LinearSolverPointerType pPostSmootherSolver)
    {
        mpPostSmootherSolver = pPostSmootherSolver;
    }
    
//    void SetRestrictionOperator(SparseMatrixType& R)
//    {
//        mpR = boost::make_shared<SparseMatrixType>(R);
//    }
//    
//    void SetProlongationOperator(SparseMatrixType& P)
//    {
//        mpP = boost::make_shared<SparseMatrixType>(P);
//    }
//    
//    void SetCoarsenMatrix(SparseMatrixType& A)
//    {
//        KRATOS_WATCH("at amg_level");
//        KRATOS_WATCH(&A);
//        mpA = boost::make_shared<SparseMatrixType>(A);
//        KRATOS_WATCH(&(*mpA));
//    }
    
//    void ComputeCoarsenMatrix(SparseMatrixPointerType& R, SparseMatrixPointerType& A, SparseMatrixPointerType& P)
//    {
////        KRATOS_WATCH("ComputeCoarsenMatrix");
////        using namespace boost::numeric::ublas;
////        SparseMatrixPointerType C;
////        noalias(*C) = prod(*R, SparseMatrixType(prod(*A, *P)));
////        SparseMatrixType T = prod(*A, *P);
////        SparseMatrixType C = prod(*R, T);

//        SparseMatrixType C = prod(*R, SparseMatrixType(prod(*A, *P)));
////        KRATOS_WATCH(C);
////        mpA->swap(C);
////        mpA = SparseMatrixPointerType(C);
////        SparseMatrixPointerType pC(C);
////        mpA.swap(pC);
//        mpA = boost::make_shared<SparseMatrixType>(C);
//    }
    
    void SetLevelDepth(IndexType lvl)
    {
        mLevelDepth = lvl;
    }
    
    ///@}
    ///@name Inquiry
    ///@{

    // this is kept to make compatible with python; should not use it to avoid copying memory
    void SetRestrictionOperator(SparseMatrixType& R)
    {
        mpR = boost::make_shared<SparseMatrixType>(R);
    }
    
    // this is kept to make compatible with python; should not use it to avoid copying memory
    void SetProlongationOperator(SparseMatrixType& P)
    {
        mpP = boost::make_shared<SparseMatrixType>(P);
    }
    
    // this is kept to make compatible with python; should not use it to avoid copying memory
    void SetCoarsenMatrix(SparseMatrixType& A)
    {
        mpA = boost::make_shared<SparseMatrixType>(A);
    }
    
    SparseMatrixPointerType GetRestrictionOperator()
    {
        return mpR;
    }

    SparseMatrixPointerType GetProlongationOperator()
    {
        return mpP;
    }

    SparseMatrixPointerType GetCoarsenMatrix()
    {
        return mpA;
    }

    SizeType GetCoarsenSize()
    {
        return TSparseSpaceType::Size1(*mpR);
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream ss;
        ss << "AMG Level depth " << mLevelDepth << std::endl;
        ss << "PreSmoother: " << (*mpPreSmootherSolver).Info() << std::endl;
        ss << "PostSmoother: " << (*mpPostSmootherSolver).Info();
        return ss.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
//        KRATOS_WATCH(*mpR);
//        KRATOS_WATCH(*mpA);
//        KRATOS_WATCH(*mpP);
//        rOStream << *mpR << *mpA << *mpP << std::endl;
//        rOStream << "mpR = " << mpR << ", mpA = " << mpA << ", mpP = " << mpP << std::endl;
//        rOStream << "*mpR = " << &(*mpR) << ", *mpA = " << &(*mpA) << ", *mpP = " << &(*mpP) << std::endl;

        rOStream << "mpR = (" << TSparseSpaceType::Size1(*mpR) << "," << TSparseSpaceType::Size2(*mpR) << "), ";
        rOStream << "mpA = (" << TSparseSpaceType::Size1(*mpA) << "," << TSparseSpaceType::Size2(*mpA) << "), ";
        rOStream << "mpP = (" << TSparseSpaceType::Size1(*mpP) << "," << TSparseSpaceType::Size2(*mpP) << "), ";
        rOStream << "nonzeros = " << (*mpA).filled2();
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
    LinearSolverPointerType mpPreSmootherSolver;
    LinearSolverPointerType mpPostSmootherSolver;

    SparseMatrixPointerType mpA;
    SparseMatrixPointerType mpP;
    SparseMatrixPointerType mpR;
    
    IndexType mLevelDepth;
    
    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{
    
    
    void Initialize()
    {
        if(mpR == NULL)
        {
            SparseMatrixPointerType pNewR = SparseMatrixPointerType(new SparseMatrixType(0, 0));
            mpR.swap(pNewR);
        }
        
        if(mpP == NULL)
        {
            SparseMatrixPointerType pNewP = SparseMatrixPointerType(new SparseMatrixType(0, 0));
            mpP.swap(pNewP);
        }
        
        if(mpA == NULL)
        {
            SparseMatrixPointerType pNewA = SparseMatrixPointerType(new SparseMatrixType(0, 0));
            mpA.swap(pNewA);
        }
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

};

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TSparseSpaceType, class TDenseSpaceType>
inline std::istream& operator >> (std::istream& IStream, AMGLevel<TSparseSpaceType, TDenseSpaceType>& rThis)
{
    return IStream;
}

/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType>
inline std::ostream& operator << (std::ostream& rOStream, const AMGLevel<TSparseSpaceType, TDenseSpaceType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

    
} // namespace Kratos.

#endif // KRATOS_AMG_LEVEL_H_INCLUDED  defined

