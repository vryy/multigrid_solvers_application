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


#if !defined(KRATOS_MULTIGRID_LEVEL_H_INCLUDED )
#define  KRATOS_MULTIGRID_LEVEL_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <sstream>
#include <cstddef>


// External includes


// Project includes
#include "includes/define.h"
#include "linear_solvers/linear_solver.h"
#include "custom_utilities/mg_projector.h"


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
 * Abstract class for a level in mutigrid hierarchy
 */
template<class TSparseSpaceType, class TDenseSpaceType>
class MGLevel
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MGLevel
    KRATOS_CLASS_POINTER_DEFINITION(MGLevel);

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::MatrixPointerType SparseMatrixPointerType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    typedef typename TDenseSpaceType::VectorType DenseVectorType;

    typedef Reorderer<TSparseSpaceType, TDenseSpaceType> ReordererType;

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, ReordererType> LinearSolverType;

    typedef typename LinearSolverType::Pointer LinearSolverPointerType;

    typedef MGProjector<TSparseSpaceType> ProjectorType;

    typedef typename ProjectorType::Pointer ProjectorPointerType;

    typedef std::size_t  SizeType;

    typedef unsigned int  IndexType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MGLevel(const IndexType& lvl) : mLevelDepth(lvl)
    {
        this->Initialize();
        KRATOS_WATCH(mpA)
    }

    /// Copy constructor. Implement copy constructor is important in order to pass the data to the container (i.e. std::vector)
    MGLevel(const MGLevel& rOther)
    : mLevelDepth(rOther.mLevelDepth)
    , mpPreSmoother(rOther.mpPreSmoother)
    , mpPostSmoother(rOther.mpPostSmoother)
    , mpRestrictor(rOther.mpRestrictor)
    , mpProlongator(rOther.mpProlongator)
    // , mpA(rOther.mpA)
    {}

    /// Destructor.
    virtual ~MGLevel()
    {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator. It's also important like the Copy constructor
    MGLevel& operator= (const MGLevel& rOther)
    {
        mLevelDepth = rOther.mLevelDepth;
        mpPreSmoother = rOther.mpPreSmoother;
        mpPostSmoother = rOther.mpPostSmoother;
        mpRestrictor = rOther.mpRestrictor;
        mpProlongator = rOther.mpProlongator;
        // mpA = rOther.mpA;
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    virtual void ApplyPreSmoother(SparseMatrixType& rA, VectorType& rX, VectorType& rB) const
    {
        if(mpPreSmoother == NULL)
        {
            std::stringstream ss;
            ss << "The pre-smoother has not been set for " << Info();
            KRATOS_THROW_ERROR(std::logic_error, ss.str(), "");
        }
        mpPreSmoother->Solve(rA, rX, rB);
    }

    virtual void ApplyPostSmoother(SparseMatrixType& rA, VectorType& rX, VectorType& rB) const
    {
        if(mpPostSmoother == NULL)
        {
            std::stringstream ss;
            ss << "The post-smoother has not been set for " << Info();
            KRATOS_THROW_ERROR(std::logic_error, ss.str(), "");
        }
        mpPostSmoother->Solve(rA, rX, rB);
    }

    virtual void ApplyRestriction(VectorType& rX, VectorType& rY) const
    {
        if(mpRestrictor == NULL)
        {
            std::stringstream ss;
            ss << "The restriction operator has not been set for " << Info();
            KRATOS_THROW_ERROR(std::logic_error, ss.str(), "");
        }
        mpRestrictor->Apply(rX, rY);
    }

    virtual void ApplyProlongation(VectorType& rX, VectorType& rY) const
    {
        if(mpProlongator == NULL)
        {
            std::stringstream ss;
            ss << "The prolongation operator has not been set for " << Info();
            KRATOS_THROW_ERROR(std::logic_error, ss.str(), "");
        }
        mpProlongator->Apply(rX, rY);
    }

    ///@}
    ///@name Access
    ///@{

    void SetPreSmoother(LinearSolverPointerType pPreSmoother)
    {
        mpPreSmoother = pPreSmoother;
    }

    void SetPostSmoother(LinearSolverPointerType pPostSmoother)
    {
        mpPostSmoother = pPostSmoother;
    }

    void SetRestrictionOperator(ProjectorPointerType pProjector)
    {
        mpRestrictor = pProjector;
    }

    void SetProlongationOperator(ProjectorPointerType pProjector)
    {
        mpProlongator = pProjector;
    }

    IndexType LevelDepth() const
    {
        return mLevelDepth;
    }

    void SetLevelDepth(const IndexType& lvl)
    {
        mLevelDepth = lvl;
    }

    // this is kept to make compatible with python; should not use it to avoid copying memory
    void SetCoarseMatrix(SparseMatrixType& A)
    {
        mpA = boost::make_shared<SparseMatrixType>(A);
    }

    SparseMatrixPointerType GetCoarseMatrix() const
    {
        return mpA;
    }

    ///@}
    ///@name Inquiry
    ///@{

    /// Get the size of the coarse matrix
    virtual SizeType GetCoarseSize() const
    {
        return TSparseSpaceType::Size1(*mpA);
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream ss;
        ss << "Multigrid Level depth: " << LevelDepth() << std::endl;
        ss << "  PreSmoother: " << mpPreSmoother->Info() << std::endl;
        ss << "  PostSmoother: " << mpPostSmoother->Info() << std::endl;
        ss << "  Restrictor: " << mpRestrictor->Info() << std::endl;
        ss << "  Prolongator: " << mpProlongator->Info() << std::endl;
        ss << "  Coarse matrix size: " << this->GetCoarseSize();
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

    LinearSolverPointerType mpPreSmoother;
    LinearSolverPointerType mpPostSmoother;

    ProjectorPointerType mpProlongator;
    ProjectorPointerType mpRestrictor;

    IndexType mLevelDepth;

    SparseMatrixPointerType mpA;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void Initialize()
    {
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
inline std::istream& operator >> (std::istream& IStream, MGLevel<TSparseSpaceType, TDenseSpaceType>& rThis)
{
    return IStream;
}

/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType>
inline std::ostream& operator << (std::ostream& rOStream, const MGLevel<TSparseSpaceType, TDenseSpaceType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


} // namespace Kratos.

#endif // KRATOS_MULTIGRID_LEVEL_H_INCLUDED  defined

