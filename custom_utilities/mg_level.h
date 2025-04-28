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
template<class TSparseSpaceType, class TDenseSpaceType, class TModelPartType>
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

    typedef typename TSparseSpaceType::VectorPointerType VectorPointerType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    typedef typename TDenseSpaceType::VectorType DenseVectorType;

    typedef Reorderer<TSparseSpaceType, TDenseSpaceType> ReordererType;

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TModelPartType, ReordererType> LinearSolverType;

    typedef typename LinearSolverType::Pointer LinearSolverPointerType;

    typedef MGProjector<TSparseSpaceType> ProjectorType;

    typedef typename ProjectorType::Pointer ProjectorPointerType;

    typedef typename TSparseSpaceType::SizeType SizeType;

    typedef typename TSparseSpaceType::IndexType IndexType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MGLevel(const IndexType lvl)
    : mLevelDepth(lvl)
    {}

    /// Copy constructor. Implement copy constructor is important in order to pass the data to the container (i.e. std::vector)
    MGLevel(const MGLevel& rOther)
    : mLevelDepth(rOther.mLevelDepth)
    , mpPreSmoother(rOther.mpPreSmoother)
    , mpPostSmoother(rOther.mpPostSmoother)
    , mpRestrictor(rOther.mpRestrictor)
    , mpProlongator(rOther.mpProlongator)
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
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /// Apply the pre-smoothing
    virtual int ApplyPreSmoother(VectorType& rX, VectorType& rB) const
    {
        KRATOS_ERROR << "Calling base class function";
        return 1;
    }

    /// Apply the post-smoothing
    virtual int ApplyPostSmoother(VectorType& rX, VectorType& rB) const
    {
        KRATOS_ERROR << "Calling base class function";
        return 1;
    }

    /// Apply the restriction operator
    virtual int ApplyRestriction(VectorType& rX, VectorType& rY) const
    {
        if(mpRestrictor == NULL)
        {
            KRATOS_ERROR << "The restriction operator has not been set for " << Info();
        }
        return mpRestrictor->Apply(rX, rY);
    }

    /// Apply the transfer operator
    virtual int ApplyTransfer(VectorType& rX, VectorType& rY) const
    {
        if(mpTransferOperator == NULL)
        {
            // use the restriction operator to transfer
            return mpRestrictor->Apply(rX, rY);
        }
        else
            return mpTransferOperator->Apply(rX, rY);
    }

    /// Apply the prolongation operator
    virtual int ApplyProlongation(VectorType& rX, VectorType& rY) const
    {
        if(mpProlongator == NULL)
        {
            KRATOS_ERROR << "The prolongation operator has not been set for " << Info();
        }
        return mpProlongator->Apply(rX, rY);
    }

    /// If the level is represented by matrix A, then rY = A*rX
    virtual int Apply(VectorType& rX, VectorType& rY) const
    {
        KRATOS_ERROR << "Calling base class function";
        return 1;
    }

    /// If the level is represented by matrix A, then rY = A^-1*rX
    virtual int Inverse(LinearSolverPointerType pCoarseSolver, VectorType& rX, VectorType& rY) const
    {
        KRATOS_ERROR << "Calling base class function";
        return 1;
    }

    ///@}
    ///@name Access
    ///@{

    LinearSolverPointerType PreSmoother() const
    {
        return mpPreSmoother;
    }

    LinearSolverPointerType PostSmoother() const
    {
        return mpPostSmoother;
    }

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

    void SetTransferOperator(ProjectorPointerType pProjector)
    {
        mpTransferOperator = pProjector;
    }

    void SetProlongationOperator(ProjectorPointerType pProjector)
    {
        mpProlongator = pProjector;
    }

    IndexType LevelDepth() const
    {
        return mLevelDepth;
    }

    ///@}
    ///@name Inquiry
    ///@{

    /// Get the size of the lower level matrix
    virtual SizeType GetCoarseSize() const
    {
        return mpRestrictor->GetProjectedSize();
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream ss;
        ss << "Multigrid Level depth: " << LevelDepth() << std::endl;

        if (mpPreSmoother != NULL)
            ss << "  PreSmoother: " << mpPreSmoother->Info() << std::endl;
        else
            ss << "  PreSmoother: null" << std::endl;

        if (mpPostSmoother != NULL)
            ss << "  PostSmoother: " << mpPostSmoother->Info() << std::endl;
        else
            ss << "  PostSmoother: null" << std::endl;

        if (mpRestrictor != NULL)
            ss << "  Restrictor: " << mpRestrictor->Info() << std::endl;
        else
            ss << "  Restrictor: null" << std::endl;

        if (mpProlongator != NULL)
            ss << "  Prolongator: " << mpProlongator->Info() << std::endl;
        else
            ss << "  Prolongator: null" << std::endl;

        if (mpRestrictor != NULL)
            ss << "  Coarse matrix size: " << this->GetCoarseSize();
        else
            ss << "  Coarse matrix size: 0";

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
    ProjectorPointerType mpTransferOperator;

    IndexType mLevelDepth;

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

};

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
template<class TSparseSpaceType, class TDenseSpaceType, class TModelPartType>
inline std::istream& operator >> (std::istream& IStream, MGLevel<TSparseSpaceType, TDenseSpaceType, TModelPartType>& rThis)
{
    return IStream;
}

/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType, class TModelPartType>
inline std::ostream& operator << (std::ostream& rOStream, const MGLevel<TSparseSpaceType, TDenseSpaceType, TModelPartType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

} // namespace Kratos.

#endif // KRATOS_MULTIGRID_LEVEL_H_INCLUDED  defined
