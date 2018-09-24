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
//   Date:                $Date: 15/7/2018 $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(KRATOS_MATRIX_BASED_MULTIGRID_LEVEL_H_INCLUDED )
#define  KRATOS_MATRIX_BASED_MULTIGRID_LEVEL_H_INCLUDED



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
class MatrixBasedMGLevel : public MGLevel<TSparseSpaceType, TDenseSpaceType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MatrixBasedMGLevel
    KRATOS_CLASS_POINTER_DEFINITION(MatrixBasedMGLevel);

    typedef MGLevel<TSparseSpaceType, TDenseSpaceType> BaseType;

    typedef typename BaseType::SparseMatrixType SparseMatrixType;

    typedef typename BaseType::SparseMatrixPointerType SparseMatrixPointerType;

    typedef typename BaseType::LinearSolverPointerType LinearSolverPointerType;

    typedef typename BaseType::VectorType VectorType;

    typedef typename BaseType::VectorPointerType VectorPointerType;

    typedef typename BaseType::IndexType IndexType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MatrixBasedMGLevel(const IndexType& lvl)
    : BaseType(lvl)
    {
        this->Initialize();
    }

    /// Copy constructor. Implement copy constructor is important in order to pass the data to the container (i.e. std::vector)
    MatrixBasedMGLevel(const MatrixBasedMGLevel& rOther)
    : BaseType()
    , mpA(rOther.mpA) // shall we make a deep copy here?
    {}

    /// Destructor.
    virtual ~MatrixBasedMGLevel()
    {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator. It's also important like the Copy constructor
    MatrixBasedMGLevel& operator= (const MatrixBasedMGLevel& rOther)
    {
        BaseType::operator=(rOther);
        mpA = rOther.mpA; // shall we make a deep copy here
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    virtual int ApplyPreSmoother(VectorType& rX, VectorType& rB) const
    {
        if(BaseType::PreSmoother() == NULL)
        {
            std::stringstream ss;
            ss << "The pre-smoother has not been set for " << Info();
            KRATOS_THROW_ERROR(std::logic_error, ss.str(), "");
        }
        int stat = !(BaseType::PreSmoother()->Solve(*mpA, rX, rB));
/*        KRATOS_WATCH(*(BaseType::PreSmoother()))*/
        return stat;
    }

    virtual int ApplyPostSmoother(VectorType& rX, VectorType& rB) const
    {
        if(BaseType::PostSmoother() == NULL)
        {
            std::stringstream ss;
            ss << "The post-smoother has not been set for " << Info();
            KRATOS_THROW_ERROR(std::logic_error, ss.str(), "");
        }
        return !(BaseType::PostSmoother()->Solve(*mpA, rX, rB));
    }

    virtual int Apply(VectorType& rX, VectorType& rY) const
    {
        if(mpA == NULL)
        {
            std::stringstream ss;
            ss << "The matrix operator has not been set for " << Info();
            KRATOS_THROW_ERROR(std::logic_error, ss.str(), "");
        }
        TSparseSpaceType::Mult(*mpA, rX, rY);
        return 0;
    }

    virtual int Inverse(LinearSolverPointerType pCoarseSolver, VectorType& rX, VectorType& rY) const
    {
        if(mpA == NULL)
        {
            std::stringstream ss;
            ss << "The matrix operator has not been set for " << Info();
            KRATOS_THROW_ERROR(std::logic_error, ss.str(), "");
        }
        bool stat = pCoarseSolver->Solve(*mpA, rX, rY);
        std::cout << *pCoarseSolver << std::endl;
        return !stat;
    }

    ///@}
    ///@name Access
    ///@{

    // this is kept to make compatible with python; should not use it to avoid copying memory
    void SetCoarseMatrix(SparseMatrixType& A)
    {
        mpA = boost::make_shared<SparseMatrixType>(A);
    }

    SparseMatrixPointerType GetCoarseMatrix()
    {
        return mpA;
    }

    VectorPointerType GetCoarseUpdateVector()
    {
        return mpDx;
    }

    VectorPointerType GetCoarseVector()
    {
        return mpb;
    }

    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream ss;
        ss << "Matrix-based " << BaseType::Info() << std::endl;
        ss << "  Fine matrix size: " << TSparseSpaceType::Size1(*mpA) << ", nonzeros = " << mpA->filled2() << std::endl;
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
        // rOStream << "Coarse Matrix: " << *mpA;
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

    SparseMatrixPointerType mpA;
    VectorPointerType mpDx;
    VectorPointerType mpb;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void Initialize()
    {
        if (mpA == NULL)
        {
            SparseMatrixPointerType pNewA = TSparseSpaceType::CreateEmptyMatrixPointer();
            mpA.swap(pNewA);
        }
        if (mpDx == NULL)
        {
            VectorPointerType pNewDx = TSparseSpaceType::CreateEmptyVectorPointer();
            mpDx.swap(pNewDx);
        }
        if (mpb == NULL)
        {
            VectorPointerType pNewb = TSparseSpaceType::CreateEmptyVectorPointer();
            mpb.swap(pNewb);
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
inline std::istream& operator >> (std::istream& IStream, MatrixBasedMGLevel<TSparseSpaceType, TDenseSpaceType>& rThis)
{
    return IStream;
}

/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType>
inline std::ostream& operator << (std::ostream& rOStream, const MatrixBasedMGLevel<TSparseSpaceType, TDenseSpaceType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


} // namespace Kratos.

#endif // KRATOS_MATRIX_BASED_MULTIGRID_LEVEL_H_INCLUDED  defined

