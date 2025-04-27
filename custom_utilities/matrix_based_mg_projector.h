/*
see multigrid_solvers_application/LICENSE.txt
*/

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 15/7/2018 $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(KRATOS_MULTIGRID_SOLVERS_APP_MATRIX_BASED_MG_PROJECTOR_H_INCLUDED )
#define  KRATOS_MULTIGRID_SOLVERS_APP_MATRIX_BASED_MG_PROJECTOR_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <sstream>


// External includes


// Project includes
#include "includes/define.h"
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
 * Class for prolongator and restrictor that uses matrix multiplication to apply the projection.
 * This can be used for both geometric multigrid and algebraic multigrid.
 */
template<class TSpaceType>
class MatrixBasedMGProjector : public MGProjector<TSpaceType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MatrixBasedMGProjector
    KRATOS_CLASS_POINTER_DEFINITION(MatrixBasedMGProjector);

    typedef MGProjector<TSpaceType> BaseType;

    typedef typename BaseType::MatrixType MatrixType;

    typedef typename BaseType::MatrixPointerType MatrixPointerType;

    typedef typename BaseType::VectorType VectorType;

    typedef typename BaseType::VectorPointerType VectorPointerType;

    typedef typename BaseType::SizeType SizeType;

    typedef typename BaseType::IndexType IndexType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    MatrixBasedMGProjector()
    {
        this->Initialize();
    }

    /// Constructor with size
    MatrixBasedMGProjector(const std::size_t nrows, const std::size_t ncols)
    {
        this->Initialize();
        TSpaceType::Resize(*mpOperator, nrows, ncols);
        TSpaceType::SetToZero(*mpOperator);
    }

    /// Destructor.
    ~MatrixBasedMGProjector() override
    {}

    /// Copy constructor
    MatrixBasedMGProjector(const MatrixBasedMGProjector& rOther)
    : BaseType(), mpOperator(rOther.mpOperator)
    {}

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator. It's also important like the Copy constructor
    MatrixBasedMGProjector& operator= (const MatrixBasedMGProjector& rOther)
    {
        this->mpOperator = rOther.mpOperator;
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /// Initialize the operator
    void Initialize() override
    {
        if(mpOperator == nullptr)
        {
            MatrixPointerType pNewP = TSpaceType::CreateEmptyMatrixPointer();
            mpOperator.swap(pNewP);
        }
    }

    /// Set the operator
    void SetOperator(MatrixPointerType pOperator)
    {
        mpOperator = pOperator;
    }

    /// populate the transformation matrix
    void AssembleOperator(const std::vector<std::size_t>& rows,
        const std::vector<std::size_t>& columns, const MatrixType& values)
    {
        // size check
        if (rows.size() != TSpaceType::Size1(values))
            KRATOS_ERROR << "The row size of the matrix is incompatible";

        if (columns.size() != TSpaceType::Size2(values))
            KRATOS_ERROR << "The column size of the matrix is incompatible";

        // populate values
        for (std::size_t i = 0; i < rows.size(); ++i)
        {
            for (std::size_t j = 0; j < columns.size(); ++j)
            {
                if (values(i, j) != 0.0)
                {
                    (*mpOperator)(rows[i], columns[j]) = values(i, j);
                }
            }
        }
    }

    /// populate the transformation matrix by block
    void AssembleOperator(const std::vector<std::size_t>& rows,
        const std::vector<std::size_t>& columns, const MatrixType& values,
        const std::size_t& block_size)
    {
        // size check
        if (rows.size() != TSpaceType::Size1(values))
            KRATOS_ERROR << "The row size of the matrix is incompatible";

        if (columns.size() != TSpaceType::Size2(values))
            KRATOS_ERROR << "The column size of the matrix is incompatible";

        // populate values
        for (std::size_t i = 0; i < rows.size(); ++i)
        {
            for (std::size_t j = 0; j < columns.size(); ++j)
            {
                if (values(i, j) != 0.0)
                {
                    for (std::size_t s = 0; s < block_size; ++s)
                        (*mpOperator)(rows[i]*block_size+s, columns[j]*block_size+s) = values(i, j);
                }
            }
        }
    }

    /// Apply the projection
    int Apply(VectorType& rX, VectorType& rY) const override
    {
        if(mpOperator == nullptr)
        {
            KRATOS_ERROR << "The matrix has not been set for " << Info();
        }

        int err = BaseType::ConsistencyCheck(rX, rY);
        if(err != 0)
            return err;

        TSpaceType::Mult(*mpOperator, rX, rY);
        return 0;
    }

    /// Apply the transpose projection
    int ApplyTranspose(VectorType& rX, VectorType& rY) const override
    {
        if(mpOperator == nullptr)
        {
            KRATOS_ERROR << "The matrix has not been set for " << Info();
        }

        int err = BaseType::ConsistencyCheck(rY, rX);
        if(err != 0) return err;

        TSpaceType::TransposeMult(*mpOperator, rX, rY);
        return 0;
    }

    ///@}
    ///@name Access
    ///@{

    MatrixPointerType GetOperator() const
    {
        return mpOperator;
    }

    ///@}
    ///@name Inquiry
    ///@{

    /// Get the size of the base space
    SizeType GetBaseSize() const override
    {
        return TSpaceType::Size2(*mpOperator);
    }

    /// Get the size of the projected space
    SizeType GetProjectedSize() const override
    {
        return TSpaceType::Size1(*mpOperator);
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream ss;
        ss << "MatrixBasedMGProjector, size = ("
           << TSpaceType::Size1(*mpOperator) << ", "
           << TSpaceType::Size2(*mpOperator) << ")";
        return ss.str();
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

    MatrixPointerType mpOperator;

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


///@}

} // namespace Kratos.

#endif // KRATOS_MULTIGRID_SOLVERS_APP_MATRIX_BASED_MG_PROJECTOR_H_INCLUDED  defined
