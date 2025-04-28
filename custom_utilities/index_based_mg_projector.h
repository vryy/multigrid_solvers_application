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
//   Date:                $Date: 24/7/2018 $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(KRATOS_MULTIGRID_SOLVERS_APP_INDEX_BASED_MG_PROJECTOR_H_INCLUDED )
#define  KRATOS_MULTIGRID_SOLVERS_APP_INDEX_BASED_MG_PROJECTOR_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <sstream>
#include <cstddef>


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
 * Projector based on indexed entries. It is useful for matrix-based refinement, e.g. isogeometric alalysis.
 */
template<class TSpaceType>
class IndexBasedMGProjector : public MGProjector<TSpaceType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of IndexBasedMGProjector
    KRATOS_CLASS_POINTER_DEFINITION(IndexBasedMGProjector);

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

    /// Default constructor.
    IndexBasedMGProjector(const SizeType& nRows, const SizeType& nCols)
    : BaseType(), m_nrows(nRows), m_ncols(nCols), m_stride(1)
    {}

    /// Destructor.
    ~IndexBasedMGProjector() override
    {}

    /// Copy constructor
    IndexBasedMGProjector(const IndexBasedMGProjector& rOther)
    : BaseType(rOther)
    , m_nrows(rOther.m_nrows), m_ncols(rOther.m_ncols)
    , m_stride(rOther.m_stride)
    {}

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator. It's also important like the Copy constructor
    IndexBasedMGProjector& operator= (const IndexBasedMGProjector& rOther)
    {
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /// Initialize the operator
    void Initialize() override
    {
        m_row_indices.resize(m_nrows);
        m_row_values.resize(m_nrows);

        m_column_indices.resize(m_ncols);
        m_column_values.resize(m_ncols);
    }

    /// Set the stride, i.e. block size
    void SetStride(const std::size_t& stride)
    {
        m_stride = stride;
    }

    /// populate the transformation matrix
    void AssembleOperator(const std::vector<std::size_t>& rows,
        const std::vector<std::size_t>& columns, const MatrixType& values)
    {
        // size check
        KRATOS_WATCH(rows.size())
        KRATOS_WATCH(columns.size())
        KRATOS_WATCH(values.size1())
        KRATOS_WATCH(values.size2())
        if (rows.size() != TSpaceType::Size1(values))
            KRATOS_THROW_ERROR(std::logic_error, "The row size of the matrix is incompatible", __FUNCTION__)

        if (columns.size() != TSpaceType::Size2(values))
            KRATOS_THROW_ERROR(std::logic_error, "The column size of the matrix is incompatible", __FUNCTION__)

        // populate values
        for (std::size_t i = 0; i < rows.size(); ++i)
        {
            for (std::size_t j = 0; j < columns.size(); ++j)
            {
                if (values(i, j) != 0.0)
                {
                    m_row_indices[rows[i]].push_back(columns[j]);
                    m_row_values[rows[i]].push_back(values(i, j));

                    m_column_indices[columns[j]].push_back(rows[i]);
                    m_column_values[columns[j]].push_back(values(i, j));
                }
            }
        }
    }

    /// Apply the projection, rX: input, rY: output
    int Apply(VectorType& rX, VectorType& rY) const override
    {
        int err = BaseType::ConsistencyCheck(rX, rY);
        if(err != 0) return err;

        TSpaceType::SetToZero(rY);

        for (std::size_t i = 0; i < m_nrows; ++i)
        {
            for (std::size_t j = 0; j < m_row_indices[i].size(); ++j)
            {
                for (std::size_t s = 0; s < m_stride; ++s)
                    rY[i*m_stride+s] += m_row_values[i][j] * rX[m_row_indices[i][j]*m_stride+s];
            }
        }

        return 0;
    }

    /// Apply the transpose of the projection, rX: input, rY: output
    /// It is noted that the GetBaseSize() and GetProjectedSize() is only applied for Apply operation. For ApplyTranspose it is reversed.
    int ApplyTranspose(VectorType& rX, VectorType& rY) const override
    {
        int err = BaseType::ConsistencyCheck(rY, rX);
        if(err != 0) return err;

        TSpaceType::SetToZero(rY);

        for (std::size_t i = 0; i < m_ncols; ++i)
        {
            for (std::size_t j = 0; j < m_column_indices[i].size(); ++j)
            {
                for (std::size_t s = 0; s < m_stride; ++s)
                    rY[i*m_stride+s] += m_column_values[i][j] * rX[m_column_indices[i][j]*m_stride+s];
            }
        }

        return 0;
    }

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{

    /// Get the size of the base space
    SizeType GetBaseSize() const override
    {
        return m_ncols*m_stride;
    }

    /// Get the size of the projected space
    SizeType GetProjectedSize() const override
    {
        return m_nrows*m_stride;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream ss;
        ss << "IndexBasedMGProjector(" << m_nrows << ", " << m_ncols << ")";
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

    // arrays to store the entries of the transformation matrix, row by column
    SizeType m_nrows;
    std::vector<std::vector<IndexType> > m_row_indices;
    std::vector<std::vector<double> > m_row_values;

    // arrays to store the entries of the transformation matrix, column by row
    SizeType m_ncols;
    std::vector<std::vector<IndexType> > m_column_indices;
    std::vector<std::vector<double> > m_column_values;

    std::size_t m_stride;

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

#endif // KRATOS_MULTIGRID_SOLVERS_APP_INDEX_BASED_MG_PROJECTOR_H_INCLUDED  defined
