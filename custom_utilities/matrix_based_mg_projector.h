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

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MatrixBasedMGProjector()
    {
        this->Initialize();
    }

    /// Destructor.
    virtual ~MatrixBasedMGProjector()
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
        this.mpOperator = rOther.mpOperator;
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    virtual void Apply(VectorType& rX, VectorType& rY) const
    {
        if(mpOperator == NULL)
        {
            std::stringstream ss;
            ss << "The matrix has not been set for " << Info();
            KRATOS_THROW_ERROR(std::logic_error, ss.str(), "");
        }
        TSpaceType::Mult(*mpOperator, rX, rY);
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


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream ss;
        ss << "MatrixBasedMGProjector";
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

    MatrixPointerType mpOperator;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void Initialize()
    {
        if(mpOperator == NULL)
        {
            MatrixPointerType pNewP = MatrixPointerType(new MatrixType(0, 0));
            mpOperator.swap(pNewP);
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
template<class TSpaceType>
inline std::istream& operator >> (std::istream& IStream, MatrixBasedMGProjector<TSpaceType>& rThis)
{
    return IStream;
}

/// output stream function
template<class TSpaceType>
inline std::ostream& operator << (std::ostream& rOStream, const MatrixBasedMGProjector<TSpaceType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


} // namespace Kratos.

#endif // KRATOS_MULTIGRID_SOLVERS_APP_MATRIX_BASED_MG_PROJECTOR_H_INCLUDED  defined

