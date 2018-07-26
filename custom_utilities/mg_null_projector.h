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
//   Date:                $Date: 25/7/2018 $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(KRATOS_MULTIGRID_SOLVERS_APP_MG_NULL_PROJECTOR_H_INCLUDED )
#define  KRATOS_MULTIGRID_SOLVERS_APP_MG_NULL_PROJECTOR_H_INCLUDED



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
 * Projector does nothing and does not return any error. Be careful to use it.
 */
template<class TSpaceType>
class MGNullProjector : public MGProjector<TSpaceType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MGNullProjector
    KRATOS_CLASS_POINTER_DEFINITION(MGNullProjector);

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
    MGNullProjector() : BaseType()
    {}

    /// Destructor.
    virtual ~MGNullProjector()
    {}

    /// Copy constructor
    MGNullProjector(const MGNullProjector& rOther)
    : BaseType(rOther)
    {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator. It's also important like the Copy constructor
    MGNullProjector& operator= (const MGNullProjector& rOther)
    {
        BaseType::operator=(rOther);
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /// Initialize the operator
    virtual void Initialize()
    {}

    /// Apply the projection, rX: input, rY: output
    virtual int Apply(VectorType& rX, VectorType& rY) const
    {
        return 0;
    }

    /// Apply the transpose of the projection, rX: input, rY: output
    /// It is noted that the GetBaseSize() and GetProjectedSize() is only applied for Apply operation. For ApplyTranspose it is reversed.
    virtual int ApplyTranspose(VectorType& rX, VectorType& rY) const
    {
        return 0;
    }

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{

    /// Get the size of the base space
    virtual SizeType GetBaseSize() const
    {
        return 0;
    }

    /// Get the size of the projected space
    virtual SizeType GetProjectedSize() const
    {
        return 0;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream ss;
        ss << "MGNullProjector";
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
template<class TSpaceType>
inline std::istream& operator >> (std::istream& IStream, MGNullProjector<TSpaceType>& rThis)
{
    return IStream;
}

/// output stream function
template<class TSpaceType>
inline std::ostream& operator << (std::ostream& rOStream, const MGNullProjector<TSpaceType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


} // namespace Kratos.

#endif // KRATOS_MULTIGRID_SOLVERS_APP_MG_NULL_PROJECTOR_H_INCLUDED  defined

