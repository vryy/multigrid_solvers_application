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


#if !defined(KRATOS_MULTIGRID_SOLVERS_APP_MG_TRANSPOSE_PROJECTOR_H_INCLUDED )
#define  KRATOS_MULTIGRID_SOLVERS_APP_MG_TRANSPOSE_PROJECTOR_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <sstream>
#include <cstddef>


// External includes


// Project includes
#include "includes/define.h"
#include "custom_utilities/multi_index.h"
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
 * Abstract class for transpose of the projector; can be used for both prolongation and restriction operator.
 */
template<class TSpaceType>
class MGTransposeProjector : public MGProjector<TSpaceType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MGTransposeProjector
    KRATOS_CLASS_POINTER_DEFINITION(MGTransposeProjector);

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
    MGTransposeProjector(typename BaseType::Pointer pProjector)
    : BaseType(), mpProjector(pProjector)
    {}

    /// Destructor.
    virtual ~MGTransposeProjector()
    {}

    /// Copy constructor
    MGTransposeProjector(const MGTransposeProjector& rOther)
    : BaseType(rOther), mpProjector(rOther.mpProjector)
    {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator. It's also important like the Copy constructor
    MGTransposeProjector& operator= (const MGTransposeProjector& rOther)
    {
        BaseType::operator=(rOther);
        mpProjector = rOther.mpProjector;
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /// Initialize the operator
    virtual void Initialize()
    {
        mpProjector->Initialize();
    }

    /// Apply the projection, rX: input, rY: output
    virtual int Apply(VectorType& rX, VectorType& rY) const
    {
        return mpProjector->ApplyTranspose(rX, rY);
    }

    /// Apply the transpose of the projection, rX: input, rY: output
    virtual int ApplyTranspose(VectorType& rX, VectorType& rY) const
    {
        return mpProjector->Apply(rX, rY);
    }

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{

    /// Get the underlying projector
    typename BaseType::Pointer pProjector() const {return mpProjector;}

    /// Get the size of the base space
    virtual SizeType GetBaseSize() const
    {
        return mpProjector->GetProjectedSize();
    }

    /// Get the size of the projected space
    virtual SizeType GetProjectedSize() const
    {
        return mpProjector->GetBaseSize();
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream ss;
        ss << "MGTransposeProjector<" << mpProjector->Info() << ">";
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
        mpProjector->PrintData(rOStream);
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

    int ConsistencyCheck(VectorType& rX, VectorType& rY) const
    {
        return mpProjector->ConsistencyCheck(rY, rX);
    }

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

    typename BaseType::Pointer mpProjector;

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
inline std::istream& operator >> (std::istream& IStream, MGTransposeProjector<TSpaceType>& rThis)
{
    return IStream;
}

/// output stream function
template<class TSpaceType>
inline std::ostream& operator << (std::ostream& rOStream, const MGTransposeProjector<TSpaceType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


} // namespace Kratos.

#endif // KRATOS_MULTIGRID_SOLVERS_APP_TRANSPOSE_PROJECTOR_H_INCLUDED  defined

