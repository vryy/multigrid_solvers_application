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


#if !defined(KRATOS_MULTIGRID_SOLVERS_APP_MG_PROJECTOR_H_INCLUDED )
#define  KRATOS_MULTIGRID_SOLVERS_APP_MG_PROJECTOR_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <sstream>
#include <cstddef>


// External includes


// Project includes
#include "includes/define.h"


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
 * Abstract class for prolongator and restrictor
 */
template<class TSpaceType>
class MGProjector
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MGProjector
    KRATOS_CLASS_POINTER_DEFINITION(MGProjector);

    typedef typename TSpaceType::MatrixType MatrixType;

    typedef typename TSpaceType::MatrixPointerType MatrixPointerType;

    typedef typename TSpaceType::VectorType VectorType;

    typedef typename TSpaceType::VectorPointerType VectorPointerType;

    typedef typename TSpaceType::SizeType SizeType;

    typedef typename TSpaceType::IndexType IndexType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MGProjector()
    {}

    /// Destructor.
    virtual ~MGProjector()
    {}

    /// Copy constructor
    MGProjector(const MGProjector& rOther)
    {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator. It's also important like the Copy constructor
    MGProjector& operator= (const MGProjector& rOther)
    {
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
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__);
        return 0;
    }

    /// Apply the transpose of the projection, rX: input, rY: output
    /// It is noted that the GetBaseSize() and GetProjectedSize() is only applied for Apply operation. For ApplyTranspose it is reversed.
    virtual int ApplyTranspose(VectorType& rX, VectorType& rY) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__);
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
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__);
    }

    /// Get the size of the projected space
    virtual SizeType GetProjectedSize() const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__);
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream ss;
        ss << "MGProjector";
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

    int ConsistencyCheck(VectorType& rX, VectorType& rY) const
    {
        if(TSpaceType::Size(rX) != this->GetBaseSize())
            return 1;

        if(TSpaceType::Size(rY) != this->GetProjectedSize())
            return 2;

        return 0;
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
inline std::istream& operator >> (std::istream& IStream, MGProjector<TSpaceType>& rThis)
{
    return IStream;
}

/// output stream function
template<class TSpaceType>
inline std::ostream& operator << (std::ostream& rOStream, const MGProjector<TSpaceType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


} // namespace Kratos.

#endif // KRATOS_MULTIGRID_SOLVERS_APP_MG_PROJECTOR_H_INCLUDED  defined

