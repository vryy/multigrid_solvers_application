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
//   Date:                $Date: 2013 Jan 11 15:06:00 $
//   Revision:            $Revision: 1.1 $
//
//


#if !defined(KRATOS_MULTIGRID_SOLVERS_APP_MULTILEVEL_SOLVER_FACTORY_H_INCLUDED )
#define  KRATOS_MULTIGRID_SOLVERS_APP_MULTILEVEL_SOLVER_FACTORY_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <sstream>
#include <cstddef>


// External includes
//#include "external_includes/pyamg/relaxation.h"
//#include "external_includes/pyamg/ruge_stuben.h"


// Project includes
#include "includes/define.h"
#include "linear_solvers/linear_solver.h"
#include "custom_utilities/parameter_list.h"
#include "custom_linear_solvers/multilevel_solver.h"


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
class MultilevelSolverFactory
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MultilevelSolverFactory
    KRATOS_CLASS_POINTER_DEFINITION(MultilevelSolverFactory);

    typedef MultilevelSolver<TSparseSpaceType, TDenseSpaceType> MultilevelSolverType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::MatrixPointerType SparseMatrixPointerType;

    typedef ParameterList<std::string> ParameterListType;

    typedef typename TSparseSpaceType::SizeType SizeType;

    typedef typename TSparseSpaceType::IndexType IndexType;

    typedef int  IntegerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MultilevelSolverFactory(ParameterListType& mg_parameter_list)
    : mmg_parameter_list(mg_parameter_list), mIsMute(false)
    {}


    /// Copy constructor.
    MultilevelSolverFactory(const MultilevelSolverFactory& rOther)
    : mmg_parameter_list(rOther.mmg_parameter_list), mIsMute(rOther.mIsMute)
    {}

    /// Destructor.
    virtual ~MultilevelSolverFactory()
    {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    MultilevelSolverFactory& operator=(const MultilevelSolverFactory& rOther)
    {
        mmg_parameter_list = rOther.mmg_parameter_list;
        mIsMute = rOther.mIsMute;
        return *this;
    }


    ///@}
    ///@name Operations
    ///@{

    void SetMute(const bool& mute) {mIsMute = mute;}

    virtual void InitializeMultilevelSolver(MultilevelSolverType& solver) const
    {
        if (!mIsMute)
            KRATOS_THROW_ERROR(std::logic_error, "Calling the base class function", __FUNCTION__);
    }

    virtual void GenerateMultilevelSolver(MultilevelSolverType& solver, SparseMatrixType& rA) const
    {
        if (!mIsMute)
            KRATOS_THROW_ERROR(std::logic_error, "Calling the base class function", __FUNCTION__);
    }

    ///@}
    ///@name Access
    ///@{



    ///@}
    ///@name Inquiry
    ///@{



    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Multilevel solver factory";
        return  buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << mmg_parameter_list;
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

    ParameterListType mmg_parameter_list;

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

    bool mIsMute;

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

}; // Class MultilevelSolverFactory

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TSparseSpaceType, class TDenseSpaceType>
inline std::istream& operator >> (std::istream& IStream, MultilevelSolverFactory<TSparseSpaceType, TDenseSpaceType>& rThis)
{
    return IStream;
}

/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType>
inline std::ostream& operator << (std::ostream& rOStream, const MultilevelSolverFactory<TSparseSpaceType, TDenseSpaceType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_MULTILEVEL_SOLVER_FACTORY_H_INCLUDED  defined


